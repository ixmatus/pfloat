//! Conversions between the IEEE 754 hardware floats (`f32`, `f64`)
//! and [`BigFloat`].
//!
//! Two directions, both pure Rust and `no_std`:
//!
//! - **Widening** ([`BigFloat::from_f32`], [`BigFloat::from_f64`]) is
//!   exact. A hardware float is a dyadic rational with at most 24
//!   (`f32`) or 53 (`f64`) significant bits, so it is represented
//!   without loss in a `BigFloat` at that precision. NaN, the
//!   infinities, and signed zero map across directly.
//!
//! - **Narrowing** ([`BigFloat::to_f32_round`],
//!   [`BigFloat::to_f64_round`]) rounds an arbitrary-precision value
//!   to the target format under an explicit [`RoundingMode`], pairing
//!   the result with a [`Status`] that carries the IEEE 754-2019
//!   sticky flags the conversion raised (`INEXACT`, `OVERFLOW`,
//!   `UNDERFLOW`, and `INVALID` for a signaling-NaN operand).
//!
//! Narrowing is the `BigFloat → float` step a libm shell performs on
//! every call, and it is a known double-rounding hazard: rounding to
//! an intermediate width and then to the format can land on the wrong
//! neighbour in the subnormal range. This implementation avoids that
//! by rounding straight to the format's grid. The mode-aware rounding
//! reuses pfloat's verified [`BigFloat::round_to_precision`] with the
//! precision the format affords at the value's magnitude (the full
//! significand for normals, the reduced significand the fixed
//! exponent floor leaves for subnormals); the result then lands
//! exactly on the format grid and its IEEE fields are read directly
//! from the mantissa limb. No decimal round trip and no second
//! rounding occur.

#![cfg(feature = "big")]

use crate::big::BigFloat;
use crate::big::Parts;
use crate::rounding::RoundingMode;
use crate::sign::Sign;
use crate::status::Status;

/// Static description of an IEEE 754 binary interchange format,
/// shared by the `f32` and `f64` conversion paths.
struct Format {
    /// Significand precision including the implicit bit (24 / 53).
    prec: u32,
    /// Stored fraction-field width in bits (23 / 52).
    mant_bits: u32,
    /// Exponent bias (127 / 1023).
    bias: i64,
    /// Minimum normal exponent `emin` (-126 / -1022).
    emin: i64,
    /// Maximum finite exponent `emax` (127 / 1023).
    emax: i64,
    /// Exponent of the smallest positive subnormal (-149 / -1074).
    sub_min_exp: i64,
    /// Total format width in bits (32 / 64).
    width: u32,
}

const F32: Format = Format {
    prec: 24,
    mant_bits: 23,
    bias: 127,
    emin: -126,
    emax: 127,
    sub_min_exp: -149,
    width: 32,
};

const F64: Format = Format {
    prec: 53,
    mant_bits: 52,
    bias: 1023,
    emin: -1022,
    emax: 1023,
    sub_min_exp: -1074,
    width: 64,
};

impl Format {
    /// Raw bit pattern of `±0`.
    fn zero(&self, neg: bool) -> u64 {
        self.sign(neg)
    }

    /// The sign bit positioned at the top of the field.
    fn sign(&self, neg: bool) -> u64 {
        if neg {
            1u64 << (self.width - 1)
        } else {
            0
        }
    }

    /// Raw bit pattern of `±∞`.
    fn infinity(&self, neg: bool) -> u64 {
        let exp_all_ones = (1u64 << (self.width - 1 - self.mant_bits)) - 1;
        self.sign(neg) | (exp_all_ones << self.mant_bits)
    }

    /// Raw bit pattern of the largest finite magnitude, `±MAX`.
    fn max_finite(&self, neg: bool) -> u64 {
        let exp_field = (self.emax + self.bias) as u64;
        let frac = (1u64 << self.mant_bits) - 1;
        self.sign(neg) | (exp_field << self.mant_bits) | frac
    }

    /// Raw bit pattern of a quiet NaN with the given sign.
    fn quiet_nan(&self, neg: bool) -> u64 {
        let exp_all_ones = (1u64 << (self.width - 1 - self.mant_bits)) - 1;
        let quiet_bit = 1u64 << (self.mant_bits - 1);
        self.sign(neg) | (exp_all_ones << self.mant_bits) | quiet_bit
    }
}

/// Build a `BigFloat` holding the exact value of an IEEE float, given
/// its decomposed bit fields. The value is dyadic, so the result is
/// exact at the format's significand precision.
fn from_ieee(sign_neg: bool, exp_field: u32, mant_field: u64, fmt: &Format) -> BigFloat {
    let sign = if sign_neg {
        Sign::Negative
    } else {
        Sign::Positive
    };
    let exp_all_ones = (1u32 << (fmt.width - 1 - fmt.mant_bits)) - 1;

    if exp_field == 0 && mant_field == 0 {
        return BigFloat::try_new_zero(sign, fmt.prec).expect("precision valid");
    }
    if exp_field == exp_all_ones {
        if mant_field == 0 {
            return BigFloat::try_new_infinity(sign, fmt.prec).expect("precision valid");
        }
        return BigFloat::try_new_quiet_nan(sign, fmt.prec, &[]).expect("precision valid");
    }

    // (integer mantissa, binary scale): value = mantissa * 2^scale.
    let (mantissa_int, scale) = if exp_field == 0 {
        // Subnormal: value = mant_field * 2^sub_min_exp.
        (mant_field as i64, fmt.sub_min_exp)
    } else {
        // Normal: value = (2^mant_bits | mant_field) * 2^(e - mant_bits)
        // with the unbiased exponent e = exp_field - bias.
        let significand = (1i64 << fmt.mant_bits) | mant_field as i64;
        let e = exp_field as i64 - fmt.bias;
        (significand, e - i64::from(fmt.mant_bits))
    };
    let signed = if sign_neg {
        -mantissa_int
    } else {
        mantissa_int
    };

    let mut bf = BigFloat::try_from_i64_exact(signed, fmt.prec).expect("significand fits in i64");

    // Apply 2^scale exactly via chained multiply or divide by powers
    // of two; each step is an exact binary-exponent shift in BigFloat
    // arithmetic, so the mantissa is preserved.
    let up = scale >= 0;
    let mut remaining = scale.unsigned_abs();
    const CHUNK: u32 = 32;
    let chunk = BigFloat::try_from_i64_exact(1i64 << CHUNK, fmt.prec).expect("2^32 fits");
    while remaining >= u64::from(CHUNK) {
        bf = step(&bf, &chunk, up);
        remaining -= u64::from(CHUNK);
    }
    if remaining > 0 {
        let rem = BigFloat::try_from_i64_exact(1i64 << remaining, fmt.prec).expect("power fits");
        bf = step(&bf, &rem, up);
    }
    bf
}

fn step(bf: &BigFloat, power: &BigFloat, up: bool) -> BigFloat {
    if up {
        bf.mul(power, RoundingMode::NearestEven).0
    } else {
        bf.div(power, RoundingMode::NearestEven).0
    }
}

/// Round `value` to the IEEE format described by `fmt` under `mode`,
/// returning the raw bit pattern and the sticky status flags.
fn to_ieee(value: &BigFloat, mode: RoundingMode, fmt: &Format) -> (u64, Status) {
    match value.parts() {
        Parts::Zero { sign } => (fmt.zero(is_neg(sign)), Status::OK),
        Parts::Infinity { sign } => (fmt.infinity(is_neg(sign)), Status::OK),
        Parts::Nan { quiet, sign, .. } => {
            // A signaling-NaN operand raises INVALID and is quieted.
            let status = if quiet { Status::OK } else { Status::INVALID };
            (fmt.quiet_nan(is_neg(sign)), status)
        }
        Parts::Normal { sign, exponent, .. } => {
            let neg = is_neg(sign);

            // Above the format's range: overflow before rounding.
            if exponent > fmt.emax {
                return overflow(neg, mode, fmt);
            }
            // Below the smallest subnormal's grid resolution: the value
            // rounds to zero or to the smallest subnormal.
            if exponent < fmt.sub_min_exp {
                return tiny(value, neg, mode, fmt);
            }

            // Significand precision the format affords at this
            // magnitude: full for normals, reduced toward the fixed
            // subnormal floor below emin.
            let p_eff = if exponent >= fmt.emin {
                fmt.prec
            } else {
                (exponent - fmt.sub_min_exp + 1) as u32
            };
            let (rounded, mut status) = value
                .round_to_precision(p_eff, mode)
                .expect("p_eff >= 1 for exponent >= sub_min_exp");

            match rounded.parts() {
                Parts::Normal {
                    exponent: er,
                    mantissa,
                    ..
                } => {
                    // A round-up across the top binade overflows.
                    if er > fmt.emax {
                        let (bits, ovf) = overflow(neg, mode, fmt);
                        return (bits, status | ovf);
                    }
                    let bits = encode(neg, er, mantissa[0], fmt);
                    // IEEE underflow: tiny (subnormal) and inexact.
                    if er < fmt.emin && status.inexact() {
                        status |= Status::UNDERFLOW;
                    }
                    (bits, status)
                }
                // round_to_precision of a non-zero value to >= 1 bit is
                // non-zero, so the remaining Parts are unreachable.
                _ => unreachable!("rounded a finite non-zero value to a non-normal result"),
            }
        }
    }
}

/// Encode an on-grid value (exponent `er`, single mantissa limb
/// `limb`, left-aligned with its leading bit at bit 63) as the raw
/// bits of the format. `er` is assumed within `[sub_min_exp, emax]`.
fn encode(neg: bool, er: i64, limb: u64, fmt: &Format) -> u64 {
    let sign = fmt.sign(neg);
    if er >= fmt.emin {
        // Normal: drop the implicit leading bit, keep `mant_bits` of
        // fraction from the top of the limb.
        let frac = (limb >> (64 - fmt.prec)) & ((1u64 << fmt.mant_bits) - 1);
        let exp_field = (er + fmt.bias) as u64;
        sign | (exp_field << fmt.mant_bits) | frac
    } else {
        // Subnormal: the leading bit sits at 2^er; align the limb so
        // that bit lands at position (er - sub_min_exp) of the field.
        let shift = 63 - er + fmt.sub_min_exp;
        let field = limb >> shift;
        sign | field
    }
}

/// Result for a magnitude above the format's largest finite value.
fn overflow(neg: bool, mode: RoundingMode, fmt: &Format) -> (u64, Status) {
    let to_inf = match mode {
        RoundingMode::NearestEven | RoundingMode::NearestAway => true,
        RoundingMode::TowardZero => false,
        RoundingMode::TowardPositive => !neg,
        RoundingMode::TowardNegative => neg,
    };
    let bits = if to_inf {
        fmt.infinity(neg)
    } else {
        fmt.max_finite(neg)
    };
    (bits, Status::OVERFLOW | Status::INEXACT)
}

/// Result for a non-zero magnitude below the smallest subnormal: the
/// value rounds either to zero or to the smallest subnormal.
fn tiny(value: &BigFloat, neg: bool, mode: RoundingMode, fmt: &Format) -> (u64, Status) {
    let smin = smallest_subnormal(fmt);
    let half = BigFloat::try_from_i64_exact(1, fmt.prec)
        .expect("precision valid")
        .div(
            &BigFloat::try_from_i64_exact(2, fmt.prec).expect("precision valid"),
            RoundingMode::NearestEven,
        )
        .0
        .mul(&smin, RoundingMode::NearestEven)
        .0;
    let magnitude = abs(value);

    // Order the magnitude against half the smallest subnormal.
    let up = match mode {
        RoundingMode::TowardZero => false,
        RoundingMode::TowardPositive => !neg,
        RoundingMode::TowardNegative => neg,
        RoundingMode::NearestEven => {
            // Above half rounds up; exactly half ties to even (zero).
            greater(&magnitude, &half)
        }
        RoundingMode::NearestAway => {
            // At or above half rounds up.
            !less(&magnitude, &half)
        }
    };
    let bits = if up { fmt.sign(neg) | 1 } else { fmt.zero(neg) };
    (bits, Status::UNDERFLOW | Status::INEXACT)
}

fn smallest_subnormal(fmt: &Format) -> BigFloat {
    from_ieee(false, 0, 1, fmt)
}

fn abs(value: &BigFloat) -> BigFloat {
    match value.parts() {
        Parts::Normal {
            sign: Sign::Negative,
            ..
        } => {
            value
                .mul(
                    &BigFloat::try_from_i64_exact(-1, value.precision()).expect("precision valid"),
                    RoundingMode::NearestEven,
                )
                .0
        }
        _ => value.clone(),
    }
}

fn greater(a: &BigFloat, b: &BigFloat) -> bool {
    matches!(a.partial_cmp(b).0, Some(core::cmp::Ordering::Greater))
}

fn less(a: &BigFloat, b: &BigFloat) -> bool {
    matches!(a.partial_cmp(b).0, Some(core::cmp::Ordering::Less))
}

fn is_neg(sign: Sign) -> bool {
    matches!(sign, Sign::Negative)
}

impl BigFloat {
    /// Widen an [`f32`] to an exact `BigFloat` at precision 24.
    ///
    /// Every `f32` is a dyadic rational with at most 24 significant
    /// bits, so the conversion is lossless. NaN, the infinities, and
    /// signed zero map across directly. To compute at a wider working
    /// precision, follow with
    /// [`round_to_precision`](BigFloat::round_to_precision) (which is
    /// exact when widening).
    pub fn from_f32(x: f32) -> BigFloat {
        let bits = x.to_bits();
        from_ieee(
            (bits >> 31) & 1 == 1,
            (bits >> 23) & 0xFF,
            u64::from(bits & 0x007F_FFFF),
            &F32,
        )
    }

    /// Widen an [`f64`] to an exact `BigFloat` at precision 53.
    ///
    /// Lossless for the same reason [`from_f32`](BigFloat::from_f32)
    /// is: an `f64` carries at most 53 significant bits.
    pub fn from_f64(x: f64) -> BigFloat {
        let bits = x.to_bits();
        from_ieee(
            (bits >> 63) & 1 == 1,
            ((bits >> 52) & 0x7FF) as u32,
            bits & 0x000F_FFFF_FFFF_FFFF,
            &F64,
        )
    }

    /// Round to [`f32`] under `mode`, returning the value paired with
    /// the IEEE 754-2019 sticky flags the conversion raised
    /// (`INEXACT`, `OVERFLOW`, `UNDERFLOW`, and `INVALID` for a
    /// signaling-NaN operand).
    ///
    /// The rounding lands directly on the `binary32` grid, so no
    /// double rounding occurs even in the subnormal range.
    pub fn to_f32_round(&self, mode: RoundingMode) -> (f32, Status) {
        let (bits, status) = to_ieee(self, mode, &F32);
        (f32::from_bits(bits as u32), status)
    }

    /// Round to [`f64`] under `mode`. The `f64` companion to
    /// [`to_f32_round`](BigFloat::to_f32_round); the same flag and
    /// grid contract holds.
    pub fn to_f64_round(&self, mode: RoundingMode) -> (f64, Status) {
        let (bits, status) = to_ieee(self, mode, &F64);
        (f64::from_bits(bits), status)
    }

    /// Round to [`f32`] under [`RoundingMode::NearestEven`], discarding
    /// the status. Convenience over
    /// [`to_f32_round`](BigFloat::to_f32_round) for the common case.
    pub fn to_f32(&self) -> f32 {
        self.to_f32_round(RoundingMode::NearestEven).0
    }

    /// Round to [`f64`] under [`RoundingMode::NearestEven`], discarding
    /// the status.
    pub fn to_f64(&self) -> f64 {
        self.to_f64_round(RoundingMode::NearestEven).0
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    fn rt32(bits: u32) {
        let x = f32::from_bits(bits);
        let y = BigFloat::from_f32(x).to_f32();
        if x.is_nan() {
            assert!(y.is_nan(), "from_f32({x}).to_f32() should be NaN");
        } else {
            assert_eq!(y.to_bits(), bits, "round-trip f32 0x{bits:08x} = {x}");
        }
    }

    fn rt64(bits: u64) {
        let x = f64::from_bits(bits);
        let y = BigFloat::from_f64(x).to_f64();
        if x.is_nan() {
            assert!(y.is_nan(), "from_f64({x}).to_f64() should be NaN");
        } else {
            assert_eq!(y.to_bits(), bits, "round-trip f64 0x{bits:016x} = {x}");
        }
    }

    #[test]
    fn f32_round_trip_boundaries() {
        for bits in [
            0x0000_0000, // +0
            0x8000_0000, // -0
            0x0000_0001, // smallest subnormal
            0x007F_FFFF, // largest subnormal
            0x0080_0000, // smallest normal
            0x7F7F_FFFF, // largest normal (f32::MAX)
            0xFF7F_FFFF, // -f32::MAX
            0x3F80_0000, // 1.0
            0x4000_0000, // 2.0
            0x7F80_0000, // +inf
            0xFF80_0000, // -inf
            0x7FC0_0000, // qNaN
        ] {
            rt32(bits);
        }
    }

    #[test]
    fn f32_round_trip_sampled() {
        let mut bits: u32 = 0;
        loop {
            rt32(bits);
            match bits.checked_add(75_013) {
                Some(n) => bits = n,
                None => break,
            }
        }
    }

    #[test]
    fn f64_round_trip_boundaries() {
        for bits in [
            0x0000_0000_0000_0000, // +0
            0x8000_0000_0000_0000, // -0
            0x0000_0000_0000_0001, // smallest subnormal
            0x000F_FFFF_FFFF_FFFF, // largest subnormal
            0x0010_0000_0000_0000, // smallest normal
            0x7FEF_FFFF_FFFF_FFFF, // f64::MAX
            0x3FF0_0000_0000_0000, // 1.0
            0x7FF0_0000_0000_0000, // +inf
            0x7FF8_0000_0000_0000, // qNaN
        ] {
            rt64(bits);
        }
    }

    #[test]
    fn f64_round_trip_sampled() {
        let mut bits: u64 = 0;
        loop {
            rt64(bits);
            match bits.checked_add(0x000A_3B91_C7D5_F00D) {
                Some(n) => bits = n,
                None => break,
            }
        }
    }

    /// The exact midpoint between 1.0 and the next f32 (1 + 2^-24) is
    /// representable in f64. It rounds per mode: to even (1.0), away
    /// (next), and per direction.
    #[test]
    fn f32_midpoint_rounding() {
        let one = 1.0f32;
        let next = f32::from_bits(one.to_bits() + 1); // 1 + 2^-23
        let mid = BigFloat::from_f64(1.0 + 2f64.powi(-24)); // exact tie

        assert_eq!(mid.to_f32_round(RoundingMode::NearestEven).0, one); // even
        assert_eq!(mid.to_f32_round(RoundingMode::NearestAway).0, next); // away
        assert_eq!(mid.to_f32_round(RoundingMode::TowardZero).0, one);
        assert_eq!(mid.to_f32_round(RoundingMode::TowardPositive).0, next);
        assert_eq!(mid.to_f32_round(RoundingMode::TowardNegative).0, one);
        assert!(mid.to_f32_round(RoundingMode::NearestEven).1.inexact());
    }

    /// A value just above `f32::MAX` overflows to `+inf` under nearest
    /// and directed-up, and saturates to `MAX` under toward-zero.
    #[test]
    fn f32_overflow() {
        let big = BigFloat::from_f64(1e40); // >> f32::MAX
        assert_eq!(big.to_f32_round(RoundingMode::NearestEven).0, f32::INFINITY);
        assert_eq!(
            big.to_f32_round(RoundingMode::TowardZero).0,
            f32::MAX,
            "toward-zero saturates to MAX"
        );
        assert_eq!(
            big.to_f32_round(RoundingMode::TowardNegative).0,
            f32::MAX,
            "positive value, toward -inf saturates to MAX"
        );
        let (v, s) = big.to_f32_round(RoundingMode::NearestEven);
        assert_eq!(v, f32::INFINITY);
        assert!(s.overflow() && s.inexact());
    }

    /// A value below half the smallest subnormal underflows to zero
    /// under nearest; directed-up reaches the smallest subnormal.
    #[test]
    fn f32_tiny_underflow() {
        let tiny = BigFloat::from_f64(2f64.powi(-151)); // < 2^-150 = half ULP
        assert_eq!(tiny.to_f32_round(RoundingMode::NearestEven).0, 0.0);
        let smallest = f32::from_bits(1);
        assert_eq!(tiny.to_f32_round(RoundingMode::TowardPositive).0, smallest);
        let (v, s) = tiny.to_f32_round(RoundingMode::NearestEven);
        assert_eq!(v, 0.0);
        assert!(s.underflow() && s.inexact());
    }
}
