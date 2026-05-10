//! [`BigFloat`]: dynamic-precision arbitrary-precision float.
//!
//! `BigFloat` carries the precision as a runtime field; the mantissa
//! lives in a heap-allocated `Vec<u64>` whose length is
//! `limbs_for(precision)`. See ADR-0001 (limb representation),
//! ADR-0002 (bit-level precision), ADR-0003 (dual API), ADR-0004
//! (storage), and ADR-0005 (`Class` tagged enum).

#[cfg(feature = "alloc")]
extern crate alloc;

use alloc::vec;
use alloc::vec::Vec;

use crate::class::Class;
use crate::mantissa::{limbs_for, storage_shift};
use crate::sign::Sign;

/// Pure-Rust correctly-rounded arbitrary-precision binary float
/// with runtime precision.
///
/// `BigFloat` matches MPFR's shape: the precision is chosen at
/// construction time and may differ between values. Use
/// [`FixedFloat<PREC>`] (slice 1g) when the precision is known at
/// compile time and stack allocation is preferable.
///
/// Slice 1a (this slice) ships construction, classification, sign
/// manipulation, comparison, and exact integer conversion. Slice 1b
/// adds the rounding pipeline; slices 1c–1f add the arithmetic
/// kernels.
///
/// [`FixedFloat<PREC>`]: ../fixed/struct.FixedFloat.html
#[cfg(feature = "big")]
#[derive(Clone, Debug, Eq, PartialEq, Hash)]
pub struct BigFloat {
    pub(crate) class: Class,
    pub(crate) precision: u32,
}

/// Failure modes for [`BigFloat`] construction.
///
/// Slice 1a uses `BuildError` for precision validation and exact
/// integer conversion. Later slices extend it with rounding-related
/// errors as needed.
#[derive(Clone, Copy, Debug, Eq, PartialEq, Hash)]
pub enum BuildError {
    /// `precision` was `0`. Precisions must be at least one bit per
    /// ADR-0002.
    PrecisionZero,

    /// The value's significant-bit count exceeds the requested
    /// precision; representing it would require rounding, which is
    /// not in scope for the `*_exact` constructors.
    ///
    /// `value_bits` is the count of bits between the value's most
    /// significant `1` and its least significant `1` (inclusive),
    /// i.e. `bit_length - trailing_zeros`. `requested` is the
    /// requested precision. The condition for `Ok` is
    /// `value_bits <= requested`.
    ValueExceedsPrecision { value_bits: u32, requested: u32 },
}

#[cfg(feature = "big")]
impl BigFloat {
    /// Constructs `±0` at the given precision.
    ///
    /// Returns [`BuildError::PrecisionZero`] when `precision == 0`.
    ///
    /// IEEE 754-2019 §6.3 distinguishes `+0` and `-0`. The sign
    /// bit propagates through arithmetic (e.g. `(-0) + (-0) == -0`,
    /// `(-0) + (+0) == +0` under round-to-nearest).
    pub fn try_new_zero(sign: Sign, precision: u32) -> Result<Self, BuildError> {
        validate_precision(precision)?;
        Ok(Self {
            class: Class::Zero { sign },
            precision,
        })
    }

    /// Constructs `±∞` at the given precision.
    ///
    /// Returns [`BuildError::PrecisionZero`] when `precision == 0`.
    pub fn try_new_infinity(sign: Sign, precision: u32) -> Result<Self, BuildError> {
        validate_precision(precision)?;
        Ok(Self {
            class: Class::Infinity { sign },
            precision,
        })
    }

    /// Constructs a quiet NaN with the given sign and (optional)
    /// payload.
    ///
    /// `payload`'s contents are copied into a `Vec<u64>` of size
    /// `limbs_for(precision)`, zero-padded if shorter and truncated
    /// (with a high-bit clear of any meaning) if longer. Most
    /// callers pass an empty slice for "no diagnostic info."
    ///
    /// IEEE 754-2019 §6.2 distinguishes quiet and signaling NaNs.
    /// Quiet NaNs propagate through arithmetic without raising
    /// [`Status::INVALID`](crate::status::Status::INVALID); use
    /// [`try_new_signaling_nan`](Self::try_new_signaling_nan) for the
    /// trapping variant.
    ///
    /// Returns [`BuildError::PrecisionZero`] when `precision == 0`.
    pub fn try_new_quiet_nan(
        sign: Sign,
        precision: u32,
        payload: &[u64],
    ) -> Result<Self, BuildError> {
        validate_precision(precision)?;
        Ok(Self {
            class: Class::Nan {
                quiet: true,
                sign,
                payload: pad_payload(payload, limbs_for(precision)),
            },
            precision,
        })
    }

    /// Constructs a signaling NaN with the given sign and payload.
    ///
    /// See [`try_new_quiet_nan`](Self::try_new_quiet_nan) for the
    /// payload semantics.
    ///
    /// Returns [`BuildError::PrecisionZero`] when `precision == 0`.
    pub fn try_new_signaling_nan(
        sign: Sign,
        precision: u32,
        payload: &[u64],
    ) -> Result<Self, BuildError> {
        validate_precision(precision)?;
        Ok(Self {
            class: Class::Nan {
                quiet: false,
                sign,
                payload: pad_payload(payload, limbs_for(precision)),
            },
            precision,
        })
    }

    /// Constructs a `BigFloat` from an `i64` exactly, at the given
    /// precision.
    ///
    /// Returns [`BuildError::ValueExceedsPrecision`] when the
    /// integer's significant-bit count (the bits between its top
    /// `1` and bottom `1`, inclusive) exceeds the precision. The
    /// rounding-required path lands in slice 1b.
    ///
    /// `n == 0` returns `+0` regardless of precision (precision
    /// must still be `>= 1`).
    pub fn try_from_i64_exact(n: i64, precision: u32) -> Result<Self, BuildError> {
        validate_precision(precision)?;

        if n == 0 {
            return Ok(Self {
                class: Class::Zero {
                    sign: Sign::Positive,
                },
                precision,
            });
        }

        let sign = if n < 0 {
            Sign::Negative
        } else {
            Sign::Positive
        };
        let magnitude: u64 = if n == i64::MIN {
            1u64 << 63
        } else {
            n.unsigned_abs()
        };

        let trailing = magnitude.trailing_zeros();
        let total_bits = u64::BITS - magnitude.leading_zeros();
        let significant = total_bits - trailing;

        if significant > precision {
            return Err(BuildError::ValueExceedsPrecision {
                value_bits: significant,
                requested: precision,
            });
        }

        let limbs = limbs_for(precision);
        let mut storage: Vec<u64> = vec![0; limbs];

        // Reduced magnitude has its top bit at position
        // `significant - 1` and zero trailing zeros (we shifted them
        // out).
        let reduced: u64 = magnitude >> trailing;

        // Place `reduced` into the storage so its top bit lands at
        // bit 63 of the most-significant limb (the storage's MSB).
        let total_shift = storage_shift(limbs, significant);
        let whole = (total_shift / 64) as usize;
        let intra = total_shift % 64;

        if intra == 0 {
            storage[whole] = reduced;
        } else {
            storage[whole] = reduced << intra;
            // The high bits of `reduced` that did not fit in the
            // first limb spill into the next one. They only
            // exist when `reduced` has more bits than `64 - intra`.
            if 64 - intra < significant && whole + 1 < limbs {
                storage[whole + 1] = reduced >> (64 - intra);
            }
        }

        // Exponent: see derivation in `mantissa.rs` storage notes
        // and `DESIGN.md`'s "Numeric representation" section. For
        // an integer `n = reduced << trailing` represented as a
        // `precision`-bit mantissa with top bit set, the exponent
        // is `trailing + significant - 1`.
        let exponent = i64::from(trailing) + i64::from(significant) - 1;

        Ok(Self {
            class: Class::Normal {
                sign,
                exponent,
                mantissa: storage,
            },
            precision,
        })
    }

    /// Returns the precision (in bits) of this value.
    #[inline]
    #[must_use]
    pub fn precision(&self) -> u32 {
        self.precision
    }
}

#[cfg(feature = "big")]
fn validate_precision(precision: u32) -> Result<(), BuildError> {
    if precision == 0 {
        Err(BuildError::PrecisionZero)
    } else {
        Ok(())
    }
}

#[cfg(feature = "big")]
fn pad_payload(payload: &[u64], target_len: usize) -> Vec<u64> {
    let mut storage: Vec<u64> = vec![0; target_len];
    let copy_len = payload.len().min(target_len);
    storage[..copy_len].copy_from_slice(&payload[..copy_len]);
    storage
}

#[cfg(test)]
#[cfg(feature = "big")]
mod tests {
    use super::*;
    use alloc::vec;

    #[test]
    fn zero_precision_is_an_error() {
        assert_eq!(
            BigFloat::try_new_zero(Sign::Positive, 0),
            Err(BuildError::PrecisionZero)
        );
        assert_eq!(
            BigFloat::try_new_infinity(Sign::Positive, 0),
            Err(BuildError::PrecisionZero)
        );
        assert_eq!(
            BigFloat::try_new_quiet_nan(Sign::Positive, 0, &[]),
            Err(BuildError::PrecisionZero)
        );
        assert_eq!(
            BigFloat::try_from_i64_exact(0, 0),
            Err(BuildError::PrecisionZero)
        );
    }

    #[test]
    fn zero_construction_round_trip() {
        let pos = BigFloat::try_new_zero(Sign::Positive, 53).unwrap();
        let neg = BigFloat::try_new_zero(Sign::Negative, 53).unwrap();
        assert_eq!(pos.precision(), 53);
        assert_eq!(neg.precision(), 53);
        assert_ne!(pos, neg); // ±0 are distinct values per IEEE.
    }

    #[test]
    fn infinity_construction_round_trip() {
        let pos = BigFloat::try_new_infinity(Sign::Positive, 113).unwrap();
        let neg = BigFloat::try_new_infinity(Sign::Negative, 113).unwrap();
        assert_eq!(pos.precision(), 113);
        assert_ne!(pos, neg);
    }

    #[test]
    fn nan_payload_is_padded_to_precision() {
        let n = BigFloat::try_new_quiet_nan(Sign::Positive, 256, &[1, 2]).unwrap();
        match &n.class {
            Class::Nan { payload, .. } => {
                assert_eq!(payload.len(), 4); // limbs_for(256)
                assert_eq!(payload[0], 1);
                assert_eq!(payload[1], 2);
                assert_eq!(payload[2], 0);
                assert_eq!(payload[3], 0);
            }
            _ => panic!("expected NaN"),
        }
    }

    #[test]
    fn nan_payload_is_truncated_when_too_long() {
        let n = BigFloat::try_new_quiet_nan(Sign::Positive, 64, &[1, 2, 3, 4]).unwrap();
        match &n.class {
            Class::Nan { payload, .. } => {
                assert_eq!(payload.len(), 1); // limbs_for(64)
                assert_eq!(payload[0], 1);
            }
            _ => panic!("expected NaN"),
        }
    }

    #[test]
    fn quiet_and_signaling_distinct() {
        let q = BigFloat::try_new_quiet_nan(Sign::Positive, 53, &[]).unwrap();
        let s = BigFloat::try_new_signaling_nan(Sign::Positive, 53, &[]).unwrap();
        assert_ne!(q, s);
    }

    #[test]
    fn from_i64_zero() {
        let z = BigFloat::try_from_i64_exact(0, 53).unwrap();
        assert_eq!(
            z.class,
            Class::Zero {
                sign: Sign::Positive
            }
        );
    }

    #[test]
    fn from_i64_one_at_53_bits() {
        let one = BigFloat::try_from_i64_exact(1, 53).unwrap();
        match one.class {
            Class::Normal {
                sign,
                exponent,
                mantissa,
            } => {
                assert_eq!(sign, Sign::Positive);
                assert_eq!(exponent, 0);
                assert_eq!(mantissa, vec![1u64 << 63]);
            }
            other => panic!("expected Normal, got {other:?}"),
        }
    }

    #[test]
    fn from_i64_two_at_53_bits() {
        let two = BigFloat::try_from_i64_exact(2, 53).unwrap();
        match two.class {
            Class::Normal {
                exponent, mantissa, ..
            } => {
                assert_eq!(exponent, 1);
                // 2 = 1 (top bit) at the same storage layout
                // (trailing zero stripped before placement). The
                // top bit of the storage is set.
                assert_eq!(mantissa, vec![1u64 << 63]);
            }
            other => panic!("expected Normal, got {other:?}"),
        }
    }

    #[test]
    fn from_i64_five_at_53_bits() {
        // 5 = 0b101, significant_bits = 3.
        // mantissa as 53-bit int = 5 << 50 = top three bits set per
        // 0b101 pattern at the top of a 53-bit field.
        // storage = (5 << 50) << 11 = 5 << 61.
        let five = BigFloat::try_from_i64_exact(5, 53).unwrap();
        match five.class {
            Class::Normal {
                exponent, mantissa, ..
            } => {
                assert_eq!(exponent, 2);
                assert_eq!(mantissa, vec![5u64 << 61]);
            }
            other => panic!("expected Normal, got {other:?}"),
        }
    }

    #[test]
    fn from_i64_negative_carries_sign() {
        let n = BigFloat::try_from_i64_exact(-7, 53).unwrap();
        match n.class {
            Class::Normal { sign, .. } => assert_eq!(sign, Sign::Negative),
            _ => panic!("expected Normal"),
        }
    }

    #[test]
    fn from_i64_min_at_64_bits() {
        // i64::MIN = -2^63. magnitude = 2^63, bit_length = 64,
        // trailing_zeros = 63, significant_bits = 1.
        let n = BigFloat::try_from_i64_exact(i64::MIN, 53).unwrap();
        match n.class {
            Class::Normal {
                sign,
                exponent,
                mantissa,
            } => {
                assert_eq!(sign, Sign::Negative);
                assert_eq!(exponent, 63);
                assert_eq!(mantissa, vec![1u64 << 63]);
            }
            _ => panic!("expected Normal"),
        }
    }

    #[test]
    fn from_i64_overflow_at_low_precision() {
        // 5 has 3 significant bits; precision 2 cannot represent it
        // exactly.
        let result = BigFloat::try_from_i64_exact(5, 2);
        assert_eq!(
            result,
            Err(BuildError::ValueExceedsPrecision {
                value_bits: 3,
                requested: 2,
            })
        );
    }

    #[test]
    fn from_i64_at_high_precision_uses_multiple_limbs() {
        // precision 128 = 2 limbs. n=1: top bit at MSL bit 63.
        let one = BigFloat::try_from_i64_exact(1, 128).unwrap();
        match one.class {
            Class::Normal {
                exponent, mantissa, ..
            } => {
                assert_eq!(exponent, 0);
                assert_eq!(mantissa.len(), 2);
                assert_eq!(mantissa[0], 0);
                assert_eq!(mantissa[1], 1u64 << 63);
            }
            _ => panic!("expected Normal"),
        }
    }
}
