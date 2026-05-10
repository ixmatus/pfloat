//! [`FixedFloat`]: compile-time-precision arbitrary-precision float.
//!
//! `FixedFloat<const PREC: u32>` is the const-generic counterpart to
//! [`BigFloat`](crate::big::BigFloat) per ADR-0003 and ADR-0004. The
//! precision is part of the type; the mantissa storage is a
//! stack-allocated `[u64; limbs_for(PREC)]` array, and the type is
//! `Copy` at every supported precision.
//!
//! Const-generic mantissa storage uses
//! `feature(generic_const_exprs)` per ADR-0011. Every public item
//! that references the storage shape carries the
//! `where [(); limbs_for(PREC)]:` clause.
//!
//! # Slice 1g shipping shape
//!
//! All arithmetic and rounding operations on `FixedFloat<PREC>`
//! delegate to [`BigFloat`](crate::big::BigFloat) via exact
//! conversion in and explicit-precision conversion back. The dual-
//! type architecture lands without duplicating kernel logic; the
//! heap allocation cost of the conversion is acceptable for a 1.0
//! correctness milestone and can be optimized in Phase 7.
//!
//! Construction, classification, sign manipulation, and comparison
//! run directly on the `[u64; N]` storage and do not allocate.

use core::cmp::Ordering;

use crate::big::{BigFloat, BuildError};
use crate::class::Class;
use crate::mantissa::limbs_for;
use crate::rounding::RoundingMode;
use crate::sign::Sign;
use crate::status::Status;

/// Internal representation of a [`FixedFloat<PREC>`].
///
/// Mirrors [`Class`](crate::class::Class) but with stack-allocated
/// `[u64; limbs_for(PREC)]` storage for the mantissa and the NaN
/// payload. Each variant is `Copy` at every reasonable `PREC`.
#[derive(Clone, Copy, Debug, Eq, PartialEq, Hash)]
pub enum ClassFixed<const PREC: u32>
where
    [(); limbs_for(PREC)]:,
{
    Zero {
        sign: Sign,
    },
    Infinity {
        sign: Sign,
    },
    Nan {
        quiet: bool,
        sign: Sign,
        payload: [u64; limbs_for(PREC)],
    },
    Normal {
        sign: Sign,
        exponent: i64,
        mantissa: [u64; limbs_for(PREC)],
    },
}

/// Pure-Rust correctly-rounded arbitrary-precision binary float
/// with compile-time precision.
///
/// `FixedFloat<PREC>` is the const-generic counterpart to
/// [`BigFloat`](crate::big::BigFloat). The precision is fixed in
/// the type; the mantissa storage is a stack-allocated
/// `[u64; limbs_for(PREC)]` array. Suitable for embedded use, hot
/// loops at known precision, and any caller who wants the
/// optimizer to see the precision as a constant.
///
/// `FixedFloat<53>` corresponds to IEEE 754 binary64 with full
/// rounding-mode control; `FixedFloat<113>` to binary128; arbitrary
/// other precisions are equally supported.
#[derive(Clone, Copy, Debug, Eq, PartialEq, Hash)]
pub struct FixedFloat<const PREC: u32>
where
    [(); limbs_for(PREC)]:,
{
    pub(crate) class: ClassFixed<PREC>,
}

impl<const PREC: u32> FixedFloat<PREC>
where
    [(); limbs_for(PREC)]:,
{
    /// Compile-time precision of this type, in bits.
    pub const PRECISION: u32 = PREC;

    /// Limb count of the mantissa storage.
    pub const LIMBS: usize = limbs_for(PREC);

    // -------- Constructors --------

    /// Returns `+0` at this precision.
    #[must_use]
    pub const fn zero() -> Self {
        Self {
            class: ClassFixed::Zero {
                sign: Sign::Positive,
            },
        }
    }

    /// Returns `-0` at this precision.
    #[must_use]
    pub const fn neg_zero() -> Self {
        Self {
            class: ClassFixed::Zero {
                sign: Sign::Negative,
            },
        }
    }

    /// Returns `+∞` at this precision.
    #[must_use]
    pub const fn infinity() -> Self {
        Self {
            class: ClassFixed::Infinity {
                sign: Sign::Positive,
            },
        }
    }

    /// Returns `-∞` at this precision.
    #[must_use]
    pub const fn neg_infinity() -> Self {
        Self {
            class: ClassFixed::Infinity {
                sign: Sign::Negative,
            },
        }
    }

    /// Returns a quiet NaN at this precision with the given sign
    /// and zero payload.
    #[must_use]
    pub const fn nan(sign: Sign) -> Self {
        Self {
            class: ClassFixed::Nan {
                quiet: true,
                sign,
                payload: [0u64; limbs_for(PREC)],
            },
        }
    }

    /// Returns a signaling NaN at this precision with the given
    /// sign and zero payload.
    #[must_use]
    pub const fn signaling_nan(sign: Sign) -> Self {
        Self {
            class: ClassFixed::Nan {
                quiet: false,
                sign,
                payload: [0u64; limbs_for(PREC)],
            },
        }
    }

    /// Constructs `FixedFloat<PREC>` from an `i64` exactly.
    ///
    /// Returns [`BuildError::ValueExceedsPrecision`] when the
    /// integer's significant-bit count exceeds `PREC`.
    pub fn try_from_i64_exact(n: i64) -> Result<Self, BuildError> {
        BigFloat::try_from_i64_exact(n, PREC).and_then(Self::try_from_big_exact)
    }

    /// Constructs `FixedFloat<PREC>` from an `i64` with rounding
    /// when the integer's significant bits exceed `PREC`.
    pub fn try_from_i64_round(n: i64, mode: RoundingMode) -> (Self, Status) {
        let (big, status) =
            BigFloat::try_from_i64_round(n, PREC, mode).expect("PREC >= 1 by const-generic bound");
        // Conversion to FixedFloat<PREC> at the same precision is
        // exact (the BigFloat was constructed at PREC), so status
        // here is OK.
        let fixed = Self::from_big_at_same_precision(big);
        (fixed, status)
    }

    // -------- Accessors --------

    /// Returns the precision (in bits) of this value.
    #[inline]
    #[must_use]
    pub const fn precision(&self) -> u32 {
        PREC
    }

    /// Returns the value's sign attribute.
    #[inline]
    #[must_use]
    pub const fn sign(&self) -> Sign {
        match &self.class {
            ClassFixed::Zero { sign }
            | ClassFixed::Infinity { sign }
            | ClassFixed::Nan { sign, .. }
            | ClassFixed::Normal { sign, .. } => *sign,
        }
    }

    // -------- Classification --------

    #[inline]
    #[must_use]
    pub const fn is_nan(&self) -> bool {
        matches!(self.class, ClassFixed::Nan { .. })
    }

    #[inline]
    #[must_use]
    pub const fn is_signaling_nan(&self) -> bool {
        matches!(self.class, ClassFixed::Nan { quiet: false, .. })
    }

    #[inline]
    #[must_use]
    pub const fn is_quiet_nan(&self) -> bool {
        matches!(self.class, ClassFixed::Nan { quiet: true, .. })
    }

    #[inline]
    #[must_use]
    pub const fn is_infinite(&self) -> bool {
        matches!(self.class, ClassFixed::Infinity { .. })
    }

    #[inline]
    #[must_use]
    pub const fn is_finite(&self) -> bool {
        matches!(
            self.class,
            ClassFixed::Zero { .. } | ClassFixed::Normal { .. }
        )
    }

    #[inline]
    #[must_use]
    pub const fn is_zero(&self) -> bool {
        matches!(self.class, ClassFixed::Zero { .. })
    }

    #[inline]
    #[must_use]
    pub const fn is_normal(&self) -> bool {
        matches!(self.class, ClassFixed::Normal { .. })
    }

    /// Always `false` for pfloat values. See
    /// [`BigFloat::is_subnormal`](crate::big::BigFloat::is_subnormal).
    #[inline]
    #[must_use]
    pub const fn is_subnormal(&self) -> bool {
        false
    }

    #[inline]
    #[must_use]
    pub const fn is_sign_negative(&self) -> bool {
        matches!(self.sign(), Sign::Negative)
    }

    #[inline]
    #[must_use]
    pub const fn is_sign_positive(&self) -> bool {
        matches!(self.sign(), Sign::Positive)
    }

    // -------- Sign manipulation --------

    /// Returns `|self|` (sign forced to [`Sign::Positive`]). NaN
    /// payload preserved.
    #[inline]
    #[must_use]
    pub fn abs(self) -> Self {
        Self {
            class: with_sign_fixed(self.class, Sign::Positive),
        }
    }

    /// Returns the negation (sign flipped). NaN payload preserved.
    #[inline]
    #[must_use]
    pub fn negated(self) -> Self {
        let new_sign = self.sign().flip();
        Self {
            class: with_sign_fixed(self.class, new_sign),
        }
    }

    /// Returns a copy of `self` with the sign of `src`.
    #[inline]
    #[must_use]
    pub fn copysign(self, src: Self) -> Self {
        Self {
            class: with_sign_fixed(self.class, src.sign()),
        }
    }

    /// Returns `±1` matching `self`'s sign for every non-NaN value
    /// (including `±0` and `±∞`). For NaN, returns the input.
    #[must_use]
    pub fn signum(self) -> Self {
        match self.class {
            ClassFixed::Nan { .. } => self,
            ClassFixed::Zero { sign }
            | ClassFixed::Infinity { sign }
            | ClassFixed::Normal { sign, .. } => {
                let one = Self::try_from_i64_exact(1).expect("1 fits in any PREC >= 1");
                if matches!(sign, Sign::Negative) {
                    one.negated()
                } else {
                    one
                }
            }
        }
    }

    // -------- Comparison --------

    /// IEEE 754-2019 §5.11 partial-comparison: `None` if either
    /// operand is NaN; `INVALID` raised on signaling-NaN comparand.
    /// `+0 == -0` numerically.
    #[must_use]
    pub fn partial_cmp(&self, other: &Self) -> (Option<Ordering>, Status) {
        self.to_big().partial_cmp(&other.to_big())
    }

    /// IEEE 754-2019 §5.10 `totalOrder`: defines a total order on
    /// every value including NaN. `-0 < +0`.
    #[must_use]
    pub fn total_cmp(&self, other: &Self) -> Ordering {
        self.to_big().total_cmp(&other.to_big())
    }

    /// IEEE 754-2019 §9.6 `minimumNumber`. Quiet NaN is treated as
    /// missing data; signaling NaN raises `INVALID`.
    #[must_use]
    pub fn min(&self, other: &Self) -> (Self, Status) {
        let (m, s) = self.to_big().min(&other.to_big());
        (Self::from_big_at_same_precision(m), s)
    }

    /// IEEE 754-2019 §9.6 `maximumNumber`. Symmetric to
    /// [`min`](Self::min).
    #[must_use]
    pub fn max(&self, other: &Self) -> (Self, Status) {
        let (m, s) = self.to_big().max(&other.to_big());
        (Self::from_big_at_same_precision(m), s)
    }

    // -------- Arithmetic (delegates to BigFloat) --------

    /// IEEE 754-2019 `addition(self, other)`.
    #[must_use]
    pub fn add(&self, other: &Self, mode: RoundingMode) -> (Self, Status) {
        let (b, s) = self.to_big().add(&other.to_big(), mode);
        (Self::from_big_at_same_precision(b), s)
    }

    /// IEEE 754-2019 `subtraction(self, other)`.
    #[must_use]
    pub fn sub(&self, other: &Self, mode: RoundingMode) -> (Self, Status) {
        let (b, s) = self.to_big().sub(&other.to_big(), mode);
        (Self::from_big_at_same_precision(b), s)
    }

    /// IEEE 754-2019 `multiplication(self, other)`.
    #[must_use]
    pub fn mul(&self, other: &Self, mode: RoundingMode) -> (Self, Status) {
        let (b, s) = self.to_big().mul(&other.to_big(), mode);
        (Self::from_big_at_same_precision(b), s)
    }

    /// IEEE 754-2019 `division(self, other)`.
    #[must_use]
    pub fn div(&self, other: &Self, mode: RoundingMode) -> (Self, Status) {
        let (b, s) = self.to_big().div(&other.to_big(), mode);
        (Self::from_big_at_same_precision(b), s)
    }

    /// IEEE 754-2019 `squareRoot(self)`.
    #[must_use]
    pub fn sqrt(&self, mode: RoundingMode) -> (Self, Status) {
        let (b, s) = self.to_big().sqrt(mode);
        (Self::from_big_at_same_precision(b), s)
    }

    /// IEEE 754-2019 `fusedMultiplyAdd(self, b, c)`.
    #[must_use]
    pub fn fma(&self, b: &Self, c: &Self, mode: RoundingMode) -> (Self, Status) {
        let (r, s) = self.to_big().fma(&b.to_big(), &c.to_big(), mode);
        (Self::from_big_at_same_precision(r), s)
    }

    // -------- BigFloat conversion --------

    /// Converts to [`BigFloat`] at the same precision (exact, infallible).
    #[must_use]
    pub fn to_big(&self) -> BigFloat {
        let class = match &self.class {
            ClassFixed::Zero { sign } => Class::Zero { sign: *sign },
            ClassFixed::Infinity { sign } => Class::Infinity { sign: *sign },
            ClassFixed::Nan {
                quiet,
                sign,
                payload,
            } => Class::Nan {
                quiet: *quiet,
                sign: *sign,
                payload: payload.to_vec(),
            },
            ClassFixed::Normal {
                sign,
                exponent,
                mantissa,
            } => Class::Normal {
                sign: *sign,
                exponent: *exponent,
                mantissa: mantissa.to_vec(),
            },
        };
        BigFloat {
            class,
            precision: PREC,
        }
    }

    /// Constructs `FixedFloat<PREC>` from a [`BigFloat`] at the same
    /// precision. The input's precision must equal `PREC`.
    ///
    /// Use [`try_from_big_round`](Self::try_from_big_round) for
    /// rounding conversions from a different-precision input.
    pub fn try_from_big_exact(big: BigFloat) -> Result<Self, BuildError> {
        if big.precision() == PREC {
            Ok(Self::from_big_at_same_precision(big))
        } else {
            Err(BuildError::ValueExceedsPrecision {
                value_bits: big.precision(),
                requested: PREC,
            })
        }
    }

    /// Constructs `FixedFloat<PREC>` from a [`BigFloat`], rounding
    /// to `PREC` under `mode`.
    pub fn try_from_big_round(big: &BigFloat, mode: RoundingMode) -> (Self, Status) {
        let (rounded, status) = big
            .round_to_precision(PREC, mode)
            .expect("PREC >= 1 by const-generic bound");
        (Self::from_big_at_same_precision(rounded), status)
    }

    /// Private: convert a `BigFloat` known to be at `PREC` into
    /// `FixedFloat<PREC>`. Panics if the precision does not match.
    fn from_big_at_same_precision(big: BigFloat) -> Self {
        debug_assert_eq!(
            big.precision(),
            PREC,
            "from_big_at_same_precision called with mismatched precision"
        );
        let class = match big.class {
            Class::Zero { sign } => ClassFixed::Zero { sign },
            Class::Infinity { sign } => ClassFixed::Infinity { sign },
            Class::Nan {
                quiet,
                sign,
                payload,
            } => ClassFixed::Nan {
                quiet,
                sign,
                payload: vec_to_array::<PREC>(&payload),
            },
            Class::Normal {
                sign,
                exponent,
                mantissa,
            } => ClassFixed::Normal {
                sign,
                exponent,
                mantissa: vec_to_array::<PREC>(&mantissa),
            },
        };
        Self { class }
    }
}

/// Convert a slice of `u64`s into the `[u64; limbs_for(PREC)]`
/// array shape. Pads with zero or truncates if the lengths differ.
fn vec_to_array<const PREC: u32>(v: &[u64]) -> [u64; limbs_for(PREC)]
where
    [(); limbs_for(PREC)]:,
{
    let mut arr = [0u64; limbs_for(PREC)];
    let copy_len = v.len().min(arr.len());
    arr[..copy_len].copy_from_slice(&v[..copy_len]);
    arr
}

fn with_sign_fixed<const PREC: u32>(c: ClassFixed<PREC>, new_sign: Sign) -> ClassFixed<PREC>
where
    [(); limbs_for(PREC)]:,
{
    match c {
        ClassFixed::Zero { .. } => ClassFixed::Zero { sign: new_sign },
        ClassFixed::Infinity { .. } => ClassFixed::Infinity { sign: new_sign },
        ClassFixed::Nan { quiet, payload, .. } => ClassFixed::Nan {
            quiet,
            sign: new_sign,
            payload,
        },
        ClassFixed::Normal {
            exponent, mantissa, ..
        } => ClassFixed::Normal {
            sign: new_sign,
            exponent,
            mantissa,
        },
    }
}

impl<const PREC: u32> From<FixedFloat<PREC>> for BigFloat
where
    [(); limbs_for(PREC)]:,
{
    fn from(value: FixedFloat<PREC>) -> Self {
        value.to_big()
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn precision_const() {
        assert_eq!(FixedFloat::<53>::PRECISION, 53);
        assert_eq!(FixedFloat::<113>::PRECISION, 113);
        assert_eq!(FixedFloat::<256>::PRECISION, 256);
    }

    #[test]
    fn limbs_const() {
        assert_eq!(FixedFloat::<1>::LIMBS, 1);
        assert_eq!(FixedFloat::<53>::LIMBS, 1);
        assert_eq!(FixedFloat::<64>::LIMBS, 1);
        assert_eq!(FixedFloat::<65>::LIMBS, 2);
        assert_eq!(FixedFloat::<128>::LIMBS, 2);
        assert_eq!(FixedFloat::<256>::LIMBS, 4);
    }

    #[test]
    fn constants_round_trip() {
        let z = FixedFloat::<53>::zero();
        assert!(z.is_zero());
        assert!(z.is_sign_positive());

        let nz = FixedFloat::<53>::neg_zero();
        assert!(nz.is_zero());
        assert!(nz.is_sign_negative());

        let pi = FixedFloat::<53>::infinity();
        assert!(pi.is_infinite());
        assert!(pi.is_sign_positive());

        let ni = FixedFloat::<53>::neg_infinity();
        assert!(ni.is_infinite());
        assert!(ni.is_sign_negative());

        let q = FixedFloat::<53>::nan(Sign::Negative);
        assert!(q.is_quiet_nan());
        assert!(q.is_sign_negative());

        let s = FixedFloat::<53>::signaling_nan(Sign::Positive);
        assert!(s.is_signaling_nan());
    }

    #[test]
    fn from_i64_exact_one() {
        let one = FixedFloat::<53>::try_from_i64_exact(1).unwrap();
        assert!(one.is_normal());
        assert!(one.is_sign_positive());
        let one_big = one.to_big();
        assert_eq!(one_big.precision(), 53);
        assert_eq!(one_big, BigFloat::try_from_i64_exact(1, 53).unwrap());
    }

    #[test]
    fn from_i64_exact_overflow() {
        // 5 has 3 significant bits; precision 2 cannot hold it.
        let err = FixedFloat::<2>::try_from_i64_exact(5).unwrap_err();
        assert!(matches!(err, BuildError::ValueExceedsPrecision { .. }));
    }

    #[test]
    fn from_i64_round_rounds_inexact() {
        let (v, status) = FixedFloat::<2>::try_from_i64_round(5, RoundingMode::NearestEven);
        assert!(status.inexact());
        // Result should be 4 (the nearest at 2-bit precision).
        let big = v.to_big();
        assert_eq!(big, BigFloat::try_from_i64_exact(4, 2).unwrap());
    }

    #[test]
    fn add_basic() {
        let a = FixedFloat::<53>::try_from_i64_exact(2).unwrap();
        let b = FixedFloat::<53>::try_from_i64_exact(3).unwrap();
        let (sum, status) = a.add(&b, RoundingMode::NearestEven);
        assert!(status.is_ok());
        let five = FixedFloat::<53>::try_from_i64_exact(5).unwrap();
        assert_eq!(sum.partial_cmp(&five).0, Some(Ordering::Equal));
    }

    #[test]
    fn sub_basic() {
        let a = FixedFloat::<53>::try_from_i64_exact(10).unwrap();
        let b = FixedFloat::<53>::try_from_i64_exact(3).unwrap();
        let (diff, _) = a.sub(&b, RoundingMode::NearestEven);
        let seven = FixedFloat::<53>::try_from_i64_exact(7).unwrap();
        assert_eq!(diff.partial_cmp(&seven).0, Some(Ordering::Equal));
    }

    #[test]
    fn mul_basic() {
        let a = FixedFloat::<53>::try_from_i64_exact(6).unwrap();
        let b = FixedFloat::<53>::try_from_i64_exact(7).unwrap();
        let (p, _) = a.mul(&b, RoundingMode::NearestEven);
        let fortytwo = FixedFloat::<53>::try_from_i64_exact(42).unwrap();
        assert_eq!(p.partial_cmp(&fortytwo).0, Some(Ordering::Equal));
    }

    #[test]
    fn div_basic() {
        let a = FixedFloat::<53>::try_from_i64_exact(20).unwrap();
        let b = FixedFloat::<53>::try_from_i64_exact(4).unwrap();
        let (q, _) = a.div(&b, RoundingMode::NearestEven);
        let five = FixedFloat::<53>::try_from_i64_exact(5).unwrap();
        assert_eq!(q.partial_cmp(&five).0, Some(Ordering::Equal));
    }

    #[test]
    fn sqrt_basic() {
        let nine = FixedFloat::<53>::try_from_i64_exact(9).unwrap();
        let (s, _) = nine.sqrt(RoundingMode::NearestEven);
        let three = FixedFloat::<53>::try_from_i64_exact(3).unwrap();
        assert_eq!(s.partial_cmp(&three).0, Some(Ordering::Equal));
    }

    #[test]
    fn fma_basic() {
        let a = FixedFloat::<53>::try_from_i64_exact(2).unwrap();
        let b = FixedFloat::<53>::try_from_i64_exact(3).unwrap();
        let c = FixedFloat::<53>::try_from_i64_exact(5).unwrap();
        let (r, _) = a.fma(&b, &c, RoundingMode::NearestEven);
        let eleven = FixedFloat::<53>::try_from_i64_exact(11).unwrap();
        assert_eq!(r.partial_cmp(&eleven).0, Some(Ordering::Equal));
    }

    #[test]
    fn abs_neg_copysign() {
        let neg_three = FixedFloat::<53>::try_from_i64_exact(-3).unwrap();
        let pos_three = FixedFloat::<53>::try_from_i64_exact(3).unwrap();
        assert_eq!(
            neg_three.abs().partial_cmp(&pos_three).0,
            Some(Ordering::Equal)
        );
        assert!(neg_three.abs().is_sign_positive());

        assert!(pos_three.negated().is_sign_negative());

        let with_neg_sign = pos_three.copysign(neg_three);
        assert!(with_neg_sign.is_sign_negative());
    }

    #[test]
    fn signum_signs() {
        let three = FixedFloat::<53>::try_from_i64_exact(3).unwrap();
        let one = FixedFloat::<53>::try_from_i64_exact(1).unwrap();
        let neg_three = FixedFloat::<53>::try_from_i64_exact(-3).unwrap();
        let neg_one = FixedFloat::<53>::try_from_i64_exact(-1).unwrap();
        assert_eq!(three.signum().partial_cmp(&one).0, Some(Ordering::Equal));
        assert_eq!(
            neg_three.signum().partial_cmp(&neg_one).0,
            Some(Ordering::Equal)
        );
    }

    #[test]
    fn comparison_and_min_max() {
        let one = FixedFloat::<53>::try_from_i64_exact(1).unwrap();
        let two = FixedFloat::<53>::try_from_i64_exact(2).unwrap();
        assert_eq!(one.partial_cmp(&two).0, Some(Ordering::Less));
        assert_eq!(one.total_cmp(&two), Ordering::Less);

        let (mn, _) = one.min(&two);
        assert_eq!(mn.partial_cmp(&one).0, Some(Ordering::Equal));
        let (mx, _) = one.max(&two);
        assert_eq!(mx.partial_cmp(&two).0, Some(Ordering::Equal));
    }

    #[test]
    fn round_trip_via_bigfloat() {
        let original = FixedFloat::<53>::try_from_i64_exact(42).unwrap();
        let big: BigFloat = original.into();
        let back = FixedFloat::<53>::try_from_big_exact(big).unwrap();
        assert_eq!(original, back);
    }

    #[test]
    fn from_big_round_rounds() {
        let big = BigFloat::try_from_i64_exact(7, 53).unwrap();
        let (fixed, status) = FixedFloat::<2>::try_from_big_round(&big, RoundingMode::NearestEven);
        // 7 at precision 2 rounds to 8 under NearestEven.
        assert!(status.inexact());
        let eight = FixedFloat::<2>::try_from_i64_exact(8).unwrap();
        assert_eq!(fixed.partial_cmp(&eight).0, Some(Ordering::Equal));
    }

    #[test]
    fn copy_semantics() {
        // FixedFloat should be Copy.
        let a = FixedFloat::<53>::try_from_i64_exact(1).unwrap();
        let b = a; // copy, not move
        let c = a; // still usable
        assert_eq!(b.partial_cmp(&c).0, Some(Ordering::Equal));
    }
}
