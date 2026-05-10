//! Tagged union for IEEE 754-2019 value kinds (per ADR-0005).
//!
//! Special values (`±0`, `±∞`, signaling and quiet NaNs) get
//! dedicated variants; finite normals carry their sign, exponent,
//! and mantissa storage. The sign rides inside each variant rather
//! than at the outer level so that IEEE-required distinctions
//! (signed zero, signed NaN) are representable by construction.
//!
//! Slice 1a defines [`Class`] for [`BigFloat`](crate::big::BigFloat).
//! `ClassFixed<const PREC: u32>` for `FixedFloat<PREC>` lands in
//! slice 1g.

#[cfg(feature = "alloc")]
extern crate alloc;

use crate::sign::Sign;

/// Internal representation of a [`BigFloat`](crate::big::BigFloat).
///
/// IEEE 754-2019 §3.2 lists five value kinds: signed zero, signed
/// infinity, quiet NaN, signaling NaN, and finite (normal /
/// subnormal). pfloat does not have subnormals at arbitrary
/// precision (no implicit minimum exponent), so [`Class::Normal`]
/// covers every finite non-zero value.
///
/// The mantissa and NaN payload both use [`Vec<u64>`] storage
/// (matching the bit-level precision via
/// [`limbs_for`](crate::mantissa::limbs_for)). The variable-width
/// payload follows ADR-0005's Update note: the diagnostic bits
/// scale with precision rather than capping at a fixed width.
#[cfg(feature = "big")]
#[derive(Clone, Debug, Eq, PartialEq, Hash)]
pub enum Class {
    /// `±0` per IEEE 754-2019 §6.3. The sign bit propagates through
    /// `copysign`, `neg`, `abs`, and the sign rule on multiplication
    /// and division.
    Zero { sign: Sign },

    /// `±∞` per IEEE 754-2019 §3.2. Produced by overflow under
    /// directed rounding and as the result of `1/0` and similar
    /// limit operations.
    Infinity { sign: Sign },

    /// NaN per IEEE 754-2019 §3.2 and §6.2.
    ///
    /// `quiet` is `true` for a quiet NaN (qNaN), `false` for a
    /// signaling NaN (sNaN). The IEEE 754-2019 §6.2.1 distinction:
    /// most operations on a sNaN raise [`INVALID`](crate::status::Status::INVALID)
    /// and return a qNaN derived from the sNaN. The
    /// [`copysign`](crate::big::BigFloat::copysign) and
    /// [`abs`](crate::big::BigFloat::abs) operations preserve
    /// signaling-ness.
    ///
    /// `payload` carries the §6.2.2 diagnostic bits, sized to
    /// `limbs_for(precision)` u64 limbs to match the mantissa
    /// shape. Most callers leave the payload all-zero.
    Nan {
        quiet: bool,
        sign: Sign,
        payload: alloc::vec::Vec<u64>,
    },

    /// Finite non-zero value.
    ///
    /// `mantissa` is interpreted as a `precision`-bit unsigned
    /// integer with the top bit set, laid out in `Vec<u64>`
    /// little-endian limbs (limb 0 = least significant 64 bits).
    /// The most-significant bit of the most-significant limb is
    /// `1`; bits below the precision are zero in canonical form.
    ///
    /// The value is `sign × mantissa × 2^(exponent - precision + 1)`
    /// where `mantissa` is the integer interpretation. ADR-0001 and
    /// ADR-0002 set the layout; ADR-0006 fixes `i64` for the
    /// exponent.
    Normal {
        sign: Sign,
        exponent: i64,
        mantissa: alloc::vec::Vec<u64>,
    },
}

#[cfg(test)]
#[cfg(feature = "big")]
mod tests {
    use super::*;
    use alloc::vec;

    #[test]
    fn zero_variants_distinct_by_sign() {
        let pos = Class::Zero {
            sign: Sign::Positive,
        };
        let neg = Class::Zero {
            sign: Sign::Negative,
        };
        assert_ne!(pos, neg);
    }

    #[test]
    fn infinity_variants_distinct_by_sign() {
        let pos = Class::Infinity {
            sign: Sign::Positive,
        };
        let neg = Class::Infinity {
            sign: Sign::Negative,
        };
        assert_ne!(pos, neg);
    }

    #[test]
    fn nan_quiet_signaling_distinct() {
        let q = Class::Nan {
            quiet: true,
            sign: Sign::Positive,
            payload: vec![0],
        };
        let s = Class::Nan {
            quiet: false,
            sign: Sign::Positive,
            payload: vec![0],
        };
        assert_ne!(q, s);
    }

    #[test]
    fn nan_payload_distinct() {
        let a = Class::Nan {
            quiet: true,
            sign: Sign::Positive,
            payload: vec![0],
        };
        let b = Class::Nan {
            quiet: true,
            sign: Sign::Positive,
            payload: vec![1],
        };
        assert_ne!(a, b);
    }

    #[test]
    fn normal_distinct_by_field() {
        let one = Class::Normal {
            sign: Sign::Positive,
            exponent: 0,
            mantissa: vec![1u64 << 63],
        };
        let other = Class::Normal {
            sign: Sign::Positive,
            exponent: 1,
            mantissa: vec![1u64 << 63],
        };
        assert_ne!(one, other);
    }
}
