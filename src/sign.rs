//! IEEE 754-2019 sign attribute (§3.4).
//!
//! pfloat carries the sign as an explicit field in the [`Class`]
//! variants rather than folding it into the mantissa, per ADR-0001.
//!
//! [`Class`]: crate::class::Class

/// Sign of a non-NaN floating-point value.
///
/// `±0`, `±∞`, and finite normals each carry a `Sign`. `NaN` also
/// carries a sign; the IEEE 754-2019 §6.3 sign-bit propagation rules
/// (e.g. through [`copysign`](crate::big::BigFloat::copysign)) require
/// the sign of a NaN to be a meaningful, propagatable bit.
///
/// `Sign::Positive` is the default per IEEE 754 convention (positive
/// zero has sign 0; the value `+1.0` is the canonical positive
/// mantissa pattern).
#[derive(Clone, Copy, Debug, Default, Eq, PartialEq, Hash)]
pub enum Sign {
    #[default]
    Positive,
    Negative,
}

impl Sign {
    /// Returns the opposite sign.
    #[inline]
    #[must_use]
    pub const fn flip(self) -> Self {
        match self {
            Sign::Positive => Sign::Negative,
            Sign::Negative => Sign::Positive,
        }
    }

    /// Sign of a product or quotient: `Positive` iff both operands
    /// share the same sign.
    ///
    /// IEEE 754-2019 §6.3 specifies the multiplicative sign rule
    /// independently of the magnitude computation; this helper
    /// implements that rule.
    #[inline]
    #[must_use]
    pub const fn xor(self, rhs: Self) -> Self {
        match (self, rhs) {
            (Sign::Positive, Sign::Positive) | (Sign::Negative, Sign::Negative) => Sign::Positive,
            _ => Sign::Negative,
        }
    }

    /// `true` for [`Sign::Negative`], `false` for [`Sign::Positive`].
    #[inline]
    #[must_use]
    pub const fn is_negative(self) -> bool {
        matches!(self, Sign::Negative)
    }

    /// `true` for [`Sign::Positive`], `false` for [`Sign::Negative`].
    #[inline]
    #[must_use]
    pub const fn is_positive(self) -> bool {
        matches!(self, Sign::Positive)
    }

    /// `Sign::Negative` if `bit` is `true`, else `Sign::Positive`.
    ///
    /// Helper for round-trips through bit patterns and IEEE 754
    /// encodings where the sign occupies a single bit.
    #[inline]
    #[must_use]
    pub const fn from_bit(bit: bool) -> Self {
        if bit {
            Sign::Negative
        } else {
            Sign::Positive
        }
    }

    /// `1` if [`Sign::Negative`], `0` if [`Sign::Positive`].
    ///
    /// Inverse of [`Sign::from_bit`].
    #[inline]
    #[must_use]
    pub const fn to_bit(self) -> bool {
        matches!(self, Sign::Negative)
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn flip_is_involution() {
        assert_eq!(Sign::Positive.flip(), Sign::Negative);
        assert_eq!(Sign::Negative.flip(), Sign::Positive);
        assert_eq!(Sign::Positive.flip().flip(), Sign::Positive);
        assert_eq!(Sign::Negative.flip().flip(), Sign::Negative);
    }

    #[test]
    fn xor_truth_table() {
        assert_eq!(Sign::Positive.xor(Sign::Positive), Sign::Positive);
        assert_eq!(Sign::Positive.xor(Sign::Negative), Sign::Negative);
        assert_eq!(Sign::Negative.xor(Sign::Positive), Sign::Negative);
        assert_eq!(Sign::Negative.xor(Sign::Negative), Sign::Positive);
    }

    #[test]
    fn xor_is_associative_and_commutative() {
        // Commutativity.
        assert_eq!(
            Sign::Positive.xor(Sign::Negative),
            Sign::Negative.xor(Sign::Positive)
        );
        // Associativity.
        let a = Sign::Negative;
        let b = Sign::Positive;
        let c = Sign::Negative;
        assert_eq!(a.xor(b).xor(c), a.xor(b.xor(c)));
    }

    #[test]
    fn predicates() {
        assert!(Sign::Negative.is_negative());
        assert!(!Sign::Positive.is_negative());
        assert!(Sign::Positive.is_positive());
        assert!(!Sign::Negative.is_positive());
    }

    #[test]
    fn bit_round_trip() {
        assert_eq!(Sign::from_bit(false), Sign::Positive);
        assert_eq!(Sign::from_bit(true), Sign::Negative);
        assert!(!Sign::Positive.to_bit());
        assert!(Sign::Negative.to_bit());
        assert_eq!(Sign::from_bit(Sign::Positive.to_bit()), Sign::Positive);
        assert_eq!(Sign::from_bit(Sign::Negative.to_bit()), Sign::Negative);
    }

    #[test]
    fn default_is_positive() {
        assert_eq!(Sign::default(), Sign::Positive);
    }
}
