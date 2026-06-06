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
use crate::rounding::{round_finite_to_precision, RoundingMode};
use crate::sign::Sign;
use crate::status::Status;

/// Pure-Rust correctly-rounded arbitrary-precision binary float
/// with runtime precision.
///
/// `BigFloat` matches MPFR's shape: the precision is chosen at
/// construction time and may differ between values. Use
/// [`FixedFloat<PREC>`] (slice 1g) when the precision is known at
/// compile time and stack allocation is preferable.
///
/// Slices 1a–1b (currently shipped) provide construction,
/// classification, sign manipulation, comparison, exact and
/// rounding-aware integer conversion, and re-rounding to a different
/// precision. Slices 1c–1f add the arithmetic kernels.
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
#[cfg_attr(feature = "serde", derive(serde::Serialize, serde::Deserialize))]
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

    /// Constructs a `BigFloat` from an `i64`, rounding under
    /// `mode` when the value's significant bits exceed
    /// `precision`. Always succeeds at any `precision >= 1`.
    ///
    /// Returns the value and a [`Status`] carrying
    /// [`Status::INEXACT`] when rounding discarded a non-zero bit.
    /// Under the `std` feature, the status is also OR-accumulated
    /// into the thread-local flag set
    /// (see [`flags`](crate::status::flags)).
    pub fn try_from_i64_round(
        n: i64,
        precision: u32,
        mode: RoundingMode,
    ) -> Result<(Self, Status), BuildError> {
        validate_precision(precision)?;

        if n == 0 {
            // Zero is exact at every precision.
            return Ok((
                Self {
                    class: Class::Zero {
                        sign: Sign::Positive,
                    },
                    precision,
                },
                Status::OK,
            ));
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

        if significant <= precision {
            // Exact path: identical to `try_from_i64_exact`.
            let value = Self::try_from_i64_exact(n, precision)
                .expect("exact-fit case already validated above");
            return Ok((value, Status::OK));
        }

        // Rounding required. Build an intermediate at `significant`-bit
        // precision (the natural width of the magnitude after stripping
        // trailing zeros), then route through the rounding pipeline to
        // shrink to the user's `precision`.
        let intermediate_precision = significant;
        let intermediate_limbs = limbs_for(intermediate_precision);
        let mut intermediate: Vec<u64> = vec![0u64; intermediate_limbs];
        let reduced: u64 = magnitude >> trailing;
        // `reduced` has `significant` bits with the top bit set.
        // Place it at the top of the intermediate storage.
        let total_shift = storage_shift(intermediate_limbs, significant);
        let whole = (total_shift / 64) as usize;
        let intra = total_shift % 64;
        if intra == 0 {
            intermediate[whole] = reduced;
        } else {
            intermediate[whole] = reduced << intra;
            if 64 - intra < significant && whole + 1 < intermediate_limbs {
                intermediate[whole + 1] = reduced >> (64 - intra);
            }
        }

        let exponent = i64::from(trailing) + i64::from(significant) - 1;

        let (value, status) = round_finite_to_precision(
            sign,
            exponent,
            &intermediate,
            intermediate_precision,
            false, // no upstream sticky
            precision,
            mode,
        );

        crate::status::auto_raise(status);
        Ok((value, status))
    }

    /// `try_from_i64_round` accumulating into a caller-supplied
    /// flag bag (`no_std`-friendly variant).
    ///
    /// Equivalent to:
    /// ```ignore
    /// let (value, status) = BigFloat::try_from_i64_round(n, precision, mode)?;
    /// *flags |= status;
    /// Ok(value)
    /// ```
    pub fn try_from_i64_round_with_flags(
        n: i64,
        precision: u32,
        mode: RoundingMode,
        flags: &mut Status,
    ) -> Result<Self, BuildError> {
        let (value, status) = Self::try_from_i64_round(n, precision, mode)?;
        *flags |= status;
        Ok(value)
    }

    /// Re-rounds `self` to a new precision under `mode`.
    ///
    /// When `new_precision >= self.precision`, the result is exact
    /// (the mantissa pads with trailing zeros). When
    /// `new_precision < self.precision`, the rounding pipeline
    /// applies and may set [`Status::INEXACT`].
    ///
    /// Special values pass through: NaN preserves payload (zero-
    /// padded or truncated to the new precision), `±0` and `±∞`
    /// re-emit at the new precision unchanged.
    ///
    /// Returns [`BuildError::PrecisionZero`] when
    /// `new_precision == 0`.
    pub fn round_to_precision(
        &self,
        new_precision: u32,
        mode: RoundingMode,
    ) -> Result<(Self, Status), BuildError> {
        validate_precision(new_precision)?;

        let (value, status) = match &self.class {
            Class::Zero { sign } => (
                Self {
                    class: Class::Zero { sign: *sign },
                    precision: new_precision,
                },
                Status::OK,
            ),
            Class::Infinity { sign } => (
                Self {
                    class: Class::Infinity { sign: *sign },
                    precision: new_precision,
                },
                Status::OK,
            ),
            Class::Nan {
                quiet,
                sign,
                payload,
            } => (
                Self {
                    class: Class::Nan {
                        quiet: *quiet,
                        sign: *sign,
                        payload: pad_payload(payload, limbs_for(new_precision)),
                    },
                    precision: new_precision,
                },
                Status::OK,
            ),
            Class::Normal {
                sign,
                exponent,
                mantissa,
            } => round_finite_to_precision(
                *sign,
                *exponent,
                mantissa,
                self.precision,
                false,
                new_precision,
                mode,
            ),
        };

        crate::status::auto_raise(status);
        Ok((value, status))
    }

    /// `round_to_precision` accumulating into a caller-supplied
    /// flag bag.
    pub fn round_to_precision_with_flags(
        &self,
        new_precision: u32,
        mode: RoundingMode,
        flags: &mut Status,
    ) -> Result<Self, BuildError> {
        let (value, status) = self.round_to_precision(new_precision, mode)?;
        *flags |= status;
        Ok(value)
    }

    /// Returns the precision (in bits) of this value.
    #[inline]
    #[must_use]
    pub fn precision(&self) -> u32 {
        self.precision
    }

    /// Borrowed read-only view of the raw IEEE-shaped structure.
    ///
    /// Pattern-match on [`Parts`] to inspect the sign, exponent,
    /// mantissa limbs, and (for NaN) payload limbs without
    /// allocating or rounding. External callers that need bit-exact
    /// access to the representation (differential testing,
    /// serialization, raw-form printers) should reach for this
    /// accessor rather than going through `Display` and a round-trip
    /// parser, which loses up to 1 ULP under non-NearestEven
    /// rounding.
    ///
    /// The accessor is `O(1)` and panic-free. The returned [`Parts`]
    /// borrows from `self` for the slice fields, so it cannot
    /// outlive this `BigFloat`.
    ///
    /// pfloat does not expose a converse constructor from raw parts.
    /// Construction goes through the validated `try_new_*` paths so
    /// the top-bit-set mantissa normalization, the
    /// `limbs_for(precision)` storage shape, and the
    /// precision-bound payload length stay invariants the type
    /// system checks.
    #[inline]
    #[must_use]
    pub fn parts(&self) -> Parts<'_> {
        match &self.class {
            Class::Zero { sign } => Parts::Zero { sign: *sign },
            Class::Infinity { sign } => Parts::Infinity { sign: *sign },
            Class::Nan {
                quiet,
                sign,
                payload,
            } => Parts::Nan {
                quiet: *quiet,
                sign: *sign,
                payload,
            },
            Class::Normal {
                sign,
                exponent,
                mantissa,
            } => Parts::Normal {
                sign: *sign,
                exponent: *exponent,
                mantissa,
                precision: self.precision,
            },
        }
    }
}

/// Read-only view of a [`BigFloat`]'s raw IEEE-shaped representation.
///
/// Public mirror of the internal tagged-union storage. The lifetime
/// `'a` ties the borrowed slice fields (`payload`, `mantissa`) to
/// the `BigFloat` that produced this `Parts`.
///
/// Each variant carries exactly the fields IEEE 754-2019 §3.2 lists
/// for the corresponding value kind. The `Normal` variant's
/// `mantissa` is a little-endian limb slice (`mantissa[0]` is the
/// least significant 64 bits) with the top bit of the most
/// significant limb set; the integer interpretation of the
/// mantissa, scaled by `2^(exponent - precision + 1)`, is the value.
/// ADR-0001 records the layout; ADR-0006 fixes the `i64` exponent;
/// ADR-0016 records the rationale for exposing this view publicly.
#[cfg(feature = "big")]
#[derive(Clone, Copy, Debug, Eq, PartialEq, Hash)]
pub enum Parts<'a> {
    /// Signed zero (IEEE 754-2019 §6.3).
    Zero { sign: Sign },
    /// Signed infinity (IEEE 754-2019 §3.2).
    Infinity { sign: Sign },
    /// NaN with payload (IEEE 754-2019 §6.2). `quiet` is `true` for
    /// quiet NaN, `false` for signaling NaN. The payload slice has
    /// length `limbs_for(precision)`; the `precision` is the
    /// originating `BigFloat`'s precision and is recoverable from
    /// `BigFloat::precision()` if needed (omitted from the variant
    /// to keep `Parts` zero-cost to construct).
    Nan {
        quiet: bool,
        sign: Sign,
        payload: &'a [u64],
    },
    /// Finite non-zero value. The integer interpretation of
    /// `mantissa` (top-bit set, `precision` bits wide) scaled by
    /// `2^(exponent - precision + 1)` and signed by `sign` is the
    /// represented value.
    Normal {
        sign: Sign,
        exponent: i64,
        mantissa: &'a [u64],
        precision: u32,
    },
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

    // ---------- Slice 1b: rounding-aware constructors ----------

    #[test]
    fn try_from_i64_round_exact_path_no_flags() {
        let (one, status) = BigFloat::try_from_i64_round(1, 53, RoundingMode::NearestEven).unwrap();
        assert!(status.is_ok());
        assert_eq!(one, BigFloat::try_from_i64_exact(1, 53).unwrap());
    }

    #[test]
    fn try_from_i64_round_inexact_sets_inexact_flag() {
        // 5 has 3 significant bits. Precision 2 cannot hold it.
        // Top 2 bits of 0b101 are 0b10 = 2, with guard=0 sticky=1.
        // NearestEven: guard=0 → no round up. Result mantissa = 0b10
        // (= 2). Exponent: top bit of 5 is at position 2, so
        // result has exponent = 2.
        let (v, status) = BigFloat::try_from_i64_round(5, 2, RoundingMode::NearestEven).unwrap();
        assert!(status.inexact());
        assert!(!status.invalid());
        assert!(!status.overflow());
        // The rounded value should be 4 (= 0b100, which is 2 << 2,
        // i.e., mantissa 0b10 with exponent 2).
        let four = BigFloat::try_from_i64_exact(4, 2).unwrap();
        assert_eq!(v, four);
    }

    #[test]
    fn try_from_i64_round_nearest_even_round_up() {
        // 7 = 0b111. Precision 2 keeps top 2 bits = 0b11.
        // Guard = bit 0 = 1, sticky = 0, lowest_kept = 1.
        // NearestEven: guard && (sticky || lowest_kept) = 1 && (0 || 1) = round up.
        // 0b11 + 1 = 0b100 → renormalize: mantissa 0b10, exponent +1.
        // Original exponent for 7 (bit_length 3, trailing 0): 0 + 3 - 1 = 2.
        // After round-up: mantissa 0b10 with exponent 3, value = 2 * 2^(3-2+1) = 2*4 = 8.
        let (v, status) = BigFloat::try_from_i64_round(7, 2, RoundingMode::NearestEven).unwrap();
        assert!(status.inexact());
        let eight = BigFloat::try_from_i64_exact(8, 2).unwrap();
        assert_eq!(v, eight);
    }

    #[test]
    fn try_from_i64_round_toward_zero_truncates() {
        // 7 = 0b111 at precision 2 with TowardZero: truncate to 0b11 = 3,
        // exponent stays at 2. Value = 3 * 2^(2-2+1) = 3*2 = 6.
        let (v, status) = BigFloat::try_from_i64_round(7, 2, RoundingMode::TowardZero).unwrap();
        assert!(status.inexact());
        let six = BigFloat::try_from_i64_exact(6, 2).unwrap();
        assert_eq!(v, six);
    }

    #[test]
    fn try_from_i64_round_toward_positive_positive_rounds_up() {
        let (v, status) = BigFloat::try_from_i64_round(5, 2, RoundingMode::TowardPositive).unwrap();
        assert!(status.inexact());
        // 5 → top 2 bits 0b10 with sticky bit set → round up to 0b11
        // → 6.
        let six = BigFloat::try_from_i64_exact(6, 2).unwrap();
        assert_eq!(v, six);
    }

    #[test]
    fn try_from_i64_round_toward_negative_negative_rounds_up() {
        // -5 magnitude 5. TowardNegative on negative sign: round up
        // toward -∞ means away from zero in magnitude. Same shape
        // as TowardPositive on a positive 5.
        let (v, status) =
            BigFloat::try_from_i64_round(-5, 2, RoundingMode::TowardNegative).unwrap();
        assert!(status.inexact());
        let neg_six = BigFloat::try_from_i64_exact(-6, 2).unwrap();
        assert_eq!(v, neg_six);
    }

    #[test]
    fn try_from_i64_round_directed_modes_truncate_wrong_sign() {
        // TowardPositive on a negative value: truncate toward zero
        // (the "ceiling" for negative is closer to zero).
        let (v, _) = BigFloat::try_from_i64_round(-5, 2, RoundingMode::TowardPositive).unwrap();
        // -5 → top 2 bits 0b10 truncated = magnitude 4 → -4.
        let neg_four = BigFloat::try_from_i64_exact(-4, 2).unwrap();
        assert_eq!(v, neg_four);
    }

    #[test]
    fn try_from_i64_round_zero_is_always_exact() {
        let (z, status) = BigFloat::try_from_i64_round(0, 53, RoundingMode::NearestEven).unwrap();
        assert!(status.is_ok());
        assert!(z.is_zero());
        assert!(z.is_sign_positive());
    }

    #[test]
    fn try_from_i64_round_with_flags_accumulates() {
        let mut flags = Status::OK;
        let _ =
            BigFloat::try_from_i64_round_with_flags(7, 2, RoundingMode::NearestEven, &mut flags)
                .unwrap();
        assert!(flags.inexact());
        // Calling again accumulates without clearing.
        let _ =
            BigFloat::try_from_i64_round_with_flags(0, 53, RoundingMode::NearestEven, &mut flags)
                .unwrap();
        assert!(flags.inexact()); // still set from before
    }

    #[test]
    fn round_to_precision_extension_is_exact() {
        let one = BigFloat::try_from_i64_exact(1, 53).unwrap();
        let (extended, status) = one
            .round_to_precision(113, RoundingMode::NearestEven)
            .unwrap();
        assert!(status.is_ok());
        assert_eq!(extended.precision(), 113);
        // Numerically equal to the original.
        assert_eq!(
            extended.partial_cmp(&one).0,
            Some(core::cmp::Ordering::Equal)
        );
    }

    #[test]
    fn round_to_precision_narrows_with_inexact() {
        // Build a 4-bit mantissa value 0b1101 (= 13) at precision 4
        // with exponent 3 (so the value is 13).
        let thirteen = BigFloat::try_from_i64_exact(13, 4).unwrap();
        // Round to precision 2: top 2 bits 0b11, guard 0, sticky 1
        // → no round up under NearestEven; result mantissa 0b11 at
        // exponent 3 → value 12.
        let (rounded, status) = thirteen
            .round_to_precision(2, RoundingMode::NearestEven)
            .unwrap();
        assert!(status.inexact());
        let twelve = BigFloat::try_from_i64_exact(12, 2).unwrap();
        assert_eq!(rounded, twelve);
    }

    #[test]
    fn round_to_precision_preserves_zero() {
        let pz = BigFloat::try_new_zero(Sign::Positive, 113).unwrap();
        let (r, _) = pz
            .round_to_precision(53, RoundingMode::NearestEven)
            .unwrap();
        assert!(r.is_zero());
        assert!(r.is_sign_positive());
        assert_eq!(r.precision(), 53);
    }

    #[test]
    fn round_to_precision_preserves_negative_zero() {
        let nz = BigFloat::try_new_zero(Sign::Negative, 113).unwrap();
        let (r, _) = nz
            .round_to_precision(53, RoundingMode::NearestEven)
            .unwrap();
        assert!(r.is_zero());
        assert!(r.is_sign_negative());
    }

    #[test]
    fn round_to_precision_preserves_infinity() {
        let pi = BigFloat::try_new_infinity(Sign::Positive, 53).unwrap();
        let (r, _) = pi
            .round_to_precision(113, RoundingMode::NearestEven)
            .unwrap();
        assert!(r.is_infinite());
        assert!(r.is_sign_positive());
    }

    #[test]
    fn round_to_precision_preserves_nan() {
        let n = BigFloat::try_new_quiet_nan(Sign::Negative, 113, &[42]).unwrap();
        let (r, _) = n.round_to_precision(53, RoundingMode::NearestEven).unwrap();
        assert!(r.is_quiet_nan());
        assert!(r.is_sign_negative());
        assert_eq!(r.precision(), 53);
        match r.class {
            Class::Nan { payload, .. } => {
                // Truncated to 1 limb (limbs_for(53) = 1).
                assert_eq!(payload.len(), 1);
                assert_eq!(payload[0], 42);
            }
            _ => panic!("expected NaN"),
        }
    }

    #[test]
    fn round_to_precision_rejects_zero_precision() {
        let one = BigFloat::try_from_i64_exact(1, 53).unwrap();
        assert_eq!(
            one.round_to_precision(0, RoundingMode::NearestEven),
            Err(BuildError::PrecisionZero)
        );
    }

    #[cfg(feature = "std")]
    #[test]
    fn rounding_ops_update_thread_local_flags() {
        crate::status::flags::clear();
        let _ = BigFloat::try_from_i64_round(7, 2, RoundingMode::NearestEven).unwrap();
        assert!(crate::status::flags::test().inexact());
        crate::status::flags::clear();
    }
}
