//! Comparison: IEEE 754-2019 §5.10 totalOrder predicate, §5.11
//! numeric comparisons, and §9.6 minimumNumber/maximumNumber.
//!
//! pfloat does **not** implement [`PartialOrd`] for [`BigFloat`].
//! NaN handling is opt-in: callers explicitly use
//! [`BigFloat::partial_cmp`], which returns
//! `(Option<Ordering>, Status)` so that signaling-NaN comparands
//! can raise [`Status::INVALID`](crate::status::Status::INVALID)
//! per IEEE 754-2019 §5.11. The convention mirrors ferrodec.

use core::cmp::Ordering;

#[cfg(feature = "big")]
use crate::big::BigFloat;
#[cfg(feature = "big")]
use crate::class::Class;
#[cfg(feature = "big")]
use crate::sign::Sign;
#[cfg(feature = "big")]
use crate::status::Status;

#[cfg(feature = "big")]
impl BigFloat {
    /// IEEE 754-2019 §5.11 `compareQuietEqual` / `compareSignaling*`
    /// fused into a single `Option<Ordering>` return.
    ///
    /// Returns `(Some(ord), Status::OK)` when the comparison is
    /// ordered, `(None, Status::OK)` when one operand is a quiet
    /// NaN, and `(None, Status::INVALID)` when one operand is a
    /// signaling NaN.
    ///
    /// `+0` and `-0` compare equal numerically (per §5.11). Use
    /// [`total_cmp`](Self::total_cmp) for the §5.10 totalOrder
    /// predicate where `-0 < +0`.
    #[must_use]
    pub fn partial_cmp(&self, other: &Self) -> (Option<Ordering>, Status) {
        let mut status = Status::OK;
        if self.is_signaling_nan() || other.is_signaling_nan() {
            status |= Status::INVALID;
            return (None, status);
        }
        if self.is_nan() || other.is_nan() {
            return (None, status);
        }
        // Neither operand is NaN; defer to the magnitude comparison
        // logic (which is the same as `total_cmp` modulo the
        // ±0 == ±0 rule).
        let ord = numeric_cmp_no_nan(self, other);
        (Some(ord), status)
    }

    /// IEEE 754-2019 §5.10 `totalOrder` predicate.
    ///
    /// Defines a total order on every pfloat value, including NaN.
    /// `-0 < +0`. Negative NaNs are less than every finite or
    /// infinite value; positive NaNs are greater. Within NaNs of
    /// the same sign, the encoding determines order (signaling
    /// before quiet on the positive side, the reverse on the
    /// negative side; payload as a tie-breaker).
    ///
    /// Cross-precision comparison is supported: two `BigFloat`s
    /// with different precisions are compared numerically, treating
    /// the smaller-precision mantissa as zero-extended.
    #[must_use]
    pub fn total_cmp(&self, other: &Self) -> Ordering {
        let lhs = total_kind_tag(&self.class);
        let rhs = total_kind_tag(&other.class);
        match lhs.cmp(&rhs) {
            Ordering::Equal => total_cmp_within_kind(&self.class, &other.class),
            ord => ord,
        }
    }

    /// IEEE 754-2019 §9.6 `minimumNumber`: smaller of two values
    /// treating quiet NaN as missing data.
    ///
    /// If exactly one operand is a quiet NaN, returns the other
    /// (NaN-as-missing). If both are quiet NaN, returns
    /// `(self.clone(), Status::OK)`. If either is signaling NaN,
    /// raises [`Status::INVALID`] and returns a quiet NaN.
    #[must_use]
    pub fn min(&self, other: &Self) -> (Self, Status) {
        if self.is_signaling_nan() || other.is_signaling_nan() {
            // For 1a we surface the INVALID flag and return a quiet
            // NaN at the higher of the two precisions. Slice 1b
            // refines the NaN payload-propagation rules.
            let prec = self.precision.max(other.precision);
            let nan = BigFloat::try_new_quiet_nan(Sign::Positive, prec, &[])
                .expect("BigFloat invariant: precision >= 1");
            return (nan, Status::INVALID);
        }
        if self.is_quiet_nan() && other.is_quiet_nan() {
            return (self.clone(), Status::OK);
        }
        if self.is_quiet_nan() {
            return (other.clone(), Status::OK);
        }
        if other.is_quiet_nan() {
            return (self.clone(), Status::OK);
        }
        match numeric_cmp_no_nan(self, other) {
            Ordering::Less | Ordering::Equal => (self.clone(), Status::OK),
            Ordering::Greater => (other.clone(), Status::OK),
        }
    }

    /// IEEE 754-2019 §9.6 `maximumNumber`: larger of two values
    /// treating quiet NaN as missing data.
    ///
    /// Symmetric to [`min`](Self::min).
    #[must_use]
    pub fn max(&self, other: &Self) -> (Self, Status) {
        if self.is_signaling_nan() || other.is_signaling_nan() {
            let prec = self.precision.max(other.precision);
            let nan = BigFloat::try_new_quiet_nan(Sign::Positive, prec, &[])
                .expect("BigFloat invariant: precision >= 1");
            return (nan, Status::INVALID);
        }
        if self.is_quiet_nan() && other.is_quiet_nan() {
            return (self.clone(), Status::OK);
        }
        if self.is_quiet_nan() {
            return (other.clone(), Status::OK);
        }
        if other.is_quiet_nan() {
            return (self.clone(), Status::OK);
        }
        match numeric_cmp_no_nan(self, other) {
            Ordering::Greater | Ordering::Equal => (self.clone(), Status::OK),
            Ordering::Less => (other.clone(), Status::OK),
        }
    }
}

/// Map a [`Class`] variant to a kind-tag for total ordering.
///
/// Tags align with IEEE 754-2019 §5.10's totalOrder rules:
/// - Negative NaNs come first; `-qNaN` precedes `-sNaN` (because
///   the negation reverses encoding-based ordering).
/// - `-∞`, then negative normals (compared by magnitude descending),
///   then `-0`.
/// - `+0`, then positive normals (compared by magnitude ascending),
///   then `+∞`.
/// - Positive NaNs come last; `+sNaN` precedes `+qNaN`.
///
/// Same-tag values are tie-broken inside `total_cmp_within_kind`.
#[cfg(feature = "big")]
fn total_kind_tag(class: &Class) -> i8 {
    match class {
        Class::Nan {
            quiet: true,
            sign: Sign::Negative,
            ..
        } => -10,
        Class::Nan {
            quiet: false,
            sign: Sign::Negative,
            ..
        } => -9,
        Class::Infinity {
            sign: Sign::Negative,
        } => -8,
        Class::Normal {
            sign: Sign::Negative,
            ..
        } => -7,
        Class::Zero {
            sign: Sign::Negative,
        } => -6,
        Class::Zero {
            sign: Sign::Positive,
        } => 6,
        Class::Normal {
            sign: Sign::Positive,
            ..
        } => 7,
        Class::Infinity {
            sign: Sign::Positive,
        } => 8,
        Class::Nan {
            quiet: false,
            sign: Sign::Positive,
            ..
        } => 9,
        Class::Nan {
            quiet: true,
            sign: Sign::Positive,
            ..
        } => 10,
    }
}

/// Tie-break two values that have the same [`total_kind_tag`].
///
/// At call time both classes share a kind, sign, and (for NaNs)
/// quiet-vs-signaling status. Compare the remaining encoding-bits
/// (mantissa for normals, payload for NaNs).
#[cfg(feature = "big")]
fn total_cmp_within_kind(a: &Class, b: &Class) -> Ordering {
    match (a, b) {
        (Class::Zero { .. }, Class::Zero { .. }) => Ordering::Equal,
        (Class::Infinity { .. }, Class::Infinity { .. }) => Ordering::Equal,
        (
            Class::Nan {
                payload: pa,
                sign: sa,
                ..
            },
            Class::Nan {
                payload: pb,
                sign: _sb,
                ..
            },
        ) => {
            let raw = limb_cmp_aligned(pa, pb);
            // Negative-signed NaNs reverse the encoding order:
            // a "more negative" NaN has the larger encoding, so we
            // flip the comparison.
            if matches!(sa, Sign::Negative) {
                raw.reverse()
            } else {
                raw
            }
        }
        (
            Class::Normal {
                sign: sa,
                exponent: ea,
                mantissa: ma,
            },
            Class::Normal {
                sign: _sb,
                exponent: eb,
                mantissa: mb,
            },
        ) => {
            let mag = match ea.cmp(eb) {
                Ordering::Equal => limb_cmp_aligned(ma, mb),
                ord => ord,
            };
            if matches!(sa, Sign::Negative) {
                mag.reverse()
            } else {
                mag
            }
        }
        // Other combinations are impossible at this site because
        // the kind tag is identical.
        _ => Ordering::Equal,
    }
}

/// Numeric magnitude comparison for two non-NaN [`BigFloat`]s.
///
/// `+0` and `-0` compare equal (IEEE 754-2019 §5.11 numeric
/// comparison rule). Cross-precision comparison falls out of the
/// limb-aligned mantissa comparison.
#[cfg(feature = "big")]
fn numeric_cmp_no_nan(a: &BigFloat, b: &BigFloat) -> Ordering {
    debug_assert!(!a.is_nan() && !b.is_nan());

    // Numeric ±0 == ±0 rule first (the special case where total_cmp
    // and partial_cmp diverge).
    if a.is_zero() && b.is_zero() {
        return Ordering::Equal;
    }

    // Different signs: negative < positive (with zero handled above).
    let a_neg = a.is_sign_negative();
    let b_neg = b.is_sign_negative();
    if a_neg && !b_neg {
        return Ordering::Less;
    }
    if !a_neg && b_neg {
        return Ordering::Greater;
    }

    // Same sign. Compare magnitudes; flip the result for negative.
    let mag = magnitude_cmp(a, b);
    if a_neg {
        mag.reverse()
    } else {
        mag
    }
}

/// Magnitude comparison ignoring sign: `|a| vs |b|`.
#[cfg(feature = "big")]
fn magnitude_cmp(a: &BigFloat, b: &BigFloat) -> Ordering {
    use Ordering as O;
    match (&a.class, &b.class) {
        (Class::Zero { .. }, Class::Zero { .. }) => O::Equal,
        (Class::Zero { .. }, _) => O::Less,
        (_, Class::Zero { .. }) => O::Greater,
        (Class::Infinity { .. }, Class::Infinity { .. }) => O::Equal,
        (Class::Infinity { .. }, _) => O::Greater,
        (_, Class::Infinity { .. }) => O::Less,
        (
            Class::Normal {
                exponent: ea,
                mantissa: ma,
                ..
            },
            Class::Normal {
                exponent: eb,
                mantissa: mb,
                ..
            },
        ) => match ea.cmp(eb) {
            O::Equal => limb_cmp_aligned(ma, mb),
            ord => ord,
        },
        // NaN cases excluded by caller's debug_assert.
        _ => O::Equal,
    }
}

/// Compare two limb arrays from MSL to LSL, treating the shorter
/// as zero-extended at the LSL end.
///
/// Both inputs are expected to be in canonical form (top-bit-set
/// for normalized non-zero mantissas; LSL-padding zeros for the
/// precision granularity invariant).
#[cfg(feature = "big")]
fn limb_cmp_aligned(a: &[u64], b: &[u64]) -> Ordering {
    let len = a.len().max(b.len());
    for offset in 1..=len {
        let av = if offset <= a.len() {
            a[a.len() - offset]
        } else {
            0
        };
        let bv = if offset <= b.len() {
            b[b.len() - offset]
        } else {
            0
        };
        let ord = av.cmp(&bv);
        if ord != Ordering::Equal {
            return ord;
        }
    }
    Ordering::Equal
}

#[cfg(test)]
#[cfg(feature = "big")]
mod tests {
    use super::*;

    fn at_each_precision<F: Fn(u32)>(f: F) {
        for &p in &[1, 53, 113, 256] {
            f(p);
        }
    }

    #[test]
    fn total_cmp_reflexive_for_constants() {
        at_each_precision(|p| {
            for v in [
                BigFloat::try_new_zero(Sign::Positive, p).unwrap(),
                BigFloat::try_new_zero(Sign::Negative, p).unwrap(),
                BigFloat::try_new_infinity(Sign::Positive, p).unwrap(),
                BigFloat::try_new_infinity(Sign::Negative, p).unwrap(),
                BigFloat::try_new_quiet_nan(Sign::Positive, p, &[]).unwrap(),
                BigFloat::try_new_quiet_nan(Sign::Negative, p, &[]).unwrap(),
                BigFloat::try_new_signaling_nan(Sign::Positive, p, &[]).unwrap(),
                BigFloat::try_new_signaling_nan(Sign::Negative, p, &[]).unwrap(),
            ] {
                assert_eq!(v.total_cmp(&v), Ordering::Equal);
            }
        });
    }

    #[test]
    fn total_cmp_orders_extremes() {
        let p = 53;
        let neg_qnan = BigFloat::try_new_quiet_nan(Sign::Negative, p, &[]).unwrap();
        let neg_snan = BigFloat::try_new_signaling_nan(Sign::Negative, p, &[]).unwrap();
        let neg_inf = BigFloat::try_new_infinity(Sign::Negative, p).unwrap();
        let neg_one = BigFloat::try_from_i64_exact(-1, p).unwrap();
        let neg_zero = BigFloat::try_new_zero(Sign::Negative, p).unwrap();
        let pos_zero = BigFloat::try_new_zero(Sign::Positive, p).unwrap();
        let pos_one = BigFloat::try_from_i64_exact(1, p).unwrap();
        let pos_inf = BigFloat::try_new_infinity(Sign::Positive, p).unwrap();
        let pos_snan = BigFloat::try_new_signaling_nan(Sign::Positive, p, &[]).unwrap();
        let pos_qnan = BigFloat::try_new_quiet_nan(Sign::Positive, p, &[]).unwrap();

        let ordered = [
            &neg_qnan, &neg_snan, &neg_inf, &neg_one, &neg_zero, &pos_zero, &pos_one, &pos_inf,
            &pos_snan, &pos_qnan,
        ];
        for win in ordered.windows(2) {
            assert_eq!(
                win[0].total_cmp(win[1]),
                Ordering::Less,
                "expected {:?} < {:?}",
                win[0].ieee_class(),
                win[1].ieee_class(),
            );
        }
    }

    #[test]
    fn neg_zero_total_lt_pos_zero() {
        let p = 53;
        let nz = BigFloat::try_new_zero(Sign::Negative, p).unwrap();
        let pz = BigFloat::try_new_zero(Sign::Positive, p).unwrap();
        assert_eq!(nz.total_cmp(&pz), Ordering::Less);
        assert_eq!(pz.total_cmp(&nz), Ordering::Greater);
    }

    #[test]
    fn partial_cmp_zero_equality() {
        let p = 53;
        let nz = BigFloat::try_new_zero(Sign::Negative, p).unwrap();
        let pz = BigFloat::try_new_zero(Sign::Positive, p).unwrap();
        // Per IEEE 754-2019 §5.11, ±0 compare equal numerically.
        let (ord, status) = nz.partial_cmp(&pz);
        assert_eq!(ord, Some(Ordering::Equal));
        assert!(status.is_ok());
    }

    #[test]
    fn partial_cmp_quiet_nan_returns_none_no_flag() {
        let p = 53;
        let q = BigFloat::try_new_quiet_nan(Sign::Positive, p, &[]).unwrap();
        let one = BigFloat::try_from_i64_exact(1, p).unwrap();
        let (ord, status) = q.partial_cmp(&one);
        assert_eq!(ord, None);
        assert!(status.is_ok()); // no INVALID for qNaN
        let (ord2, status2) = one.partial_cmp(&q);
        assert_eq!(ord2, None);
        assert!(status2.is_ok());
    }

    #[test]
    fn partial_cmp_signaling_nan_raises_invalid() {
        let p = 53;
        let s = BigFloat::try_new_signaling_nan(Sign::Positive, p, &[]).unwrap();
        let one = BigFloat::try_from_i64_exact(1, p).unwrap();
        let (ord, status) = s.partial_cmp(&one);
        assert_eq!(ord, None);
        assert!(status.invalid());
        let (ord2, status2) = one.partial_cmp(&s);
        assert_eq!(ord2, None);
        assert!(status2.invalid());
    }

    #[test]
    fn partial_cmp_orders_finite_values() {
        let p = 53;
        let one = BigFloat::try_from_i64_exact(1, p).unwrap();
        let two = BigFloat::try_from_i64_exact(2, p).unwrap();
        let neg_two = BigFloat::try_from_i64_exact(-2, p).unwrap();
        assert_eq!(one.partial_cmp(&two).0, Some(Ordering::Less));
        assert_eq!(two.partial_cmp(&one).0, Some(Ordering::Greater));
        assert_eq!(one.partial_cmp(&one).0, Some(Ordering::Equal));
        assert_eq!(neg_two.partial_cmp(&one).0, Some(Ordering::Less));
        assert_eq!(one.partial_cmp(&neg_two).0, Some(Ordering::Greater));
    }

    #[test]
    fn partial_cmp_infinities() {
        let p = 53;
        let pi = BigFloat::try_new_infinity(Sign::Positive, p).unwrap();
        let ni = BigFloat::try_new_infinity(Sign::Negative, p).unwrap();
        let one = BigFloat::try_from_i64_exact(1, p).unwrap();
        assert_eq!(pi.partial_cmp(&one).0, Some(Ordering::Greater));
        assert_eq!(one.partial_cmp(&pi).0, Some(Ordering::Less));
        assert_eq!(ni.partial_cmp(&one).0, Some(Ordering::Less));
        assert_eq!(pi.partial_cmp(&pi).0, Some(Ordering::Equal));
        assert_eq!(pi.partial_cmp(&ni).0, Some(Ordering::Greater));
    }

    #[test]
    fn cross_precision_comparison() {
        // Same numeric value at different precisions compares Equal.
        let one_53 = BigFloat::try_from_i64_exact(1, 53).unwrap();
        let one_128 = BigFloat::try_from_i64_exact(1, 128).unwrap();
        assert_eq!(one_53.partial_cmp(&one_128).0, Some(Ordering::Equal));
        assert_eq!(one_53.total_cmp(&one_128), Ordering::Equal);

        // Different values across precisions still compare correctly.
        let two_128 = BigFloat::try_from_i64_exact(2, 128).unwrap();
        assert_eq!(one_53.partial_cmp(&two_128).0, Some(Ordering::Less));
        assert_eq!(one_53.total_cmp(&two_128), Ordering::Less);
    }

    #[test]
    fn min_max_basic() {
        let p = 53;
        let one = BigFloat::try_from_i64_exact(1, p).unwrap();
        let two = BigFloat::try_from_i64_exact(2, p).unwrap();
        let (mn, ms) = one.min(&two);
        assert_eq!(mn, one);
        assert!(ms.is_ok());
        let (mx, _) = one.max(&two);
        assert_eq!(mx, two);
    }

    #[test]
    fn min_max_nan_handling() {
        let p = 53;
        let q = BigFloat::try_new_quiet_nan(Sign::Positive, p, &[]).unwrap();
        let one = BigFloat::try_from_i64_exact(1, p).unwrap();
        // qNaN treated as missing data: min(qNaN, 1) == 1.
        let (m, s) = q.min(&one);
        assert_eq!(m, one);
        assert!(s.is_ok());
        let (m2, s2) = one.min(&q);
        assert_eq!(m2, one);
        assert!(s2.is_ok());
        // Both qNaN: min returns self.
        let (m3, s3) = q.min(&q);
        assert!(m3.is_quiet_nan());
        assert!(s3.is_ok());
        // sNaN: raises INVALID and returns qNaN.
        let s_nan = BigFloat::try_new_signaling_nan(Sign::Positive, p, &[]).unwrap();
        let (m4, s4) = s_nan.min(&one);
        assert!(m4.is_quiet_nan());
        assert!(s4.invalid());
    }

    #[test]
    fn total_cmp_is_total_order() {
        // Reflexive, antisymmetric, transitive on three values.
        let p = 53;
        let a = BigFloat::try_from_i64_exact(-3, p).unwrap();
        let b = BigFloat::try_new_zero(Sign::Negative, p).unwrap();
        let c = BigFloat::try_from_i64_exact(7, p).unwrap();

        // Reflexivity.
        assert_eq!(a.total_cmp(&a), Ordering::Equal);
        assert_eq!(b.total_cmp(&b), Ordering::Equal);
        assert_eq!(c.total_cmp(&c), Ordering::Equal);

        // Antisymmetry.
        assert_eq!(a.total_cmp(&c), Ordering::Less);
        assert_eq!(c.total_cmp(&a), Ordering::Greater);

        // Transitivity (a < b < c implies a < c).
        assert_eq!(a.total_cmp(&b), Ordering::Less);
        assert_eq!(b.total_cmp(&c), Ordering::Less);
        assert_eq!(a.total_cmp(&c), Ordering::Less);
    }
}
