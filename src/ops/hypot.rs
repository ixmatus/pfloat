//! IEEE 754-2019 §9.2.1 `hypot(x, y) = sqrt(x² + y²)` for [`BigFloat`]
//! (ADR-0032, ADR-0056).
//!
//! `hypot` is a direct primary kernel, not the naive `sqrt(x*x + y*y)`
//! composition: squaring at the target precision loses half the input
//! precision, so the composed result is not correctly rounded for
//! hard-to-round inputs. The direct kernel computes the squares, sum, and
//! square root at an inflated Ziv working precision inside one
//! `ziv_round` closure and rounds once at the target.
//!
//! No scaling by `max(|x|, |y|)` (Moler's trick) is needed. That guards a
//! fixed-width exponent field against `x²` overflowing before the `sqrt`
//! pulls it back; `BigFloat`'s exponent is an `i64` with saturating
//! arithmetic, so `x²` never overflows the field for any finite operand,
//! and the only overflow that can occur is the final result saturating to
//! `±∞`, which the rounding pipeline already handles.
//!
//! Special cases per IEEE 754-2019 §9.2.1 (infinity is checked before
//! NaN — the counterintuitive rule):
//!
//! - `hypot(±∞, y)` and `hypot(x, ±∞)` are `+∞`, *even when the other
//!   operand is NaN*. A signaling-NaN operand still raises `INVALID` (the
//!   ∞ override fixes the value, not the §7.2 signal).
//! - With no infinity present: sNaN → qNaN + `INVALID`; qNaN propagates.
//! - Both finite: `sqrt(x² + y²)`, always non-negative. Symmetric and
//!   sign-independent: `hypot(x, y) = hypot(y, x) = hypot(|x|, |y|)`.

use crate::big::{BigFloat, BuildError};
use crate::class::Class;
use crate::rounding::RoundingMode;
use crate::sign::Sign;
use crate::status::{auto_raise, Status};

#[cfg(feature = "fixed")]
use crate::fixed::FixedFloat;
#[cfg(feature = "fixed")]
use crate::mantissa::limbs_for;

use crate::math::ziv::ziv_round_with_depth;
use crate::math::ziv_calibration::HYPOT_ERROR_GUARD;

impl BigFloat {
    /// IEEE 754-2019 `hypot(self, other)` = `sqrt(self² + other²)`,
    /// rounded under `mode` to `self.precision`.
    #[must_use]
    pub fn hypot(&self, other: &Self, mode: RoundingMode) -> (Self, Status) {
        self.hypot_round(other, self.precision, mode)
            .expect("self.precision >= 1 by invariant")
    }

    /// `hypot(self, other)` with an explicit result precision.
    pub fn hypot_round(
        &self,
        other: &Self,
        target_precision: u32,
        mode: RoundingMode,
    ) -> Result<(Self, Status), BuildError> {
        if target_precision == 0 {
            return Err(BuildError::PrecisionZero);
        }
        Ok(hypot_kernel(self, other, target_precision, mode))
    }

    /// `hypot` accumulating into a caller-supplied flag bag.
    #[must_use]
    pub fn hypot_with_flags(&self, other: &Self, mode: RoundingMode, flags: &mut Status) -> Self {
        let (value, status) = self.hypot(other, mode);
        *flags |= status;
        value
    }

    /// `hypot_round` accumulating into a caller-supplied flag bag.
    pub fn hypot_round_with_flags(
        &self,
        other: &Self,
        target_precision: u32,
        mode: RoundingMode,
        flags: &mut Status,
    ) -> Result<Self, BuildError> {
        let (value, status) = self.hypot_round(other, target_precision, mode)?;
        *flags |= status;
        Ok(value)
    }
}

#[cfg(feature = "fixed")]
impl<const PREC: u32> FixedFloat<PREC>
where
    [(); limbs_for(PREC)]:,
{
    /// `hypot(self, other)` for `FixedFloat`. Delegates to
    /// [`BigFloat::hypot`].
    #[must_use]
    pub fn hypot(&self, other: &Self, mode: RoundingMode) -> (Self, Status) {
        let (big, status) = self.to_big().hypot(&other.to_big(), mode);
        (
            Self::try_from_big_exact(big).expect("precision matches"),
            status,
        )
    }
}

fn hypot_kernel(x: &BigFloat, y: &BigFloat, target: u32, mode: RoundingMode) -> (BigFloat, Status) {
    // Infinity dominates the value, even against a NaN operand (§9.2.1).
    if x.is_infinite() || y.is_infinite() {
        let inf = BigFloat::try_new_infinity(Sign::Positive, target).expect("precision >= 1");
        if x.is_signaling_nan() || y.is_signaling_nan() {
            // The ∞ override fixes the value; the sNaN still signals.
            auto_raise(Status::INVALID);
            return (inf, Status::INVALID);
        }
        return (inf, Status::OK);
    }

    // No infinity: standard two-operand NaN propagation (sNaN → INVALID +
    // qNaN; qNaN propagates).
    if let Some((nan, status)) = super::propagate_nan2(x, y, target) {
        return (nan, status);
    }

    // Deep exponent gap (pf-71u2, ADR-0102): with gap g = e_big −
    // e_small, the truth is |big| + δ with δ = small²/(hypot + |big|)
    // ∈ (2^(e_big − 2g − 3), 2^(e_big − 2g + 1)]. The Ziv eval's sum
    // absorbs small² whenever 2g ≥ working, collapsing onto |big| —
    // exactly on-grid — so the interval test never converges and the
    // exhausted fall-through certified the collapsed value: falsely
    // EXACT whenever |big| rounds exactly at the target. Past the
    // representable band (2g ≥ max(p_big, target) + 6, two bits of
    // slack over the derived bound 2g > max + 4) δ sits strictly
    // inside the boundary-free zone above |big| (width at least
    // 2^(e_big − max − 2)), so the grow-direction infinitesimal
    // rounding is exact in every mode and honestly INEXACT (an exact
    // hypot at the target is impossible there: the true square's
    // bit-span exceeds the target). Inside the band the driver still
    // certifies (its cap clears the depth) and stays untouched.
    if let (Class::Normal { exponent: ex, .. }, Class::Normal { exponent: ey, .. }) =
        (&x.class, &y.class)
    {
        let (e_b, e_s, big) = if ex >= ey {
            (*ex, *ey, x)
        } else {
            (*ey, *ex, y)
        };
        let two_gap = e_b.saturating_sub(e_s).saturating_mul(2);
        let max_pt = i64::from(big.precision.max(target));
        // Rim guard (ADR-0102 verifier finding): round_with_infinitesimal
        // places its residue at e_b − (max_pt + 3) + 1; within that reach
        // of i64::MIN the placement saturates, base + ε becomes exactly
        // representable, and the rounding certifies a wrong value with
        // Status OK. Refuse the dispatch there (the driver fall-through
        // keeps the pre-existing rim behavior) until pf-a77o fixes the
        // residue placement at its root.
        if two_gap >= max_pt.saturating_add(6)
            && e_b >= i64::MIN.saturating_add(max_pt).saturating_add(5)
        {
            return crate::rounding::round_with_infinitesimal(
                &big.abs(),
                Sign::Positive,
                false, // magnitude grows: δ > 0 strictly (small ≠ 0)
                target,
                mode,
            );
        }
    }

    // Both finite: sqrt(x² + y²) at the Ziv working precision, rounded
    // once at the target.
    let (v, status) = ziv_round_with_depth(
        |w| {
            let xw = x
                .round_to_precision(w, RoundingMode::NearestEven)
                .expect("precision >= 1")
                .0;
            let yw = y
                .round_to_precision(w, RoundingMode::NearestEven)
                .expect("precision >= 1")
                .0;
            let x2 = xw.mul(&xw, RoundingMode::NearestEven).0;
            let y2 = yw.mul(&yw, RoundingMode::NearestEven).0;
            let s = x2.add(&y2, RoundingMode::NearestEven).0;
            s.sqrt(RoundingMode::NearestEven).0
        },
        target,
        mode,
        HYPOT_ERROR_GUARD,
        // Hole-residue certification depth (pf-fbjn, ADR-0104): a
        // big operand whose precision exceeds the dispatch arm with
        // a gap past the legacy cap resolves at the deep rung. The
        // truth's distance to any rounding-change point B is bounded
        // below through the exact dyadic truth² − B², whose bottom
        // bit the operands' precisions and the gap bound, so the
        // rung is certain outside exact-square coincidences. Lazy:
        // free unless the schedule exhausts.
        || match (&x.class, &y.class) {
            (Class::Normal { exponent: ex, .. }, Class::Normal { exponent: ey, .. }) => {
                let (big_p, small_p) = if ex >= ey {
                    (x.precision, y.precision)
                } else {
                    (y.precision, x.precision)
                };
                let two_gap =
                    u32::try_from(ex.saturating_sub(*ey).saturating_abs().saturating_mul(2))
                        .unwrap_or(u32::MAX);
                two_gap
                    .saturating_add(small_p.saturating_mul(2))
                    .max(big_p.saturating_mul(2))
                    .saturating_add(target)
                    .saturating_add(64)
            }
            _ => 0,
        },
    );

    // Exactness audit on claimed-exact results (pf-fbjn, ADR-0104).
    // The interval test certifies VALUES, not statuses: a collapsed
    // eval (small² absorbed below the working precision) can round
    // exactly, certify under the nearest modes — the interval sits
    // strictly inside the half-ulp — and claim OK while the truth is
    // inexact. hypot's exact set is decidable by construct-and-check:
    // v is exact iff v² = x² + y², all dyadic. The squares are exact
    // at twice their operands' precisions and the audit width covers
    // every genuinely-exact case's span (the difference x² − v² spans
    // at most ~2p + 2t bits around the shared binade, and a true
    // exact has x² − v² = −y² with span ≤ 2·p_small); a nonzero
    // residue, or any audit step rounding (unreachable for true
    // exacts), forces the honest INEXACT. Runs only when the driver
    // claims OK — the rare exact-result path.
    if status.is_ok() {
        let p_audit = x
            .precision
            .saturating_mul(2)
            .saturating_add(y.precision.saturating_mul(2))
            .saturating_add(target.saturating_mul(2))
            .saturating_add(64);
        let mut audit = Status::OK;
        let x2 = x
            .mul_round(x, p_audit, RoundingMode::NearestEven)
            .map(|(r, s)| {
                audit |= s;
                r
            })
            .expect("p_audit >= 1");
        let y2 = y
            .mul_round(y, p_audit, RoundingMode::NearestEven)
            .map(|(r, s)| {
                audit |= s;
                r
            })
            .expect("p_audit >= 1");
        let v2 = v
            .mul_round(&v, p_audit, RoundingMode::NearestEven)
            .map(|(r, s)| {
                audit |= s;
                r
            })
            .expect("p_audit >= 1");
        let (d1, s1) = x2
            .sub_round(&v2, p_audit, RoundingMode::NearestEven)
            .expect("p_audit >= 1");
        let (d2, s2) = d1
            .add_round(&y2, p_audit, RoundingMode::NearestEven)
            .expect("p_audit >= 1");
        audit |= s1 | s2;
        if !(audit.is_ok() && d2.is_zero()) {
            auto_raise(Status::INEXACT);
            return (v, status | Status::INEXACT);
        }
    }
    (v, status)
}

#[cfg(test)]
mod tests {
    use super::*;
    use core::cmp::Ordering;

    fn from_i64(n: i64, p: u32) -> BigFloat {
        BigFloat::try_from_i64_exact(n, p).unwrap()
    }

    fn eq(a: &BigFloat, b: &BigFloat) -> bool {
        matches!(a.partial_cmp(b).0, Some(Ordering::Equal))
    }

    #[test]
    fn hypot_three_four_is_five() {
        let (r, s) = from_i64(3, 53).hypot(&from_i64(4, 53), RoundingMode::NearestEven);
        assert!(s.is_ok());
        assert!(eq(&r, &from_i64(5, 53)));
    }

    #[test]
    fn hypot_five_twelve_is_thirteen() {
        let (r, _) = from_i64(5, 53).hypot(&from_i64(12, 53), RoundingMode::NearestEven);
        assert!(eq(&r, &from_i64(13, 53)));
    }

    #[test]
    fn hypot_with_zero_is_abs() {
        let zero = BigFloat::try_new_zero(Sign::Positive, 53).unwrap();
        let (r, _) = from_i64(7, 53).hypot(&zero, RoundingMode::NearestEven);
        assert!(eq(&r, &from_i64(7, 53)));
        // hypot(-7, 0) = 7 (sign-independent).
        let (r2, _) = from_i64(-7, 53).hypot(&zero, RoundingMode::NearestEven);
        assert!(eq(&r2, &from_i64(7, 53)));
    }

    #[test]
    fn hypot_zero_zero_is_pos_zero() {
        let pz = BigFloat::try_new_zero(Sign::Positive, 53).unwrap();
        let nz = BigFloat::try_new_zero(Sign::Negative, 53).unwrap();
        let (r, _) = nz.hypot(&nz, RoundingMode::NearestEven);
        assert!(r.is_zero() && r.is_sign_positive());
        let (r2, _) = pz.hypot(&nz, RoundingMode::NearestEven);
        assert!(r2.is_zero() && r2.is_sign_positive());
    }

    #[test]
    fn hypot_is_symmetric_and_sign_independent() {
        let a = from_i64(-3, 80);
        let b = from_i64(4, 80);
        let (ab, _) = a.hypot(&b, RoundingMode::NearestEven);
        let (ba, _) = b.hypot(&a, RoundingMode::NearestEven);
        assert!(eq(&ab, &ba));
        // |a|,|b| gives the same magnitude.
        let (pp, _) = from_i64(3, 80).hypot(&from_i64(4, 80), RoundingMode::NearestEven);
        assert!(eq(&ab, &pp));
        assert!(eq(&ab, &from_i64(5, 80)));
    }

    #[test]
    fn hypot_one_one_is_sqrt_two_inexact() {
        let one = from_i64(1, 53);
        let (r, s) = one.hypot(&one, RoundingMode::NearestEven);
        assert!(s.inexact());
        // sqrt(2) == hypot(1,1) bit-for-bit (both are sqrt of 2).
        let two = from_i64(2, 53);
        let (sqrt2, _) = two.sqrt(RoundingMode::NearestEven);
        assert!(eq(&r, &sqrt2));
    }

    #[test]
    fn hypot_infinity_dominates() {
        let pinf = BigFloat::try_new_infinity(Sign::Positive, 53).unwrap();
        let ninf = BigFloat::try_new_infinity(Sign::Negative, 53).unwrap();
        for inf in [&pinf, &ninf] {
            let (r, s) = inf.hypot(&from_i64(5, 53), RoundingMode::NearestEven);
            assert!(r.is_infinite() && r.is_sign_positive() && s.is_ok());
            let (r2, _) = from_i64(5, 53).hypot(inf, RoundingMode::NearestEven);
            assert!(r2.is_infinite() && r2.is_sign_positive());
        }
    }

    #[test]
    fn hypot_infinity_beats_qnan() {
        let pinf = BigFloat::try_new_infinity(Sign::Positive, 53).unwrap();
        let q = BigFloat::try_new_quiet_nan(Sign::Positive, 53, &[]).unwrap();
        let (r, s) = pinf.hypot(&q, RoundingMode::NearestEven);
        assert!(r.is_infinite() && r.is_sign_positive() && !s.invalid());
    }

    #[test]
    fn hypot_infinity_with_snan_still_signals() {
        let pinf = BigFloat::try_new_infinity(Sign::Positive, 53).unwrap();
        let sn = BigFloat::try_new_signaling_nan(Sign::Positive, 53, &[]).unwrap();
        let (r, s) = pinf.hypot(&sn, RoundingMode::NearestEven);
        assert!(r.is_infinite() && s.invalid());
    }

    #[test]
    fn hypot_qnan_propagates() {
        let q = BigFloat::try_new_quiet_nan(Sign::Positive, 53, &[]).unwrap();
        let (r, s) = q.hypot(&from_i64(3, 53), RoundingMode::NearestEven);
        assert!(r.is_quiet_nan() && !s.invalid());
    }

    #[test]
    fn hypot_snan_raises_invalid() {
        let sn = BigFloat::try_new_signaling_nan(Sign::Positive, 53, &[]).unwrap();
        let (r, s) = sn.hypot(&from_i64(3, 53), RoundingMode::NearestEven);
        assert!(r.is_quiet_nan() && s.invalid());
    }

    #[test]
    fn hypot_round_rejects_zero_precision() {
        assert_eq!(
            from_i64(3, 53).hypot_round(&from_i64(4, 53), 0, RoundingMode::NearestEven),
            Err(BuildError::PrecisionZero)
        );
    }

    #[test]
    fn hypot_higher_precision() {
        // 8² + 15² = 289 = 17².
        let (r, _) = from_i64(8, 200).hypot(&from_i64(15, 200), RoundingMode::NearestEven);
        assert_eq!(r.precision(), 200);
        assert!(eq(&r, &from_i64(17, 200)));
    }

    #[test]
    fn hypot_near_equal_operands() {
        // x ≈ y exercises the sum-of-squares path with no cancellation.
        let a = BigFloat::parse_str("1.0000001", 113, RoundingMode::NearestEven)
            .unwrap()
            .0;
        let b = from_i64(1, 113);
        let (r, _) = a.hypot(&b, RoundingMode::NearestEven);
        // result is between max(a,b) and a+b.
        assert_eq!(r.partial_cmp(&a).0, Some(Ordering::Greater));
        let (sum, _) = a.add(&b, RoundingMode::NearestEven);
        assert_eq!(r.partial_cmp(&sum).0, Some(Ordering::Less));
    }
}
