//! `cot(x)`, `sec(x)`, `csc(x)`: the reciprocal circular functions
//! (DLMF §4.14; ADR-0032, ADR-0056).
//!
//! These are direct primary kernels, not `1 / tan(x)`, `1 / cos(x)`,
//! `1 / sin(x)` applied to an already-rounded sibling: composing two
//! correctly-rounded operations double-rounds and can be 1 ULP off for
//! hard-to-round inputs (ADR-0032). The blessed pattern, proven by
//! [`super::tan`], computes `sin` and `cos` at one inflated Ziv working
//! precision inside the eval closure, forms the reciprocal/ratio at that
//! working precision, and lets [`super::ziv::ziv_round`] perform the
//! single rounding to the target.
//!
//! All three share [`super::tan`]'s machinery: the Payne-Hanek reduction
//! [`super::trig_reduce::reduce`] (which already inflates its internal
//! precision with the argument magnitude, so large arguments are handled),
//! the [`super::sin::sin_taylor`] / [`super::sin::cos_taylor`] series, the
//! range-cap pre-check at `target + ZIV_BASE_GUARD`, and the post-Ziv
//! NaN→INVALID surfacing. No `cancellation_boosted` path is needed: `cot`
//! at its zeros (the odd multiples of π/2) reduces to a small accurately
//! computed quotient, and near the poles the reduced reciprocal is a
//! large-but-finite value, exactly as [`super::tan`] behaves near its
//! poles.
//!
//! Quadrant identities (`s = sin(r)`, `c = cos(r)`, `r ∈ [−π/4, π/4]`):
//!
//! | q | sin(x) | cos(x) | cot=cos/sin | sec=1/cos | csc=1/sin |
//! |---|--------|--------|-------------|-----------|-----------|
//! | 0 | +s     | +c     | c/s         | 1/c       | 1/s       |
//! | 1 | +c     | −s     | −s/c        | −1/s      | 1/c       |
//! | 2 | −s     | −c     | c/s         | −1/c      | −1/s      |
//! | 3 | −c     | +s     | −s/c        | 1/s       | −1/c      |
//!
//! Special cases (DLMF §4.14; the pole conventions match pfloat's
//! existing `Ci(+0) = −∞` / `K_n(+0) = +∞` poles):
//!
//! - `cot(±0)` and `csc(±0)` are poles: `±∞` (sign of the zero, both are
//!   odd) with `DIV_BY_ZERO`. `sec(±0) = 1` exactly (`sec` is even).
//! - `cot(±∞)`, `sec(±∞)`, `csc(±∞)` = qNaN + `INVALID`.
//! - NaN propagates; sNaN raises `INVALID`.
//! - `|x|` past the reduction table budget: qNaN + `INVALID`.

use crate::big::{BigFloat, BuildError};
use crate::class::Class;
use crate::rounding::RoundingMode;
use crate::sign::Sign;
use crate::status::{auto_raise, Status};

#[cfg(feature = "fixed")]
use crate::fixed::FixedFloat;
#[cfg(feature = "fixed")]
use crate::mantissa::limbs_for;

use super::sin::{cos_taylor, sin_taylor};
use super::trig_reduce::{reduce, Reduction};
use super::ziv::{ziv_round, ZIV_BASE_GUARD};
use super::ziv_calibration::{COT_ERROR_GUARD, CSC_ERROR_GUARD, SEC_ERROR_GUARD};

macro_rules! reciprocal_api {
    ($name:ident, $round:ident, $with:ident, $round_with:ident, $kernel:ident, $disp:literal) => {
        impl BigFloat {
            #[doc = concat!("`", $disp, "(self)` rounded under `mode` to `self.precision`.")]
            #[must_use]
            pub fn $name(&self, mode: RoundingMode) -> (Self, Status) {
                self.$round(self.precision, mode)
                    .expect("self.precision >= 1 by invariant")
            }

            #[doc = concat!("`", $disp, "(self)` with explicit result precision.")]
            pub fn $round(
                &self,
                target_precision: u32,
                mode: RoundingMode,
            ) -> Result<(Self, Status), BuildError> {
                if target_precision == 0 {
                    return Err(BuildError::PrecisionZero);
                }
                Ok($kernel(self, target_precision, mode))
            }

            #[doc = concat!("`", $disp, "` accumulating into a caller-supplied flag bag.")]
            #[must_use]
            pub fn $with(&self, mode: RoundingMode, flags: &mut Status) -> Self {
                let (value, status) = self.$name(mode);
                *flags |= status;
                value
            }

            #[doc = concat!("`", $disp, "_round` accumulating into a caller-supplied flag bag.")]
            pub fn $round_with(
                &self,
                target_precision: u32,
                mode: RoundingMode,
                flags: &mut Status,
            ) -> Result<Self, BuildError> {
                let (value, status) = self.$round(target_precision, mode)?;
                *flags |= status;
                Ok(value)
            }
        }

        #[cfg(feature = "fixed")]
        impl<const PREC: u32> FixedFloat<PREC>
        where
            [(); limbs_for(PREC)]:,
        {
            #[doc = concat!("`", $disp, "(self)` for `FixedFloat`. Delegates to `BigFloat`.")]
            #[must_use]
            pub fn $name(&self, mode: RoundingMode) -> (Self, Status) {
                let (big, status) = self.to_big().$name(mode);
                (
                    Self::try_from_big_exact(big).expect("precision matches"),
                    status,
                )
            }
        }
    };
}

reciprocal_api!(
    cot,
    cot_round,
    cot_with_flags,
    cot_round_with_flags,
    cot_kernel,
    "cot"
);
reciprocal_api!(
    sec,
    sec_round,
    sec_with_flags,
    sec_round_with_flags,
    sec_kernel,
    "sec"
);
reciprocal_api!(
    csc,
    csc_round,
    csc_with_flags,
    csc_round_with_flags,
    csc_kernel,
    "csc"
);

/// NaN-input result shared by all three kernels: a quiet NaN propagates;
/// a signaling NaN raises `INVALID` and returns a canonical quiet NaN.
fn nan_result(quiet: bool, sign: Sign, payload: &[u64], target: u32) -> (BigFloat, Status) {
    if quiet {
        (
            BigFloat::try_new_quiet_nan(sign, target, payload).expect("precision >= 1"),
            Status::OK,
        )
    } else {
        auto_raise(Status::INVALID);
        (
            BigFloat::try_new_quiet_nan(Sign::Positive, target, &[]).expect("precision >= 1"),
            Status::INVALID,
        )
    }
}

/// `qNaN + INVALID` shared by the `±∞` arms.
fn invalid_nan(target: u32) -> (BigFloat, Status) {
    auto_raise(Status::INVALID);
    (
        BigFloat::try_new_quiet_nan(Sign::Positive, target, &[]).expect("precision >= 1"),
        Status::INVALID,
    )
}

/// The shared reduce-precheck + `ziv_round` + post-Ziv NaN→INVALID
/// plumbing (mirrors [`super::tan`]). `compose(quadrant, r, w)` builds the
/// reciprocal/ratio at the Ziv working precision `w`.
fn reciprocal_via_ziv(
    x: &BigFloat,
    target: u32,
    mode: RoundingMode,
    error_guard: u32,
    compose: impl Fn(u8, &BigFloat, u32) -> BigFloat,
) -> (BigFloat, Status) {
    // Range-cap at the first Ziv working precision (pf-1axr; see sin.rs).
    let ziv_first_working = target.saturating_add(ZIV_BASE_GUARD);
    if reduce(x, ziv_first_working).is_none() {
        return invalid_nan(target);
    }

    let (result, status) = ziv_round(
        |w| match reduce(x, w) {
            Some(Reduction { quadrant, r }) => compose(quadrant, &r, w),
            None => BigFloat::try_new_quiet_nan(Sign::Positive, w, &[]).expect("precision >= 1"),
        },
        target,
        mode,
        error_guard,
    );

    // A Ziv iteration's reduce hitting the table cap propagated NaN through
    // the driver; surface as INVALID, matching the pre-check path.
    if matches!(result.class, Class::Nan { .. }) && !status.invalid() {
        let merged = status.merge(Status::INVALID);
        auto_raise(Status::INVALID);
        return (result, merged);
    }
    // cot/sec/csc of a finite normal x are transcendental (Lindemann–
    // Weierstrass), hence irrational, hence INEXACT even where the result
    // rounds onto a grid value (pf-uqd1, ADR-0063). The only exact input,
    // sec(±0) = 1, and the cot/csc poles at 0 are dispatched in the
    // kernels above before this helper runs.
    let status = if matches!(result.class, Class::Normal { .. }) {
        status | Status::INEXACT
    } else {
        status
    };
    auto_raise(status);
    (result, status)
}

fn cot_kernel(x: &BigFloat, target: u32, mode: RoundingMode) -> (BigFloat, Status) {
    match &x.class {
        Class::Nan {
            quiet,
            sign,
            payload,
        } => nan_result(*quiet, *sign, payload, target),
        Class::Zero { sign } => {
            // cot(±0) is a pole: ±∞ (cot is odd) with DIV_BY_ZERO.
            auto_raise(Status::DIV_BY_ZERO);
            (
                BigFloat::try_new_infinity(*sign, target).expect("precision >= 1"),
                Status::DIV_BY_ZERO,
            )
        }
        Class::Infinity { .. } => invalid_nan(target),
        Class::Normal { .. } => {
            reciprocal_via_ziv(x, target, mode, COT_ERROR_GUARD, |quadrant, r, w| {
                let s = sin_taylor(r, w);
                let c = cos_taylor(r, w);
                match quadrant {
                    0 | 2 => c.div(&s, RoundingMode::NearestEven).0,
                    // q1: −s/c ; q3: s/(−c) = −s/c.
                    _ => s.negated().div(&c, RoundingMode::NearestEven).0,
                }
            })
        }
    }
}

fn sec_kernel(x: &BigFloat, target: u32, mode: RoundingMode) -> (BigFloat, Status) {
    match &x.class {
        Class::Nan {
            quiet,
            sign,
            payload,
        } => nan_result(*quiet, *sign, payload, target),
        Class::Zero { .. } => {
            // sec(±0) = 1 exactly (sec is even).
            (
                BigFloat::try_from_i64_exact(1, target).expect("precision >= 1"),
                Status::OK,
            )
        }
        Class::Infinity { .. } => invalid_nan(target),
        Class::Normal { .. } => {
            reciprocal_via_ziv(x, target, mode, SEC_ERROR_GUARD, |quadrant, r, w| {
                // sec = 1 / cos(x); cos(x) per quadrant.
                let cos_x = match quadrant {
                    0 => cos_taylor(r, w),
                    1 => sin_taylor(r, w).negated(),
                    2 => cos_taylor(r, w).negated(),
                    _ => sin_taylor(r, w),
                };
                let one = BigFloat::try_from_i64_exact(1, w).expect("precision >= 1");
                one.div(&cos_x, RoundingMode::NearestEven).0
            })
        }
    }
}

fn csc_kernel(x: &BigFloat, target: u32, mode: RoundingMode) -> (BigFloat, Status) {
    match &x.class {
        Class::Nan {
            quiet,
            sign,
            payload,
        } => nan_result(*quiet, *sign, payload, target),
        Class::Zero { sign } => {
            // csc(±0) is a pole: ±∞ (csc is odd) with DIV_BY_ZERO.
            auto_raise(Status::DIV_BY_ZERO);
            (
                BigFloat::try_new_infinity(*sign, target).expect("precision >= 1"),
                Status::DIV_BY_ZERO,
            )
        }
        Class::Infinity { .. } => invalid_nan(target),
        Class::Normal { .. } => {
            reciprocal_via_ziv(x, target, mode, CSC_ERROR_GUARD, |quadrant, r, w| {
                // csc = 1 / sin(x); sin(x) per quadrant.
                let sin_x = match quadrant {
                    0 => sin_taylor(r, w),
                    1 => cos_taylor(r, w),
                    2 => sin_taylor(r, w).negated(),
                    _ => cos_taylor(r, w).negated(),
                };
                let one = BigFloat::try_from_i64_exact(1, w).expect("precision >= 1");
                one.div(&sin_x, RoundingMode::NearestEven).0
            })
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use core::cmp::Ordering;

    fn close_at(v: &BigFloat, expected: &BigFloat, bits: u32) -> bool {
        let (diff, _) = v.sub(expected, RoundingMode::NearestEven);
        let abs_diff = diff.abs();
        if abs_diff.is_zero() {
            return true;
        }
        let p = v.precision().max(expected.precision());
        let one = BigFloat::try_from_i64_exact(1, p).unwrap();
        let two = BigFloat::try_from_i64_exact(2, p).unwrap();
        let abs_b = expected.abs();
        let mut bound = if abs_b.is_zero() { one } else { abs_b };
        for _ in 0..bits {
            bound = bound.div(&two, RoundingMode::NearestEven).0;
        }
        matches!(
            abs_diff.partial_cmp(&bound).0,
            Some(Ordering::Less | Ordering::Equal)
        )
    }

    fn from_i64(n: i64, p: u32) -> BigFloat {
        BigFloat::try_from_i64_exact(n, p).unwrap()
    }

    #[test]
    fn cot_pi_over_4_is_one() {
        let pi_2 = super::super::pi_over_2_at(160);
        let (pi_4, _) = pi_2.div(&from_i64(2, 160), RoundingMode::NearestEven);
        let (r, _) = pi_4.cot(RoundingMode::NearestEven);
        assert!(close_at(&r, &from_i64(1, 160), 140), "cot(pi/4) = {r}");
    }

    #[test]
    fn sec_zero_is_one_exact() {
        for s in [Sign::Positive, Sign::Negative] {
            let z = BigFloat::try_new_zero(s, 53).unwrap();
            let (r, st) = z.sec(RoundingMode::NearestEven);
            assert!(st.is_ok() && !st.inexact());
            assert_eq!(r.partial_cmp(&from_i64(1, 53)).0, Some(Ordering::Equal));
        }
    }

    #[test]
    fn csc_pi_over_2_is_one() {
        let pi_2 = super::super::pi_over_2_at(160);
        let (r, _) = pi_2.csc(RoundingMode::NearestEven);
        assert!(close_at(&r, &from_i64(1, 160), 140), "csc(pi/2) = {r}");
    }

    #[test]
    fn cot_is_reciprocal_of_tan() {
        // tan(x)·cot(x) = 1 away from the poles.
        let x = from_i64(1, 160);
        let (t, _) = x.tan(RoundingMode::NearestEven);
        let (c, _) = x.cot(RoundingMode::NearestEven);
        let (prod, _) = t.mul(&c, RoundingMode::NearestEven);
        assert!(close_at(&prod, &from_i64(1, 160), 150), "tan·cot = {prod}");
    }

    #[test]
    fn sec_matches_one_over_cos() {
        let x = from_i64(1, 160);
        let (sx, _) = x.sec(RoundingMode::NearestEven);
        let (cx, _) = x.cos(RoundingMode::NearestEven);
        let (recip, _) = from_i64(1, 160).div(&cx, RoundingMode::NearestEven);
        assert!(close_at(&sx, &recip, 150), "sec vs 1/cos");
    }

    #[test]
    fn csc_matches_one_over_sin() {
        let x = from_i64(1, 160);
        let (cx, _) = x.csc(RoundingMode::NearestEven);
        let (sx, _) = x.sin(RoundingMode::NearestEven);
        let (recip, _) = from_i64(1, 160).div(&sx, RoundingMode::NearestEven);
        assert!(close_at(&cx, &recip, 150), "csc vs 1/sin");
    }

    #[test]
    fn cot_csc_odd_sec_even() {
        let x = from_i64(2, 160);
        let nx = from_i64(-2, 160);
        let (c, _) = x.cot(RoundingMode::NearestEven);
        let (cn, _) = nx.cot(RoundingMode::NearestEven);
        assert!(close_at(&cn, &c.negated(), 150), "cot odd");
        let (k, _) = x.csc(RoundingMode::NearestEven);
        let (kn, _) = nx.csc(RoundingMode::NearestEven);
        assert!(close_at(&kn, &k.negated(), 150), "csc odd");
        let (s, _) = x.sec(RoundingMode::NearestEven);
        let (sn, _) = nx.sec(RoundingMode::NearestEven);
        assert!(close_at(&sn, &s, 150), "sec even");
    }

    #[test]
    fn cot_csc_zero_are_poles() {
        for s in [Sign::Positive, Sign::Negative] {
            let z = BigFloat::try_new_zero(s, 53).unwrap();
            let (rc, sc) = z.cot(RoundingMode::NearestEven);
            assert!(rc.is_infinite() && sc.div_by_zero());
            assert_eq!(rc.is_sign_negative(), matches!(s, Sign::Negative));
            let (rk, sk) = z.csc(RoundingMode::NearestEven);
            assert!(rk.is_infinite() && sk.div_by_zero());
            assert_eq!(rk.is_sign_negative(), matches!(s, Sign::Negative));
        }
    }

    #[test]
    fn reciprocals_infinity_is_invalid() {
        for s in [Sign::Positive, Sign::Negative] {
            let inf = BigFloat::try_new_infinity(s, 53).unwrap();
            for (r, st) in [
                inf.cot(RoundingMode::NearestEven),
                inf.sec(RoundingMode::NearestEven),
                inf.csc(RoundingMode::NearestEven),
            ] {
                assert!(r.is_quiet_nan() && st.invalid());
            }
        }
    }

    #[test]
    fn reciprocals_nan_propagation() {
        let q = BigFloat::try_new_quiet_nan(Sign::Positive, 53, &[]).unwrap();
        assert!(q.cot(RoundingMode::NearestEven).0.is_quiet_nan());
        assert!(q.sec(RoundingMode::NearestEven).0.is_quiet_nan());
        assert!(q.csc(RoundingMode::NearestEven).0.is_quiet_nan());
        let sn = BigFloat::try_new_signaling_nan(Sign::Positive, 53, &[]).unwrap();
        for (r, st) in [
            sn.cot(RoundingMode::NearestEven),
            sn.sec(RoundingMode::NearestEven),
            sn.csc(RoundingMode::NearestEven),
        ] {
            assert!(r.is_quiet_nan() && st.invalid());
        }
    }

    #[test]
    fn reciprocals_out_of_range_is_invalid() {
        // 2^5000 exceeds the reduction table budget.
        let one = from_i64(1, 53);
        let two = from_i64(2, 53);
        let mut big = one;
        for _ in 0..5000 {
            big = big.mul(&two, RoundingMode::NearestEven).0;
        }
        for (r, st) in [
            big.cot(RoundingMode::NearestEven),
            big.sec(RoundingMode::NearestEven),
            big.csc(RoundingMode::NearestEven),
        ] {
            assert!(r.is_quiet_nan() && st.invalid());
        }
    }

    #[test]
    fn reciprocals_round_rejects_zero_precision() {
        assert_eq!(
            from_i64(1, 53).cot_round(0, RoundingMode::NearestEven),
            Err(BuildError::PrecisionZero)
        );
        assert_eq!(
            from_i64(1, 53).sec_round(0, RoundingMode::NearestEven),
            Err(BuildError::PrecisionZero)
        );
        assert_eq!(
            from_i64(1, 53).csc_round(0, RoundingMode::NearestEven),
            Err(BuildError::PrecisionZero)
        );
    }

    #[test]
    fn sec_csc_at_least_one_in_magnitude() {
        // |sec| >= 1 and |csc| >= 1 everywhere they are finite.
        let x = from_i64(1, 113);
        let one = from_i64(1, 113);
        let (s, _) = x.sec(RoundingMode::NearestEven);
        let (k, _) = x.csc(RoundingMode::NearestEven);
        assert!(matches!(
            s.abs().partial_cmp(&one).0,
            Some(Ordering::Greater | Ordering::Equal)
        ));
        assert!(matches!(
            k.abs().partial_cmp(&one).0,
            Some(Ordering::Greater | Ordering::Equal)
        ));
    }
}
