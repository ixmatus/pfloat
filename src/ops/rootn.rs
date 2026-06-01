//! IEEE 754-2019 §9.2.1 `rootn(x, n)`: the real `n`-th root, `x^(1/n)`,
//! for an integer order `n` (ADR-0032, ADR-0056).
//!
//! `rootn` is a direct primary kernel, not `pow(x, 1/n)`: `1/n` is not
//! representable, so the composition can never be correctly rounded, and
//! the even/odd domain rules below are not expressible through `pow`.
//!
//! Algorithm (no `exp`/`ln`, which would be the forbidden composition and
//! sit behind a feature `rootn` does not require): reduce to the positive
//! `m`-th root with `m = |n|`, compute it by Newton's iteration at the
//! Ziv working precision, then apply the sign (odd roots preserve it) and,
//! for `n < 0`, the reciprocal. One `ziv_round` rounds the assembled value
//! once at the target. The Newton step
//! `y' = ((m-1)·y + a/y^(m-1)) / m` raises `y^(m-1)` by
//! square-and-multiply, so the per-step cost is `O(log m)`, never `O(m)`:
//! a caller-supplied `i32` order must not be able to drive a linear loop
//! (the same denial-of-service discipline as the parse and `pow`
//! integer-exponent caps).
//! Newton is self-correcting, so the accumulated error is dominated by the
//! final step's handful of `NearestEven` operations, comfortably inside
//! [`ROOTN_ERROR_GUARD`].
//!
//! Special cases per IEEE 754-2019 §9.2.1:
//!
//! - `rootn(x, 0)` = qNaN + `INVALID` for every `x` (including NaN/∞).
//! - NaN propagates; sNaN raises `INVALID`.
//! - `rootn(±0, n>0)` = `+0` for even `n`, sign-preserving for odd `n`.
//! - `rootn(±0, n<0)` is a pole: `+∞` for even `n`, sign-preserving `∞`
//!   for odd `n`, with `DIV_BY_ZERO`.
//! - `rootn(+∞, n>0)` = `+∞`; `rootn(+∞, n<0)` = `+0`.
//! - `rootn(−∞, n)` for odd `n`: `−∞` (`n>0`) / `−0` (`n<0`); for even
//!   `n`: qNaN + `INVALID`.
//! - `rootn(negative finite, n)` for odd `n` is the negative real root;
//!   for even `n` it is qNaN + `INVALID`.

use crate::big::{BigFloat, BuildError};
use crate::class::Class;
use crate::rounding::RoundingMode;
use crate::sign::Sign;
use crate::status::{auto_raise, Status};

#[cfg(feature = "fixed")]
use crate::fixed::FixedFloat;
#[cfg(feature = "fixed")]
use crate::mantissa::limbs_for;

use crate::math::ziv::ziv_round;
use crate::math::ziv_calibration::ROOTN_ERROR_GUARD;

/// Newton-iteration cap. Convergence from the power-of-two seed is
/// quadratic, reaching the `≤ target + 1024` working precision in roughly
/// `log2(working)` steps (~12 at the cap); the bound is a generous
/// backstop, and the Ziv guard absorbs any residual sub-ulp wobble.
const ROOTN_NEWTON_MAX_ITERS: u32 = 100;

impl BigFloat {
    /// IEEE 754-2019 `rootn(self, n)`: the real `n`-th root.
    ///
    /// Returns `self^(1/n)` rounded under `mode` to `self.precision`.
    #[must_use]
    pub fn rootn(&self, n: i32, mode: RoundingMode) -> (Self, Status) {
        self.rootn_round(n, self.precision, mode)
            .expect("self.precision >= 1 by invariant")
    }

    /// `rootn(self, n)` with an explicit result precision.
    pub fn rootn_round(
        &self,
        n: i32,
        target_precision: u32,
        mode: RoundingMode,
    ) -> Result<(Self, Status), BuildError> {
        if target_precision == 0 {
            return Err(BuildError::PrecisionZero);
        }
        Ok(rootn_kernel(self, n, target_precision, mode))
    }

    /// `rootn` accumulating into a caller-supplied flag bag.
    #[must_use]
    pub fn rootn_with_flags(&self, n: i32, mode: RoundingMode, flags: &mut Status) -> Self {
        let (value, status) = self.rootn(n, mode);
        *flags |= status;
        value
    }

    /// `rootn_round` accumulating into a caller-supplied flag bag.
    pub fn rootn_round_with_flags(
        &self,
        n: i32,
        target_precision: u32,
        mode: RoundingMode,
        flags: &mut Status,
    ) -> Result<Self, BuildError> {
        let (value, status) = self.rootn_round(n, target_precision, mode)?;
        *flags |= status;
        Ok(value)
    }
}

#[cfg(feature = "fixed")]
impl<const PREC: u32> FixedFloat<PREC>
where
    [(); limbs_for(PREC)]:,
{
    /// `rootn(self, n)` for `FixedFloat`. Delegates to
    /// [`BigFloat::rootn`].
    #[must_use]
    pub fn rootn(&self, n: i32, mode: RoundingMode) -> (Self, Status) {
        let (big, status) = self.to_big().rootn(n, mode);
        (
            Self::try_from_big_exact(big).expect("precision matches"),
            status,
        )
    }
}

fn qnan(target: u32) -> BigFloat {
    BigFloat::try_new_quiet_nan(Sign::Positive, target, &[]).expect("precision >= 1")
}

fn rootn_kernel(x: &BigFloat, n: i32, target: u32, mode: RoundingMode) -> (BigFloat, Status) {
    // n == 0 is a domain error for every operand, including NaN / ∞.
    if n == 0 {
        auto_raise(Status::INVALID);
        return (qnan(target), Status::INVALID);
    }

    let odd = n & 1 != 0;

    match &x.class {
        Class::Nan {
            quiet,
            sign,
            payload,
        } => {
            if *quiet {
                (
                    BigFloat::try_new_quiet_nan(*sign, target, payload).expect("precision >= 1"),
                    Status::OK,
                )
            } else {
                auto_raise(Status::INVALID);
                (qnan(target), Status::INVALID)
            }
        }
        Class::Zero { sign } => {
            // Odd n keeps the zero's sign; even n collapses to +0 / +∞.
            let result_sign = if odd { *sign } else { Sign::Positive };
            if n > 0 {
                (
                    BigFloat::try_new_zero(result_sign, target).expect("precision >= 1"),
                    Status::OK,
                )
            } else {
                // Pole: 0^(negative) = ∞ from a finite operand → DIV_BY_ZERO.
                auto_raise(Status::DIV_BY_ZERO);
                (
                    BigFloat::try_new_infinity(result_sign, target).expect("precision >= 1"),
                    Status::DIV_BY_ZERO,
                )
            }
        }
        Class::Infinity { sign } => {
            match sign {
                Sign::Positive => {
                    let inf_or_zero = if n > 0 {
                        BigFloat::try_new_infinity(Sign::Positive, target)
                    } else {
                        BigFloat::try_new_zero(Sign::Positive, target)
                    };
                    (inf_or_zero.expect("precision >= 1"), Status::OK)
                }
                Sign::Negative => {
                    if odd {
                        // n>0 odd: −∞; n<0 odd: −0.
                        let v = if n > 0 {
                            BigFloat::try_new_infinity(Sign::Negative, target)
                        } else {
                            BigFloat::try_new_zero(Sign::Negative, target)
                        };
                        (v.expect("precision >= 1"), Status::OK)
                    } else {
                        // Even root of −∞ is not real.
                        auto_raise(Status::INVALID);
                        (qnan(target), Status::INVALID)
                    }
                }
            }
        }
        Class::Normal { sign, .. } => {
            let x_neg = matches!(sign, Sign::Negative);
            if x_neg && !odd {
                // Even root of a negative finite is not real.
                auto_raise(Status::INVALID);
                return (qnan(target), Status::INVALID);
            }

            let m = n.unsigned_abs();
            let neg_order = n < 0;
            ziv_round(
                |w| {
                    let a = x
                        .abs()
                        .round_to_precision(w, RoundingMode::NearestEven)
                        .expect("precision >= 1")
                        .0;
                    let root = mth_root_at(&a, m, w);
                    let mag = if neg_order {
                        let one = BigFloat::try_from_i64_exact(1, w).expect("precision >= 1");
                        one.div(&root, RoundingMode::NearestEven).0
                    } else {
                        root
                    };
                    if x_neg {
                        mag.negated()
                    } else {
                        mag
                    }
                },
                target,
                mode,
                ROOTN_ERROR_GUARD,
            )
        }
    }
}

/// `a^(1/m)` for a positive finite `a`, `m >= 1`, evaluated at working
/// precision `w` under `NearestEven`. `m == 1` returns `a`; otherwise
/// Newton's iteration `y' = ((m-1)·y + a/y^(m-1)) / m` from the
/// power-of-two seed `2^(e_a / m)`.
fn mth_root_at(a: &BigFloat, m: u32, w: u32) -> BigFloat {
    if m == 1 {
        return a.clone();
    }
    let e_a = match &a.class {
        Class::Normal { exponent, .. } => *exponent,
        _ => 0,
    };
    let mut y = pow2_at(e_a / i64::from(m), w);
    let m_big = BigFloat::try_from_i64_exact(i64::from(m), w).expect("precision >= 1");
    let m1_big = BigFloat::try_from_i64_exact(i64::from(m - 1), w).expect("precision >= 1");

    let mut prev: Option<BigFloat> = None;
    for _ in 0..ROOTN_NEWTON_MAX_ITERS {
        let y_pow = int_pow_at(&y, m - 1, w);
        let (a_over, _) = a.div(&y_pow, RoundingMode::NearestEven);
        let (m1y, _) = m1_big.mul(&y, RoundingMode::NearestEven);
        let (num, _) = m1y.add(&a_over, RoundingMode::NearestEven);
        let (next, _) = num.div(&m_big, RoundingMode::NearestEven);

        // Fixed point, or a 1-ulp 2-cycle around the root — either member
        // is within 1 ulp, far inside the Ziv guard.
        if matches!(next.partial_cmp(&y).0, Some(core::cmp::Ordering::Equal)) {
            return next;
        }
        if let Some(p) = &prev {
            if matches!(next.partial_cmp(p).0, Some(core::cmp::Ordering::Equal)) {
                return next;
            }
        }
        prev = Some(y);
        y = next;
    }
    y
}

/// `base^k` (`k >= 1`) at working precision `w` under `NearestEven`, via
/// square-and-multiply (`O(log k)` multiplies).
fn int_pow_at(base: &BigFloat, k: u32, w: u32) -> BigFloat {
    let mut result = BigFloat::try_from_i64_exact(1, w).expect("precision >= 1");
    let mut b = base.clone();
    let mut e = k;
    while e > 0 {
        if e & 1 == 1 {
            result = result.mul(&b, RoundingMode::NearestEven).0;
        }
        e >>= 1;
        if e > 0 {
            b = b.mul(&b, RoundingMode::NearestEven).0;
        }
    }
    result
}

/// `2^k` at precision `w`: `1.0` with its exponent shifted by `k`.
fn pow2_at(k: i64, w: u32) -> BigFloat {
    let mut v = BigFloat::try_from_i64_exact(1, w).expect("precision >= 1");
    if let Class::Normal { exponent, .. } = &mut v.class {
        *exponent = exponent.saturating_add(k);
    }
    v
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

    fn close(a: &BigFloat, b: &BigFloat, bits: u32) -> bool {
        let (diff, _) = a.sub(b, RoundingMode::NearestEven);
        let d = diff.abs();
        if d.is_zero() {
            return true;
        }
        let p = a.precision().max(b.precision());
        let two = from_i64(2, p);
        let mut bound = b.abs();
        if bound.is_zero() {
            bound = from_i64(1, p);
        }
        for _ in 0..bits {
            bound = bound.div(&two, RoundingMode::NearestEven).0;
        }
        matches!(
            d.partial_cmp(&bound).0,
            Some(Ordering::Less | Ordering::Equal)
        )
    }

    #[test]
    fn rootn_zero_order_is_invalid() {
        for x in [
            from_i64(8, 53),
            BigFloat::try_new_zero(Sign::Positive, 53).unwrap(),
            BigFloat::try_new_infinity(Sign::Positive, 53).unwrap(),
            BigFloat::try_new_quiet_nan(Sign::Positive, 53, &[]).unwrap(),
        ] {
            let (r, s) = x.rootn(0, RoundingMode::NearestEven);
            assert!(r.is_quiet_nan() && s.invalid(), "rootn(_, 0)");
        }
    }

    #[test]
    fn rootn_perfect_powers() {
        assert!(eq(
            &from_i64(8, 53).rootn(3, RoundingMode::NearestEven).0,
            &from_i64(2, 53)
        ));
        assert!(eq(
            &from_i64(27, 53).rootn(3, RoundingMode::NearestEven).0,
            &from_i64(3, 53)
        ));
        assert!(eq(
            &from_i64(16, 53).rootn(4, RoundingMode::NearestEven).0,
            &from_i64(2, 53)
        ));
        assert!(eq(
            &from_i64(32, 53).rootn(5, RoundingMode::NearestEven).0,
            &from_i64(2, 53)
        ));
        assert!(eq(
            &from_i64(1000, 53).rootn(3, RoundingMode::NearestEven).0,
            &from_i64(10, 53)
        ));
    }

    #[test]
    fn rootn_negative_base_odd_order() {
        // Odd root of a negative is the negative real root.
        let (r, s) = from_i64(-8, 53).rootn(3, RoundingMode::NearestEven);
        assert!(!s.invalid());
        assert!(eq(&r, &from_i64(-2, 53)));
        let (r2, _) = from_i64(-32, 53).rootn(5, RoundingMode::NearestEven);
        assert!(eq(&r2, &from_i64(-2, 53)));
    }

    #[test]
    fn rootn_negative_base_even_order_is_invalid() {
        let (r, s) = from_i64(-16, 53).rootn(4, RoundingMode::NearestEven);
        assert!(r.is_quiet_nan() && s.invalid());
        let (r2, s2) = from_i64(-4, 53).rootn(2, RoundingMode::NearestEven);
        assert!(r2.is_quiet_nan() && s2.invalid());
    }

    #[test]
    fn rootn_negative_order_reciprocates() {
        // rootn(8, -3) = 1/2.
        let one = from_i64(1, 60);
        let two = from_i64(2, 60);
        let (half, _) = one.div(&two, RoundingMode::NearestEven);
        let (r, _) = from_i64(8, 60).rootn(-3, RoundingMode::NearestEven);
        assert!(close(&r, &half, 50), "rootn(8,-3) = {r}");
    }

    #[test]
    fn rootn_order_one_is_identity() {
        let x = BigFloat::parse_str("3.5", 53, RoundingMode::NearestEven)
            .unwrap()
            .0;
        assert!(eq(&x.rootn(1, RoundingMode::NearestEven).0, &x));
        // rootn(x, -1) = 1/x.
        let one = from_i64(1, 53);
        let (recip, _) = one.div(&x, RoundingMode::NearestEven);
        assert!(close(&x.rootn(-1, RoundingMode::NearestEven).0, &recip, 50));
    }

    #[test]
    fn rootn_order_two_matches_sqrt() {
        for v in [2i64, 3, 5, 10, 99] {
            let (a, _) = from_i64(v, 113).rootn(2, RoundingMode::NearestEven);
            let (b, _) = from_i64(v, 113).sqrt(RoundingMode::NearestEven);
            assert!(eq(&a, &b), "rootn({v},2) != sqrt({v})");
        }
    }

    #[test]
    fn rootn_round_trip_pow() {
        // rootn(x, m)^m ≈ x for a non-perfect root.
        let x = from_i64(7, 200);
        let (root, _) = x.rootn(3, RoundingMode::NearestEven);
        let back = int_pow_at(&root, 3, 200);
        assert!(close(&back, &x, 180), "rootn(7,3)^3 = {back}");
    }

    #[test]
    fn rootn_zero_signed() {
        let pz = BigFloat::try_new_zero(Sign::Positive, 53).unwrap();
        let nz = BigFloat::try_new_zero(Sign::Negative, 53).unwrap();
        // n>0: +0 / sign-preserving.
        assert!(pz.rootn(2, RoundingMode::NearestEven).0.is_sign_positive());
        assert!(nz.rootn(3, RoundingMode::NearestEven).0.is_sign_negative()); // odd keeps sign
        assert!(nz.rootn(2, RoundingMode::NearestEven).0.is_sign_positive()); // even -> +0
        for z in [&pz, &nz] {
            let (r, _) = z.rootn(2, RoundingMode::NearestEven);
            assert!(r.is_zero());
        }
    }

    #[test]
    fn rootn_zero_negative_order_is_pole() {
        let pz = BigFloat::try_new_zero(Sign::Positive, 53).unwrap();
        let nz = BigFloat::try_new_zero(Sign::Negative, 53).unwrap();
        let (r, s) = pz.rootn(-2, RoundingMode::NearestEven);
        assert!(r.is_infinite() && r.is_sign_positive() && s.div_by_zero());
        let (r2, s2) = nz.rootn(-3, RoundingMode::NearestEven);
        assert!(r2.is_infinite() && r2.is_sign_negative() && s2.div_by_zero());
        let (r3, s3) = nz.rootn(-2, RoundingMode::NearestEven);
        assert!(r3.is_infinite() && r3.is_sign_positive() && s3.div_by_zero()); // even -> +inf
    }

    #[test]
    fn rootn_infinity() {
        let pinf = BigFloat::try_new_infinity(Sign::Positive, 53).unwrap();
        let ninf = BigFloat::try_new_infinity(Sign::Negative, 53).unwrap();
        assert!(pinf.rootn(2, RoundingMode::NearestEven).0.is_infinite());
        let (r, _) = pinf.rootn(-2, RoundingMode::NearestEven);
        assert!(r.is_zero() && r.is_sign_positive()); // +inf, n<0 -> +0
        let (r2, _) = ninf.rootn(3, RoundingMode::NearestEven);
        assert!(r2.is_infinite() && r2.is_sign_negative()); // -inf odd -> -inf
        let (r3, _) = ninf.rootn(-3, RoundingMode::NearestEven);
        assert!(r3.is_zero() && r3.is_sign_negative()); // -inf, n<0 odd -> -0
        let (r4, s4) = ninf.rootn(2, RoundingMode::NearestEven);
        assert!(r4.is_quiet_nan() && s4.invalid()); // even root of -inf
    }

    #[test]
    fn rootn_nan_propagation() {
        let q = BigFloat::try_new_quiet_nan(Sign::Positive, 53, &[]).unwrap();
        assert!(q.rootn(3, RoundingMode::NearestEven).0.is_quiet_nan());
        let sn = BigFloat::try_new_signaling_nan(Sign::Positive, 53, &[]).unwrap();
        let (r, s) = sn.rootn(3, RoundingMode::NearestEven);
        assert!(r.is_quiet_nan() && s.invalid());
    }

    #[test]
    fn rootn_round_rejects_zero_precision() {
        assert_eq!(
            from_i64(8, 53).rootn_round(3, 0, RoundingMode::NearestEven),
            Err(BuildError::PrecisionZero)
        );
    }

    #[test]
    fn rootn_large_order_terminates() {
        // A large order must not drive a linear loop; the square-and-multiply
        // power keeps this fast and the result near 1.
        let (r, _) = from_i64(2, 64).rootn(1_000_000, RoundingMode::NearestEven);
        // 2^(1/1e6) is just above 1.
        let one = from_i64(1, 64);
        assert_eq!(r.partial_cmp(&one).0, Some(Ordering::Greater));
        let two = from_i64(2, 64);
        assert_eq!(r.partial_cmp(&two).0, Some(Ordering::Less));
    }
}
