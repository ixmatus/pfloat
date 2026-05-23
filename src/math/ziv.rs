//! Ziv correct-rounding driver.
//!
//! The interval test, factored out of [`crate::math::pow`] at slice
//! p1.2 so the rest of the elementary-transcendental surface can move
//! off the slice-3a fixed-64-bit-guard convention without each kernel
//! reimplementing the loop. The pow kernel was the first caller
//! (ADR-0022, slice 7c); subsequent slices wire `exp`, `ln`, `tanh`,
//! and `lgamma` through [`ziv_round`] and inherit the same correctness
//! argument.
//!
//! The driver is the realization of DESIGN.md §"Ziv's strategy": the
//! caller supplies a closure that evaluates the function at variable
//! working precision; the driver rounds the value to the target under
//! the caller's mode, then asks whether both ends of the bounded
//! evaluation-error interval round to the same value. If they do, the
//! true value rounds there too and the result is correctly rounded;
//! otherwise a rounding boundary lies inside the uncertainty, the
//! guard doubles, and the loop retries, bounded by
//! [`ZIV_MAX_ITERS`]. The cap is the honest measure-zero caveat MPFR
//! also documents.
//!
//! The interval test is the sound termination criterion. Comparing
//! two adjacent guards' rounded values would false-converge on a
//! hard-to-round input (both insufficient guards agree on the same
//! wrong value); the interval test does not.

use core::cmp::Ordering;

use crate::big::BigFloat;
use crate::class::Class;
use crate::rounding::RoundingMode;
use crate::sign::Sign;
use crate::status::{auto_raise, Status};

/// First Ziv guard: the initial evaluation uses
/// `target + ZIV_BASE_GUARD` extra bits.
const ZIV_BASE_GUARD: u32 = 64;

/// Maximum extra guard bits above the target precision. The doubling
/// schedule (64, 128, 256, 512, 1024) reaches this at the last
/// iteration.
const ZIV_GUARD_CAP: u32 = 1024;

/// Maximum guard-doubling iterations. On the measure-zero exact-tie
/// inputs that exhaust this many iterations the result may be 1 ULP
/// off in directed modes — the honest caveat MPFR also documents
/// (DESIGN.md §"Ziv's strategy", lines 287-299).
pub(crate) const ZIV_MAX_ITERS: u32 = 5;

/// Slack, in bits below the working precision, charged to `eval`'s
/// accumulated `NearestEven` rounding error. `pow_int`'s
/// square-and-multiply uses ≤ 64 multiplies (≤ 2⁶ ULP) and the
/// `exp·ln` path a handful of operations, all far under 2²⁴ ULP, so
/// the half-width `|y|·2^-(w-24)` is a sound upper bound on
/// `|eval(w) − f(x)|` for the domain the elementary kernels serve
/// (exp/ln series at ~4·w iterations, lgamma Stirling+product at ~30
/// ops, all comfortably under 2²⁴ ULP at the working precisions this
/// driver runs).
const ZIV_ERROR_GUARD: u32 = 24;

/// `|y| · 2^-shift`, formed by decrementing the binary exponent (the
/// exact power-of-two scaling used elsewhere in the crate, e.g.
/// `math::mod::pi_over_2_at`). Non-normal `y` has no boundary
/// uncertainty, so a zero half-width is returned.
fn half_width(y: &BigFloat, shift: i64) -> BigFloat {
    match &y.class {
        Class::Normal {
            exponent, mantissa, ..
        } => BigFloat {
            class: Class::Normal {
                sign: Sign::Positive,
                exponent: exponent - shift,
                mantissa: mantissa.clone(),
            },
            precision: y.precision,
        },
        _ => BigFloat::try_new_zero(Sign::Positive, y.precision).expect("precision >= 1"),
    }
}

/// Correctly round `eval`'s value to `target` precision under `mode`
/// by the Ziv interval test (DESIGN.md §"Ziv's strategy").
///
/// `eval(working)` returns the kernel's value computed at the working
/// precision with `NearestEven` internal rounding; its error against
/// the true value is bounded by the half-width `|y|·2^-(working −
/// ZIV_ERROR_GUARD)`. If both ends of that uncertainty interval round
/// to the same `target`-precision value under `mode`, every point in
/// the interval — including the true value — rounds there too, so
/// that value is correctly rounded. Otherwise a rounding boundary
/// lies within the uncertainty: the guard doubles (capped at
/// [`ZIV_GUARD_CAP`]) and the loop retries, bounded by
/// [`ZIV_MAX_ITERS`]. Comparing two *adjacent* guards would falsely
/// converge on a hard-to-round input (both insufficient guards agree
/// on the wrong value); the interval test does not.
pub(crate) fn ziv_round(
    eval: impl Fn(u32) -> BigFloat,
    target: u32,
    mode: RoundingMode,
) -> (BigFloat, Status) {
    let mut guard = ZIV_BASE_GUARD;
    let mut fallback: Option<(BigFloat, Status)> = None;
    for _ in 0..ZIV_MAX_ITERS {
        let working = target.saturating_add(guard);
        let y = eval(working);
        let (cand, status) = y
            .round_to_precision(target, mode)
            .expect("target precision >= 1");

        let shift = i64::from(working) - i64::from(ZIV_ERROR_GUARD);
        let d = half_width(&y, shift);
        let lo = y.sub(&d, RoundingMode::NearestEven).0;
        let hi = y.add(&d, RoundingMode::NearestEven).0;
        let lo_r = lo.round_to_precision(target, mode).expect("target >= 1").0;
        let hi_r = hi.round_to_precision(target, mode).expect("target >= 1").0;
        if matches!(lo_r.partial_cmp(&hi_r).0, Some(Ordering::Equal)) {
            // The whole uncertainty interval rounds to one value:
            // correct rounding is settled.
            auto_raise(status);
            return (cand, status);
        }

        fallback = Some((cand, status));
        guard = guard.saturating_mul(2).min(ZIV_GUARD_CAP);
    }
    // Cap reached on a pathologically hard input: best effort.
    let (cand, status) = fallback.expect("ZIV_MAX_ITERS >= 1");
    auto_raise(status);
    (cand, status)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn ziv_round_exact_value_is_bit_exact() {
        // An exact value (8) is identical at every working precision,
        // so the first two guards agree immediately and the result is
        // returned bit-exactly — the MPFR integer-fast-path parity
        // the pow integer branch relies on.
        for &mode in &[
            RoundingMode::NearestEven,
            RoundingMode::TowardPositive,
            RoundingMode::TowardNegative,
            RoundingMode::TowardZero,
            RoundingMode::NearestAway,
        ] {
            let (r, _) = ziv_round(
                |w| BigFloat::try_from_i64_exact(8, w).expect("w >= 1"),
                53,
                mode,
            );
            let eight = BigFloat::try_from_i64_exact(8, 53).unwrap();
            assert_eq!(
                r.partial_cmp(&eight).0,
                Some(Ordering::Equal),
                "ziv_round exact 8 under {mode:?} = {r}"
            );
            assert_eq!(r.precision, 53);
        }
    }

    #[test]
    fn ziv_round_converges_to_correctly_rounded() {
        // A transcendental-shaped constant: as the guard grows the
        // target-rounding stabilizes, and ziv_round must land on the
        // value obtained by rounding the constant directly to the
        // target precision under the same mode (correct rounding).
        let digits = "1.41421356237309504880168872420969807856967187537694\
                       8073176679737990732478462107038850387534327641572735";
        for &mode in &[
            RoundingMode::NearestEven,
            RoundingMode::TowardPositive,
            RoundingMode::TowardNegative,
            RoundingMode::TowardZero,
            RoundingMode::NearestAway,
        ] {
            let target = 113u32;
            let (r, _) = ziv_round(
                |w| {
                    BigFloat::parse_str(digits, w, RoundingMode::NearestEven)
                        .expect("parse")
                        .0
                },
                target,
                mode,
            );
            let direct = BigFloat::parse_str(digits, target, mode).expect("parse").0;
            assert_eq!(
                r.partial_cmp(&direct).0,
                Some(Ordering::Equal),
                "ziv_round vs direct round under {mode:?}: ziv={r}, direct={direct}"
            );
        }
    }

    #[test]
    fn ziv_round_is_idempotent_under_recompute() {
        // Stable input → stable output across an independent re-run:
        // the driver does not introduce nondeterminism.
        let eval = |w: u32| {
            BigFloat::parse_str(
                "2.7182818284590452353602874713527",
                w,
                RoundingMode::NearestEven,
            )
            .expect("parse")
            .0
        };
        let (a, _) = ziv_round(eval, 80, RoundingMode::NearestEven);
        let (b, _) = ziv_round(eval, 80, RoundingMode::NearestEven);
        assert_eq!(a.partial_cmp(&b).0, Some(Ordering::Equal));
    }
}
