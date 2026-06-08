//! Shared directed-pair enclosure machinery for the complex elementary
//! kernels (`csqrt`, `cexp`, `clog`).
//!
//! Each kernel computes a component's true value bracketed by a directed pair
//! `[lo, hi]` (`lo ≤ true ≤ hi`) at a growing working precision, then rounds
//! both ends to the output precision under the caller's mode and accepts when
//! they agree. This is the same correctness argument the C3 divide uses
//! (`div.rs`): when both ends of the enclosure round to one value, the true
//! value rounds there too, so that value is correctly rounded. The convergence
//! test requires the ends to agree in value AND sign, which separates `±0`.
//!
//! `INEXACT` is computed from whether the bracket collapsed to a single
//! representable value, never forced, so an exact algebraic output (a Gaussian
//! integer square root, `csqrt(3 + 4i) = 2 + i`) reports `OK`.

use core::cmp::Ordering;

use pfloat::{BigFloat, RoundingMode, Status};

/// Working-precision guards, following pfloat's Ziv schedule (`ziv.rs`) and the
/// C3 divide (`div.rs`): iteration `k` uses `working = p + GUARDS[k]`, capped
/// at five (the MPFR measure-zero caveat for the rare hard-to-round input in a
/// directed mode).
pub(crate) const GUARDS: [u32; 5] = [64, 128, 256, 512, 1024];

/// The outcome of rounding one component's directed bracket to the output
/// precision: the rounded value, whether the bracket converged, and the IEEE
/// status (only `INEXACT`/`OK` here; the kernels add `INVALID`/`DIV_BY_ZERO`).
pub(crate) struct Resolved {
    pub value: BigFloat,
    pub converged: bool,
    pub status: Status,
}

/// Round a directed bracket `[lo, hi]` (with `lo ≤ true ≤ hi`) to precision `p`
/// under `mode` and decide convergence. Convergence needs the two rounded ends
/// to agree in value AND sign (the sign clause separates `+0` from `−0`, which
/// IEEE comparison treats as equal). A bracket that straddles a rounding
/// boundary (or `±0`) does not converge, so the caller grows the guard.
pub(crate) fn resolve_bracket(
    lo: &BigFloat,
    hi: &BigFloat,
    p: u32,
    mode: RoundingMode,
) -> Resolved {
    let (lo_r, _) = lo.round_to_precision(p, mode).expect("p >= 1");
    let (hi_r, _) = hi.round_to_precision(p, mode).expect("p >= 1");

    // A NaN end means the bracket itself degenerated (a NaN propagated through
    // the enclosure); converged only when both ends are NaN. The caller's
    // special-value dispatch owns the genuine NaN rows, so this is a backstop.
    if lo_r.is_nan() || hi_r.is_nan() {
        return Resolved {
            converged: lo_r.is_nan() && hi_r.is_nan(),
            status: Status::OK,
            value: lo_r,
        };
    }

    let value_equal = matches!(lo_r.partial_cmp(&hi_r).0, Some(Ordering::Equal));
    let sign_equal = lo_r.is_sign_negative() == hi_r.is_sign_negative();
    if value_equal && sign_equal {
        let status = if bracket_is_exact(lo, hi, &lo_r) {
            Status::OK
        } else {
            Status::INEXACT
        };
        Resolved {
            value: lo_r,
            converged: true,
            status,
        }
    } else {
        Resolved {
            value: lo_r,
            converged: false,
            status: Status::INEXACT,
        }
    }
}

/// The component is exact iff the directed bracket collapsed to one signed
/// value and rounding to `p` did not move it.
pub(crate) fn bracket_is_exact(lo: &BigFloat, hi: &BigFloat, rounded: &BigFloat) -> bool {
    let collapsed = matches!(lo.partial_cmp(hi).0, Some(Ordering::Equal))
        && lo.is_sign_negative() == hi.is_sign_negative();
    let representable = matches!(rounded.partial_cmp(lo).0, Some(Ordering::Equal))
        && rounded.is_sign_negative() == lo.is_sign_negative();
    collapsed && representable
}
