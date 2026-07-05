//! `exp(x)`: natural exponential function.
//!
//! Algorithm:
//!
//! 1. Special cases: NaN propagates; `exp(±0) = 1`;
//!    `exp(+∞) = +∞`; `exp(-∞) = +0`. Past the exponent rim
//!    (ADR-0096): certain overflow gives `+∞` (`MaxFinite` under the
//!    inward modes) with `OVERFLOW|INEXACT`; results below
//!    representability give `+0` (`MinPos` under `TowardPositive`, and
//!    in the `[MinPos/2, MinPos)` sliver also under the nearest
//!    modes) with `UNDERFLOW|INEXACT`. Inputs whose result exponent
//!    lands inside the rim compute exactly like any other.
//! 2. Range reduce: choose integer `k = round(x / ln(2))`, then
//!    `r = x − k · ln(2)`. After reduction `|r| ≤ ln(2)/2 ≈ 0.347`.
//!    `ln(2)` is hardcoded at 1024-bit precision (see
//!    `super::LN2_LIMBS_1024`).
//! 3. Taylor series: `exp(r) = 1 + r + r²/2! + r³/3! + …`. With
//!    `|r| < 0.5`, the series converges geometrically faster than
//!    one bit per term once `n > |r| · target_precision`.
//! 4. Compose: `exp(x) = exp(r) · 2^k`. The `2^k` factor is a free
//!    exponent shift on the `BigFloat` (no arithmetic needed).
//!
//! Correctly rounded under every IEEE rounding mode via the shared
//! [`crate::math::ziv::ziv_round`] driver (slice p1.2, ADR-0022): the
//! kernel function [`exp_at_w`] evaluates the algorithm above at the
//! working precision the driver supplies, and the driver grows the
//! guard until the Ziv interval test certifies that the rounded
//! target-precision value lies in the bracket of all working-
//! precision evaluations. CORE-MATH's worst-case-rounding corpus
//! (sourced from the Lefèvre–Muller ARITH-15 2001 table) including
//! the leading underflow block exercises the boundary cases in
//! `tests/differential_lefevre_muller.rs`. Slice p1.2 closed slice
//! 8b's documented exp underflow defect: the corpus's underflow
//! stress block now lands as the regression guard.

use crate::big::BigFloat;
use crate::class::Class;
use crate::ops::limbs::extract_as_integer;
use crate::rounding::RoundingMode;
use crate::sign::Sign;
use crate::status::{auto_raise, Status};

#[cfg(feature = "fixed")]
use crate::fixed::FixedFloat;
#[cfg(feature = "fixed")]
use crate::mantissa::limbs_for;

use super::ln_2_at;
use super::ziv::ziv_round;
use super::ziv_calibration::EXP_ERROR_GUARD;

impl BigFloat {
    /// `exp(self)`: returns `e^x` rounded under `mode` to
    /// `self.precision`.
    ///
    /// Correctly rounded under every IEEE rounding mode via the
    /// shared [`crate::math::ziv::ziv_round`] driver (slice p1.2,
    /// ADR-0022). See [the module docs](self) for the algorithm.
    #[must_use]
    pub fn exp(&self, mode: RoundingMode) -> (Self, Status) {
        let target = self.precision;
        self.exp_round(target, mode)
    }

    /// `exp(self)` with an explicit result precision.
    #[must_use]
    pub fn exp_round(&self, target_precision: u32, mode: RoundingMode) -> (Self, Status) {
        debug_assert!(target_precision >= 1);
        exp_kernel(self, target_precision, mode)
    }
}

#[cfg(feature = "fixed")]
impl<const PREC: u32> FixedFloat<PREC>
where
    [(); limbs_for(PREC)]:,
{
    /// `exp(self)` for `FixedFloat`. Delegates to
    /// [`BigFloat::exp`].
    #[must_use]
    pub fn exp(&self, mode: RoundingMode) -> (Self, Status) {
        let (big, status) = self.to_big().exp(mode);
        (
            Self::try_from_big_exact(big).expect("precision matches"),
            status,
        )
    }
}

fn exp_kernel(x: &BigFloat, target_precision: u32, mode: RoundingMode) -> (BigFloat, Status) {
    // Special cases per IEEE 754-2019 §9.
    match &x.class {
        Class::Nan {
            quiet,
            sign,
            payload,
        } => {
            if !*quiet {
                let nan = BigFloat::try_new_quiet_nan(Sign::Positive, target_precision, &[])
                    .expect("precision >= 1");
                auto_raise(Status::INVALID);
                return (nan, Status::INVALID);
            }
            let nan = BigFloat::try_new_quiet_nan(*sign, target_precision, payload)
                .expect("precision >= 1");
            return (nan, Status::OK);
        }
        Class::Zero { .. } => {
            // exp(±0) = 1.
            let one = BigFloat::try_from_i64_exact(1, target_precision).expect("precision >= 1");
            return (one, Status::OK);
        }
        Class::Infinity {
            sign: Sign::Positive,
        } => {
            return (
                BigFloat::try_new_infinity(Sign::Positive, target_precision)
                    .expect("precision >= 1"),
                Status::OK,
            );
        }
        Class::Infinity {
            sign: Sign::Negative,
        } => {
            // exp(-∞) = +0.
            return (
                BigFloat::try_new_zero(Sign::Positive, target_precision).expect("precision >= 1"),
                Status::OK,
            );
        }
        Class::Normal { .. } => {}
    }

    // Exponent-ceiling triage (pf-7z66, ADR-0096). For |x| ≥ 2^62,
    // exp(x) sits outside or at the rim of the representable
    // exponent range, and the generic reduction below would wrap or
    // saturate its k = round(x/ln2): the review confirmed garbage
    // Normals (one with Status OK certified) and a spurious +inf on
    // representable results. e_x ≤ 61 keeps the established path
    // bit-for-bit untouched (|x|/ln2 < 2^62.6 there, so k and the
    // result exponent stay far from the i64 rim).
    let e_x = match &x.class {
        Class::Normal { exponent, .. } => *exponent,
        _ => unreachable!("specials handled above"),
    };
    // Near-1 tiny x: exp(x) = 1 + x + … perturbs 1 by the LINEAR term x
    // (grows for x > 0, shrinks for x < 0). For e_x ≤ -(target+2) the
    // perturbation is below the target boundary near 1, so exp(x) rounds
    // to 1 (nearest) or the neighbour (directed); past the Ziv guard cap
    // the reduced series collapses to exactly 1 and the directed modes
    // returned 1 where the neighbour is due (pf-767j, ADR-0127; the near-1
    // analogue of the expm1 tiny-x family, ADR-0059). subtracts_magnitude
    // tracks the sign of x (exp(x) < 1 iff x < 0). The `e_x ≥ 62` ceiling
    // below is the opposite extreme and disjoint.
    if e_x <= -(i64::from(target_precision) + 2) {
        let one = BigFloat::try_from_i64_exact(1, target_precision).expect("precision >= 1");
        return crate::rounding::round_with_infinitesimal(
            &one,
            Sign::Positive,
            x.is_sign_negative(),
            target_precision,
            mode,
        );
    }
    if e_x >= 62 {
        return exp_extreme(x, target_precision, mode);
    }

    // Correctly rounded under `mode` via the Ziv interval test
    // (ADR-0022). The driver supplies a working precision, the
    // closure evaluates exp at that precision, and the driver grows
    // the guard until correct rounding is certified.
    let (result, status) = ziv_round(
        |working_prec| exp_at_w(x, working_prec),
        target_precision,
        mode,
        EXP_ERROR_GUARD,
    );
    // exp(x) for finite normal x ≠ 0 is transcendental (Lindemann–
    // Weierstrass: e^α is transcendental for nonzero algebraic α, and
    // a dyadic x is algebraic), hence irrational, hence never exactly
    // representable. The result is therefore INEXACT even when it
    // rounds onto a grid value because the residual fell below the
    // kernel's working precision — e.g. exp(2^-1074) → 1.0, where the
    // true 1 + 2^-1074 ≠ 1 but the working-precision evaluation never
    // observes the rounding (pf-njs5 under-report, ADR-0060). x = 0 is
    // the only exact input and is special-cased above to exp(0) = 1.
    let status = status | Status::INEXACT;
    auto_raise(status);
    (result, status)
}

/// Evaluate `exp(x)` at the supplied working precision via
/// range-reduction + Taylor series + power-of-two scale. Returns the
/// unrounded value; the Ziv driver handles rounding to the caller's
/// target precision and mode.
fn exp_at_w(x: &BigFloat, working_prec: u32) -> BigFloat {
    // The reduction r = x − k·ln2 cancels the leading ~e_x bits of
    // the operands (both have magnitude ~2^e_x while r is O(1)), so
    // computing it AT the working precision leaves r with only
    // working − e_x good bits — while the Ziv driver charges this
    // closure a flat 24-bit guard. For e_x ≳ 25 that mismatch
    // certified 1-ulp NE misroundings at percent-level density
    // (pf-t6ht's probed band, confirmed by run at e_x = 61 during
    // the ADR-0101 verification when cosh/sinh began forwarding
    // arguments here). Carry the reduction at working + e_x + 8
    // bits — the ADR-0097 deterministic pre-boost, costing extra
    // bits only on the multiply/subtract and only for large
    // arguments — then run the Taylor series at the working
    // precision as before.
    let e_x = match &x.class {
        Class::Normal { exponent, .. } => *exponent,
        _ => 0,
    };
    let reduce_prec =
        working_prec.saturating_add(u32::try_from(e_x.max(0)).unwrap_or(0).saturating_add(8));
    let x_w = x
        .round_to_precision(reduce_prec, RoundingMode::NearestEven)
        .expect("precision >= 1")
        .0;

    let ln_2 = ln_2_at(reduce_prec);

    // k = round(x / ln(2)) as i64.
    let (x_over_ln2, _) = x_w.div(&ln_2, RoundingMode::NearestEven);
    let k = round_bigfloat_to_i64(&x_over_ln2);

    // r = x - k * ln(2), carried wide then narrowed for the series.
    let k_big = BigFloat::try_from_i64_exact(k, reduce_prec)
        .or_else(|_| {
            BigFloat::try_from_i64_round(k, reduce_prec, RoundingMode::NearestEven).map(|(v, _)| v)
        })
        .expect("i64 fits in the reduction precision");
    let (k_ln2, _) = k_big.mul(&ln_2, RoundingMode::NearestEven);
    let (r_wide, _) = x_w.sub(&k_ln2, RoundingMode::NearestEven);
    let r = r_wide
        .round_to_precision(working_prec, RoundingMode::NearestEven)
        .expect("precision >= 1")
        .0;

    // exp(x) = exp(r) × 2^k. Apply k as a free exponent shift.
    shift_exponent(exp_taylor(&r, working_prec), k)
}

/// Taylor series `exp(r) = 1 + r + r²/2! + r³/3! + …` at the
/// supplied working precision. `|r| ≲ 1.4` (the reduced argument of
/// either reduction path); the series converges geometrically once
/// `n > |r|`, and the iteration cap covers pathological stalls.
fn exp_taylor(r: &BigFloat, working_prec: u32) -> BigFloat {
    let one = BigFloat::try_from_i64_exact(1, working_prec).expect("precision >= 1");
    let mut sum = one.clone();
    let mut term = one;
    let max_iter = 4u32.saturating_mul(working_prec).max(256);
    for n in 1u32..=max_iter {
        let (new_numer, _) = term.mul(r, RoundingMode::NearestEven);
        let n_big =
            BigFloat::try_from_i64_exact(i64::from(n), working_prec).expect("precision >= 1");
        let (new_term, _) = new_numer.div(&n_big, RoundingMode::NearestEven);
        term = new_term;
        let (new_sum, _) = sum.add(&term, RoundingMode::NearestEven);
        sum = new_sum;
        // Termination: term below 2^(-working_prec) compared to sum (~1).
        if let Class::Normal { exponent, .. } = &term.class {
            if *exponent < -i64::from(working_prec) - 4 {
                break;
            }
        } else {
            break;
        }
    }
    sum
}

/// `exp(x)` for `|x| ≥ 2^62` (pf-7z66, ADR-0096): the exponent-rim
/// band. Classifies `n = floor(x/ln2)` — the result's binary
/// exponent — against the representable i64 range and dispatches:
///
/// - `n ≥ 2^63`: the truth is at least `2^(i64::MAX + 1)`, past
///   `MaxFinite` by more than an ulp: certain overflow, mode-aware
///   (+∞ for nearest/upward, `MaxFinite` inward), `OVERFLOW|INEXACT`.
/// - `n ≤ i64::MIN − 2`: the truth is below half of `MinPos`: +0 for
///   every mode except `TowardPositive` (`MinPos`), `UNDERFLOW|INEXACT`.
/// - `n = i64::MIN − 1` (the sliver): the truth lies in
///   `[MinPos/2, MinPos)`, strictly above the to-nearest midpoint
///   (its mantissa `2^{x/ln2 − n} ∈ (1, 2)` strictly, since `x/ln2`
///   is irrational): `MinPos` for nearest/upward, +0 inward,
///   `UNDERFLOW|INEXACT`.
/// - otherwise the truth is representable: Ziv-certify the unscaled
///   `s = exp(x − k·ln2)` and compose with the exact
///   `scale_by_pow2(k)`. Exact power-of-two scaling commutes with
///   rounding, so rounding `s` at the target under `mode` rounds
///   the result; only a genuine carry across the top binade can
///   saturate, and `scale_by_pow2` flags that honestly.
///
/// For `e_x ≥ 63`, `|x|/ln2 > 2^63.5` clears both rims with margin
/// and the sign decides directly; only `e_x = 62` (where `|x|/ln2`
/// spans `[2^62.5, 2^63.5)`) needs the certified floor.
fn exp_extreme(x: &BigFloat, target_precision: u32, mode: RoundingMode) -> (BigFloat, Status) {
    let positive = matches!(x.sign(), Sign::Positive);
    let e_x = match &x.class {
        Class::Normal { exponent, .. } => *exponent,
        _ => unreachable!("caller dispatches on Normal"),
    };
    if e_x >= 63 {
        return if positive {
            exp_certain_overflow(target_precision, mode)
        } else {
            exp_deep_underflow(target_precision, mode)
        };
    }

    let n = certified_floor_x_over_ln2(x);
    if n >= 1i128 << 63 {
        return exp_certain_overflow(target_precision, mode);
    }
    if n < i128::from(i64::MIN) - 1 {
        return exp_deep_underflow(target_precision, mode);
    }
    if n == i128::from(i64::MIN) - 1 {
        return exp_underflow_sliver(target_precision, mode);
    }

    // Representable window. Clamp k one step inside the rim so the
    // unscaled value s = exp(x − k·ln2) ∈ (0.5, 4) composes without
    // saturating for in-range truths: at k = n the bottom-window
    // s ∈ [1, 2) can round down across the binade and the top-window
    // s can carry to 2.0, either of which would push the exact
    // compose onto the rim. With the clamp, exponent(s) ∈ {−1, 0, 1}
    // and k + exponent(s) ∈ [i64::MIN, i64::MAX] always; only a
    // target-rounding carry past the true top binade reaches
    // scale_by_pow2's saturation contract, which flags it.
    let k = n
        .max(i128::from(i64::MIN) + 1)
        .min(i128::from(i64::MAX) - 1) as i64;
    let (s, ziv_status) = ziv_round(
        |w| exp_reduced_pinned(x, k, w),
        target_precision,
        mode,
        EXP_ERROR_GUARD,
    );
    let (result, scale_status) = s.scale_by_pow2(k);
    // A carry in the certified rounding of s (e.g. s → 4.0 at
    // k = i64::MAX − 1) means the infinitely-precise rounding of
    // exp(x) under `mode` lands at 2^(i64::MAX + 1): a genuine IEEE
    // §7.4 overflow. scale_by_pow2 clamps the exponent AFTER the
    // carry has already replaced the mantissa with 1.0, which would
    // return a non-monotone 2^(i64::MAX) (about half the truth and
    // below the same input's TowardZero answer) — route to the
    // mode-aware overflow result instead. The inward modes cannot
    // carry upward (their certified rounding never exceeds the true
    // s < 4), and the bottom compose cannot saturate at all
    // (k ≥ i64::MIN + 1 and exponent(s) ≥ −1), so no underflow
    // analogue exists. Caught by the slice's adversarial
    // verification (x = 2^63·RD_130(ln2), truth (2 − 2^-65)·2^MAX).
    if scale_status.overflow() {
        return exp_certain_overflow(target_precision, mode);
    }
    // exp of a nonzero finite x is transcendental (the ADR-0060
    // posture of the main path), so INEXACT is unconditional.
    let status = ziv_status | scale_status | Status::INEXACT;
    auto_raise(status);
    (result, status)
}

/// Mode-aware certain-overflow result: IEEE 754-2019 §7.4 shape.
/// The truth exceeds `2^(i64::MAX+1)`, i.e. more than an ulp past
/// `MaxFinite`, so the to-nearest modes give +∞ and the inward modes
/// give `MaxFinite`. This is a deliberate divergence from the ops'
/// saturate-to-finite exponent contract (mul/div/fma): their
/// saturated mantissa still carries the true leading bits, while a
/// deep exp overflow has no representable approximation within any
/// bounded relative error, so +∞ + OVERFLOW is the honest §7.4
/// answer and matches this module's documented promise.
pub(super) fn exp_certain_overflow(
    target_precision: u32,
    mode: RoundingMode,
) -> (BigFloat, Status) {
    let value = match mode {
        RoundingMode::TowardZero | RoundingMode::TowardNegative => {
            BigFloat::try_new_infinity(Sign::Positive, target_precision)
                .expect("precision >= 1")
                .next_down()
                .0
        }
        _ => BigFloat::try_new_infinity(Sign::Positive, target_precision).expect("precision >= 1"),
    };
    let status = Status::OVERFLOW | Status::INEXACT;
    auto_raise(status);
    (value, status)
}

/// Mode-aware deep-underflow result: the truth is strictly below
/// `MinPos/2` (= `2^(i64::MIN − 1)`), so every mode except
/// `TowardPositive` rounds to +0, and `TowardPositive` rounds up to
/// `MinPos`.
pub(super) fn exp_deep_underflow(target_precision: u32, mode: RoundingMode) -> (BigFloat, Status) {
    let value = if matches!(mode, RoundingMode::TowardPositive) {
        BigFloat::try_new_zero(Sign::Positive, target_precision)
            .expect("precision >= 1")
            .next_up()
            .0
    } else {
        BigFloat::try_new_zero(Sign::Positive, target_precision).expect("precision >= 1")
    };
    let status = Status::UNDERFLOW | Status::INEXACT;
    auto_raise(status);
    (value, status)
}

/// Mode-aware sliver result for `floor(x/ln2) = i64::MIN − 1`: the
/// truth lies in `(MinPos/2, MinPos)` strictly, so nearest modes and
/// `TowardPositive` give `MinPos` while the inward modes give +0.
pub(super) fn exp_underflow_sliver(
    target_precision: u32,
    mode: RoundingMode,
) -> (BigFloat, Status) {
    let value = match mode {
        RoundingMode::TowardZero | RoundingMode::TowardNegative => {
            BigFloat::try_new_zero(Sign::Positive, target_precision).expect("precision >= 1")
        }
        _ => {
            BigFloat::try_new_zero(Sign::Positive, target_precision)
                .expect("precision >= 1")
                .next_up()
                .0
        }
    };
    let status = Status::UNDERFLOW | Status::INEXACT;
    auto_raise(status);
    (value, status)
}

/// `floor(x/ln2)` certified by interval arithmetic, for
/// `|x| ∈ [2^62, 2^63)` (so the result magnitude is below `2^64`
/// and fits i128 with room). Brackets `ln2` by the neighbours of
/// its correctly rounded value, brackets `x` by directed rounding,
/// takes the min/max over all directed quotients, and accepts when
/// both ends floor to the same integer; otherwise the precision
/// doubles up to a cap that SCALES WITH x's PRECISION: by the
/// irrationality measure of `ln 2` (`μ(ln 2) ≤ 3.57455`,
/// Marcovecchio 2009; `docs/references/`), a dyadic `x` of
/// precision `px` cannot place `x/ln2` closer to an integer than
/// `~2^(64 − μ·(px + 64))`, so a cap of `4·(px + 64) + 1024`
/// bracket bits always certifies (`4 > μ` with slack — the same
/// derivation as `exp_reduced_pinned`'s retry cap). A fixed
/// `q = 1024` cap was REFUTED by this slice's adversarial
/// verification: an 1100-bit `x` one part in `2^1037` from
/// `i64::MIN·ln2` defeated the bracket, and the fall-through
/// crossed a DISPATCH rim (sliver instead of window), returning a
/// wrong `TowardZero` value and a spurious UNDERFLOW — unlike a
/// one-off `k` in the window path, which is self-correcting. The
/// fall-through return of the lower floor remains as the defensive
/// total fallback, now genuinely unreachable by the measure bound.
fn certified_floor_x_over_ln2(x: &BigFloat) -> i128 {
    use core::cmp::Ordering;
    let q_cap = 4u32
        .saturating_mul(x.precision().saturating_add(64))
        .saturating_add(1024);
    let mut q = 128u32;
    loop {
        let ln_2 = ln_2_at(q);
        let ln2_ends = [ln_2.next_down().0, ln_2.next_up().0];
        let x_ends = [
            x.round_to_precision(q, RoundingMode::TowardNegative)
                .expect("q >= 1")
                .0,
            x.round_to_precision(q, RoundingMode::TowardPositive)
                .expect("q >= 1")
                .0,
        ];
        let mut t_lo: Option<BigFloat> = None;
        let mut t_hi: Option<BigFloat> = None;
        for xb in &x_ends {
            for lb in &ln2_ends {
                for m in [RoundingMode::TowardNegative, RoundingMode::TowardPositive] {
                    let t = xb.div_round(lb, q, m).expect("q >= 1").0;
                    let lower = t_lo
                        .as_ref()
                        .is_none_or(|c| matches!(t.partial_cmp(c).0, Some(Ordering::Less)));
                    if lower {
                        t_lo = Some(t.clone());
                    }
                    let higher = t_hi
                        .as_ref()
                        .is_none_or(|c| matches!(t.partial_cmp(c).0, Some(Ordering::Greater)));
                    if higher {
                        t_hi = Some(t);
                    }
                }
            }
        }
        let f_lo = floor_to_i128(&t_lo.expect("loop ran"));
        let f_hi = floor_to_i128(&t_hi.expect("loop ran"));
        if f_lo == f_hi || q >= q_cap {
            return f_lo;
        }
        q = q.saturating_mul(2).min(q_cap);
    }
}

/// `floor(v)` as i128 for `|v| < 2^100` (caller domain: quotients of
/// magnitude below `2^64`). Reads the integer part straight off the
/// mantissa limbs.
fn floor_to_i128(v: &BigFloat) -> i128 {
    let (sign, e, mantissa) = match &v.class {
        Class::Normal {
            sign,
            exponent,
            mantissa,
        } => (*sign, *exponent, mantissa),
        _ => return 0,
    };
    if e < 0 {
        // |v| < 1: floor is 0 (positive) or −1 (negative non-zero).
        return if matches!(sign, Sign::Positive) {
            0
        } else {
            -1
        };
    }
    debug_assert!(e < 100, "floor_to_i128 caller domain");
    let p = v.precision;
    let m_int = extract_as_integer(mantissa, p);
    // value = m · 2^(e − p + 1); integer part = m >> (p − 1 − e).
    // p (≥ 128 in the caller) exceeds e + 1, so the shift is
    // non-negative.
    let s = p - 1 - (e as u32);
    let int_part = shifted_low_u128(&m_int, s);
    let frac_nonzero = low_bits_nonzero(&m_int, s);
    if matches!(sign, Sign::Positive) {
        int_part as i128
    } else {
        -(int_part as i128) - i128::from(frac_nonzero)
    }
}

/// Bits `s..s+128` of a little-endian limb integer, as u128. The
/// caller guarantees the value `m >> s` fits (its top set bit is
/// below position 100).
fn shifted_low_u128(m: &[u64], s: u32) -> u128 {
    let limb = (s / 64) as usize;
    let bit = s % 64;
    let g = |i: usize| u128::from(m.get(i).copied().unwrap_or(0));
    if bit == 0 {
        g(limb) | (g(limb + 1) << 64)
    } else {
        (g(limb) >> bit) | (g(limb + 1) << (64 - bit)) | (g(limb + 2) << (128 - bit))
    }
}

/// Any set bit strictly below position `s` of a little-endian limb
/// integer.
fn low_bits_nonzero(m: &[u64], s: u32) -> bool {
    let limb = (s / 64) as usize;
    let bit = s % 64;
    if m.iter().take(limb.min(m.len())).any(|&l| l != 0) {
        return true;
    }
    bit > 0 && limb < m.len() && (m[limb] & ((1u64 << bit) - 1)) != 0
}

/// The unscaled `exp(x − k·ln2)` at working precision `w`, with the
/// reduction carried at a precision that absorbs the realized
/// cancellation. `x` agrees with `k·ln2` in its leading ~63 bits by
/// construction; an adversarial `x` (one built as a high-precision
/// rounding of `k·ln2`) can agree far deeper, bounded by the
/// irrationality measure of `ln 2` at roughly `μ·(precision(x) +
/// 64)` bits (`μ(ln 2) ≤ 3.57455`, Marcovecchio 2009;
/// `docs/references/`). Start with a 256-bit allowance (covers
/// every `x` up to ~70 bits of precision outright) and grow on
/// realized collapse up to the measure-derived cap; past the cap —
/// unreachable by the bound — the collapsed reduction reproduces
/// the Ziv-cap measure-zero caveat (a possible final-ulp error in
/// directed modes) rather than a certified-garbage value.
fn exp_reduced_pinned(x: &BigFloat, k: i64, w: u32) -> BigFloat {
    let px = x.precision();
    let cap = w
        .saturating_add(4u32.saturating_mul(px.saturating_add(64)))
        .saturating_add(1024);
    let mut wr = w.saturating_add(256);
    loop {
        let x_w = x
            .round_to_precision(wr, RoundingMode::NearestEven)
            .expect("precision >= 1")
            .0;
        let ln_2 = ln_2_at(wr);
        let k_big = BigFloat::try_from_i64_exact(k, wr).expect("wr >= 64 holds any i64");
        let (k_ln2, _) = k_big.mul(&ln_2, RoundingMode::NearestEven);
        let (r_wide, _) = x_w.sub(&k_ln2, RoundingMode::NearestEven);

        // Resolved when r clears the reduction's absolute noise
        // floor 2^(e(k·ln2) + 1 − wr) by w + 8 bits.
        let noise = match &k_ln2.class {
            Class::Normal { exponent, .. } => exponent.saturating_sub(i64::from(wr)) + 1,
            _ => i64::MIN,
        };
        let resolved = match &r_wide.class {
            Class::Normal { exponent, .. } => *exponent >= noise.saturating_add(i64::from(w) + 8),
            _ => false,
        };
        if resolved || wr >= cap {
            let r = r_wide
                .round_to_precision(w, RoundingMode::NearestEven)
                .expect("precision >= 1")
                .0;
            return exp_taylor(&r, w);
        }
        wr = wr.saturating_mul(2).min(cap);
    }
}

/// Round a `BigFloat` to the nearest `i64` (banker's rounding for
/// ties).
///
/// Saturates to `i64::MAX` / `i64::MIN` for out-of-range inputs.
/// Returns `0` for `NaN` and other non-Normal/Zero classes.
pub(super) fn round_bigfloat_to_i64(v: &BigFloat) -> i64 {
    let (sign, exponent, mantissa) = match &v.class {
        Class::Normal {
            sign,
            exponent,
            mantissa,
        } => (*sign, *exponent, mantissa),
        Class::Zero { .. } => return 0,
        _ => return 0,
    };

    if exponent < -1 {
        // |v| < 0.5, rounds to 0.
        return 0;
    }
    if exponent > 62 {
        return if matches!(sign, Sign::Positive) {
            i64::MAX
        } else {
            i64::MIN
        };
    }

    let p = v.precision;
    let m_int = extract_as_integer(mantissa, p);
    let scale = exponent - i64::from(p) + 1;

    let magnitude: u64 = if scale >= 0 {
        shift_left_to_u64(&m_int, scale as u32)
    } else {
        shift_right_round_to_u64(&m_int, (-scale) as u32)
    };

    // Saturate instead of `as i64`: a magnitude that rounded up to
    // exactly 2^63 wrapped to i64::MIN here (pf-7z66 failure (b),
    // the k-wrap window). The exponent-ceiling triage now keeps
    // such magnitudes out of this path, but the helper's contract
    // is saturating regardless.
    if matches!(sign, Sign::Negative) {
        if magnitude >= 1u64 << 63 {
            i64::MIN
        } else {
            -(magnitude as i64)
        }
    } else if magnitude > i64::MAX as u64 {
        i64::MAX
    } else {
        magnitude as i64
    }
}

/// Shift a multi-limb integer left by `s` bits and return the low
/// 64 bits as u64. Used only via `round_bigfloat_to_i64` for inputs
/// whose result is guaranteed to fit in u64 (exp ≤ 62, scale ≥ 0
/// implies the high limbs are all zero).
fn shift_left_to_u64(m: &[u64], s: u32) -> u64 {
    let m0 = m.first().copied().unwrap_or(0);
    if s >= 64 {
        0
    } else {
        m0.checked_shl(s).unwrap_or(0)
    }
}

/// Shift a multi-limb integer right by `s` bits with
/// round-to-nearest-even on the discarded bits, returning the
/// result as a u64 (saturating).
fn shift_right_round_to_u64(m: &[u64], s: u32) -> u64 {
    if s == 0 {
        return m.first().copied().unwrap_or(0);
    }

    let limb_shift = (s / 64) as usize;
    let bit_shift = s % 64;

    let int_part: u64 = if bit_shift == 0 {
        m.get(limb_shift).copied().unwrap_or(0)
    } else {
        let lo = m.get(limb_shift).copied().unwrap_or(0) >> bit_shift;
        let hi = m
            .get(limb_shift + 1)
            .copied()
            .unwrap_or(0)
            .checked_shl(64 - bit_shift)
            .unwrap_or(0);
        lo | hi
    };

    // Guard bit at position s - 1 of m.
    let guard_pos = s - 1;
    let guard_limb = (guard_pos / 64) as usize;
    let guard_bit = guard_pos % 64;
    let guard = if guard_limb < m.len() {
        (m[guard_limb] >> guard_bit) & 1
    } else {
        0
    };

    // Sticky: any non-zero bit below position s - 1.
    let mut sticky = false;
    for &limb in m.iter().take(guard_limb.min(m.len())) {
        if limb != 0 {
            sticky = true;
            break;
        }
    }
    if !sticky && guard_limb < m.len() && guard_bit > 0 {
        let mask = (1u64 << guard_bit) - 1;
        if m[guard_limb] & mask != 0 {
            sticky = true;
        }
    }

    let round_up = guard == 1 && (sticky || (int_part & 1) == 1);
    if round_up {
        int_part.saturating_add(1)
    } else {
        int_part
    }
}

/// Multiply a finite `BigFloat` by `2^k` by adding `k` to its
/// exponent. Non-normal classes pass through.
fn shift_exponent(mut v: BigFloat, k: i64) -> BigFloat {
    if let Class::Normal {
        exponent,
        mantissa: _,
        sign: _,
    } = &mut v.class
    {
        *exponent = exponent.saturating_add(k);
    }
    v
}

#[cfg(test)]
mod tests {
    use super::*;
    use core::cmp::Ordering;

    fn parse(s: &str, p: u32) -> BigFloat {
        BigFloat::parse_str(s, p, RoundingMode::NearestEven)
            .expect("parse should succeed")
            .0
    }

    fn close_at(v: &BigFloat, expected: &BigFloat, prec: u32) -> bool {
        // |v - expected| <= 4 ULPs at `prec`. ULP at a value of
        // magnitude ~1 at precision `prec` is 2^(-prec).
        let (diff, _) = v.sub(expected, RoundingMode::NearestEven);
        let abs_diff = diff.abs();
        if abs_diff.is_zero() {
            return true;
        }
        if !abs_diff.is_normal() {
            return false;
        }
        let exp = match &abs_diff.class {
            Class::Normal { exponent, .. } => *exponent,
            _ => return false,
        };
        // 4 ULP at expected's magnitude is 2^(expected_exp - prec + 3).
        let expected_exp = match &expected.class {
            Class::Normal { exponent, .. } => *exponent,
            _ => 0,
        };
        let tolerance_exp = expected_exp - i64::from(prec) + 3;
        exp <= tolerance_exp
    }

    #[test]
    fn exp_zero_is_one() {
        let zero = BigFloat::try_new_zero(Sign::Positive, 53).unwrap();
        let (r, status) = zero.exp(RoundingMode::NearestEven);
        assert!(status.is_ok());
        let one = BigFloat::try_from_i64_exact(1, 53).unwrap();
        assert_eq!(r.partial_cmp(&one).0, Some(Ordering::Equal));
    }

    #[test]
    fn exp_neg_zero_is_one() {
        let nz = BigFloat::try_new_zero(Sign::Negative, 53).unwrap();
        let (r, _) = nz.exp(RoundingMode::NearestEven);
        let one = BigFloat::try_from_i64_exact(1, 53).unwrap();
        assert_eq!(r.partial_cmp(&one).0, Some(Ordering::Equal));
    }

    #[test]
    fn exp_pos_inf() {
        let pi = BigFloat::try_new_infinity(Sign::Positive, 53).unwrap();
        let (r, _) = pi.exp(RoundingMode::NearestEven);
        assert!(r.is_infinite());
        assert!(r.is_sign_positive());
    }

    #[test]
    fn exp_neg_inf() {
        let ni = BigFloat::try_new_infinity(Sign::Negative, 53).unwrap();
        let (r, _) = ni.exp(RoundingMode::NearestEven);
        assert!(r.is_zero());
        assert!(r.is_sign_positive());
    }

    #[test]
    fn exp_qnan_propagates() {
        let q = BigFloat::try_new_quiet_nan(Sign::Positive, 53, &[]).unwrap();
        let (r, _) = q.exp(RoundingMode::NearestEven);
        assert!(r.is_quiet_nan());
    }

    #[test]
    fn exp_snan_raises_invalid() {
        let sn = BigFloat::try_new_signaling_nan(Sign::Positive, 53, &[]).unwrap();
        let (r, status) = sn.exp(RoundingMode::NearestEven);
        assert!(r.is_quiet_nan());
        assert!(status.invalid());
    }

    #[test]
    fn exp_one_is_e() {
        let one = BigFloat::try_from_i64_exact(1, 53).unwrap();
        let (r, _) = one.exp(RoundingMode::NearestEven);
        // e ≈ 2.718281828459045
        let e = parse("2.718281828459045", 53);
        assert!(close_at(&r, &e, 53), "exp(1) = {r:?}, expected ≈ {e:?}");
    }

    #[test]
    fn exp_neg_one() {
        let neg_one = BigFloat::try_from_i64_exact(-1, 53).unwrap();
        let (r, _) = neg_one.exp(RoundingMode::NearestEven);
        // e^-1 ≈ 0.36787944117144233
        let expected = parse("0.36787944117144233", 53);
        assert!(close_at(&r, &expected, 53));
    }

    #[test]
    fn exp_ln2_is_two() {
        // exp(ln(2)) = 2
        let ln2 = parse("0.6931471805599453", 53);
        let (r, _) = ln2.exp(RoundingMode::NearestEven);
        let two = BigFloat::try_from_i64_exact(2, 53).unwrap();
        assert!(close_at(&r, &two, 53), "exp(ln(2)) should ≈ 2, got {r:?}");
    }

    #[test]
    fn exp_ten_is_about_22026() {
        // e^10 ≈ 22026.465794806718
        let ten = BigFloat::try_from_i64_exact(10, 53).unwrap();
        let (r, _) = ten.exp(RoundingMode::NearestEven);
        let expected = parse("22026.465794806718", 53);
        assert!(
            close_at(&r, &expected, 53),
            "exp(10) = {r}, expected ≈ {expected}"
        );
    }

    #[test]
    fn exp_large_negative() {
        // e^-30 ≈ 9.357622968840175e-14
        let neg_thirty = BigFloat::try_from_i64_exact(-30, 53).unwrap();
        let (r, _) = neg_thirty.exp(RoundingMode::NearestEven);
        let expected = parse("9.357622968840175e-14", 53);
        assert!(close_at(&r, &expected, 53), "exp(-30) = {r}");
    }

    #[test]
    fn exp_at_higher_precision() {
        let one = BigFloat::try_from_i64_exact(1, 113).unwrap();
        let (r, _) = one.exp(RoundingMode::NearestEven);
        assert_eq!(r.precision(), 113);
        let e = parse("2.718281828459045235360287471352662", 113);
        assert!(close_at(&r, &e, 113));
    }

    #[test]
    fn exp_with_explicit_round() {
        let one = BigFloat::try_from_i64_exact(1, 113).unwrap();
        let (r, _) = one.exp_round(53, RoundingMode::NearestEven);
        assert_eq!(r.precision(), 53);
        let e = parse("2.718281828459045", 53);
        assert!(close_at(&r, &e, 53));
    }

    #[cfg(feature = "fixed")]
    #[test]
    fn fixed_exp() {
        let one = FixedFloat::<53>::try_from_i64_exact(1).unwrap();
        let (r, _) = one.exp(RoundingMode::NearestEven);
        let e = parse("2.718281828459045", 53);
        assert!(close_at(&r.to_big(), &e, 53));
    }

    #[test]
    fn ln_2_constant_top_bits() {
        let ln2 = ln_2_at(53);
        // Top 8 bits of mantissa = 0xB1 (= 0b10110001), matching
        // ln(2) ≈ 0.1011000101110...
        match &ln2.class {
            Class::Normal { mantissa, .. } => {
                // Top byte of MSL.
                let top_byte = (*mantissa.last().unwrap()) >> 56;
                assert_eq!(top_byte, 0xB1);
            }
            _ => panic!("expected Normal"),
        }
    }
}
