//! Riemann zeta function `ζ(s)` for real argument `s` (DLMF
//! Chapter 25). Real-valued, single argument (no order parameter,
//! unlike Bessel).
//!
//! Two evaluation regimes (filled by slices 6r.2 / 6r.3):
//!
//! - `s > 0`, `s ≠ 1`: the Borwein / Cohen–Villegas–Zagier
//!   alternating-series acceleration ([`zeta_borwein`]). DLMF
//!   25.2.3 `ζ(s) = (1−2^{1−s})⁻¹ Σ_{k≥0} (−1)^k/(k+1)^s` holds for
//!   `ℜ s > 0`; `a_k = 1/(k+1)^s` is the moment sequence of a
//!   positive measure for `s > 0`, so CVZ Proposition 1's
//!   `2·(3+√8)^{−n}` relative-error bound applies on the whole open
//!   right half-line. **The roadmap/DESIGN.md originally specified
//!   Euler–Maclaurin reusing the 17-pair `gamma_stirling` Bernoulli
//!   table; that was found insufficient — `|B_{2k}/(2k)!| ≈
//!   2/(2π)^{2k}` caps a 17-term correction at ≈ 90 bits regardless
//!   of the sum length, well short of the bit-exact p = 1024
//!   differential lane. The algorithm was changed to Borwein; the
//!   deviation and its rationale are recorded in ADR-0026.**
//! - `s < 0`: the functional equation DLMF 25.4.2
//!   `ζ(s) = 2·(2π)^{s−1}·sin(πs/2)·Γ(1−s)·ζ(1−s)`, reflecting into
//!   `1−s > 1` where the Borwein core is well-conditioned (slice
//!   6r.3). Routes through the in-crate `π`, `pow`, `sin`, `Γ`.
//!
//! Special values are handled directly per this domain table
//! (DLMF 25.2 / 25.6, derived not recalled; ADR-0026):
//!
//! - `ζ(NaN) = NaN`; `sNaN` raises `INVALID`.
//! - `ζ(1) = +∞`, raising `DIV_BY_ZERO`: the only singularity is a
//!   simple pole at `s = 1` with residue `+1` (DLMF 25.2). The
//!   `s → 1⁺` side is `+∞`; `+∞` is the documented pole convention,
//!   the [`super::ci`] / [`super::li`] / [`super::bessel_k`]
//!   precedent.
//! - `ζ(0) = −1/2` exact (DLMF 25.6.1), for `±0`.
//! - `ζ(−2n) = +0` exact for `n ≥ 1` (DLMF 25.6.4): the trivial
//!   zeros at the negative even integers. Special-cased here so the
//!   functional-equation path's `sin(πs/2) = 0` cancellation does
//!   not have to produce an exact zero.
//! - `ζ(+∞) = 1`, `Status::OK`: a genuine limit (the Dirichlet
//!   series DLMF 25.2.1 collapses to its first term).
//! - `ζ(−∞) = NaN`, raising `INVALID`. Via the functional equation
//!   `Γ(1−s) → +∞` super-exponentially while `sin(πs/2)` oscillates
//!   in `[−1, 1]`, so `|ζ(s)|` grows without bound *and* does not
//!   converge. This is an **unbounded non-converging oscillation**,
//!   explicitly **not** the `J`/`Y`/Airy decaying-envelope
//!   convention (a *bounded* non-converging oscillation, where `+0`
//!   is a total-keeping choice). With no limit and no finite
//!   total-keeping value, the honest convention is `NaN` +
//!   `INVALID`. ADR-0026 records the distinction (the K-vs-Y `+∞`
//!   precedent, ADR-0025).

use super::lgamma::is_integer_test;
use super::pi_at;

use super::ziv_calibration::ZETA_ERROR_GUARD;
use crate::big::{BigFloat, BuildError};
use crate::class::Class;
use crate::rounding::RoundingMode;
use crate::sign::Sign;
use crate::status::{auto_raise, Status};

use core::cmp::Ordering;

#[cfg(feature = "fixed")]
use crate::fixed::FixedFloat;
#[cfg(feature = "fixed")]
use crate::mantissa::limbs_for;

impl BigFloat {
    /// `ζ(self)` rounded under `mode` to `self.precision`.
    #[must_use]
    pub fn zeta(&self, mode: RoundingMode) -> (Self, Status) {
        self.zeta_round(self.precision, mode)
            .expect("self.precision >= 1 by invariant")
    }

    /// `ζ(self)` with explicit result precision.
    ///
    /// Correctly rounded under every IEEE 754-2019 rounding mode
    /// via the shared `crate::math::ziv::ziv_round` driver (slice
    /// p1.34, ADR-0038).
    pub fn zeta_round(
        &self,
        target_precision: u32,
        mode: RoundingMode,
    ) -> Result<(Self, Status), BuildError> {
        if target_precision == 0 {
            return Err(BuildError::PrecisionZero);
        }
        Ok(zeta_kernel(self, target_precision, mode))
    }
}

#[cfg(feature = "fixed")]
impl<const PREC: u32> FixedFloat<PREC>
where
    [(); limbs_for(PREC)]:,
{
    /// `ζ(self)` for `FixedFloat`. Delegates to [`BigFloat::zeta`].
    #[must_use]
    pub fn zeta(&self, mode: RoundingMode) -> (Self, Status) {
        let (big, status) = self.to_big().zeta(mode);
        (
            Self::try_from_big_exact(big).expect("precision matches"),
            status,
        )
    }
}

/// Integer `v` as a `BigFloat` at precision `p` (exact for the small
/// integers this kernel forms).
fn ci(v: i64, p: u32) -> BigFloat {
    BigFloat::try_from_i64_exact(v, p).expect("precision >= 1")
}

/// `true` if `x` is a negative even integer (`−2, −4, −6, …`), the
/// trivial zeros `ζ(−2n) = 0`, `n ≥ 1`. `x` is integral and `x/2`
/// is also integral (division by two is exact, a pure exponent
/// shift, so the test is total and exact). `x = 0` is excluded by
/// the caller (the `Class::Zero` arm); odd negative integers fail
/// the `x/2` integrality test and route to the functional equation.
fn is_negative_even_integer(x: &BigFloat) -> bool {
    if !matches!(x.sign(), Sign::Negative) || !is_integer_test(x) {
        return false;
    }
    let two = ci(2, x.precision());
    let (half, _) = x.div(&two, RoundingMode::NearestEven);
    is_integer_test(&half)
}

/// `ζ(s)` for real `s`.
///
/// Special values are handled directly per the module-level domain
/// table; a finite non-special argument routes to [`zeta_finite`]
/// (the Borwein core for `s > 0`, the functional equation for
/// `s < 0`), then the single final round is applied.
fn zeta_kernel(x: &BigFloat, target_precision: u32, mode: RoundingMode) -> (BigFloat, Status) {
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
            (nan, Status::OK)
        }
        Class::Zero { .. } => {
            // ζ(0) = −1/2 exact (DLMF 25.6.1), for ±0. −1 ÷ 2 is an
            // exact exponent shift, so −1/2 is exact at any precision.
            let two = ci(2, target_precision);
            let (half, _) = ci(-1, target_precision).div(&two, RoundingMode::NearestEven);
            (half, Status::OK)
        }
        Class::Infinity { sign } => {
            if matches!(sign, Sign::Negative) {
                // ζ(−∞): unbounded non-converging oscillation via the
                // functional equation (Γ(1−s) → ∞, sin oscillates).
                // Not the decaying-envelope convention; no limit.
                let nan = BigFloat::try_new_quiet_nan(Sign::Positive, target_precision, &[])
                    .expect("precision >= 1");
                auto_raise(Status::INVALID);
                return (nan, Status::INVALID);
            }
            // ζ(+∞) = 1, the genuine Dirichlet-series limit (DLMF
            // 25.2.1: every term past the first vanishes).
            (ci(1, target_precision), Status::OK)
        }
        Class::Normal { .. } => {
            // Pole: ζ(1) = +∞ + DIV_BY_ZERO (simple pole, residue +1,
            // DLMF 25.2; +∞ is the s → 1⁺ side, the Ci/li/K pole
            // convention).
            let one = ci(1, x.precision());
            if matches!(x.partial_cmp(&one).0, Some(Ordering::Equal)) {
                let pinf = BigFloat::try_new_infinity(Sign::Positive, target_precision)
                    .expect("precision >= 1");
                auto_raise(Status::DIV_BY_ZERO);
                return (pinf, Status::DIV_BY_ZERO);
            }

            // Trivial zeros: ζ(−2n) = +0 exact, n ≥ 1 (DLMF 25.6.4).
            if is_negative_even_integer(x) {
                let zero = BigFloat::try_new_zero(Sign::Positive, target_precision)
                    .expect("precision >= 1");
                return (zero, Status::OK);
            }

            // Large s: ζ(s) = 1 + Σ_{k≥2} k^{−s} approaches 1 from above
            // (the 2^{−s} term dominates the positive residual), never
            // reaching it. Past the Ziv guard cap the residual underflows
            // every working precision, the interval test never converges,
            // and the fallback returns the on-grid limit — correct under
            // nearest and TowardNegative, but wrong under TowardPositive
            // (which must round up to 1 + ulp) and TowardZero (the capped
            // composition lands a hair the wrong side of 1). Short-circuit
            // to the mode-aware rounding of 1 with a magnitude-growing
            // infinitesimal (the saturation analogue of the constant
            // special cases; see `super::saturation_threshold_exponent`).
            if matches!(x.sign(), Sign::Positive)
                && exponent_of(x) >= super::saturation_threshold_exponent(target_precision)
            {
                let one = ci(1, target_precision);
                return crate::rounding::round_with_infinitesimal(
                    &one,
                    Sign::Positive,
                    false,
                    target_precision,
                    mode,
                );
            }

            // Ziv-driven correct rounding under every IEEE mode. The
            // eval closure dispatches on sign(s) (binary, mode- and
            // precision-independent) to zeta_borwein (s > 0) or
            // zeta_fe (s < 0). Both helpers compose through gamma
            // (Ziv-driven, p1.29), sin (p1.26), pow (already
            // five-mode), and a recursive zeta(1-s) call that routes
            // through the Borwein branch (1 - s > 1 when s < 0), so
            // there is no infinite recursion. The FE branch's
            // composition is now correct under every mode because
            // every constituent is.
            // The certification depth near the pole is input-encoded
            // (pf-jl35, ADR-0103): ζ(1 ± ε) ≈ ±1/ε + γ sits |e(ε)|
            // bits of relative distance from the on-grid ±2^k, past
            // the driver's fixed cap for deep ε. The lazy hint
            // recomputes the ADR-0098 conditioning probe (Sterbenz-
            // exact s − 1 at the input's own precision) only when
            // the legacy schedule exhausts.
            let (result, status) = super::ziv::ziv_round_with_depth(
                |w| zeta_finite(x, w),
                target_precision,
                mode,
                ZETA_ERROR_GUARD,
                || {
                    let probe = x.precision().max(target_precision.saturating_add(8));
                    let one_probe = ci(1, probe);
                    let (s_minus_1, _) = x
                        .round_to_precision(probe, RoundingMode::NearestEven)
                        .expect("precision >= 1")
                        .0
                        .sub(&one_probe, RoundingMode::NearestEven);
                    u32::try_from(exponent_of(&s_minus_1).min(0).unsigned_abs()).unwrap_or(u32::MAX)
                },
            );
            // Defensive INEXACT guard (pf-umlm, ADR-0066). The dispatched
            // dyadic outputs (ζ(0) = −1/2, the trivial zeros ζ(−2n) = 0,
            // ζ(+∞) = 1) return above; the negative-odd-integer rationals
            // (ζ(−1) = −1/12, …) are non-dyadic and already flag INEXACT,
            // and ζ(s) > 1 strictly for s > 1 makes the large-s collapse
            // to 1.0 a true INEXACT. The ADR-0065 sweep found this path
            // flags INEXACT everywhere, so the force is a no-op hardening
            // against regression; its worst-case soundness rests on the
            // irrationality of ζ at the non-special set (ζ(5) is open).
            let status = if matches!(result.class, Class::Normal { .. }) {
                status | Status::INEXACT
            } else {
                status
            };
            auto_raise(status);
            (result, status)
        }
    }
}

/// `ζ(s)` for a finite, non-special real `s` (not `0`, not `1`, not
/// a negative even integer).
///
/// - `s > 0`: the Borwein / Cohen–Villegas–Zagier alternating-series
///   acceleration ([`zeta_borwein`]). DLMF 25.2.3
///   `ζ(s) = (1−2^{1−s})⁻¹ Σ_{k≥0} (−1)^k/(k+1)^s` is valid for
///   `ℜ s > 0`, and `a_k = 1/(k+1)^s = Γ(s)⁻¹∫₀¹ xᵏ(−ln x)^{s−1}dx`
///   is the moment sequence of a positive measure for `s > 0`, so
///   the CVZ Proposition 1 error bound applies on the whole open
///   right half-line.
/// - `s < 0`: the functional equation DLMF 25.4.2 ([`zeta_fe`])
///   reflecting into `1−s > 1`, where the Borwein core is
///   well-conditioned.
fn zeta_finite(x: &BigFloat, target_precision: u32) -> BigFloat {
    if matches!(x.sign(), Sign::Positive) {
        zeta_borwein(x, target_precision)
    } else {
        zeta_fe(x, target_precision)
    }
}

/// `ζ(s)` for real `s < 0`, `s` not a negative even integer (those
/// are the trivial zeros, special-cased exactly upstream), via the
/// functional equation DLMF 25.4.2
///
/// ```text
/// ζ(s) = 2·(2π)^{s−1}·sin(πs/2)·Γ(1−s)·ζ(1−s)
/// ```
///
/// quoted verbatim from the DLMF primary source, **not recalled**
/// (the 25.4.1 `cos(πs/2)·Γ(s)` companion is the wrong branch for
/// the negative axis; the sin-vs-cos / Γ(1−s)-vs-Γ(s) choice is the
/// derive-don't-recall catch, the ADR-0025 K-recurrence-sign
/// precedent). For `s < 0` the reflected argument `1−s > 1` lands
/// strictly inside the Borwein moment region, so `ζ(1−s)` is the
/// well-conditioned [`zeta_borwein`] path. Numerically pinned by
/// `ζ(−1) = 2·(2π)⁻²·(−1)·1·(π²/6) = −1/12`. The factors multiply
/// with no cancellation (`Γ(1−s)` grows, `(2π)^{s−1}` decays, the
/// product is the genuine — possibly large — value), so a single
/// `+96`-bit working boost over the gamma/sin/pow composition
/// suffices. Returns the unrounded working-precision value.
fn zeta_fe(x: &BigFloat, target_precision: u32) -> BigFloat {
    let working = target_precision
        .saturating_add(96)
        .min(target_precision.saturating_add(4096));

    let sw = x
        .round_to_precision(working, RoundingMode::NearestEven)
        .expect("precision >= 1")
        .0;
    let one = ci(1, working);
    let two = ci(2, working);
    let pi = pi_at(working);

    // 2·(2π)^{s−1}.
    let (two_pi, _) = two.mul(&pi, RoundingMode::NearestEven);
    let (s_minus_1, _) = sw.sub(&one, RoundingMode::NearestEven);
    let (tp_pow, _) = two_pi.pow(&s_minus_1, RoundingMode::NearestEven);
    let (coeff, _) = two.mul(&tp_pow, RoundingMode::NearestEven);

    // sin(πs/2).
    let (pi_s, _) = pi.mul(&sw, RoundingMode::NearestEven);
    let (arg, _) = pi_s.div(&two, RoundingMode::NearestEven);
    let (sin_term, _) = arg.sin(RoundingMode::NearestEven);

    // Γ(1−s), with 1−s > 1.
    let (one_minus_s, _) = one.sub(&sw, RoundingMode::NearestEven);
    let (gamma_term, _) = one_minus_s.gamma(RoundingMode::NearestEven);

    // ζ(1−s): 1−s > 1 > 0, the Borwein moment region.
    let zeta_reflected = zeta_borwein(&one_minus_s, working);

    let (p1, _) = coeff.mul(&sin_term, RoundingMode::NearestEven);
    let (p2, _) = p1.mul(&gamma_term, RoundingMode::NearestEven);
    let (result, _) = p2.mul(&zeta_reflected, RoundingMode::NearestEven);
    result
}

/// Binary exponent of `v`, or `0` for zero / non-finite (the
/// conditioning-estimate idiom; only used to size the working
/// precision, never for a value decision).
fn exponent_of(v: &BigFloat) -> i64 {
    match &v.class {
        Class::Normal { exponent, .. } => *exponent,
        _ => 0,
    }
}

/// `ζ(s)` for real `s > 0`, `s ≠ 1`, via the Cohen–Villegas–Zagier
/// "Algorithm 1" acceleration of the Dirichlet eta series.
///
/// DLMF 25.2.3: `ζ(s) = (1−2^{1−s})⁻¹ Σ_{k=0}^∞ (−1)^k/(k+1)^s`
/// for `ℜ s > 0`. With `a_k = 1/(k+1)^s` (the moments of the
/// positive measure `Γ(s)⁻¹(−ln x)^{s−1}dx` on `(0,1)` for
/// `s > 0`), CVZ Algorithm 1 (Cohen, Rodriguez-Villegas, Zagier,
/// *Experiment. Math.* 9 (2000), the Borwein algorithm) computes
/// `S = Σ (−1)^k a_k` to relative error `≤ 2·(3+√8)^{−n}` from the
/// first `n` terms (CVZ Proposition 1). The algorithm, **verified
/// against the paper's own worked `n = 1, 2` examples** (`2a₀/3`,
/// `(16a₀−8a₁)/17`) rather than recalled (the derive-don't-recall
/// reflex; ADR-0026):
///
/// ```text
/// d = (3+√8)^n;  d = (d + 1/d)/2;   // = d_n, an integer
/// b = −1;  c = −d;  s = 0;
/// for k = 0 .. n−1:
///     c = b − c
///     s = s + c · a_k                // a_k = (k+1)^{−s}
///     b = (k+n)(k−n)·b / ((k+½)(k+1))
/// return s / d
/// ```
///
/// Then `ζ(s) = (s/d) / (1−2^{1−s})`. The intermediate `s` reaches
/// magnitude `≈ d_n·ζ(s)`, and `d_n ≈ ½(3+√8)^n ≈ 2^{2.543 n}`; with
/// `n ≈ p·ln2/ln(3+√8) ≈ 0.3933 p` terms for `p` target bits this is
/// `≈ 2^p`, so recovering `p` bits of `ζ` after the `/d_n` needs a
/// working precision `≈ 2p` (derived, not guessed — the
/// `gamma_stirling` boost analog). The `(1−2^{1−s})⁻¹` factor blows
/// up as `s → 1`; the working precision absorbs that with an
/// additional `−log₂|s−1|` bits. Returns the unrounded
/// working-precision value; [`zeta_kernel`] does the single round.
fn zeta_borwein(x: &BigFloat, target_precision: u32) -> BigFloat {
    // n: CVZ Proposition 1 relative error ≤ 2·(3+√8)^{−n}. To reach
    // 2^{−p}: n ≥ p·ln2/ln(3+√8) = p·0.39321…; 3933/10000 is a
    // safe rational over-estimate, plus a guard.
    let n: i64 = (u64::from(target_precision).saturating_mul(3933) / 10_000) as i64 + 16;

    // Conditioning: 1/(1−2^{1−s}) loses ≈ −log₂|s−1| bits near
    // s = 1. |s−1| ≥ 2^{exponent_of(s−1)}, so boost by that many
    // bits (plus slack). s = 1 itself is special-cased upstream.
    // The probe must see the INPUT's full precision (pf-gg96,
    // ADR-0098): rounding s to target+8 bits first collapsed
    // s = 1 + 2^-5000 (p5001) to exactly 1, the boost came out 0,
    // the working round then made 1 − 2^{1−s} exactly 0, and the
    // discarded DIV_BY_ZERO became a certified +Inf with Status OK.
    // s − 1 at max(input precision, probe) is Sterbenz-exact near 1
    // and its exponent is what the boost needs everywhere else.
    let probe = x.precision().max(target_precision.saturating_add(8));
    let one_probe = ci(1, probe);
    let (s_minus_1, _) = x
        .round_to_precision(probe, RoundingMode::NearestEven)
        .expect("precision >= 1")
        .0
        .sub(&one_probe, RoundingMode::NearestEven);
    let cond_boost = exponent_of(&s_minus_1).min(0).unsigned_abs() as u32;

    // d_n ≈ 2^{2.543 n}; intermediate ≈ d_n·ζ, so absolute bits
    // needed ≈ p + log₂ d_n ≈ 2p. 2·target + 96 + the conditioning
    // boost is the derived working precision. The cap scales with
    // the input precision: the conditioning boost is bounded by the
    // input-encoded proximity to 1, itself bounded by the input
    // precision, so the cost stays proportional to what the caller
    // already supplied (a fixed cap re-collapses deep inputs).
    let working = target_precision
        .saturating_mul(2)
        .saturating_add(96)
        .saturating_add(cond_boost)
        .min(
            target_precision
                .saturating_mul(2)
                .saturating_add(8192)
                .saturating_add(x.precision()),
        );

    let sw = x
        .round_to_precision(working, RoundingMode::NearestEven)
        .expect("precision >= 1")
        .0;
    let neg_s = sw.negated();
    let one = ci(1, working);
    let two = ci(2, working);

    // d = (3+√8)^n; d = (d + 1/d)/2  →  d_n = ½((3+√8)^n+(3−√8)^n).
    let (sqrt8, _) = ci(8, working).sqrt(RoundingMode::NearestEven);
    let (base, _) = ci(3, working).add(&sqrt8, RoundingMode::NearestEven);
    let mut d = one.clone();
    for _ in 0..n {
        d = d.mul(&base, RoundingMode::NearestEven).0;
    }
    let (inv_d, _) = one.div(&d, RoundingMode::NearestEven);
    let (d_sum, _) = d.add(&inv_d, RoundingMode::NearestEven);
    let (d, _) = d_sum.div(&two, RoundingMode::NearestEven);

    // b = −1; c = −d; s = 0.
    let mut b = one.negated();
    let mut c = d.negated();
    let mut sum = BigFloat::try_new_zero(Sign::Positive, working).expect("precision >= 1");

    for k in 0..n {
        // c = b − c.
        let (c_new, _) = b.sub(&c, RoundingMode::NearestEven);
        c = c_new;

        // a_k = (k+1)^{−s}.
        let base_k = ci(k + 1, working);
        let (a_k, _) = base_k.pow(&neg_s, RoundingMode::NearestEven);
        let (term, _) = c.mul(&a_k, RoundingMode::NearestEven);
        let (s_new, _) = sum.add(&term, RoundingMode::NearestEven);
        sum = s_new;

        // b = (k+n)(k−n)·b / ((k+½)(k+1)).
        let num = (k + n) * (k - n); // negative (0 ≤ k < n)
        let (b1, _) = ci(num, working).mul(&b, RoundingMode::NearestEven);
        // k + ½ = (2k+1)/2.
        let (half, _) = ci(2 * k + 1, working).div(&two, RoundingMode::NearestEven);
        let (den, _) = half.mul(&ci(k + 1, working), RoundingMode::NearestEven);
        let (b_new, _) = b1.div(&den, RoundingMode::NearestEven);
        b = b_new;
    }

    // η(s) = s / d_n.
    let (eta, _) = sum.div(&d, RoundingMode::NearestEven);

    // ζ(s) = η(s) / (1 − 2^{1−s}).
    let (one_minus_s, _) = one.sub(&sw, RoundingMode::NearestEven);
    let (two_pow, _) = two.pow(&one_minus_s, RoundingMode::NearestEven);
    let (factor, _) = one.sub(&two_pow, RoundingMode::NearestEven);
    // Defensive belt (pf-gg96, ADR-0098): a factor of exactly 0
    // would make eta/factor an infinity whose half_width(inf) = 0
    // the Ziv driver certifies silently — the collapsed-special
    // trap. The input-precision conditioning boost above makes this
    // unreachable (s ≠ 1 is dispatched upstream and the working
    // precision now resolves 2^{1−s} ≠ 1), so refuse loudly with a
    // NaN the driver cannot certify rather than divide.
    if factor.is_zero() {
        return BigFloat::try_new_quiet_nan(Sign::Positive, working, &[]).expect("precision >= 1");
    }
    let (zeta, _) = eta.div(&factor, RoundingMode::NearestEven);
    zeta
}

#[cfg(test)]
mod tests {
    use super::*;
    use core::cmp::Ordering as Ord2;

    /// `|v − expected| ≤ 2^(−bits)·|expected|` (the bessel test
    /// helper). References: `mpmath.zeta` at `mp.dps = 340`
    /// (`nix-shell -p 'python3.withPackages(ps:[ps.mpmath])'`),
    /// treated as a fact.
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
            Some(Ord2::Less | Ord2::Equal)
        )
    }

    fn py(s: &str, p: u32) -> BigFloat {
        BigFloat::parse_str(s, p, RoundingMode::NearestEven)
            .unwrap()
            .0
    }

    fn at(n: i64, d: i64, p: u32) -> BigFloat {
        let nn = BigFloat::try_from_i64_exact(n, p).unwrap();
        if d == 1 {
            nn
        } else {
            nn.div(
                &BigFloat::try_from_i64_exact(d, p).unwrap(),
                RoundingMode::NearestEven,
            )
            .0
        }
    }

    // mpmath.zeta, dps = 340 (truncated to fit p).
    const Z2: &str = "1.64493406684822643647241516664602518921894990120679843773555822937000747040320087383362890061975870530400431896233719067962872468700500778793510294633086627683173330936776260509525100687214005479681155879489036082327776191984075645587696323563670971009694890208593200805163647887833884604444518405982514525068";
    const Z3: &str = "1.20205690315959428539973816151144999076498629234049888179227155534183820578631309018645587360933525814619915779526071941849199599867328321377639683720790016145394178294936006671919157552224249424396156390966410329115909578096551465127991840510571525598801543710978110203982753256678760352233698494166181105701";
    const Z4: &str = "1.08232323371113819151600369654116790277475095191872690768297621544412061618696884655690963594169991723299081390804274241458407157457004534928200351471621920708778348091083702932618873482617527360423550621937375061711174534929686775073307606686934118905862833795279512033449589046886262694822083503298363214902";
    const Z15: &str = "2.61237534868548834334856756792407163057080065240006340757332824881492776768827286099624386812631195238297635877214975569815763296843445913443832056180833600833933396280548054166294852684829798168645847550187899242552790919645625985746620957819178983247798052614814070472260846524069586856423142070771015331232";
    const Z_HALF: &str = "-1.46035450880958681288949915251529801246722933101258149054288608782553052947450062527641937546335681951449637467986952958389234371035889426181923283975376292518263335864916412789122939415410119791731044810824194092788169842885717682395579918451788361465548665937991689152316352160424275374940796571353042261007";
    const Z10: &str = "1.00099457512781808533714595890031901700601953156447751725778899463629146515191295439704196861038565275400689206320530767736809020353629380731906959498428739536216033347223525967320521789323288320665440138759279913286048883976147693647789769806971192063361022944054388731501219022076400989382492087774683640358";
    // ζ(2) and ζ(3) at 1024-bit precision (dps = 340).
    const Z2_BIG: &str = Z2;
    const Z3_BIG: &str = Z3;
    // ζ near the pole: s = 1.0009765625 (= 1 + 1/1024). The
    // 1/(1−2^{1−s}) factor amplifies to ≈ 1024, exercising the
    // derived conditioning boost.
    const Z_NEAR1: &str = "1024.5772867695045940578681624248887776501597556226467113160352190702981219581341444863800913012818856950855998236564103910760853485875372916860380523003939719972891239530972602843354325731818616935853605173448667011466739027602092030574417526225201845179601125220008848142816634721046490224598777861367925203549";

    /// Borwein/CVZ core on the moment region `s > 0`: closed-form
    /// pins ζ(2)=π²/6, ζ(4)=π⁴/90, the irrational ζ(3), the
    /// critical-line value ζ(1/2) (negative), ζ(3/2), ζ(10), each
    /// vs the 340-digit mpmath reference at p = 240. Spans
    /// 0 < s < 1 (where the Dirichlet series itself diverges but
    /// the eta acceleration converges) and s > 1.
    #[test]
    fn borwein_matches_reference() {
        let p = 240;
        for (num, den, want) in [
            (2i64, 1i64, Z2),
            (4, 1, Z4),
            (3, 1, Z3),
            (1, 2, Z_HALF),
            (3, 2, Z15),
            (10, 1, Z10),
        ] {
            let s = at(num, den, p);
            let (r, st) = s.zeta(RoundingMode::NearestEven);
            assert!(!st.invalid() && !st.div_by_zero());
            assert!(close_at(&r, &py(want, p), p - 16), "ζ({num}/{den})");
        }
    }

    /// Large `s`: ζ(50) is 1 to within 2^{−48} (every eta term past
    /// the first is negligible). Sanity that the accelerator does
    /// not mis-handle the near-constant regime.
    #[test]
    fn borwein_large_argument_is_one() {
        let p = 113;
        let s = at(50, 1, p);
        let (r, _) = s.zeta(RoundingMode::NearestEven);
        let one = BigFloat::try_from_i64_exact(1, p).unwrap();
        assert!(close_at(&r, &one, 48), "ζ(50) ≈ 1");
    }

    /// Conditioning near the pole: s = 1 + 1/1024. ζ ≈ 1024.577;
    /// the 1/(1−2^{1−s}) factor loses ≈ 10 bits, absorbed by the
    /// derived `cond_boost`. A wrong boost (or sign) fails here.
    #[test]
    fn borwein_near_pole_conditioning() {
        let p = 160;
        // 1 + 1/1024 exactly (dyadic).
        let s = {
            let one = BigFloat::try_from_i64_exact(1, p).unwrap();
            let d = BigFloat::try_from_i64_exact(1024, p).unwrap();
            one.add(
                &BigFloat::try_from_i64_exact(1, p)
                    .unwrap()
                    .div(&d, RoundingMode::NearestEven)
                    .0,
                RoundingMode::NearestEven,
            )
            .0
        };
        let (r, st) = s.zeta(RoundingMode::NearestEven);
        assert!(!st.invalid());
        assert!(close_at(&r, &py(Z_NEAR1, p), p - 24), "ζ(1+1/1024)");
    }

    /// Second-term-matters pin at p = 1024 (the derive-don't-recall
    /// reflex): ζ(2)=π²/6 and ζ(3) bit-accurate to p−4 against the
    /// 340-digit references. A `d_k` recurrence error, a wrong sign
    /// in `c = b − c` / the `b` update, or an under-sized working
    /// precision is invisible at low p and fails catastrophically
    /// here.
    #[test]
    fn borwein_high_precision_pin() {
        let s = at(2, 1, 1024);
        let (r, _) = s.zeta(RoundingMode::NearestEven);
        assert!(close_at(&r, &py(Z2_BIG, 1024), 1020), "ζ(2) p=1024");
        let s = at(3, 1, 1024);
        let (r, _) = s.zeta(RoundingMode::NearestEven);
        assert!(close_at(&r, &py(Z3_BIG, 1024), 1020), "ζ(3) p=1024");
    }

    // mpmath.zeta on the negative axis, dps = 340. ζ(−1)=−1/12,
    // ζ(−3)=1/120, ζ(−5)=−1/252, ζ(−7)=1/240 (exact rationals,
    // routed entirely through the functional equation since the
    // negative odd integers are NOT special-cased).
    const ZN1: &str = "-0.0833333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333";
    const ZN3: &str = "0.00833333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333";
    const ZN5: &str = "-0.00396825396825396825396825396825396825396825396825396825396825396825396825396825396825396825396825396825396825396825396825396825396825396825396825396825396825396825396825396825396825396825396825396825396825396825396825396825396825396825396825396825397";
    const ZN_HALF: &str = "-0.207886224977354566017306725397049302226268531287672537610113557106147291932292340487543266940733215643109975614128689565661326914694458311965705623294109531061640017807007041375078320755666248787786920661504691428291233832569371613677729383610945938789";
    const ZN15: &str = "-0.0254852018898330359495429869107047454690249846009729968346454983492493771883392785970925189475243606083956786090522338379231706922081480226595534733619982128781578133045264443099251947318696331135537227085510583845600615587030880946621048463536127324092";
    const ZN25: &str = "0.00851692877785033054235856702834448693627599022007447776588885495191457755599180493669488160134326161961109140844500462979255917397088304193018048008544285079380916536592281145828175781513758799235867264721263349427536024814776909920531663911328610074";
    const ZN105: &str = "0.0111461224739428141361386754700669760847882533145070950020015456357735060427323428518682330397325401032561182052740841597292413751057605763603111924201930408938517933090311436848472127517446387995128075773617523581154995570298259266384590180139117838251";

    /// Functional-equation negative axis (DLMF 25.4.2): exact
    /// rationals ζ(−1)=−1/12, ζ(−3)=1/120, ζ(−5)=−1/252, ζ(−7)
    /// =1/240 (each computed via the full FE, the negative odd
    /// integers deliberately NOT special-cased), plus the
    /// non-integer ζ(−1/2), ζ(−3/2), ζ(−5/2), ζ(−21/2) vs the
    /// 340-digit mpmath reference at p = 128. (The FE composes
    /// Γ+sin+pow+Borwein, so it is costly per
    /// `feedback_differential_lane_cost`; p = 128 keeps the
    /// lib-test tier fast — exhaustive precision coverage is the
    /// 6r.4 differential lane's job, and `fe_high_precision_pin`
    /// carries the p = 512 second-term pin.)
    #[test]
    fn fe_matches_reference() {
        let p = 128;
        for (num, den, want) in [
            (-1i64, 1i64, ZN1),
            (-3, 1, ZN3),
            (-5, 1, ZN5),
            (-7, 1, "0.004166666666666666666666666666666666666666666666666666666666666666666666666666666666666666666666666666666666666666666666666666666666666666666666666666666666666666666666666666666666666666666666666666666666666666666666666666666666666666666666666666667"),
            (-1, 2, ZN_HALF),
            (-3, 2, ZN15),
            (-5, 2, ZN25),
            (-21, 2, ZN105),
        ] {
            let s = at(num, den, p);
            let (r, st) = s.zeta(RoundingMode::NearestEven);
            assert!(!st.invalid() && !st.div_by_zero());
            assert!(close_at(&r, &py(want, p), p - 24), "ζ({num}/{den})");
        }
    }

    /// Functional-equation cross-tie / second-term pin at p = 512:
    /// ζ(−1) = −1/12 to p−4 bits. Binds the FE composition
    /// (2·(2π)^{s−1}·sin·Γ·ζ(1−s), Borwein on the s>0 side) against
    /// the independent exact rational. A sin-vs-cos or
    /// Γ(1−s)-vs-Γ(s) error fails here.
    #[test]
    fn fe_high_precision_pin() {
        let s = at(-1, 1, 512);
        let (r, _) = s.zeta(RoundingMode::NearestEven);
        assert!(close_at(&r, &py(ZN1, 512), 508), "ζ(−1) = −1/12 p=512");
    }

    #[test]
    fn zeta_quiet_nan_propagates() {
        let q = BigFloat::try_new_quiet_nan(Sign::Positive, 53, &[]).unwrap();
        let (r, s) = q.zeta(RoundingMode::NearestEven);
        assert!(r.is_quiet_nan());
        assert!(!s.invalid());
    }

    #[test]
    fn zeta_signaling_nan_raises_invalid() {
        let sn = BigFloat::try_new_signaling_nan(Sign::Positive, 53, &[]).unwrap();
        let (r, st) = sn.zeta(RoundingMode::NearestEven);
        assert!(r.is_quiet_nan());
        assert!(st.invalid());
    }

    #[test]
    fn zeta_pole_at_one() {
        let one = BigFloat::try_from_i64_exact(1, 113).unwrap();
        let (r, s) = one.zeta(RoundingMode::NearestEven);
        assert!(r.is_infinite() && r.is_sign_positive(), "ζ(1) = +∞");
        assert!(s.div_by_zero());
    }

    #[test]
    fn zeta_at_zero_is_minus_half() {
        // ζ(0) = −1/2 exact, for both +0 and −0.
        for sign in [Sign::Positive, Sign::Negative] {
            let z = BigFloat::try_new_zero(sign, 113).unwrap();
            let (r, s) = z.zeta(RoundingMode::NearestEven);
            assert!(!s.invalid() && !s.div_by_zero());
            let expected = {
                let two = BigFloat::try_from_i64_exact(2, 113).unwrap();
                BigFloat::try_from_i64_exact(-1, 113)
                    .unwrap()
                    .div(&two, RoundingMode::NearestEven)
                    .0
            };
            assert!(
                matches!(r.partial_cmp(&expected).0, Some(Ordering::Equal)),
                "ζ(0) = −1/2"
            );
        }
    }

    #[test]
    fn zeta_trivial_zeros_at_negative_even_integers() {
        for k in [-2i64, -4, -6, -10, -42] {
            let s = BigFloat::try_from_i64_exact(k, 113).unwrap();
            let (r, st) = s.zeta(RoundingMode::NearestEven);
            assert!(
                r.is_zero() && r.is_sign_positive(),
                "ζ({k}) = +0 (trivial zero)"
            );
            assert!(!st.invalid());
        }
    }

    #[test]
    fn zeta_negative_odd_integers_are_not_trivial_zeros() {
        // ζ(−1), ζ(−3), … are nonzero rationals; they must NOT be
        // special-cased to zero (they route to the functional
        // equation in 6r.3, here the stubbed finite path → NaN).
        for k in [-1i64, -3, -5] {
            let s = BigFloat::try_from_i64_exact(k, 53).unwrap();
            let (r, _) = s.zeta(RoundingMode::NearestEven);
            assert!(!r.is_zero(), "ζ({k}) is not the zero special-case");
        }
    }

    #[test]
    fn zeta_plus_infinity_is_one() {
        let inf = BigFloat::try_new_infinity(Sign::Positive, 53).unwrap();
        let (r, st) = inf.zeta(RoundingMode::NearestEven);
        let one = BigFloat::try_from_i64_exact(1, 53).unwrap();
        assert!(
            matches!(r.partial_cmp(&one).0, Some(Ordering::Equal)),
            "ζ(+∞) = 1"
        );
        assert!(!st.invalid());
    }

    #[test]
    fn zeta_negative_infinity_is_invalid() {
        let ninf = BigFloat::try_new_infinity(Sign::Negative, 53).unwrap();
        let (r, s) = ninf.zeta(RoundingMode::NearestEven);
        assert!(r.is_quiet_nan(), "ζ(−∞) = NaN (no limit)");
        assert!(s.invalid());
    }

    #[test]
    fn precision_zero_is_rejected() {
        let x = BigFloat::try_from_i64_exact(2, 53).unwrap();
        assert!(x.zeta_round(0, RoundingMode::NearestEven).is_err());
    }

    // ---- FE-branch constituent diagnostic probe ------------------
    //
    // Investigation gate for the `differential_zeta` p=1024 failure
    // on ζ(-1/2). The differential lane reports a ~100-bit residual
    // between pfloat and MPFR; the prior multi-precision probe
    // showed bit-exact agreement at p ∈ {53, 113, 256, 512}. The
    // boundary is the FE branch's working precision (target + 96 =
    // 1120 at target=1024) crossing the 1024-bit `PI_LIMBS_1024`
    // table cap, forcing fallback to `agm_constants::pi_via_agm`.
    // This probe pins the constituents at working precision and
    // emits each per-step comparison.

    /// Print bit-position of disagreement (or "agree") for each FE
    /// constituent at p=1120. Compile-gated on unix because the
    /// `rug` dev-dependency is unix-only in `Cargo.toml`.
    #[cfg(unix)]
    #[test]
    #[ignore = "diagnostic only; un-ignore to investigate the FE branch divergence"]
    fn diag_zeta_neg_half_fe_constituents_p1024() {
        use core::cmp::Ordering;

        use rug::float::Round;
        use rug::ops::Pow as _;
        use rug::Float;

        use crate::math::pi_at;

        let target: u32 = 1024;
        let working: u32 = target + 96;

        fn report(label: &str, working: u32, pf: &BigFloat, rg: &Float) {
            // Bit-exact construction of rug::Float from pfloat
            // BigFloat. Mirrors `tests/differential/mod.rs::bigfloat_to_rug`
            // inline because that helper lives in the dev-test
            // crate tree, not in src/.
            let pf_as_rug: Float = match &pf.class {
                Class::Zero { sign } => {
                    let mut f = Float::with_val(pf.precision, 0);
                    if matches!(sign, Sign::Negative) {
                        f = -f;
                    }
                    f
                }
                Class::Infinity { sign } => {
                    let mut f = Float::with_val(pf.precision, rug::float::Special::Infinity);
                    if matches!(sign, Sign::Negative) {
                        f = -f;
                    }
                    f
                }
                Class::Nan { .. } => Float::with_val(pf.precision, rug::float::Special::Nan),
                Class::Normal {
                    sign,
                    exponent,
                    mantissa,
                } => {
                    let int = rug::Integer::from_digits(mantissa, rug::integer::Order::Lsf);
                    let mut f = Float::with_val(pf.precision, &int);
                    let stored_bits = (mantissa.len() as i64) * 64;
                    let shift: i64 = exponent + 1 - stored_bits;
                    let mut remaining = shift;
                    while remaining != 0 {
                        let step = if remaining >= 0 {
                            remaining.min(i64::from(i32::MAX)) as i32
                        } else {
                            remaining.max(i64::from(i32::MIN)) as i32
                        };
                        f <<= step;
                        remaining -= i64::from(step);
                    }
                    if matches!(sign, Sign::Negative) {
                        f = -f;
                    }
                    f
                }
            };

            let agree = pf_as_rug == *rg;
            if agree {
                eprintln!("{label:36} p={working:5}  AGREE");
            } else {
                let diff = Float::with_val(working + 64, &pf_as_rug - rg);
                let value_exp: i64 = rg.get_exp().unwrap_or(0).into();
                let diff_exp: i64 = diff.get_exp().unwrap_or(0).into();
                let bits_agree = i64::from(working) - 1 - (diff_exp - value_exp);
                eprintln!("{label:36} p={working:5}  DIFFER bits-of-agreement={bits_agree}");
            }
        }

        // π at working precision: pfloat via pi_at(1120) (AGM path
        // since 1120 > 1024), MPFR via Constant::Pi.
        let pi_pf = pi_at(working);
        let pi_rg = Float::with_val(working, rug::float::Constant::Pi);
        report("pi_at(working)", working, &pi_pf, &pi_rg);

        // Also: the leading 1024 bits of pi_pf should match the
        // static PI_LIMBS_1024 table (which is multiply-pinned).
        let pi_at_table_prec = pi_at(1024);
        let pi_rg_at_table_prec = Float::with_val(1024, rug::float::Constant::Pi);
        report(
            "pi_at(1024) — table path",
            1024,
            &pi_at_table_prec,
            &pi_rg_at_table_prec,
        );

        // sin(-π/4) at working.
        let neg_one = BigFloat::try_from_i64_exact(-1, working).unwrap();
        let two = BigFloat::try_from_i64_exact(2, working).unwrap();
        let (sw, _) = neg_one.div(&two, RoundingMode::NearestEven);
        let (pi_s, _) = pi_pf.mul(&sw, RoundingMode::NearestEven);
        let (arg, _) = pi_s.div(&two, RoundingMode::NearestEven);
        let (sin_pf, _) = arg.sin(RoundingMode::NearestEven);

        let sw_rg = Float::with_val(working, -1) / Float::with_val(working, 2);
        let pi_s_rg = Float::with_val(working, &pi_rg * &sw_rg);
        let arg_rg = pi_s_rg / Float::with_val(working, 2);
        let sin_rg = Float::with_val_round(working, arg_rg.sin_ref(), Round::Nearest).0;
        report("sin(-π/4)", working, &sin_pf, &sin_rg);

        // Γ(3/2) at working.
        let one = BigFloat::try_from_i64_exact(1, working).unwrap();
        let (one_minus_s, _) = one.sub(&sw, RoundingMode::NearestEven);
        let (gamma_pf, _) = one_minus_s.gamma(RoundingMode::NearestEven);
        let one_minus_s_rg = Float::with_val(working, 1) - sw_rg.clone();
        let gamma_rg = Float::with_val_round(working, one_minus_s_rg.gamma_ref(), Round::Nearest).0;
        report("Γ(3/2)", working, &gamma_pf, &gamma_rg);

        // (2π)^(-3/2) at working.
        let (two_pi_pf, _) = two.mul(&pi_pf, RoundingMode::NearestEven);
        let (s_minus_1, _) = sw.sub(&one, RoundingMode::NearestEven);
        let (pow_pf, _) = two_pi_pf.pow(&s_minus_1, RoundingMode::NearestEven);
        let two_pi_rg = Float::with_val(working, 2) * &pi_rg;
        let s_minus_1_rg = sw_rg - Float::with_val(working, 1);
        let pow_rg =
            Float::with_val_round(working, (&two_pi_rg).pow(&s_minus_1_rg), Round::Nearest).0;
        report("(2π)^(-3/2)", working, &pow_pf, &pow_rg);

        // ζ(3/2) at working (Borwein branch — recursive zeta call).
        let (zeta_reflected_pf, _) = one_minus_s.zeta(RoundingMode::NearestEven);
        let zeta_reflected_rg =
            Float::with_val_round(working, one_minus_s_rg.zeta_ref(), Round::Nearest).0;
        report("ζ(3/2)", working, &zeta_reflected_pf, &zeta_reflected_rg);

        // Suppress the unused-import warnings if all comparisons agree.
        let _ = Ordering::Equal;
    }
}
