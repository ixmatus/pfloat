//! Bessel functions of the second kind `Y0`, `Y1`, `Yn` (DLMF
//! Chapter 10): ordinary Bessel of integer order, real argument.
//!
//! Unlike [`super::bessel_j`], `Y` is real-valued only for `x > 0`:
//! `Y_n` has a logarithmic branch point at the origin and is complex
//! for `x < 0`. The domain convention follows the [`super::ci`] /
//! [`super::li`] precedent (cosine / logarithmic integral, same
//! "real-only, complex off the positive axis" shape):
//!
//! - `Y_n(±0)` is a pole for BOTH zero signs (`log(±0) = −∞` groups
//!   `−0` with the pole, IEEE 754-2019 §9.2 / C11 F.10.3.7; pf-k8ax,
//!   ADR-0123), raising `DIV_BY_ZERO`. The sign follows the order
//!   parity `Y_{−n} = (−1)^n Y_n` (DLMF 10.4.1) and POSIX `yn`:
//!   `Y_n(±0) = −∞` unless `n` is negative AND odd, then `+∞`.
//! - `x < 0` (and `−∞`) ⇒ `NaN` + `INVALID` (`Y` is complex in
//!   the reals there).
//! - `Y_n(+∞) = +0` for every order, by the decaying-envelope
//!   convention (ADR-0021/0023, the [`super::airy`] / J precedent):
//!   the true behaviour at `+∞` is a bounded decaying oscillation
//!   with no limit; the conservative total result is `+0`,
//!   `Status::OK`.
//! - `Y_n(NaN) = NaN`; `sNaN` raises `INVALID`.
//!
//! Negative order reduces before evaluation: `Y₋ₙ(x) = (−1)ⁿ Yₙ(x)`
//! (DLMF 10.4.1, the same parity as `J`), so the kernel evaluates
//! `Y_m(x)` for `m = |n| ≥ 0` and applies one parity sign. There is
//! no argument-parity reduction (the domain is `x > 0`; a negative
//! argument is `INVALID`, not folded to `|x|`).
//!
//! `Y` is the *dominant* solution of the Bessel three-term
//! recurrence, so the kernel computes `Y₀` and `Y₁` directly and
//! climbs to `Yₙ` by **upward** recurrence
//! `Y_{k+1}(x) = (2k/x)·Y_k(x) − Y_{k−1}(x)` (DLMF 10.6.1), which is
//! stable for the dominant solution. This is the opposite of `J`'s
//! Miller backward descent; [`super::bessel_j::bessel_j_miller`] is
//! not reused (there is no recessive solution to renormalise). The
//! base pair `Y₀`/`Y₁` uses two regimes on the binary exponent of
//! `x`, sharing [`super::bessel_j::bessel_j_threshold`] with `J`:
//!
//! - Below threshold: the DLMF 10.8.1 logarithmic series, with
//!   working precision boosted `≈ x·log₂e` for the alternating
//!   cancellation (the [`super::ci`] guard idiom). `Y` has no
//!   recessive-normalisation trick, so unlike `J` there is no cheap
//!   middle "moderate" regime; the log series carries everything
//!   below the asymptotic cut.
//! - At/above threshold: the DLMF 10.17.4 Hankel asymptotic, reusing
//!   the J `a_k(ν)` coefficients (DLMF 10.17.1, pinned in ADR-0023)
//!   with `Y`'s trig combination.
//!
//! ADR-0024 records the design and the coefficient provenance.

use super::ziv::ziv_round;
use super::ziv_calibration::BESSEL_Y_ERROR_GUARD;
use super::{euler_gamma_at, pi_at, pi_over_2_at};
use crate::big::{BigFloat, BuildError};
use crate::class::Class;
use crate::rounding::RoundingMode;
use crate::sign::Sign;
use crate::status::{auto_raise, Status};

#[cfg(feature = "fixed")]
use crate::fixed::FixedFloat;
#[cfg(feature = "fixed")]
use crate::mantissa::limbs_for;

impl BigFloat {
    /// `Y₀(self)` rounded under `mode` to `self.precision`.
    #[must_use]
    pub fn y0(&self, mode: RoundingMode) -> (Self, Status) {
        self.y0_round(self.precision, mode)
            .expect("self.precision >= 1 by invariant")
    }

    /// `Y₀(self)` with explicit result precision.
    ///
    /// Correctly rounded under every IEEE 754-2019 rounding mode
    /// via the shared `crate::math::ziv::ziv_round` driver (slice
    /// p1.32, ADR-0038).
    pub fn y0_round(
        &self,
        target_precision: u32,
        mode: RoundingMode,
    ) -> Result<(Self, Status), BuildError> {
        if target_precision == 0 {
            return Err(BuildError::PrecisionZero);
        }
        Ok(bessel_y_kernel(0, self, target_precision, mode))
    }

    /// `Y₁(self)` rounded under `mode` to `self.precision`.
    #[must_use]
    pub fn y1(&self, mode: RoundingMode) -> (Self, Status) {
        self.y1_round(self.precision, mode)
            .expect("self.precision >= 1 by invariant")
    }

    /// `Y₁(self)` with explicit result precision.
    ///
    /// Correctly rounded under every IEEE 754-2019 rounding mode
    /// via the shared `crate::math::ziv::ziv_round` driver (slice
    /// p1.32, ADR-0038).
    pub fn y1_round(
        &self,
        target_precision: u32,
        mode: RoundingMode,
    ) -> Result<(Self, Status), BuildError> {
        if target_precision == 0 {
            return Err(BuildError::PrecisionZero);
        }
        Ok(bessel_y_kernel(1, self, target_precision, mode))
    }

    /// `Yₙ(self)` for integer order `n` (any `i32`, including
    /// negative: `Y₋ₙ = (−1)ⁿ Yₙ`), rounded under `mode` to
    /// `self.precision`.
    #[must_use]
    pub fn yn(&self, n: i32, mode: RoundingMode) -> (Self, Status) {
        self.yn_round(n, self.precision, mode)
            .expect("self.precision >= 1 by invariant")
    }

    /// `Yₙ(self)` with explicit result precision.
    ///
    /// Correctly rounded under every IEEE 754-2019 rounding mode
    /// via the shared `crate::math::ziv::ziv_round` driver (slice
    /// p1.32, ADR-0038).
    pub fn yn_round(
        &self,
        n: i32,
        target_precision: u32,
        mode: RoundingMode,
    ) -> Result<(Self, Status), BuildError> {
        if target_precision == 0 {
            return Err(BuildError::PrecisionZero);
        }
        Ok(bessel_y_kernel(n, self, target_precision, mode))
    }
}

#[cfg(feature = "fixed")]
impl<const PREC: u32> FixedFloat<PREC>
where
    [(); limbs_for(PREC)]:,
{
    /// `Y₀(self)` for `FixedFloat`. Delegates to [`BigFloat::y0`].
    #[must_use]
    pub fn y0(&self, mode: RoundingMode) -> (Self, Status) {
        let (big, status) = self.to_big().y0(mode);
        (
            Self::try_from_big_exact(big).expect("precision matches"),
            status,
        )
    }

    /// `Y₁(self)` for `FixedFloat`. Delegates to [`BigFloat::y1`].
    #[must_use]
    pub fn y1(&self, mode: RoundingMode) -> (Self, Status) {
        let (big, status) = self.to_big().y1(mode);
        (
            Self::try_from_big_exact(big).expect("precision matches"),
            status,
        )
    }

    /// `Yₙ(self)` for `FixedFloat`. Delegates to [`BigFloat::yn`].
    #[must_use]
    pub fn yn(&self, n: i32, mode: RoundingMode) -> (Self, Status) {
        let (big, status) = self.to_big().yn(n, mode);
        (
            Self::try_from_big_exact(big).expect("precision matches"),
            status,
        )
    }
}

/// `Yₙ(x)` for integer order `n` and real `x`.
///
/// Special values are handled directly per the module-level domain
/// table; for a normal positive argument the order is reduced to
/// `Y_m(x)`, `m = |n|`, with one parity sign, then the regime
/// evaluator runs.
fn bessel_y_kernel(
    n: i32,
    x: &BigFloat,
    target_precision: u32,
    mode: RoundingMode,
) -> (BigFloat, Status) {
    let m = n.unsigned_abs();

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
            // Y_n(±0) is the POLE, not the negative axis: `log(±0) = −∞`
            // groups −0 with the pole (IEEE 754-2019 §9.2, C11 F.10.3.7;
            // `Y_0 ~ (2/π) ln x` inherits `ln`'s zero convention), for
            // BOTH zero signs — −0 was wrongly returning NaN + INVALID
            // (pf-k8ax, ADR-0123). The pole's sign follows the order
            // parity `Y_{−n} = (−1)^n Y_n` (DLMF 10.4.1) and POSIX yn:
            // `Y_n(0) = −∞` unless `n` is negative AND odd, then `+∞`
            // (the prior code was unconditionally −∞, a latent
            // negative-odd-order sign bug the audit surfaced).
            let pole_sign = if n < 0 && m % 2 == 1 {
                Sign::Positive
            } else {
                Sign::Negative
            };
            let inf =
                BigFloat::try_new_infinity(pole_sign, target_precision).expect("precision >= 1");
            auto_raise(Status::DIV_BY_ZERO);
            (inf, Status::DIV_BY_ZERO)
        }
        Class::Infinity { sign } => {
            if matches!(sign, Sign::Negative) {
                // −∞: complex (groups with x < 0).
                let nan = BigFloat::try_new_quiet_nan(Sign::Positive, target_precision, &[])
                    .expect("precision >= 1");
                auto_raise(Status::INVALID);
                return (nan, Status::INVALID);
            }
            // Y_n(+∞) = +0 (decaying-envelope convention,
            // ADR-0021/0023).
            let zero =
                BigFloat::try_new_zero(Sign::Positive, target_precision).expect("precision >= 1");
            (zero, Status::OK)
        }
        Class::Normal { sign, .. } => {
            if matches!(sign, Sign::Negative) {
                // x < 0: Y_n(−x) is complex in the reals.
                let nan = BigFloat::try_new_quiet_nan(Sign::Positive, target_precision, &[])
                    .expect("precision >= 1");
                auto_raise(Status::INVALID);
                return (nan, Status::INVALID);
            }

            // Resource budget (pf-ap01, ADR-0124): the upward recurrence
            // runs O(m) steps at ~4·m-bit working precision, so an
            // attacker-controlled order is an unbounded DoS. Refuse an
            // order past the feasible cap — representable but uncomputable
            // within budget — with NaN + INVALID.
            if m > super::bessel_j::MAX_BESSEL_ORDER {
                let nan = BigFloat::try_new_quiet_nan(Sign::Positive, target_precision, &[])
                    .expect("precision >= 1");
                auto_raise(Status::INVALID);
                return (nan, Status::INVALID);
            }

            // Y₋ₙ(x) = (−1)ⁿ Yₙ(x) (DLMF 10.4.1): order parity is
            // binary and pinned outside the Ziv envelope.
            let negate = (m % 2 == 1) && (n < 0);

            // Pin the base-pair regime decision from target_precision
            // so it does not flip across Ziv retries (slice p1.4 erf
            // precedent). Sub-slice 2b.2.a: Y uses the tightened
            // `bessel_yik_threshold` (not the conservative
            // `bessel_j_threshold` J keeps) because the log-series
            // base path is dramatically more expensive than the
            // asymptotic at the dispatch boundary — the bench
            // measured 89%–96% reductions at the flip cells
            // (`Y0_p256_x256` 57.97ms → 6.15ms; `Y0_p1024_x1024`
            // 2.05s → 85.08ms). See `bessel_yik_threshold`'s doc for
            // the four-part risk enumeration.
            let e_x = match &x.class {
                Class::Normal { exponent, .. } => *exponent,
                _ => 0,
            };
            let use_asymptotic = e_x >= super::bessel_j::bessel_yik_threshold(target_precision);
            // Pin the near-zero fallback for the direct Y₀/Y₁ asymptotic:
            // an input near a Y₀/Y₁ zero falls below the divergent series'
            // truncation floor; resolve with the convergent log series,
            // which cancellation_boosted drives to any depth (pf-1vzg,
            // ADR-0125). n ≥ 2 climbs by recurrence from the base pair and
            // is out of this hand-off.
            let maclaurin_fallback = m <= 1
                && use_asymptotic
                && !super::ziv::asymptotic_reliable(target_precision, BESSEL_Y_ERROR_GUARD, |w| {
                    bessel_y_asymptotic(m, x, w)
                });

            let (result, status) = ziv_round(
                |w| {
                    let value = if maclaurin_fallback {
                        super::ziv::cancellation_boosted(w, |ww| bessel_y_series(m, x, ww))
                    } else {
                        bessel_y_eval_normal_at_w(m, x, w, use_asymptotic)
                    };
                    if negate {
                        value.negated()
                    } else {
                        value
                    }
                },
                target_precision,
                mode,
                BESSEL_Y_ERROR_GUARD,
            );
            // Yₙ(x) for finite-normal x > 0 is transcendental at nonzero
            // algebraic argument (ADR-0064; Y is a second-kind solution
            // with a logarithmic term at the origin, so the guarantee
            // comes from the transcendence theory of the Bessel
            // differential system rather than the E-function theorem, and
            // no named open problem obstructs it). A finite-normal result
            // is therefore INEXACT even where the working-precision
            // evaluation lands on a grid value. The +0 pole and the ±∞
            // limits are dispatched above.
            let status = super::force_transcendental_inexact(&result, status);
            auto_raise(status);
            (result, status)
        }
    }
}

/// `Y_m(x)` for `m ≥ 0`, normal `x > 0`: the base pair plus upward
/// recurrence. Returns the unrounded working-precision value;
/// [`bessel_y_kernel`] does the single final round.
///
/// The base pair `Y₀`/`Y₁` go through the two-regime dispatch
/// ([`bessel_y01`]); `Yₙ (n ≥ 2)` climbs from that pair by the
/// **upward** recurrence (DLMF 10.6.1)
///
/// ```text
/// Y_{k+1}(x) = (2k/x)·Y_k(x) − Y_{k−1}(x)
/// ```
///
/// which is numerically stable because `Y` is the *dominant*
/// solution of the three-term relation (forward recurrence grows the
/// dominant solution; the `−Y_{k−1}` subtraction is never
/// catastrophic since `(2k/x)·Y_k` dominates it). This is the
/// opposite of `J`'s Miller backward descent; there is no recessive
/// solution to renormalise, so no sum rule and no
/// [`super::bessel_j::bessel_j_miller`] reuse. The climb is a
/// rolling pair, `O(m)` time and `O(1)` space, no `Vec`.
fn bessel_y_eval_normal_at_w(
    m: u32,
    x: &BigFloat,
    target_precision: u32,
    use_asymptotic_base: bool,
) -> BigFloat {
    if m <= 1 {
        return bessel_y01(m, x, target_precision, use_asymptotic_base);
    }

    // pf-1axr: the pre-fix `extra = mag·23/16` budget (matching
    // `bessel_y_series`'s alternating-series-cancellation boost)
    // was wildly over-provisioned for the upward recurrence. The
    // recurrence step `Y_{k+1} = (2k/x)·Y_k − Y_{k−1}` cancels only
    // when `(2k/x)·Y_k` and `Y_{k−1}` have the same sign and similar
    // magnitude; for `x ≫ k` the `(2k/x)` factor is small and
    // `|Y_{k+1}| ≈ |Y_{k−1}|` with no cancellation. For `x ≲ k` the
    // amplification per step is bounded by `1 + x/(2k)` (≤ 1 bit per
    // step). The original `|x|`-scaled budget pushed working past
    // the 4096-bit `2/π` reduction table's range for `x ≳ 2^11`,
    // breaking `bessel_y_asymptotic`'s internal `cos`/`sin` calls
    // (`Y2(2049)` at `p=53` NaN). Use a constant boost matching the
    // worst-case `m`-step error amplification: `32` baseline plus
    // `4·m` for the per-step bit-loss bound.
    let recurrence_extra = 32u32.saturating_add(m.saturating_mul(4));
    let working = target_precision
        .saturating_add(64)
        .saturating_add(recurrence_extra);

    let xw = x
        .round_to_precision(working, RoundingMode::NearestEven)
        .expect("precision >= 1")
        .0;
    let (inv_x, _) = ci(1, working).div(&xw, RoundingMode::NearestEven);

    // Seeds Y₀, Y₁ at the recurrence working precision, using the
    // base-pair regime flag pinned by the caller from target_precision.
    let mut y_prev = bessel_y01(0, x, working, use_asymptotic_base)
        .round_to_precision(working, RoundingMode::NearestEven)
        .expect("precision >= 1")
        .0;
    let mut y_cur = bessel_y01(1, x, working, use_asymptotic_base)
        .round_to_precision(working, RoundingMode::NearestEven)
        .expect("precision >= 1")
        .0;

    // Y_{k+1} = (2k/x)·Y_k − Y_{k−1}, k = 1 … m−1.
    for k in 1..i64::from(m) {
        let two_k = ci(2 * k, working);
        let (c1, _) = two_k.mul(&inv_x, RoundingMode::NearestEven);
        let (c2, _) = c1.mul(&y_cur, RoundingMode::NearestEven);
        let (y_next, _) = c2.sub(&y_prev, RoundingMode::NearestEven);
        y_prev = y_cur;
        y_cur = y_next;
    }
    y_cur
}

/// `Y₀(x)` (`which = 0`) or `Y₁(x)` (`which = 1`) for `x > 0`, the
/// two-regime base dispatch on the binary exponent of `x`. Shares
/// [`super::bessel_j::bessel_j_threshold`] with `J`: the DLMF 10.17.4
/// `Y` asymptotic has the same `e^{−2|x|}` error order as the DLMF
/// 10.17.3 `J` asymptotic, so the same conservative cut applies.
/// `Y` has no recessive-normalisation analog, so unlike `J` there is
/// no cheap middle regime; the log series carries everything below
/// the asymptotic cut. Returns the unrounded working-precision value.
fn bessel_y01(which: u32, x: &BigFloat, target_precision: u32, use_asymptotic: bool) -> BigFloat {
    if use_asymptotic {
        // Reliable asymptotic (above the floor), boosted so a near-zero
        // cancellation past the Ziv guard cap still resolves (pf-1vzg).
        super::ziv::cancellation_boosted(target_precision, |w| {
            let (v, op_scale, _floor) = bessel_y_asymptotic(which, x, w);
            (v, op_scale)
        })
    } else {
        // The DLMF 10.8.1 log series' head/log/tail composition cancels
        // (Y_n small toward the boundary); charge the realised depth via
        // bessel_y_series's exposed op_scale instead of the old fixed cap,
        // which undershot once target exceeded it (pf-6naq, ADR-0126). The
        // asymptotic branch already boosts for its near-zero case (pf-1vzg).
        super::ziv::cancellation_boosted(target_precision, |w| bessel_y_series(which, x, w))
    }
}

/// Binary exponent of `v`, or `i64::MIN`/`i64::MAX` for zero /
/// non-finite (the [`super::si`] / [`super::bessel_j`] `magnitude`
/// idiom; `bessel_j`'s copy is private, so `Y` carries its own).
fn magnitude(v: &BigFloat) -> i64 {
    match &v.class {
        Class::Normal { exponent, .. } => *exponent,
        Class::Zero { .. } => i64::MIN,
        _ => i64::MAX,
    }
}

/// `true` once `term` has fallen `working + 8` bits below the running
/// `sum`, so further terms cannot perturb the rounded result (the
/// [`super::si`] / [`super::bessel_j`] `negligible` idiom).
fn negligible(term: &BigFloat, sum: &BigFloat, working: u32) -> bool {
    match &term.class {
        Class::Zero { .. } => true,
        Class::Normal { exponent, .. } => *exponent < magnitude(sum) - i64::from(working) - 8,
        _ => false,
    }
}

/// Integer `v` as a `BigFloat` at precision `p` (exact for the small
/// integers this kernel forms).
fn ci(v: i64, p: u32) -> BigFloat {
    BigFloat::try_from_i64_exact(v, p).expect("precision >= 1")
}

/// `Y_n(x)` for integer `n ≥ 0`, `x > 0`, via the DLMF 10.8.1
/// logarithmic series
///
/// ```text
/// Y_n(x) = −(x/2)^{−n}/π · Σ_{k=0}^{n−1} (n−k−1)!/k! (x²/4)^k
///        + (2/π) ln(x/2) · J_n(x)
///        − (x/2)^n/π · Σ_{k≥0} (ψ(k+1)+ψ(n+k+1)) (−x²/4)^k/(k!(n+k)!)
/// ```
///
/// The finite first ("head") sum is empty for `n = 0`. The digamma
/// terms reduce to elementary running sums: `ψ(k+1) = −γ + H_k` and
/// `ψ(n+k+1) = −γ + H_{n+k}` (`H_j = Σ_{i=1}^{j} 1/i`, `H_0 = 0`), so
/// `ψ(k+1)+ψ(n+k+1) = −2γ + H_k + H_{n+k}` — **no digamma kernel
/// needed**, only harmonic partial sums and the in-tree
/// Euler–Mascheroni `γ` (slice 6m0, [`super::euler_gamma_at`]). The
/// `J_n(x)` piece is **6o's `bessel_j` kernel** (`BigFloat::jn_round`
/// at working precision).
///
/// Cross-check (the `derive, don't recall` reflex): for `n = 0` the
/// finite sum is empty and `−2γ + 2H_k` substituted into the tail,
/// with `Σ (−x²/4)^k/(k!)² = J_0(x)`, folds the `−2γ J_0` into a
/// `+γ` and (since `H_0 = 0`) kills the `k = 0` tail term, exactly
/// reproducing DLMF 10.8.2
/// `Y_0 = (2/π)(ln(x/2)+γ)J_0 + (2/π)[(x²/4)/(1!)² − …]`. Hand-check
/// `k = 1`: `−(1/π)·2H_1·(−x²/4)/(1!)² = +(2/π)(x²/4)`, the first
/// 10.8.2 bracket term. The two DLMF forms agree, pinning the
/// reduction.
///
/// Returns `(value, op_scale)`. A fixed head-start covers the moderate
/// `≈ x·log₂ e` cancellation for direct callers; production callers wrap
/// this in [`super::ziv::cancellation_boosted`], which raises the working
/// precision past the cap for the deep cases the cap alone undershot (near
/// a `Y` zero, or a large moderate `x` at high target — pf-6naq, ADR-0126).
fn bessel_y_series(n: u32, x: &BigFloat, target_precision: u32) -> (BigFloat, i64) {
    // A self-sufficiency head-start so DIRECT callers (the cross-check
    // tests, which call this outside `cancellation_boosted`) resolve the
    // moderate `≈ |x|·log₂ e` cancellation without a boost. Production
    // callers wrap this in `cancellation_boosted` (via the returned
    // op_scale), which grows the working precision past the cap for the
    // deep cases the cap alone undershot (pf-6naq, ADR-0126) — so the cap
    // bounds only the head-start, not the ceiling.
    let e_x = match &x.class {
        Class::Normal { exponent, .. } => *exponent,
        _ => 0,
    };
    let extra = if e_x <= 0 {
        64
    } else {
        let shift = (e_x + 1).min(20) as u32;
        let mag: u64 = 1u64 << shift;
        (mag.saturating_mul(23) / 16).min(4096) as u32
    };
    let working = target_precision
        .saturating_add(64)
        .saturating_add(extra)
        .min(target_precision.saturating_add(4096));

    let nn = i64::from(n);
    let xw = x
        .round_to_precision(working, RoundingMode::NearestEven)
        .expect("precision >= 1")
        .0;
    let two = ci(2, working);
    let (half, _) = xw.div(&two, RoundingMode::NearestEven); // x/2
    let (x_sq, _) = xw.mul(&xw, RoundingMode::NearestEven);
    let (qxs, _) = x_sq.div(&ci(4, working), RoundingMode::NearestEven); // x²/4
    let neg_qxs = qxs.negated(); // −x²/4
    let pi = pi_at(working);
    let gamma = euler_gamma_at(working);
    let one = ci(1, working);

    // (x/2)^n by repeated multiply (n is small: the base orders are
    // 0 and 1; an in-module recurrence-from-a-base-term, not
    // pow.rs::pow_int).
    let mut half_pow_n = one.clone();
    for _ in 0..n {
        half_pow_n = half_pow_n.mul(&half, RoundingMode::NearestEven).0;
    }

    // log term: (2/π) ln(x/2) · J_n(x).
    let (ln_half_x, _) = half.ln(RoundingMode::NearestEven);
    let (jn, _) = xw
        .jn_round(n as i32, working, RoundingMode::NearestEven)
        .expect("precision >= 1");
    let (two_over_pi, _) = two.div(&pi, RoundingMode::NearestEven);
    let (lt0, _) = two_over_pi.mul(&ln_half_x, RoundingMode::NearestEven);
    let (log_term, _) = lt0.mul(&jn, RoundingMode::NearestEven);

    // Finite head −(x/2)^{−n}/π · Σ_{k=0}^{n−1} (n−k−1)!/k! (x²/4)^k.
    // h_0 = (n−1)!; h_k = h_{k−1}·(x²/4)/(k·(n−k)). Empty for n = 0.
    let head_term = if n == 0 {
        BigFloat::try_new_zero(Sign::Positive, working).expect("precision >= 1")
    } else {
        let mut fact = one.clone(); // (n−1)!
        for j in 1..nn {
            fact = fact.mul(&ci(j, working), RoundingMode::NearestEven).0;
        }
        let mut term = fact;
        let mut hd = term.clone();
        for k in 1..nn {
            term = term.mul(&qxs, RoundingMode::NearestEven).0;
            term = term
                .div(&ci(k * (nn - k), working), RoundingMode::NearestEven)
                .0;
            hd = hd.add(&term, RoundingMode::NearestEven).0;
        }
        let (inv_hp, _) = one.div(&half_pow_n, RoundingMode::NearestEven); // (x/2)^{−n}
        let (a, _) = hd.mul(&inv_hp, RoundingMode::NearestEven);
        let (b, _) = a.div(&pi, RoundingMode::NearestEven);
        b.negated()
    };

    // Tail −(x/2)^n/π · Σ_{k≥0} (−2γ + H_k + H_{n+k}) c_k,
    // c_0 = 1/n!, c_k = c_{k−1}·(−x²/4)/(k·(n+k)).
    let mut n_fact = one.clone();
    for j in 2..=nn {
        n_fact = n_fact.mul(&ci(j, working), RoundingMode::NearestEven).0;
    }
    let (mut c, _) = one.div(&n_fact, RoundingMode::NearestEven); // c_0 = 1/n!
    let mut h_k = BigFloat::try_new_zero(Sign::Positive, working).expect("precision >= 1"); // H_0
    let mut h_nk = {
        let mut s = BigFloat::try_new_zero(Sign::Positive, working).expect("precision >= 1");
        for i in 1..=nn {
            let (r, _) = one.div(&ci(i, working), RoundingMode::NearestEven);
            s = s.add(&r, RoundingMode::NearestEven).0;
        }
        s // H_n
    };
    let (two_gamma, _) = gamma.mul(&two, RoundingMode::NearestEven);
    let coef0 = {
        let (s, _) = h_k.add(&h_nk, RoundingMode::NearestEven);
        s.sub(&two_gamma, RoundingMode::NearestEven).0
    };
    let mut tail = c.mul(&coef0, RoundingMode::NearestEven).0;
    // Largest tail partial term: the convergent log series' terms peak at
    // ≈ 2^{x·log₂ e} before cancelling to Y_n; that peak (lifted to the
    // result scale) is the operand scale that charges the deep near-zero
    // cancellation to cancellation_boosted (pf-1vzg, ADR-0125).
    let mut max_tail_term = magnitude(&tail);
    let max_iter: i64 = 1 << 22;
    for k in 1..=max_iter {
        c = c.mul(&neg_qxs, RoundingMode::NearestEven).0;
        c = c
            .div(&ci(k * (k + nn), working), RoundingMode::NearestEven)
            .0;
        let (rk, _) = one.div(&ci(k, working), RoundingMode::NearestEven);
        h_k = h_k.add(&rk, RoundingMode::NearestEven).0;
        let (rnk, _) = one.div(&ci(k + nn, working), RoundingMode::NearestEven);
        h_nk = h_nk.add(&rnk, RoundingMode::NearestEven).0;
        let coef = {
            let (s, _) = h_k.add(&h_nk, RoundingMode::NearestEven);
            s.sub(&two_gamma, RoundingMode::NearestEven).0
        };
        let t = c.mul(&coef, RoundingMode::NearestEven).0;
        max_tail_term = max_tail_term.max(magnitude(&t));
        tail = tail.add(&t, RoundingMode::NearestEven).0;
        if negligible(&t, &tail, working) {
            break;
        }
    }
    let (tp0, _) = half_pow_n.mul(&tail, RoundingMode::NearestEven);
    let (tp1, _) = tp0.div(&pi, RoundingMode::NearestEven);
    let tail_term = tp1.negated();

    let (s1, _) = head_term.add(&log_term, RoundingMode::NearestEven);
    let (y, _) = s1.add(&tail_term, RoundingMode::NearestEven);
    // op_scale: the tail's internal peak lifted to the result scale by
    // the (x/2)^n/π prefactor, plus the head/log pieces (already at the
    // result scale). The largest is what cancellation_boosted charges.
    let op_scale = max_tail_term
        .saturating_add(magnitude(&half_pow_n))
        .saturating_sub(magnitude(&pi))
        .max(magnitude(&head_term))
        .max(magnitude(&log_term));
    (y, op_scale)
}

/// `Y_n(x)` for `n ≥ 0`, large `x > 0`, via the DLMF 10.17.4
/// Hankel-form asymptotic
///
/// ```text
/// Y_n(x) ∼ √(2/(πx))·[sin ω · Σ_{k≥0} (−1)^k a_{2k}(n)/x^{2k}
///                    + cos ω · Σ_{k≥0} (−1)^k a_{2k+1}(n)/x^{2k+1}]
/// ω = x − n·π/2 − π/4
/// ```
///
/// summed to its smallest term (the [`super::bessel_j`] /
/// [`super::si`] optimal-truncation idiom). The coefficients
/// `a_k(n)` and the phase `ω` are **identical to `J`'s** (DLMF
/// 10.17.1/10.17.2): `a_0(n)=1`,
/// `a_k(n) = a_{k−1}(n)·(4n²−(2k−1)²)/(8k)`, already derived and
/// cross-pinned in ADR-0023 (the closed-form ratio and the
/// Pochhammer form, agreeing at `k=1,2`); 6p reuses that recurrence
/// verbatim rather than re-deriving it. `Y` differs from `J` only
/// in the trig combination: DLMF 10.17.4 is `sinω·ΣP + cosω·ΣQ`
/// where `J`'s 10.17.3 is `cosω·ΣP − sinω·ΣQ`
/// (`ΣP = Σ(−1)^k a_{2k}/x^{2k}`, `ΣQ = Σ(−1)^k a_{2k+1}/x^{2k+1}`).
/// Folding the explicit `(−1)^k` into the trig assignment, the
/// per-`j` factor on `a_j(n)/x^j` is the period-4 cycle
/// `[+sinω, +cosω, −sinω, −cosω]` for `j ≡ 0,1,2,3 (mod 4)` (vs
/// `J`'s `[+cosω, −sinω, −cosω, +sinω]`). Hand-check: `j=0` →
/// `+sinω·a_0` (`k=0` of `ΣP`); `j=1` → `+cosω·a_1` (`k=0` of
/// `ΣQ`); `j=2` → `−sinω·a_2` (`k=1` of `ΣP`); `j=3` → `−cosω·a_3`
/// (`k=1` of `ΣQ`). Returns `(value, op_scale, floor_exp)`: near a `Y_n`
/// zero the bracket cancels below its largest term (`op_scale`, lifted
/// through the `√(2/πx)` prefactor), but this is a DIVERGENT asymptotic
/// with an irreducible truncation floor (`floor_exp`). Below the floor
/// the value is uncomputable at any working precision, so the caller
/// detects it and hands off to the convergent log series (pf-1vzg,
/// ADR-0125).
fn bessel_y_asymptotic(n: u32, x: &BigFloat, target_precision: u32) -> (BigFloat, i64, i64) {
    // ω = x − n·(π/2) − π/4 has magnitude ~2^e_x; cos/sin(ω) depend on
    // ω mod 2π, so the working precision must exceed x's integer width
    // to keep `target_precision` accurate FRACTIONAL bits. The prior
    // target+64 hard-capped at target+512 ignored |x|, so the phase was
    // garbage (often wrong sign) for large x and the Ziv loop could
    // never grow past the cap (review 2026-05-29, root cause 1).
    let e_x = match &x.class {
        Class::Normal { exponent, .. } => *exponent,
        _ => 0,
    };
    let x_int_bits = u32::try_from(e_x.max(0)).unwrap_or(u32::MAX);
    let working = target_precision
        .saturating_add(64)
        .saturating_add(x_int_bits);
    let xw = x
        .round_to_precision(working, RoundingMode::NearestEven)
        .expect("precision >= 1")
        .0;

    // ω = x − n·(π/2) − π/4.
    let half_pi = pi_over_2_at(working);
    let two = ci(2, working);
    let (quarter_pi, _) = half_pi.div(&two, RoundingMode::NearestEven);
    let n_big = ci(i64::from(n), working);
    let (n_half_pi, _) = n_big.mul(&half_pi, RoundingMode::NearestEven);
    let (w0, _) = xw.sub(&n_half_pi, RoundingMode::NearestEven);
    let (omega, _) = w0.sub(&quarter_pi, RoundingMode::NearestEven);
    let (cw, _) = omega.cos(RoundingMode::NearestEven);
    let (sw, _) = omega.sin(RoundingMode::NearestEven);

    // prefactor √(2/(πx)).
    let pi = pi_at(working);
    let (pi_x, _) = pi.mul(&xw, RoundingMode::NearestEven);
    let (ratio, _) = two.div(&pi_x, RoundingMode::NearestEven);
    let (prefac, _) = ratio.sqrt(RoundingMode::NearestEven);

    let (inv_x, _) = ci(1, working).div(&xw, RoundingMode::NearestEven);
    let four_n2: i64 = 4 * i64::from(n) * i64::from(n);

    // j = 0: g_0 = a_0/x^0 = 1, factor +sin ω (Y's j≡0 term).
    let mut g = ci(1, working);
    let (mut bracket, _) = sw.mul(&g, RoundingMode::NearestEven);
    // Largest bracket term seen (the j = 0 term +sin ω·g₀ to start); its
    // exponent lifted by the prefactor is the operand scale that charges
    // the near-zero cancellation to cancellation_boosted.
    let mut max_term_exp = magnitude(&bracket);
    let mut prev_mag = magnitude(&g);
    let max_iter: i64 = 1 << 22;
    for j in 1..=max_iter {
        // a_j/a_{j−1} = (4n²−(2j−1)²)/(8j); g_j = g_{j−1}·that·(1/x).
        let odd = 2 * j - 1;
        let num = four_n2 - odd * odd;
        let num_b = ci(num, working);
        let den = ci(8 * j, working);
        let (t1, _) = g.mul(&num_b, RoundingMode::NearestEven);
        let (t2, _) = t1.div(&den, RoundingMode::NearestEven);
        let (cand, _) = t2.mul(&inv_x, RoundingMode::NearestEven);
        let mag = magnitude(&cand);
        if mag > prev_mag {
            break; // smallest term passed: optimal truncation.
        }
        prev_mag = mag;
        g = cand;
        // Y's period-4 trig/sign cycle on a_j/x^j (DLMF 10.17.4).
        let contribution = match j % 4 {
            0 => sw.mul(&g, RoundingMode::NearestEven).0,
            1 => cw.mul(&g, RoundingMode::NearestEven).0,
            2 => sw.mul(&g, RoundingMode::NearestEven).0.negated(),
            _ => cw.mul(&g, RoundingMode::NearestEven).0.negated(),
        };
        max_term_exp = max_term_exp.max(magnitude(&contribution));
        let (b, _) = bracket.add(&contribution, RoundingMode::NearestEven);
        bracket = b;
    }
    let (result, _) = prefac.mul(&bracket, RoundingMode::NearestEven);
    let op_scale = max_term_exp.saturating_add(magnitude(&prefac));
    // Truncation floor: the smallest retained coefficient `prev_mag`
    // lifted by the prefactor. Below it the divergent asymptotic cannot
    // reach the true value, so the caller routes a near-zero result to
    // the convergent log series (pf-1vzg, ADR-0125).
    let floor_exp = prev_mag.saturating_add(magnitude(&prefac));
    (result, op_scale, floor_exp)
}

#[cfg(test)]
mod tests {
    use super::*;
    use core::cmp::Ordering;

    /// pf-ap01 (ADR-0124): an integer order past `MAX_BESSEL_ORDER` is a
    /// resource blowup (O(|n|) steps at ~4·|n|-bit precision), so
    /// `yn`/`jn`/`in` refuse it with NaN + INVALID (a budget), fast; an
    /// in-cap order and
    /// the exact `x = 0` cases still work.
    #[test]
    fn bessel_order_resource_budget_refuses_exotic_orders() {
        let huge = (super::super::bessel_j::MAX_BESSEL_ORDER + 1) as i32;
        let x = BigFloat::try_from_i64_exact(1, 53).unwrap();
        let (ry, sy) = x.yn(huge, RoundingMode::NearestEven);
        assert!(ry.is_nan() && sy.invalid(), "yn(huge) refused");
        let (rj, sj) = x.jn(huge, RoundingMode::NearestEven);
        assert!(rj.is_nan() && sj.invalid(), "jn(huge) refused");
        let (ri, si) = x.in_(huge, RoundingMode::NearestEven);
        assert!(ri.is_nan() && si.invalid(), "in(huge) refused");
        // In-cap order still computes.
        let (rc, _) = x.jn(10, RoundingMode::NearestEven);
        assert!(!rc.is_nan(), "jn(10) still computes");
        // x = 0 exact cases are not capped.
        let z = BigFloat::try_new_zero(Sign::Positive, 53).unwrap();
        let (rz, _) = z.jn(huge, RoundingMode::NearestEven);
        assert!(rz.is_zero(), "J_n(0) = 0 exact, not capped");
    }

    /// pf-k8ax (ADR-0123): `Y_n(±0)` is the pole for BOTH zero signs —
    /// −0 wrongly returned NaN + INVALID where −∞ + `DIV_BY_ZERO` is due
    /// (log(±0) = −∞ groups −0 with the pole). Order 0 → −∞.
    #[test]
    fn y0_negative_zero_is_the_pole_not_nan() {
        let nz = BigFloat::try_new_zero(Sign::Negative, 53).unwrap();
        let (r, st) = nz.y0_round(53, RoundingMode::NearestEven).unwrap();
        assert!(
            r.is_infinite() && r.is_sign_negative(),
            "y0(-0) = -inf (the pole), got {r}"
        );
        assert!(!st.invalid(), "y0(-0) is a pole (DIV_BY_ZERO), not INVALID");
    }

    /// `|v − expected| ≤ 2^(−bits)·|expected|` (the `bessel_j`/`erf`
    /// test helper). Reference decimals: `mpmath` `bessely(n, x)` at
    /// `mp.dps = 330`
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
            Some(Ordering::Less | Ordering::Equal)
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

    // mpmath bessely(n, x), dps = 330 (truncated to fit p = 113).
    const Y0_HALF: &str = "-0.44451873350670655714839847506833191037356512440151102041489117938823968793141728600040643311534111986534166279374081164055681338354759796819979290414603607710322924252297191321517536861397042466137749158319359986535079";
    const Y1_HALF: &str = "-1.4714723926702430691885846353232974532410880554357483229559223834066940909523570678868634705980974507427173395904931472603822957011919081229168205907226988262044905682586324008798759255110477030702664113187480772586178";
    const Y2_HALF: &str = "-5.4413708371742657196059400662248579025907870973414822714087983542385366758780109855470474492770486831055276955682317774009723694212200345234674894587447592277147330305115576903043283334302203876196881536917987091691206";
    const Y0_52: &str = "0.49807035961523188782747235036208980611506253265681530429974629407406321557566995812355070589537398192385672943001934612775337251591109454305007402401299973085144637759449327183580472195082322152061855168008649846205288";
    const Y1_52: &str = "0.14591813796678579887875994053587757127608019654670099984510336876847982786986827339694694495998904089288127698548695420739592301285341299984688991547700642348292700296345898227891129232942949834606441993069787616256066";
    const Y2_52: &str = "-0.38133584924180324872446439793338774909419837541945450442366359905927935327977533940599314992738274920955170784162978276183663410562836414317256209163139459206510477522372608601267568808727962284376701573552819753200435";
    const Y0_7: &str = "-0.025949743967209264884284963135722970186930836172718446607766370553193513857582489525922712522537272092128691592370516582921565613566622223247970896269009905783095511803775145797385878119661716988114304058052984665018584";
    const Y1_7: &str = "-0.30266723702418487006076816955839496834131089203393953780158416495055339422290387639178721348655609063443268492972537685185103628707041049719364672402943351702708751523804853131680988618472817659828972258775307059424644";
    const Y2_7: &str = "-0.060526609468272126561648799595247020767729418694121421335543390861250313063247189443159348473621610946280646958979591089035873325596352204521642453453685384796072349692810148864559803647403476325682759538447892647623256";
    const Y0_M3: &str = "-4.4714166113759232689802886934264955747044811557836578355293457089934050985223905404716180228246204738260124225710113837270898795280105746607399842023055143947979883646042194261470348839760198332670650930447884411058852";
    const Y1_M3: &str = "-636.62216723113942807437320601957569637320004286220162988363010289166468020567331606275382078239146840454146131313959710498920403749964474385892562171738689411877766793201649541296641816020550741373175777008581457495406030";

    /// DLMF 10.8.1 log-series path vs mpmath at `p = 113`: orders
    /// `0,1,2` over small (`0.5`), moderate (`2.5`, `7`) `x`, plus
    /// the pole approach (`1e-3`, where the `(2/x)^n` head dominates
    /// for `n ≥ 1`). Large enough orders/arguments that the second
    /// `c_k` term materially contributes (the `derive, don't recall`
    /// reflex: a low-`x`-only check would hide a coefficient error).
    #[test]
    fn series_matches_mpmath() {
        let p = 113;
        let x = at(1, 2, p); // 0.5
        for (n, want) in [(0i32, Y0_HALF), (1, Y1_HALF), (2, Y2_HALF)] {
            let (r, _) = x.yn(n, RoundingMode::NearestEven);
            assert!(close_at(&r, &py(want, p), p - 8), "Y_{n}(0.5)");
        }
        let x = at(5, 2, p); // 2.5
        for (n, want) in [(0i32, Y0_52), (1, Y1_52), (2, Y2_52)] {
            let (r, _) = x.yn(n, RoundingMode::NearestEven);
            assert!(close_at(&r, &py(want, p), p - 10), "Y_{n}(2.5)");
        }
        let x = at(7, 1, p);
        for (n, want) in [(0i32, Y0_7), (1, Y1_7), (2, Y2_7)] {
            let (r, _) = x.yn(n, RoundingMode::NearestEven);
            assert!(close_at(&r, &py(want, p), p - 12), "Y_{n}(7)");
        }
        let x = at(1, 1000, p); // 1e-3, pole approach
        for (n, want) in [(0i32, Y0_M3), (1, Y1_M3)] {
            let (r, _) = x.yn(n, RoundingMode::NearestEven);
            assert!(close_at(&r, &py(want, p), p - 12), "Y_{n}(1e-3)");
        }
    }

    // mpmath bessely at large |x| (dps = 90).
    const Y0_200: &str =
        "-0.054265775249817910693500115141763046843260653011054063148190169321587553804947092";
    const Y1_200: &str =
        "0.01530182458038998921966780774391601230422488536804796783556856489053218931466921";
    const Y2_200: &str =
        "0.054418793495621810585696793219202206966302901864734542826545854970492875698093784";
    const Y3_200: &str =
        "-0.014213448710477553007953871879531968164898827330753276979037647791122331800707334";
    const Y5_200: &str =
        "0.012019640832200107520916455504508441524448832663366991333600968414953269974622734";
    const Y0_1000: &str =
        "0.0047159179776228133997732614656652550098590048968019771852800026829151397602999922";
    const Y1_1000: &str =
        "-0.024784331292351778914862356097141290938631854864870528758349019940178780233621799";
    const Y2_1000: &str =
        "-0.0047654866402075169576029861778595375917362686065317182427967007227954973207672358";

    /// DLMF 10.17.4 Hankel-asymptotic path, called directly so the
    /// reused `a_k(n)` recurrence and `Y`'s period-4 trig cycle are
    /// pinned independently of the regime dispatch (added in 6p.4):
    /// `Y_n` at `x = 200` (`p = 53`) and `x = 1000` (`p = 113`) vs
    /// mpmath. Large enough that the second and later `a_k` terms
    /// materially contribute, so a recalled coefficient or a wrong
    /// trig cycle fails here.
    #[test]
    fn asymptotic_matches_mpmath() {
        let p = 53;
        let x = at(200, 1, p);
        for (n, want) in [
            (0u32, Y0_200),
            (1, Y1_200),
            (2, Y2_200),
            (3, Y3_200),
            (5, Y5_200),
        ] {
            let (r, _, _) = bessel_y_asymptotic(n, &x, p);
            assert!(close_at(&r, &py(want, p), p - 8), "Y_{n}(200)");
        }
        let p = 113;
        let x = at(1000, 1, p);
        for (n, want) in [(0u32, Y0_1000), (1, Y1_1000), (2, Y2_1000)] {
            let (r, _, _) = bessel_y_asymptotic(n, &x, p);
            assert!(close_at(&r, &py(want, p), p - 12), "Y_{n}(1000)");
        }

        // Public path at large x must route through bessel_y01 into
        // the asymptotic (pins the bessel_j_threshold dispatch).
        let x = at(200, 1, 53);
        let (r0, _) = x.y0(RoundingMode::NearestEven);
        assert!(close_at(&r0, &py(Y0_200, 53), 45), "y0(200) via dispatch");
        let (r1, _) = x.y1(RoundingMode::NearestEven);
        assert!(close_at(&r1, &py(Y1_200, 53), 45), "y1(200) via dispatch");
    }

    /// Upward-recurrence cross-tie `Y_{n−1}(x) + Y_{n+1}(x) =
    /// (2n/x)·Y_n(x)` (DLMF 10.6.1), binding three orders that all
    /// climb from the same `(Y₀, Y₁)` base pair.
    #[test]
    fn recurrence_spot_check() {
        let p = 160;
        let x = at(5, 2, p); // 2.5
        let (y2, _) = x.yn(2, RoundingMode::NearestEven);
        let (y3, _) = x.yn(3, RoundingMode::NearestEven);
        let (y4, _) = x.yn(4, RoundingMode::NearestEven);
        let (lhs, _) = y2.add(&y4, RoundingMode::NearestEven);
        let six = BigFloat::try_from_i64_exact(6, p).unwrap();
        let (r1, _) = six.mul(&y3, RoundingMode::NearestEven);
        let (rhs, _) = r1.div(&x, RoundingMode::NearestEven);
        assert!(close_at(&lhs, &rhs, p - 10), "Y_2+Y_4 = (6/x)Y_3");
    }

    /// Recurrence agrees with the direct DLMF 10.8.1 series for
    /// `n ≥ 2` (independent paths), at moderate `x`.
    #[test]
    fn recurrence_matches_series() {
        let p = 160;
        let x = at(7, 2, p); // 3.5
        for n in [2u32, 3, 5] {
            let recur = bessel_y_eval_normal_at_w(n, &x, p, false);
            let (series, _) = bessel_y_series(n, &x, p);
            assert!(close_at(&recur, &series, p - 12), "Y_{n}(3.5) recur=series");
        }
    }

    /// Negative-order parity `Y₋ₙ(x) = (−1)ⁿ Yₙ(x)` (DLMF 10.4.1):
    /// the kernel reduces on `m = |n|` and flips the sign, so this
    /// is bit-exact.
    #[test]
    fn negative_order_parity() {
        let p = 160;
        let x = at(9, 4, p); // 2.25
        for n in [1i32, 2, 3, 4] {
            let (pos, _) = x.yn(n, RoundingMode::NearestEven);
            let (neg, _) = x.yn(-n, RoundingMode::NearestEven);
            let expected = if n % 2 == 0 {
                pos.clone()
            } else {
                pos.negated()
            };
            assert_eq!(
                neg.partial_cmp(&expected).0,
                Some(Ordering::Equal),
                "Y_(-{n}) = (-1)^{n} Y_{n}"
            );
        }
    }

    // mpmath bessely (dps = 340) for the regime-boundary check.
    const Y0_256: &str = "-0.0338129017179245490904802001579105142989436828054872732411660094390095268631908096790145694330311129900788922935589008234";
    const Y1_256: &str = "0.0365875273992917400808021213386381993171223828654781838181095392316295493832100467605787409863810864302865546691367985571";

    /// Cross-regime continuity: at `x = 256` the DLMF 10.17.4
    /// asymptotic and the DLMF 10.8.1 series (called directly) agree
    /// and both match mpmath. Pins the `bessel_j_threshold`
    /// crossover and the reused `a_k(n)` recurrence against the
    /// independent log-series path (the 6o
    /// `asymptotic_miller_continuity` analog).
    #[test]
    fn asymptotic_series_continuity() {
        let p = 113;
        let x = at(256, 1, p);
        let (a0, _, _) = bessel_y_asymptotic(0, &x, p);
        let (s0, _) = bessel_y_series(0, &x, p);
        assert!(close_at(&a0, &s0, p - 14), "asymp vs series Y_0(256)");
        assert!(close_at(&a0, &py(Y0_256, p), p - 12), "asymp Y_0(256)");
        assert!(close_at(&s0, &py(Y0_256, p), p - 12), "series Y_0(256)");
        let (a1, _, _) = bessel_y_asymptotic(1, &x, p);
        let (s1, _) = bessel_y_series(1, &x, p);
        assert!(close_at(&a1, &s1, p - 14), "asymp vs series Y_1(256)");
        assert!(close_at(&a1, &py(Y1_256, p), p - 12), "asymp Y_1(256)");
    }

    // mpmath bessely (dps = 330) for the p = 1024 second-term pin.
    const Y0_52_BIG: &str = "0.498070359615231887827472350362089806115062532656815304299746294074063215575669958123550705895373981923856729430019346127753372515911094543050074024012999730851446377594493271835804721950823221520618551680086498462052887794347122834260248470497165541676516778255376763944618616628570103413678087256903882209128454569424527607827318";
    const Y2_256_BIG: &str = "0.0340987417757315158098614667308686252311087014216238215522449902142566327177471381693315908469872152278155060019115320621472548227534413709121545010974754454341942206357994610878177648119695691004956102002044082882671327695788506702585311791097572449883224212162085878637866065222838728564465367755616196261472949877860202203482873";

    /// Second-term-matters pin (the `derive, don't recall` reflex):
    /// the series path `Y₀(2.5)` and the asymptotic+recurrence path
    /// `Y₂(256)` at `p = 1024`, validated to `p − 2` against the
    /// 330-digit references. A coefficient or harmonic-reduction
    /// error invisible at low precision fails here.
    #[test]
    fn high_precision_pin() {
        let x = at(5, 2, 1024);
        let (r, _) = x.y0(RoundingMode::NearestEven);
        assert!(close_at(&r, &py(Y0_52_BIG, 1024), 1022), "Y_0(2.5) p=1024");

        let x = at(256, 1, 1024);
        let (r, _) = x.yn(2, RoundingMode::NearestEven);
        assert!(close_at(&r, &py(Y2_256_BIG, 1024), 1022), "Y_2(256) p=1024");
    }

    /// J/Y Wronskian `J_{n+1}(x)·Y_n(x) − J_n(x)·Y_{n+1}(x) =
    /// 2/(πx)` (DLMF 10.5.2): the load-bearing cross-tie binding the
    /// 6o `J` kernel to the new `Y` kernel, the deliverable 6o could
    /// not produce.
    #[test]
    fn jy_wronskian() {
        let p = 200;
        for &(num, den) in &[(5i64, 2i64), (7, 1), (9, 4)] {
            let x = at(num, den, p);
            for n in [0i32, 1, 2, 3] {
                let (jn, _) = x.jn(n, RoundingMode::NearestEven);
                let (jn1, _) = x.jn(n + 1, RoundingMode::NearestEven);
                let (yn, _) = x.yn(n, RoundingMode::NearestEven);
                let (yn1, _) = x.yn(n + 1, RoundingMode::NearestEven);
                let (a, _) = jn1.mul(&yn, RoundingMode::NearestEven);
                let (b, _) = jn.mul(&yn1, RoundingMode::NearestEven);
                let (lhs, _) = a.sub(&b, RoundingMode::NearestEven);
                let two = BigFloat::try_from_i64_exact(2, p).unwrap();
                let (pi_x, _) = pi_at(p).mul(&x, RoundingMode::NearestEven);
                let (rhs, _) = two.div(&pi_x, RoundingMode::NearestEven);
                assert!(
                    close_at(&lhs, &rhs, p - 12),
                    "Wronskian n={n} x={num}/{den}"
                );
            }
        }
    }

    #[test]
    fn y_positive_zero_is_pole() {
        let z = BigFloat::try_new_zero(Sign::Positive, 53).unwrap();
        let (r, s) = z.y0(RoundingMode::NearestEven);
        assert!(r.is_infinite() && r.is_sign_negative(), "Y0(+0) = −∞");
        assert!(s.div_by_zero());
        let (r1, s1) = z.y1(RoundingMode::NearestEven);
        assert!(r1.is_infinite() && r1.is_sign_negative(), "Y1(+0) = −∞");
        assert!(s1.div_by_zero());
        let (rn, sn) = z.yn(3, RoundingMode::NearestEven);
        assert!(rn.is_infinite() && rn.is_sign_negative(), "Y3(+0) = −∞");
        assert!(sn.div_by_zero());
    }

    #[test]
    fn y_negative_zero_is_the_pole() {
        // pf-k8ax (ADR-0123): Y0(±0) is the pole −∞ + DIV_BY_ZERO for
        // BOTH zero signs (−0 groups with the pole, not the negative
        // axis); it used to wrongly return NaN + INVALID.
        let z = BigFloat::try_new_zero(Sign::Negative, 53).unwrap();
        let (r, s) = z.y0(RoundingMode::NearestEven);
        assert!(r.is_infinite() && r.is_sign_negative(), "Y0(−0) = −∞");
        assert!(s.div_by_zero() && !s.invalid());
    }

    #[test]
    fn y_negative_argument_is_invalid() {
        let x = BigFloat::try_from_i64_exact(-3, 53).unwrap();
        for n in [0i32, 1, 2, -2] {
            let (r, s) = x.yn(n, RoundingMode::NearestEven);
            assert!(r.is_quiet_nan(), "Y_{n}(−3) = NaN (complex)");
            assert!(s.invalid());
        }
    }

    #[test]
    fn y_positive_infinity_is_zero() {
        let inf = BigFloat::try_new_infinity(Sign::Positive, 53).unwrap();
        for n in [0i32, 1, 5, -3] {
            let (r, st) = inf.yn(n, RoundingMode::NearestEven);
            assert!(r.is_zero() && r.is_sign_positive(), "Y_{n}(+∞) = +0");
            assert!(!st.invalid());
        }
    }

    #[test]
    fn y_negative_infinity_is_invalid() {
        let ninf = BigFloat::try_new_infinity(Sign::Negative, 53).unwrap();
        let (r, s) = ninf.y1(RoundingMode::NearestEven);
        assert!(r.is_quiet_nan(), "Y1(−∞) = NaN (complex)");
        assert!(s.invalid());
    }

    #[test]
    fn y_quiet_nan_propagates() {
        let q = BigFloat::try_new_quiet_nan(Sign::Positive, 53, &[]).unwrap();
        let (r, _) = q.yn(2, RoundingMode::NearestEven);
        assert!(r.is_quiet_nan());
    }

    #[test]
    fn y_signaling_nan_raises_invalid() {
        let sn = BigFloat::try_new_signaling_nan(Sign::Positive, 53, &[]).unwrap();
        let (r, st) = sn.y0(RoundingMode::NearestEven);
        assert!(r.is_quiet_nan());
        assert!(st.invalid());
    }

    #[test]
    fn precision_zero_is_rejected() {
        let x = BigFloat::try_from_i64_exact(2, 53).unwrap();
        assert!(x.y0_round(0, RoundingMode::NearestEven).is_err());
        assert!(x.yn_round(3, 0, RoundingMode::NearestEven).is_err());
    }
}
