//! Airy functions `Ai`, `Bi` and their derivatives `Ai′`, `Bi′`
//! (DLMF Chapter 9). All four are entire on the real line: no
//! poles, no domain restriction.
//!
//! Three regimes, dispatched on the binary exponent of `|x|` and
//! its sign (the [`super::erf`] / [`super::si`] dispatch idiom
//! extended with a sign-aware asymptotic split):
//!
//! - Small `|x|`: the Maclaurin series (DLMF 9.4.1–9.4.6) in the
//!   two entire solutions `f`, `g`, combined with the boundary
//!   constants `Ai(0)`, `Ai′(0)` (DLMF 9.2.3–9.2.6).
//! - Large positive `x`: the exponential asymptotic (DLMF
//!   9.7.5–9.7.8) in `ζ = (2/3)·x^{3/2}`, summed to its smallest
//!   term.
//! - Large negative `x`: the oscillatory asymptotic (DLMF
//!   9.7.9–9.7.12) in `ζ = (2/3)·|x|^{3/2}` with the phase
//!   `φ = ζ − π/4`.
//!
//! Special cases:
//!
//! - `Ai(±0)`, `Bi(±0)`, `Ai′(±0)`, `Bi′(±0)`: the exact boundary
//!   constants (finite, normal).
//! - `Ai(+∞) = +0`, `Ai′(+∞) = −0`, `Bi(+∞) = +∞`, `Bi′(+∞) = +∞`
//!   (the exact limits at an infinite argument, `Status::OK`, the
//!   `exp(+∞)`/`gamma(+∞)` convention).
//! - `Ai(−∞) = Bi(−∞) = Ai′(−∞) = Bi′(−∞) = +0` by the
//!   decaying-envelope convention: the true behaviour at `−∞` is a
//!   bounded oscillation with no limit; the conservative total
//!   result is `+0` with `Status::OK`. ADR-0021 records this.
//! - `f(NaN) = NaN`; `sNaN` raises `INVALID`.

use crate::big::{BigFloat, BuildError};
use crate::class::Class;
use crate::rounding::RoundingMode;
use crate::sign::Sign;
use crate::status::{auto_raise, Status};

#[cfg(feature = "fixed")]
use crate::fixed::FixedFloat;
#[cfg(feature = "fixed")]
use crate::mantissa::limbs_for;

use super::pi_at;
use super::ziv::ziv_round;
use super::ziv_calibration::AIRY_ERROR_GUARD;

/// Which Airy function a kernel invocation evaluates. The four share
/// the boundary constants, the `f`/`g` Maclaurin series, the
/// `u_k`/`v_k` asymptotic coefficient recurrence, and `ζ`/`x^{1/4}`,
/// so the kernel is parameterised rather than duplicated.
#[derive(Clone, Copy, PartialEq, Eq, Debug)]
pub(super) enum AiryFn {
    Ai,
    Bi,
    AiPrime,
    BiPrime,
}

impl BigFloat {
    /// `Ai(self)` rounded under `mode` to `self.precision`.
    #[must_use]
    pub fn ai(&self, mode: RoundingMode) -> (Self, Status) {
        self.ai_round(self.precision, mode)
            .expect("self.precision >= 1 by invariant")
    }

    /// `Ai(self)` with explicit result precision.
    ///
    /// Correctly rounded under every IEEE 754-2019 rounding mode
    /// via the shared `crate::math::ziv::ziv_round` driver (slice
    /// p1.31, ADR-0038).
    pub fn ai_round(
        &self,
        target_precision: u32,
        mode: RoundingMode,
    ) -> Result<(Self, Status), BuildError> {
        if target_precision == 0 {
            return Err(BuildError::PrecisionZero);
        }
        Ok(airy_kernel(AiryFn::Ai, self, target_precision, mode))
    }

    /// `Bi(self)` rounded under `mode` to `self.precision`.
    #[must_use]
    pub fn bi(&self, mode: RoundingMode) -> (Self, Status) {
        self.bi_round(self.precision, mode)
            .expect("self.precision >= 1 by invariant")
    }

    /// `Bi(self)` with explicit result precision.
    ///
    /// Correctly rounded under every IEEE 754-2019 rounding mode
    /// via the shared `crate::math::ziv::ziv_round` driver (slice
    /// p1.31, ADR-0038).
    pub fn bi_round(
        &self,
        target_precision: u32,
        mode: RoundingMode,
    ) -> Result<(Self, Status), BuildError> {
        if target_precision == 0 {
            return Err(BuildError::PrecisionZero);
        }
        Ok(airy_kernel(AiryFn::Bi, self, target_precision, mode))
    }

    /// `Ai′(self)` (derivative of `Ai`) rounded under `mode` to
    /// `self.precision`.
    #[must_use]
    pub fn ai_prime(&self, mode: RoundingMode) -> (Self, Status) {
        self.ai_prime_round(self.precision, mode)
            .expect("self.precision >= 1 by invariant")
    }

    /// `Ai′(self)` with explicit result precision.
    ///
    /// Correctly rounded under every IEEE 754-2019 rounding mode
    /// via the shared `crate::math::ziv::ziv_round` driver (slice
    /// p1.31, ADR-0038).
    pub fn ai_prime_round(
        &self,
        target_precision: u32,
        mode: RoundingMode,
    ) -> Result<(Self, Status), BuildError> {
        if target_precision == 0 {
            return Err(BuildError::PrecisionZero);
        }
        Ok(airy_kernel(AiryFn::AiPrime, self, target_precision, mode))
    }

    /// `Bi′(self)` (derivative of `Bi`) rounded under `mode` to
    /// `self.precision`.
    #[must_use]
    pub fn bi_prime(&self, mode: RoundingMode) -> (Self, Status) {
        self.bi_prime_round(self.precision, mode)
            .expect("self.precision >= 1 by invariant")
    }

    /// `Bi′(self)` with explicit result precision.
    ///
    /// Correctly rounded under every IEEE 754-2019 rounding mode
    /// via the shared `crate::math::ziv::ziv_round` driver (slice
    /// p1.31, ADR-0038).
    pub fn bi_prime_round(
        &self,
        target_precision: u32,
        mode: RoundingMode,
    ) -> Result<(Self, Status), BuildError> {
        if target_precision == 0 {
            return Err(BuildError::PrecisionZero);
        }
        Ok(airy_kernel(AiryFn::BiPrime, self, target_precision, mode))
    }
}

#[cfg(feature = "fixed")]
impl<const PREC: u32> FixedFloat<PREC>
where
    [(); limbs_for(PREC)]:,
{
    /// `Ai(self)` for `FixedFloat`. Delegates to [`BigFloat::ai`].
    #[must_use]
    pub fn ai(&self, mode: RoundingMode) -> (Self, Status) {
        let (big, status) = self.to_big().ai(mode);
        (
            Self::try_from_big_exact(big).expect("precision matches"),
            status,
        )
    }

    /// `Bi(self)` for `FixedFloat`. Delegates to [`BigFloat::bi`].
    #[must_use]
    pub fn bi(&self, mode: RoundingMode) -> (Self, Status) {
        let (big, status) = self.to_big().bi(mode);
        (
            Self::try_from_big_exact(big).expect("precision matches"),
            status,
        )
    }

    /// `Ai′(self)` for `FixedFloat`. Delegates to
    /// [`BigFloat::ai_prime`].
    #[must_use]
    pub fn ai_prime(&self, mode: RoundingMode) -> (Self, Status) {
        let (big, status) = self.to_big().ai_prime(mode);
        (
            Self::try_from_big_exact(big).expect("precision matches"),
            status,
        )
    }

    /// `Bi′(self)` for `FixedFloat`. Delegates to
    /// [`BigFloat::bi_prime`].
    #[must_use]
    pub fn bi_prime(&self, mode: RoundingMode) -> (Self, Status) {
        let (big, status) = self.to_big().bi_prime(mode);
        (
            Self::try_from_big_exact(big).expect("precision matches"),
            status,
        )
    }
}

fn airy_kernel(
    which: AiryFn,
    x: &BigFloat,
    target_precision: u32,
    mode: RoundingMode,
) -> (BigFloat, Status) {
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
            // f(±0) = the exact boundary constant (finite, normal).
            let working = target_precision.saturating_add(64);
            let value = airy_zero_value(which, working);
            let (rounded, status) = value
                .round_to_precision(target_precision, mode)
                .expect("precision >= 1");
            auto_raise(status);
            return (rounded, status);
        }
        Class::Infinity { sign } => {
            let result = airy_at_infinity(which, *sign, target_precision);
            return (result, Status::OK);
        }
        Class::Normal { .. } => {}
    }

    // Regime decision pinned from target_precision so it does not
    // flip across Ziv retries (slice p1.4 erf precedent). The
    // three-regime helpers (Maclaurin, +x asymptotic, -x oscillatory
    // asymptotic) each take their working precision and apply their
    // internal +64 + cancellation-tax guard on top.
    let e_x = match &x.abs().class {
        Class::Normal { exponent, .. } => *exponent,
        _ => 0,
    };
    let threshold = airy_threshold_exponent(target_precision);
    let use_asymptotic = e_x >= threshold;
    let asymptotic_neg = use_asymptotic && matches!(x.sign(), Sign::Negative);

    let (result, status) = ziv_round(
        |w| {
            if use_asymptotic {
                if asymptotic_neg {
                    airy_asymptotic_neg(which, &x.abs(), w)
                } else {
                    airy_asymptotic_pos(which, x, w)
                }
            } else {
                airy_series(which, x, w)
            }
        },
        target_precision,
        mode,
        AIRY_ERROR_GUARD,
    );
    auto_raise(status);
    (result, status)
}

/// The exact limit of an Airy function at `±∞` (DLMF 9.7). At `+∞`:
/// `Ai → +0`, `Ai′ → −0`, `Bi → +∞`, `Bi′ → +∞` (true limits at an
/// infinite argument, mirroring `exp(+∞)`). At `−∞` all four → `+0`
/// by the decaying-envelope convention (ADR-0021): the true
/// behaviour is a bounded oscillation with no limit; the
/// conservative total result is `+0`.
fn airy_at_infinity(which: AiryFn, sign: Sign, target_precision: u32) -> BigFloat {
    if matches!(sign, Sign::Negative) {
        return BigFloat::try_new_zero(Sign::Positive, target_precision).expect("precision >= 1");
    }
    match which {
        AiryFn::Ai => {
            BigFloat::try_new_zero(Sign::Positive, target_precision).expect("precision >= 1")
        }
        AiryFn::AiPrime => {
            BigFloat::try_new_zero(Sign::Negative, target_precision).expect("precision >= 1")
        }
        AiryFn::Bi | AiryFn::BiPrime => {
            BigFloat::try_new_infinity(Sign::Positive, target_precision).expect("precision >= 1")
        }
    }
}

/// `Γ(num/den)` at `working` precision via the in-crate gamma
/// kernel. `num`, `den` are small integers (1, 2, 3) so the rational
/// is formed exactly enough at working precision; gamma carries its
/// own internal guard.
fn gamma_of_ratio(num: i64, den: i64, working: u32) -> BigFloat {
    let n = BigFloat::try_from_i64_exact(num, working).expect("precision >= 1");
    let d = BigFloat::try_from_i64_exact(den, working).expect("precision >= 1");
    let (ratio, _) = n.div(&d, RoundingMode::NearestEven);
    let (g, _) = ratio.gamma(RoundingMode::NearestEven);
    g
}

/// `3^(num/den) = exp((num/den)·ln 3)` at `working` precision.
/// Built from `ln`/`exp` directly rather than `pow`: `pow` composes
/// the same exp/ln and slice 7c (not yet shipped) is what tightens
/// it to 1 ULP, so the direct form keeps the error budget explicit
/// (ADR-0021, risk 2).
fn three_pow_ratio(num: i64, den: i64, working: u32) -> BigFloat {
    let three = BigFloat::try_from_i64_exact(3, working).expect("precision >= 1");
    let (ln3, _) = three.ln(RoundingMode::NearestEven);
    let n = BigFloat::try_from_i64_exact(num, working).expect("precision >= 1");
    let d = BigFloat::try_from_i64_exact(den, working).expect("precision >= 1");
    let (a, _) = ln3.mul(&n, RoundingMode::NearestEven);
    let (b, _) = a.div(&d, RoundingMode::NearestEven);
    let (r, _) = b.exp(RoundingMode::NearestEven);
    r
}

/// The boundary constant `f(0)` for `f ∈ {Ai, Bi, Ai′, Bi′}`,
/// DLMF 9.2.3–9.2.6:
///
/// ```text
/// Ai(0)  =  1 / (3^(2/3)·Γ(2/3))
/// Ai′(0) = −1 / (3^(1/3)·Γ(1/3))
/// Bi(0)  =  1 / (3^(1/6)·Γ(2/3))
/// Bi′(0) =  3^(1/6) / Γ(1/3)
/// ```
///
/// Composed at `working` precision from the in-crate gamma kernel and
/// `exp`/`ln`; no hardcoded table (ADR-0021). Memoization is
/// deferred pending a bench.
fn airy_zero_value(which: AiryFn, working_prec: u32) -> BigFloat {
    let one = BigFloat::try_from_i64_exact(1, working_prec).expect("precision >= 1");
    match which {
        AiryFn::Ai => {
            let p = three_pow_ratio(2, 3, working_prec);
            let g = gamma_of_ratio(2, 3, working_prec);
            let (denom, _) = p.mul(&g, RoundingMode::NearestEven);
            let (r, _) = one.div(&denom, RoundingMode::NearestEven);
            r
        }
        AiryFn::AiPrime => {
            let p = three_pow_ratio(1, 3, working_prec);
            let g = gamma_of_ratio(1, 3, working_prec);
            let (denom, _) = p.mul(&g, RoundingMode::NearestEven);
            let (r, _) = one.div(&denom, RoundingMode::NearestEven);
            r.negated()
        }
        AiryFn::Bi => {
            let p = three_pow_ratio(1, 6, working_prec);
            let g = gamma_of_ratio(2, 3, working_prec);
            let (denom, _) = p.mul(&g, RoundingMode::NearestEven);
            let (r, _) = one.div(&denom, RoundingMode::NearestEven);
            r
        }
        AiryFn::BiPrime => {
            let p = three_pow_ratio(1, 6, working_prec);
            let g = gamma_of_ratio(1, 3, working_prec);
            let (r, _) = p.div(&g, RoundingMode::NearestEven);
            r
        }
    }
}

/// Smallest binary exponent of `|x|` at which the Airy asymptotic
/// expansions already give `target_precision + 32` bits.
///
/// The Airy coefficient ratio is `u_k/u_{k−1} ≈ k²` for large `k`,
/// so the asymptotic series `Σ u_k/ζ^k` minimises near `k ≈ √ζ` and
/// its optimally-truncated error is `≈ e^{−2√ζ}` (NOT `e^{−2ζ}`;
/// the smallest term, not the late tail). The requirement is
/// therefore `2√ζ·log₂e ≥ p+32` with `ζ = (2/3)|x|^{3/2}`, i.e.
///
/// ```text
/// |x|³ ≥ (9/4)·((p+32)/(2·log₂e))⁴
/// ```
///
/// Solved in integer arithmetic on the cube of the `|x|` lower
/// bound `2^{e_x}` (so the helper stays no_std-clean), with the
/// conservative rational `log₂e ≈ 23/16`. The fourth-power growth
/// means the asymptotic only takes over for genuinely large `|x|`;
/// everything below routes through the (uncapped-guard) Maclaurin
/// series, which is correct at any precision. ADR-0021 records the
/// `e^{−2√ζ}` accuracy law and this regime boundary.
pub(super) fn airy_threshold_exponent(target_precision: u32) -> i64 {
    let bits: u128 = u128::from(target_precision) + 32;
    // K = (p+32)/(2·log₂e) ≈ (p+32)·8/23.
    // need = ⌈(9/4)·K⁴⌉ = ⌈9·bits⁴·8⁴ / (4·23⁴)⌉.
    let bits2 = bits.saturating_mul(bits);
    let bits4 = bits2.saturating_mul(bits2);
    let need: u128 = bits4.saturating_mul(9 * 4096).div_ceil(4 * 279_841);
    let mut e: i64 = 1;
    let mut pow_8: u128 = 8; // 8^1 = 2^{3·1}, the lower bound on |x|³
    while pow_8 < need && e < 90 {
        e += 1;
        pow_8 = pow_8.saturating_mul(8);
    }
    e
}

/// Accumulate the asymptotic coefficient sums. `u_k` follows the
/// DLMF 9.7.2 recurrence `u_k = (6k−5)(6k−3)(6k−1)/(216k)·u_{k−1}`,
/// `u_0 = 1`; `v_k = −(6k+1)/(6k−1)·u_k`. With `inv_zeta = 1/ζ`, the
/// running term is `t_k = u_k/ζ^k`. Returns, summed to the smallest
/// term (optimal truncation, the [`super::si::si_ci_f`] idiom):
///
/// - `s_u_alt = Σ (−1)^k u_k/ζ^k`, `s_u = Σ u_k/ζ^k`
/// - `s_v_alt = Σ (−1)^k v_k/ζ^k`, `s_v = Σ v_k/ζ^k`
fn airy_uv_sums(inv_zeta: &BigFloat, working: u32) -> (BigFloat, BigFloat, BigFloat, BigFloat) {
    let one = BigFloat::try_from_i64_exact(1, working).expect("precision >= 1");
    let mut t_u = one.clone(); // u_0/ζ^0
    let mut s_u = one.clone();
    let mut s_u_alt = one.clone();
    let mut s_v = one.clone(); // v_0 = 1
    let mut s_v_alt = one.clone();
    let mut prev_mag = series_magnitude(&t_u);
    let max_iter: i64 = 1 << 20;
    for k in 1..=max_iter {
        // DLMF 9.7.2: u_k = (6k−5)(6k−3)(6k−1) / [(2k−1)·216·k] · u_{k−1}.
        let num1 = BigFloat::try_from_i64_exact(6 * k - 5, working).expect("precision >= 1");
        let num2 = BigFloat::try_from_i64_exact(6 * k - 3, working).expect("precision >= 1");
        let num3 = BigFloat::try_from_i64_exact(6 * k - 1, working).expect("precision >= 1");
        let den1 = BigFloat::try_from_i64_exact(216 * k, working).expect("precision >= 1");
        let den2 = BigFloat::try_from_i64_exact(2 * k - 1, working).expect("precision >= 1");
        let (a, _) = t_u.mul(&num1, RoundingMode::NearestEven);
        let (a, _) = a.mul(&num2, RoundingMode::NearestEven);
        let (a, _) = a.mul(&num3, RoundingMode::NearestEven);
        let (a, _) = a.div(&den1, RoundingMode::NearestEven);
        let (a, _) = a.div(&den2, RoundingMode::NearestEven);
        let (cand, _) = a.mul(inv_zeta, RoundingMode::NearestEven);
        let mag = series_magnitude(&cand);
        if mag > prev_mag {
            break; // smallest term passed: optimal truncation.
        }
        prev_mag = mag;
        t_u = cand;

        // v_k/ζ^k = (t_u) · −(6k+1)/(6k−1).
        let vn = BigFloat::try_from_i64_exact(6 * k + 1, working).expect("precision >= 1");
        let vd = BigFloat::try_from_i64_exact(6 * k - 1, working).expect("precision >= 1");
        let (tv, _) = t_u.mul(&vn, RoundingMode::NearestEven);
        let (tv, _) = tv.div(&vd, RoundingMode::NearestEven);
        let t_v = tv.negated();

        let (s, _) = s_u.add(&t_u, RoundingMode::NearestEven);
        s_u = s;
        let (s, _) = s_v.add(&t_v, RoundingMode::NearestEven);
        s_v = s;
        if k % 2 == 0 {
            let (s, _) = s_u_alt.add(&t_u, RoundingMode::NearestEven);
            s_u_alt = s;
            let (s, _) = s_v_alt.add(&t_v, RoundingMode::NearestEven);
            s_v_alt = s;
        } else {
            let (s, _) = s_u_alt.sub(&t_u, RoundingMode::NearestEven);
            s_u_alt = s;
            let (s, _) = s_v_alt.sub(&t_v, RoundingMode::NearestEven);
            s_v_alt = s;
        }
    }
    (s_u_alt, s_u, s_v_alt, s_v)
}

/// `x^{3/2} = x·√x` and `x^{1/4} = √√x` at `working` precision,
/// built from `sqrt` (never `pow`; ADR-0021 risk 2).
fn x_three_halves_and_quarter(xw: &BigFloat) -> (BigFloat, BigFloat) {
    let (sqrt_x, _) = xw.sqrt(RoundingMode::NearestEven);
    let (x32, _) = xw.mul(&sqrt_x, RoundingMode::NearestEven);
    let (x_q, _) = sqrt_x.sqrt(RoundingMode::NearestEven);
    (x32, x_q)
}

fn two_thirds(v: &BigFloat, working: u32) -> BigFloat {
    let two = BigFloat::try_from_i64_exact(2, working).expect("precision >= 1");
    let three = BigFloat::try_from_i64_exact(3, working).expect("precision >= 1");
    let (a, _) = v.mul(&two, RoundingMode::NearestEven);
    let (b, _) = a.div(&three, RoundingMode::NearestEven);
    b
}

/// Large positive `x`: the exponential asymptotic, DLMF 9.7.5–9.7.8,
/// with `ζ = (2/3)x^{3/2}`:
///
/// ```text
/// Ai (x) ~  e^{−ζ}/(2√π x^{1/4}) · Σ (−1)^k u_k/ζ^k
/// Bi (x) ~  e^{+ζ}/( √π x^{1/4}) · Σ        u_k/ζ^k
/// Ai′(x) ~ −x^{1/4} e^{−ζ}/(2√π) · Σ (−1)^k v_k/ζ^k
/// Bi′(x) ~  x^{1/4} e^{+ζ}/( √π) · Σ        v_k/ζ^k
/// ```
fn airy_asymptotic_pos(which: AiryFn, x: &BigFloat, target_precision: u32) -> BigFloat {
    let working = target_precision.saturating_add(64);
    let xw = x
        .round_to_precision(working, RoundingMode::NearestEven)
        .expect("precision >= 1")
        .0;
    let (x32, x_q) = x_three_halves_and_quarter(&xw);
    let zeta = two_thirds(&x32, working);
    let one = BigFloat::try_from_i64_exact(1, working).expect("precision >= 1");
    let (inv_zeta, _) = one.div(&zeta, RoundingMode::NearestEven);
    let (s_u_alt, s_u, s_v_alt, s_v) = airy_uv_sums(&inv_zeta, working);

    let (sqrt_pi, _) = pi_at(working).sqrt(RoundingMode::NearestEven);
    let two_sqrt_pi = {
        let two = BigFloat::try_from_i64_exact(2, working).expect("precision >= 1");
        two.mul(&sqrt_pi, RoundingMode::NearestEven).0
    };
    let (e_pos, _) = zeta.exp(RoundingMode::NearestEven);
    let (e_neg, _) = zeta.negated().exp(RoundingMode::NearestEven);

    match which {
        AiryFn::Ai => {
            // e^{−ζ}·s_u_alt / (2√π·x^{1/4})
            let (n, _) = e_neg.mul(&s_u_alt, RoundingMode::NearestEven);
            let (d, _) = two_sqrt_pi.mul(&x_q, RoundingMode::NearestEven);
            n.div(&d, RoundingMode::NearestEven).0
        }
        AiryFn::Bi => {
            let (n, _) = e_pos.mul(&s_u, RoundingMode::NearestEven);
            let (d, _) = sqrt_pi.mul(&x_q, RoundingMode::NearestEven);
            n.div(&d, RoundingMode::NearestEven).0
        }
        AiryFn::AiPrime => {
            // −x^{1/4} e^{−ζ} s_v_alt / (2√π)
            let (n, _) = x_q.mul(&e_neg, RoundingMode::NearestEven);
            let (n, _) = n.mul(&s_v_alt, RoundingMode::NearestEven);
            let (r, _) = n.div(&two_sqrt_pi, RoundingMode::NearestEven);
            r.negated()
        }
        AiryFn::BiPrime => {
            let (n, _) = x_q.mul(&e_pos, RoundingMode::NearestEven);
            let (n, _) = n.mul(&s_v, RoundingMode::NearestEven);
            n.div(&sqrt_pi, RoundingMode::NearestEven).0
        }
    }
}

/// Large negative `x` (`t = |x| > 0`): the oscillatory asymptotic,
/// DLMF 9.7.9–9.7.12, with `ζ = (2/3)t^{3/2}` and phase
/// `φ = ζ − π/4`. With the even/odd index splits
/// `Pu = Σ (−1)^k u_{2k}/ζ^{2k}`, `Qu = Σ (−1)^k u_{2k+1}/ζ^{2k+1}`
/// and `Pv`, `Qv` the `v` analogues:
///
/// ```text
/// Ai (−t) =  π^{−1/2} t^{−1/4} ( cos φ · Pu + sin φ · Qu )
/// Bi (−t) =  π^{−1/2} t^{−1/4} (−sin φ · Pu + cos φ · Qu )
/// Ai′(−t) =  π^{−1/2} t^{ 1/4} ( sin φ · Pv − cos φ · Qv )
/// Bi′(−t) =  π^{−1/2} t^{ 1/4} ( cos φ · Pv + sin φ · Qv )
/// ```
fn airy_asymptotic_neg(which: AiryFn, t: &BigFloat, target_precision: u32) -> BigFloat {
    let working = target_precision.saturating_add(64);
    let tw = t
        .round_to_precision(working, RoundingMode::NearestEven)
        .expect("precision >= 1")
        .0;
    let (t32, t_q) = x_three_halves_and_quarter(&tw);
    let zeta = two_thirds(&t32, working);
    let one = BigFloat::try_from_i64_exact(1, working).expect("precision >= 1");
    let (inv_zeta, _) = one.div(&zeta, RoundingMode::NearestEven);

    let (pu, qu, pv, qv) = airy_uv_sums_split(&inv_zeta, working);

    // Phase φ = ζ − π/4. sin/cos range-reduce internally.
    let pi_over_4 = {
        let four = BigFloat::try_from_i64_exact(4, working).expect("precision >= 1");
        pi_at(working).div(&four, RoundingMode::NearestEven).0
    };
    let (phi, _) = zeta.sub(&pi_over_4, RoundingMode::NearestEven);
    let (cos_phi, _) = phi.cos(RoundingMode::NearestEven);
    let (sin_phi, _) = phi.sin(RoundingMode::NearestEven);

    let (sqrt_pi, _) = pi_at(working).sqrt(RoundingMode::NearestEven);
    let (inv_sqrt_pi, _) = one.div(&sqrt_pi, RoundingMode::NearestEven);

    let mul = |a: &BigFloat, b: &BigFloat| a.mul(b, RoundingMode::NearestEven).0;
    let add = |a: &BigFloat, b: &BigFloat| a.add(b, RoundingMode::NearestEven).0;
    let sub = |a: &BigFloat, b: &BigFloat| a.sub(b, RoundingMode::NearestEven).0;

    match which {
        AiryFn::Ai => {
            // π^{−1/2} t^{−1/4} ( cos φ·Pu + sin φ·Qu )
            let inner = add(&mul(&cos_phi, &pu), &mul(&sin_phi, &qu));
            let (r, _) = mul(&inv_sqrt_pi, &inner).div(&t_q, RoundingMode::NearestEven);
            r
        }
        AiryFn::Bi => {
            // π^{−1/2} t^{−1/4} ( −sin φ·Pu + cos φ·Qu )
            let inner = sub(&mul(&cos_phi, &qu), &mul(&sin_phi, &pu));
            let (r, _) = mul(&inv_sqrt_pi, &inner).div(&t_q, RoundingMode::NearestEven);
            r
        }
        AiryFn::AiPrime => {
            // DLMF 9.7.10: +π^{−1/2} t^{1/4} ( sin φ·Pv − cos φ·Qv )
            let inner = sub(&mul(&sin_phi, &pv), &mul(&cos_phi, &qv));
            mul(&mul(&inv_sqrt_pi, &t_q), &inner)
        }
        AiryFn::BiPrime => {
            // π^{−1/2} t^{1/4} ( cos φ·Pv + sin φ·Qv )
            let inner = add(&mul(&cos_phi, &pv), &mul(&sin_phi, &qv));
            mul(&mul(&inv_sqrt_pi, &t_q), &inner)
        }
    }
}

/// The even/odd-index split of the u- and v-coefficient asymptotic
/// series (DLMF 9.7.9–9.7.12), summed to the smallest term:
///
/// ```text
/// Pu = Σ_{k≥0} (−1)^k u_{2k}  /ζ^{2k}     Qu = Σ (−1)^k u_{2k+1}/ζ^{2k+1}
/// Pv = Σ_{k≥0} (−1)^k v_{2k}  /ζ^{2k}     Qv = Σ (−1)^k v_{2k+1}/ζ^{2k+1}
/// ```
///
/// Iterating the single sequence `t_j = u_j/ζ^j` and routing term
/// `j` by parity: `j = 2k → P` with sign `(−1)^k`, `j = 2k+1 → Q`
/// with sign `(−1)^k`. `v_j = −(6j+1)/(6j−1)·u_j`.
fn airy_uv_sums_split(
    inv_zeta: &BigFloat,
    working: u32,
) -> (BigFloat, BigFloat, BigFloat, BigFloat) {
    let zero = BigFloat::try_new_zero(Sign::Positive, working).expect("precision >= 1");
    let one = BigFloat::try_from_i64_exact(1, working).expect("precision >= 1");
    let mut t_u = one.clone(); // j = 0
    let mut pu = one.clone();
    let mut qu = zero.clone();
    let mut pv = one.clone(); // v_0 = 1
    let mut qv = zero;
    let mut prev_mag = series_magnitude(&t_u);
    let max_iter: i64 = 1 << 20;
    for j in 1..=max_iter {
        // DLMF 9.7.2: u_j = (6j−5)(6j−3)(6j−1) / [(2j−1)·216·j] · u_{j−1}.
        let num1 = BigFloat::try_from_i64_exact(6 * j - 5, working).expect("precision >= 1");
        let num2 = BigFloat::try_from_i64_exact(6 * j - 3, working).expect("precision >= 1");
        let num3 = BigFloat::try_from_i64_exact(6 * j - 1, working).expect("precision >= 1");
        let den1 = BigFloat::try_from_i64_exact(216 * j, working).expect("precision >= 1");
        let den2 = BigFloat::try_from_i64_exact(2 * j - 1, working).expect("precision >= 1");
        let (a, _) = t_u.mul(&num1, RoundingMode::NearestEven);
        let (a, _) = a.mul(&num2, RoundingMode::NearestEven);
        let (a, _) = a.mul(&num3, RoundingMode::NearestEven);
        let (a, _) = a.div(&den1, RoundingMode::NearestEven);
        let (a, _) = a.div(&den2, RoundingMode::NearestEven);
        let (cand, _) = a.mul(inv_zeta, RoundingMode::NearestEven);
        let mag = series_magnitude(&cand);
        if mag > prev_mag {
            break; // optimal truncation.
        }
        prev_mag = mag;
        t_u = cand;

        let vn = BigFloat::try_from_i64_exact(6 * j + 1, working).expect("precision >= 1");
        let vd = BigFloat::try_from_i64_exact(6 * j - 1, working).expect("precision >= 1");
        let (tv, _) = t_u.mul(&vn, RoundingMode::NearestEven);
        let (tv, _) = tv.div(&vd, RoundingMode::NearestEven);
        let t_v = tv.negated();

        // Sign (−1)^{⌊j/2⌋}: negative when ⌊j/2⌋ is odd.
        let negate = (j / 2) % 2 == 1;
        if j % 2 == 0 {
            if negate {
                pu = pu.sub(&t_u, RoundingMode::NearestEven).0;
                pv = pv.sub(&t_v, RoundingMode::NearestEven).0;
            } else {
                pu = pu.add(&t_u, RoundingMode::NearestEven).0;
                pv = pv.add(&t_v, RoundingMode::NearestEven).0;
            }
        } else if negate {
            qu = qu.sub(&t_u, RoundingMode::NearestEven).0;
            qv = qv.sub(&t_v, RoundingMode::NearestEven).0;
        } else {
            qu = qu.add(&t_u, RoundingMode::NearestEven).0;
            qv = qv.add(&t_v, RoundingMode::NearestEven).0;
        }
    }
    (pu, qu, pv, qv)
}

fn series_magnitude(v: &BigFloat) -> i64 {
    match &v.class {
        Class::Normal { exponent, .. } => *exponent,
        Class::Zero { .. } => i64::MIN,
        _ => i64::MAX,
    }
}

/// `true` once `term` has fallen `working_prec + 8` bits below the
/// running `sum`, so further terms cannot perturb the rounded
/// result (the [`super::si`] `negligible` idiom).
fn series_negligible(term: &BigFloat, sum: &BigFloat, working_prec: u32) -> bool {
    match &term.class {
        Class::Zero { .. } => true,
        Class::Normal { exponent, .. } => {
            *exponent < series_magnitude(sum) - i64::from(working_prec) - 8
        }
        _ => false,
    }
}

/// Maclaurin series for an Airy function (DLMF 9.4.1–9.4.6),
/// evaluated at the signed argument directly (Airy has no parity).
///
/// The two entire solutions and their derivatives, with
/// `z3 = z³`:
///
/// ```text
/// f (z) = Σ_{k≥0} c_k z^{3k},      c_0 = 1, c_k = c_{k-1}/((3k)(3k−1))
/// g (z) = Σ_{k≥0} d_k z^{3k+1},    d_0 = 1, d_k = d_{k-1}/((3k+1)(3k))
/// g′(z) = Σ_{k≥0} (3k+1) d_k z^{3k}   [term ratio z³/((3k−2)(3k))]
/// f′(z) = Σ_{k≥1} 3k c_k z^{3k−1}     [first term z²/2,
///                                      ratio z³/(3(k−1)(3k−1))]
/// ```
///
/// combined with the boundary constants `c1 = Ai(0)`,
/// `c2 = −Ai′(0)` and `√3` (DLMF 9.4):
///
/// ```text
/// Ai  = c1·f  − c2·g       Bi  = √3·(c1·f  + c2·g)
/// Ai′ = c1·f′ − c2·g′      Bi′ = √3·(c1·f′ + c2·g′)
/// ```
///
/// Working precision is boosted by `≈ (2/3)|x|^{3/2}·log₂e` so the
/// peak term and the `c1·f − c2·g` cancellation do not exhaust the
/// budget (the [`super::erf::erf_maclaurin`] guard idiom scaled for
/// the `z^{3/2}` growth of `f`, `g`).
fn airy_series(which: AiryFn, x: &BigFloat, target_precision: u32) -> BigFloat {
    let e_x = match &x.abs().class {
        Class::Normal { exponent, .. } => *exponent,
        _ => 0,
    };
    // The peak term of f/g sits near 3k ≈ |x|^{3/2}; the c1·f − c2·g
    // combination then cancels down to an O(|x|^{−1/4}) result. The
    // guard must absorb ≈ (2/3)|x|^{3/2}·log₂e bits. It is NOT
    // capped: the Maclaurin path is the correctness backstop valid at
    // any precision and any |x| (the library's no-caps ethos), and
    // `airy_threshold_exponent` hands large |x| to the asymptotic, so
    // the series is only ever used where this guard is bounded.
    let extra: u32 = if e_x <= 0 {
        64
    } else {
        // |x| < 2^{e_x+1} ⇒ |x|^{3/2} < 2^{⌈(3/2)(e_x+1)⌉}; scale by
        // (2/3)·log₂e ≈ 23/24. Saturate (a precision wider than u32
        // is not representable; in practice the threshold keeps this
        // far below that).
        let shift = (((3 * (e_x + 1)) / 2) as u32).min(62);
        let mag: u64 = 1u64 << shift;
        let bits = mag.saturating_mul(23) / 24;
        u32::try_from(bits.saturating_add(64)).unwrap_or(u32::MAX)
    };
    let working_prec = target_precision.saturating_add(64).saturating_add(extra);

    let z = x
        .round_to_precision(working_prec, RoundingMode::NearestEven)
        .expect("precision >= 1")
        .0;
    let (z2, _) = z.mul(&z, RoundingMode::NearestEven);
    let (z3, _) = z2.mul(&z, RoundingMode::NearestEven);

    // k = 0 contributions, and the k = 1 term of f′.
    let mut tf = BigFloat::try_from_i64_exact(1, working_prec).expect("precision >= 1");
    let mut sum_f = tf.clone();
    let mut tg = z.clone();
    let mut sum_g = z.clone();
    let mut tgp = BigFloat::try_from_i64_exact(1, working_prec).expect("precision >= 1");
    let mut sum_gp = tgp.clone();
    let two = BigFloat::try_from_i64_exact(2, working_prec).expect("precision >= 1");
    let (mut tfp, _) = z2.div(&two, RoundingMode::NearestEven); // z²/2 (k = 1)
    let mut sum_fp = tfp.clone();

    let max_iter: i64 = 1 << 22;
    for k in 1..=max_iter {
        let three_k = 3 * k;
        let d_f1 = BigFloat::try_from_i64_exact(three_k, working_prec).expect("precision >= 1");
        let d_f2 = BigFloat::try_from_i64_exact(three_k - 1, working_prec).expect("precision >= 1");
        let (a, _) = tf.mul(&z3, RoundingMode::NearestEven);
        let (a, _) = a.div(&d_f1, RoundingMode::NearestEven);
        let (a, _) = a.div(&d_f2, RoundingMode::NearestEven);
        tf = a;
        let (s, _) = sum_f.add(&tf, RoundingMode::NearestEven);
        sum_f = s;

        let d_g1 = BigFloat::try_from_i64_exact(three_k + 1, working_prec).expect("precision >= 1");
        let d_g2 = BigFloat::try_from_i64_exact(three_k, working_prec).expect("precision >= 1");
        let (b, _) = tg.mul(&z3, RoundingMode::NearestEven);
        let (b, _) = b.div(&d_g1, RoundingMode::NearestEven);
        let (b, _) = b.div(&d_g2, RoundingMode::NearestEven);
        tg = b;
        let (s, _) = sum_g.add(&tg, RoundingMode::NearestEven);
        sum_g = s;

        let d_gp1 =
            BigFloat::try_from_i64_exact(three_k - 2, working_prec).expect("precision >= 1");
        let d_gp2 = BigFloat::try_from_i64_exact(three_k, working_prec).expect("precision >= 1");
        let (c, _) = tgp.mul(&z3, RoundingMode::NearestEven);
        let (c, _) = c.div(&d_gp1, RoundingMode::NearestEven);
        let (c, _) = c.div(&d_gp2, RoundingMode::NearestEven);
        tgp = c;
        let (s, _) = sum_gp.add(&tgp, RoundingMode::NearestEven);
        sum_gp = s;

        if k >= 2 {
            let d_fp1 =
                BigFloat::try_from_i64_exact(3 * (k - 1), working_prec).expect("precision >= 1");
            let d_fp2 =
                BigFloat::try_from_i64_exact(three_k - 1, working_prec).expect("precision >= 1");
            let (e, _) = tfp.mul(&z3, RoundingMode::NearestEven);
            let (e, _) = e.div(&d_fp1, RoundingMode::NearestEven);
            let (e, _) = e.div(&d_fp2, RoundingMode::NearestEven);
            tfp = e;
            let (s, _) = sum_fp.add(&tfp, RoundingMode::NearestEven);
            sum_fp = s;
        }

        if k > (1i64 << e_x.max(0)).saturating_add(4)
            && series_negligible(&tf, &sum_f, working_prec)
            && series_negligible(&tg, &sum_g, working_prec)
            && series_negligible(&tgp, &sum_gp, working_prec)
            && series_negligible(&tfp, &sum_fp, working_prec)
        {
            break;
        }
    }

    let c1 = airy_zero_value(AiryFn::Ai, working_prec);
    // c2 = −Ai′(0) = |Ai′(0)| (airy_zero_value(AiPrime) is negative).
    let c2 = airy_zero_value(AiryFn::AiPrime, working_prec).negated();

    match which {
        AiryFn::Ai => {
            let (p, _) = c1.mul(&sum_f, RoundingMode::NearestEven);
            let (q, _) = c2.mul(&sum_g, RoundingMode::NearestEven);
            let (r, _) = p.sub(&q, RoundingMode::NearestEven);
            r
        }
        AiryFn::AiPrime => {
            let (p, _) = c1.mul(&sum_fp, RoundingMode::NearestEven);
            let (q, _) = c2.mul(&sum_gp, RoundingMode::NearestEven);
            let (r, _) = p.sub(&q, RoundingMode::NearestEven);
            r
        }
        AiryFn::Bi => {
            let sqrt3 = {
                let three = BigFloat::try_from_i64_exact(3, working_prec).expect("precision >= 1");
                three.sqrt(RoundingMode::NearestEven).0
            };
            let (p, _) = c1.mul(&sum_f, RoundingMode::NearestEven);
            let (q, _) = c2.mul(&sum_g, RoundingMode::NearestEven);
            let (s, _) = p.add(&q, RoundingMode::NearestEven);
            let (r, _) = sqrt3.mul(&s, RoundingMode::NearestEven);
            r
        }
        AiryFn::BiPrime => {
            let sqrt3 = {
                let three = BigFloat::try_from_i64_exact(3, working_prec).expect("precision >= 1");
                three.sqrt(RoundingMode::NearestEven).0
            };
            let (p, _) = c1.mul(&sum_fp, RoundingMode::NearestEven);
            let (q, _) = c2.mul(&sum_gp, RoundingMode::NearestEven);
            let (s, _) = p.add(&q, RoundingMode::NearestEven);
            let (r, _) = sqrt3.mul(&s, RoundingMode::NearestEven);
            r
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use core::cmp::Ordering;

    /// `|v − expected| ≤ 2^(−bits)·|expected|` (the erf.rs test
    /// helper). Source of the reference decimals: `mpmath`
    /// `airyai/airybi(0[,1])` at 60 digits; treated as a fact.
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

    fn zero_at(p: u32) -> BigFloat {
        BigFloat::try_new_zero(Sign::Positive, p).unwrap()
    }

    #[test]
    fn ai_zero_boundary_constant() {
        // Ai(0) = 1/(3^(2/3)·Γ(2/3)) ≈ 0.3550280538878172392600631860…
        let (r, _) = zero_at(160).ai(RoundingMode::NearestEven);
        let want = BigFloat::parse_str(
            "0.355028053887817239260063186004183176397979174199177240583327",
            160,
            RoundingMode::NearestEven,
        )
        .unwrap()
        .0;
        assert!(close_at(&r, &want, 160 - 12));
    }

    #[test]
    fn ai_prime_zero_boundary_constant() {
        // Ai′(0) = −1/(3^(1/3)·Γ(1/3)) ≈ −0.2588194037928067984051835…
        let (r, _) = zero_at(160).ai_prime(RoundingMode::NearestEven);
        let want = BigFloat::parse_str(
            "-0.258819403792806798405183560189203963479091138354934582210002",
            160,
            RoundingMode::NearestEven,
        )
        .unwrap()
        .0;
        assert!(close_at(&r, &want, 160 - 12));
    }

    #[test]
    fn bi_zero_boundary_constant() {
        // Bi(0) = 1/(3^(1/6)·Γ(2/3)) ≈ 0.6149266274460007351509223690…
        let (r, _) = zero_at(160).bi(RoundingMode::NearestEven);
        let want = BigFloat::parse_str(
            "0.614926627446000735150922369093613553594728188648596505040879",
            160,
            RoundingMode::NearestEven,
        )
        .unwrap()
        .0;
        assert!(close_at(&r, &want, 160 - 12));
    }

    #[test]
    fn bi_prime_zero_boundary_constant() {
        // Bi′(0) = 3^(1/6)/Γ(1/3) ≈ 0.4482883573538263579148237103…
        let (r, _) = zero_at(160).bi_prime(RoundingMode::NearestEven);
        let want = BigFloat::parse_str(
            "0.448288357353826357914823710398828390866226799212262061082809",
            160,
            RoundingMode::NearestEven,
        )
        .unwrap()
        .0;
        assert!(close_at(&r, &want, 160 - 12));
    }

    #[test]
    fn airy_negative_zero_same_as_positive_zero() {
        // Airy functions are entire: f(−0) = f(+0) = f(0).
        let neg0 = BigFloat::try_new_zero(Sign::Negative, 113).unwrap();
        let (a, _) = neg0.ai(RoundingMode::NearestEven);
        let (b, _) = zero_at(113).ai(RoundingMode::NearestEven);
        assert_eq!(a.partial_cmp(&b).0, Some(Ordering::Equal));
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

    #[test]
    fn ai_small_positive() {
        // Ai(1) ≈ 0.135292416312881415524147423515466306…
        let (r, _) = at(1, 1, 160).ai(RoundingMode::NearestEven);
        let want = BigFloat::parse_str(
            "0.135292416312881415524147423515466306174944142988330706009102",
            160,
            RoundingMode::NearestEven,
        )
        .unwrap()
        .0;
        assert!(close_at(&r, &want, 160 - 16));
    }

    #[test]
    fn bi_small_positive() {
        // Bi(1) ≈ 1.207423594952871259436378817028286995…
        let (r, _) = at(1, 1, 160).bi(RoundingMode::NearestEven);
        let want = BigFloat::parse_str(
            "1.20742359495287125943637881702828699538534894464444253753862",
            160,
            RoundingMode::NearestEven,
        )
        .unwrap()
        .0;
        assert!(close_at(&r, &want, 160 - 16));
    }

    #[test]
    fn ai_prime_small_positive() {
        // Ai′(1) ≈ −0.159147441296793212787500252497229686…
        let (r, _) = at(1, 1, 160).ai_prime(RoundingMode::NearestEven);
        let want = BigFloat::parse_str(
            "-0.1591474412967932127875002524972296865738892015116109694",
            160,
            RoundingMode::NearestEven,
        )
        .unwrap()
        .0;
        assert!(close_at(&r, &want, 160 - 16));
    }

    #[test]
    fn bi_prime_small_positive() {
        // Bi′(1) ≈ 0.932435933392775632959451453674435344…
        let (r, _) = at(1, 1, 160).bi_prime(RoundingMode::NearestEven);
        let want = BigFloat::parse_str(
            "0.9324359333927756329594514536744353442695653752386283955",
            160,
            RoundingMode::NearestEven,
        )
        .unwrap()
        .0;
        assert!(close_at(&r, &want, 160 - 16));
    }

    #[test]
    fn ai_negative_argument() {
        // Ai(−2) ≈ 0.227407428201685575991924436037873799…
        let (r, _) = at(-2, 1, 160).ai(RoundingMode::NearestEven);
        let want = BigFloat::parse_str(
            "0.22740742820168557599192443603787379946077222541709671649579",
            160,
            RoundingMode::NearestEven,
        )
        .unwrap()
        .0;
        assert!(close_at(&r, &want, 160 - 16));
    }

    #[test]
    fn bi_negative_argument() {
        // Bi(−2) ≈ −0.412302587956398488083234054611461042…
        let (r, _) = at(-2, 1, 160).bi(RoundingMode::NearestEven);
        let want = BigFloat::parse_str(
            "-0.4123025879563984880832340546114610420345348344724047288",
            160,
            RoundingMode::NearestEven,
        )
        .unwrap()
        .0;
        assert!(close_at(&r, &want, 160 - 16));
    }

    #[test]
    fn wronskian_identity_small_x() {
        // Ai(x)·Bi′(x) − Ai′(x)·Bi(x) = 1/π for every x (DLMF 9.2.7).
        let inv_pi = BigFloat::parse_str(
            "0.3183098861837906715377675267450287240689192914809128975",
            200,
            RoundingMode::NearestEven,
        )
        .unwrap()
        .0;
        for (n, d) in [(1, 2), (1, 1), (-2, 1), (3, 2)] {
            let x = at(n, d, 200);
            let (ai, _) = x.ai(RoundingMode::NearestEven);
            let (bi, _) = x.bi(RoundingMode::NearestEven);
            let (aip, _) = x.ai_prime(RoundingMode::NearestEven);
            let (bip, _) = x.bi_prime(RoundingMode::NearestEven);
            let (t1, _) = ai.mul(&bip, RoundingMode::NearestEven);
            let (t2, _) = aip.mul(&bi, RoundingMode::NearestEven);
            let (w, _) = t1.sub(&t2, RoundingMode::NearestEven);
            assert!(
                close_at(&w, &inv_pi, 200 - 24),
                "Wronskian at {n}/{d}: got {w}"
            );
        }
    }

    /// Exercises the asymptotic regime (DLMF 9.7.5–9.7.12) at
    /// `|x| = 300`, where `ζ ≈ 3464` makes the optimally-truncated
    /// asymptotic accurate well past `p = 150`. References: `mpmath`
    /// `airyai/airybi(±300[,1])` at 45 digits; treated as a fact.
    #[test]
    fn asymptotic_large_argument() {
        let p = 150u32;
        let pos = [
            (
                AiryFn::Ai,
                "2.45974362033695840226898797290286351741382337e-1506",
            ),
            (
                AiryFn::Bi,
                "3.73567997123772793952321696938542273282941398e+1503",
            ),
            (
                AiryFn::AiPrime,
                "-4.26060587800406700838701857147799339670857865e-1505",
            ),
            (
                AiryFn::BiPrime,
                "6.47007616688173094631461580802880020007629144e+1504",
            ),
        ];
        let x = BigFloat::try_from_i64_exact(300, p).unwrap();
        for (which, refstr) in pos {
            let got = airy_asymptotic_pos(which, &x, p);
            let want = BigFloat::parse_str(refstr, p, RoundingMode::NearestEven)
                .unwrap()
                .0;
            assert!(close_at(&got, &want, 130), "asymptotic+ {which:?}: {got}");
        }
        let neg = [
            (
                AiryFn::Ai,
                "0.0387263629051379071866597753986827190288953645",
            ),
            (
                AiryFn::Bi,
                "-0.12991496664041682548622062250765880525502632",
            ),
            (
                AiryFn::AiPrime,
                "2.25022551383809411125438701102366844058468647",
            ),
            (
                AiryFn::BiPrime,
                "0.670652022853768111569637871724183669881112175",
            ),
        ];
        let t = BigFloat::try_from_i64_exact(300, p).unwrap();
        for (which, refstr) in neg {
            let got = airy_asymptotic_neg(which, &t, p);
            let want = BigFloat::parse_str(refstr, p, RoundingMode::NearestEven)
                .unwrap()
                .0;
            assert!(close_at(&got, &want, 130), "asymptotic− {which:?}: {got}");
        }
    }

    /// `airy_threshold_exponent` only hands `|x|` to the asymptotic
    /// once that path is accurate: the threshold grows with the
    /// fourth power of the precision (the `e^{−2√ζ}` accuracy law),
    /// so it is small at low precision and large at high precision.
    #[test]
    fn threshold_grows_with_precision() {
        let t53 = airy_threshold_exponent(53);
        let t1024 = airy_threshold_exponent(1024);
        assert!(t53 >= 1, "threshold must be a positive exponent");
        assert!(
            t1024 > t53,
            "higher precision needs larger |x| before the asymptotic \
             is accurate: p=53 → e_x≥{t53}, p=1024 → e_x≥{t1024}"
        );
    }

    #[test]
    fn airy_nan_propagates() {
        for which in [AiryFn::Ai, AiryFn::Bi, AiryFn::AiPrime, AiryFn::BiPrime] {
            let q = BigFloat::try_new_quiet_nan(Sign::Positive, 53, &[]).unwrap();
            let (r, status) = airy_kernel(which, &q, 53, RoundingMode::NearestEven);
            assert!(r.is_quiet_nan(), "{which:?} qNaN");
            assert!(status.is_ok(), "{which:?} qNaN status");
        }
    }

    #[test]
    fn airy_snan_raises_invalid() {
        for which in [AiryFn::Ai, AiryFn::Bi, AiryFn::AiPrime, AiryFn::BiPrime] {
            let sn = BigFloat::try_new_signaling_nan(Sign::Positive, 53, &[]).unwrap();
            let (r, status) = airy_kernel(which, &sn, 53, RoundingMode::NearestEven);
            assert!(r.is_quiet_nan(), "{which:?} sNaN");
            assert!(status.invalid(), "{which:?} sNaN INVALID");
        }
    }

    #[test]
    fn airy_pos_inf_limits() {
        let pi = BigFloat::try_new_infinity(Sign::Positive, 53).unwrap();
        let (ai, _) = pi.ai(RoundingMode::NearestEven);
        assert!(ai.is_zero() && !ai.is_sign_negative(), "Ai(+∞) = +0");
        let (aip, _) = pi.ai_prime(RoundingMode::NearestEven);
        assert!(aip.is_zero() && aip.is_sign_negative(), "Ai′(+∞) = −0");
        let (bi, _) = pi.bi(RoundingMode::NearestEven);
        assert!(bi.is_infinite() && !bi.is_sign_negative(), "Bi(+∞) = +∞");
        let (bip, _) = pi.bi_prime(RoundingMode::NearestEven);
        assert!(bip.is_infinite() && !bip.is_sign_negative(), "Bi′(+∞) = +∞");
    }

    #[test]
    fn airy_neg_inf_is_pos_zero_by_convention() {
        let ni = BigFloat::try_new_infinity(Sign::Negative, 53).unwrap();
        for (got, name) in [
            (ni.ai(RoundingMode::NearestEven).0, "Ai"),
            (ni.bi(RoundingMode::NearestEven).0, "Bi"),
            (ni.ai_prime(RoundingMode::NearestEven).0, "Ai′"),
            (ni.bi_prime(RoundingMode::NearestEven).0, "Bi′"),
        ] {
            assert!(
                got.is_zero() && !got.is_sign_negative(),
                "{name}(−∞) = +0 by the decaying-envelope convention"
            );
        }
    }

    #[test]
    fn airy_fn_enum_is_copy() {
        // Guards the parameterised-kernel design: AiryFn must stay a
        // trivial Copy tag so the four entry points share one kernel.
        let a = AiryFn::Ai;
        let b = a;
        assert_eq!(a, b);
        assert_ne!(AiryFn::Ai, AiryFn::Bi);
    }

    #[test]
    fn public_api_round_and_precision_zero() {
        let x = at(1, 1, 53);
        // The convenience method equals `*_round(self.precision, …)`.
        for (conv, rnd) in [
            (
                x.ai(RoundingMode::NearestEven).0,
                x.ai_round(53, RoundingMode::NearestEven),
            ),
            (
                x.bi(RoundingMode::NearestEven).0,
                x.bi_round(53, RoundingMode::NearestEven),
            ),
            (
                x.ai_prime(RoundingMode::NearestEven).0,
                x.ai_prime_round(53, RoundingMode::NearestEven),
            ),
            (
                x.bi_prime(RoundingMode::NearestEven).0,
                x.bi_prime_round(53, RoundingMode::NearestEven),
            ),
        ] {
            assert_eq!(conv.partial_cmp(&rnd.unwrap().0).0, Some(Ordering::Equal));
        }
        // target_precision == 0 is a typed error on every entry point.
        assert!(x.ai_round(0, RoundingMode::NearestEven).is_err());
        assert!(x.bi_round(0, RoundingMode::NearestEven).is_err());
        assert!(x.ai_prime_round(0, RoundingMode::NearestEven).is_err());
        assert!(x.bi_prime_round(0, RoundingMode::NearestEven).is_err());
    }

    #[cfg(feature = "fixed")]
    #[test]
    fn fixed_float_delegation() {
        use crate::fixed::FixedFloat;
        // FixedFloat delegates to BigFloat and round-trips precision.
        let one = FixedFloat::<160>::try_from_i64_exact(1).unwrap();
        let (ai, _) = one.ai(RoundingMode::NearestEven);
        let want = BigFloat::parse_str(
            "0.135292416312881415524147423515466306174944142988330706009102",
            160,
            RoundingMode::NearestEven,
        )
        .unwrap()
        .0;
        let (bip, _) = one.bi_prime(RoundingMode::NearestEven);
        assert!(close_at(&ai.to_big(), &want, 160 - 16));
        assert!(!bip.to_big().is_zero());
    }
}
