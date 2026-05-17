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
//!   `ξ = ζ + π/4`.
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

    let value = airy_eval_normal(which, x, target_precision);
    let (rounded, status) = value
        .round_to_precision(target_precision, mode)
        .expect("precision >= 1");
    auto_raise(status);
    (rounded, status)
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

/// Evaluate an Airy function at a finite non-zero argument. The
/// three-regime sign-aware dispatch (Maclaurin / exponential
/// asymptotic / oscillatory asymptotic) lands in beads pf-nfb; until
/// then every finite argument routes through the convergent
/// Maclaurin series, which is valid for all `z` (the entire
/// solutions `f`, `g`).
fn airy_eval_normal(which: AiryFn, x: &BigFloat, target_precision: u32) -> BigFloat {
    airy_series(which, x, target_precision)
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
    let extra = if e_x <= 0 {
        64
    } else {
        // |x| < 2^{e_x+1} ⇒ (2/3)|x|^{3/2} < (2/3)·2^{(3/2)(e_x+1)};
        // bound the shift and scale by log₂e ≈ 23/16.
        let shift = ((3 * (e_x + 1)) / 2).min(20) as u32;
        let mag: u64 = 1u64 << shift;
        (mag.saturating_mul(23) / 24).min(4096) as u32
    };
    let working_prec = target_precision
        .saturating_add(64)
        .saturating_add(extra)
        .min(target_precision.saturating_add(4096));

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
}
