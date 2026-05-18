//! Bessel functions of the second kind `Y0`, `Y1`, `Yn` (DLMF
//! Chapter 10): ordinary Bessel of integer order, real argument.
//!
//! Unlike [`super::bessel_j`], `Y` is real-valued only for `x > 0`:
//! `Y_n` has a logarithmic branch point at the origin and is complex
//! for `x < 0`. The domain convention follows the [`super::ci`] /
//! [`super::li`] precedent (cosine / logarithmic integral, same
//! "real-only, complex off the positive axis" shape):
//!
//! - `Y_n(+0) = −∞`, raising `DIV_BY_ZERO` (a pole: the DLMF 10.8.1
//!   `−(½x)^{−n}/π` head, and `(2/π) ln(½x) J_0` for `n = 0`, both
//!   diverge to `−∞` as `x → 0⁺`).
//! - `x < 0` (and `−0`, `−∞`) ⇒ `NaN` + `INVALID` (`Y` is complex in
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

use super::{euler_gamma_at, pi_at};
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
        Class::Zero { sign } => {
            if matches!(sign, Sign::Negative) {
                // −0: Y is complex off the positive axis (the Ci/li
                // convention; −0 groups with x < 0).
                let nan = BigFloat::try_new_quiet_nan(Sign::Positive, target_precision, &[])
                    .expect("precision >= 1");
                auto_raise(Status::INVALID);
                return (nan, Status::INVALID);
            }
            // Y_n(+0) = −∞ + DIV_BY_ZERO (a pole, DLMF 10.8.1).
            let ninf = BigFloat::try_new_infinity(Sign::Negative, target_precision)
                .expect("precision >= 1");
            auto_raise(Status::DIV_BY_ZERO);
            (ninf, Status::DIV_BY_ZERO)
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

            // Y₋ₙ(x) = (−1)ⁿ Yₙ(x) (DLMF 10.4.1): order parity only,
            // negate when m is odd and n < 0.
            let negate = (m % 2 == 1) && (n < 0);
            let value = bessel_y_eval_normal(m, x, target_precision);
            let value = if negate { value.negated() } else { value };

            let (rounded, status) = value
                .round_to_precision(target_precision, mode)
                .expect("precision >= 1");
            auto_raise(status);
            (rounded, status)
        }
    }
}

/// `Y_m(x)` for `m ≥ 0`, normal `x > 0`: the base pair plus upward
/// recurrence. Returns the unrounded working-precision value;
/// [`bessel_y_kernel`] does the single final round.
///
/// Slice 6p.2: the DLMF 10.8.1 logarithmic series carries every
/// order directly (valid for all `n ≥ 0`, all `x > 0`). Slice 6p.3
/// adds the large-`x` Hankel asymptotic and 6p.4 the regime dispatch
/// plus the upward recurrence that makes `Yₙ (n ≥ 2)` cheap; until
/// then `eval_normal` is the series alone.
fn bessel_y_eval_normal(m: u32, x: &BigFloat, target_precision: u32) -> BigFloat {
    bessel_y_series(m, x, target_precision)
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
/// Working precision is boosted `≈ x·log₂e` for the alternating
/// cancellation (the [`super::ci`] / [`super::bessel_j`] capped
/// guard). Returns the unrounded working-precision value.
fn bessel_y_series(n: u32, x: &BigFloat, target_precision: u32) -> BigFloat {
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
    y
}

#[cfg(test)]
mod tests {
    use super::*;
    use core::cmp::Ordering;

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
    fn y_negative_zero_is_invalid() {
        let z = BigFloat::try_new_zero(Sign::Negative, 53).unwrap();
        let (r, s) = z.y0(RoundingMode::NearestEven);
        assert!(r.is_quiet_nan(), "Y0(−0) = NaN");
        assert!(s.invalid());
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
