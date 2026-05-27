//! Bessel functions of the first kind `J0`, `J1`, `Jn` (DLMF
//! Chapter 10): ordinary Bessel, integer order, real argument.
//! Entire on the real line: no poles, no domain restriction.
//!
//! Order and sign are reduced before evaluation. `Jₙ(−x) =
//! (−1)ⁿ Jₙ(x)` (DLMF 10.11.1) and `J₋ₙ(x) = (−1)ⁿ Jₙ(x)` (DLMF
//! 10.4.1), so the kernel evaluates `J_m(|x|)` for `m = |n| ≥ 0`
//! and applies one parity sign.
//!
//! Three regimes dispatch on the binary exponent of `|x|` (the
//! [`super::airy`] / [`super::si`] integer exponent selector idiom):
//!
//! - Tiny `|x| < 1`: the convergent Maclaurin series (DLMF 10.2.2),
//!   the small argument backstop that keeps the `2k/x` recurrence
//!   away from `x → 0`.
//! - Moderate `|x|`: Miller backward recurrence
//!   `f_{k−1} = (2k/x)·f_k − f_{k+1}` from a seed index derived from
//!   the DLMF 10.19.1 large order decay, normalised by the sum rule
//!   `J₀ + 2·Σ J_{2k} = 1` (DLMF 10.6.1, 10.12.4).
//! - Large `|x|`: the Hankel asymptotic (DLMF 10.17.3) summed to
//!   its smallest term, with coefficients `a_k(m)` derived from
//!   DLMF 10.17.1 (cross checked against the Pochhammer form).
//!
//! ADR-0023 records the design and the coefficient provenance.
//!
//! Correctly rounded under every IEEE rounding mode via the shared
//! [`crate::math::ziv::ziv_round`] driver (slice p1.4, ADR-0022).
//! The regime decision (tiny / Miller / asymptotic) is fixed at
//! kernel entry from the input's binary exponent; the Ziv driver's
//! working precision flows into the chosen evaluator. This brings
//! `J0`, `J1`, `Jn` onto the same correctness scaffolding the slice
//! p1.2 pattern wires around `exp`, `ln`, `tanh`, `lgamma`, and
//! `erf` (the latter via slice p1.4.3).
//!
//! Special cases:
//!
//! - `J₀(±0) = 1`, `Jₙ(±0) = 0` for `n ≠ 0` (exact, DLMF 10.2.2).
//! - `Jₙ(±∞) = +0` for every order, by the decaying-envelope
//!   convention (ADR-0021, the [`super::airy`] precedent): the true
//!   behaviour at `±∞` is a bounded decaying oscillation with no
//!   limit; the conservative total result is `+0`, `Status::OK`.
//! - `Jₙ(NaN) = NaN`; `sNaN` raises `INVALID`.

use super::ziv::ziv_round;
use super::ziv_calibration::BESSEL_J_ERROR_GUARD;
use super::{pi_at, pi_over_2_at};
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
    /// `J₀(self)` rounded under `mode` to `self.precision`.
    #[must_use]
    pub fn j0(&self, mode: RoundingMode) -> (Self, Status) {
        self.j0_round(self.precision, mode)
            .expect("self.precision >= 1 by invariant")
    }

    /// `J₀(self)` with explicit result precision.
    pub fn j0_round(
        &self,
        target_precision: u32,
        mode: RoundingMode,
    ) -> Result<(Self, Status), BuildError> {
        if target_precision == 0 {
            return Err(BuildError::PrecisionZero);
        }
        Ok(bessel_j_kernel(0, self, target_precision, mode))
    }

    /// `J₁(self)` rounded under `mode` to `self.precision`.
    #[must_use]
    pub fn j1(&self, mode: RoundingMode) -> (Self, Status) {
        self.j1_round(self.precision, mode)
            .expect("self.precision >= 1 by invariant")
    }

    /// `J₁(self)` with explicit result precision.
    pub fn j1_round(
        &self,
        target_precision: u32,
        mode: RoundingMode,
    ) -> Result<(Self, Status), BuildError> {
        if target_precision == 0 {
            return Err(BuildError::PrecisionZero);
        }
        Ok(bessel_j_kernel(1, self, target_precision, mode))
    }

    /// `Jₙ(self)` for integer order `n` (any `i32`, including
    /// negative: `J₋ₙ = (−1)ⁿ Jₙ`), rounded under `mode` to
    /// `self.precision`.
    #[must_use]
    pub fn jn(&self, n: i32, mode: RoundingMode) -> (Self, Status) {
        self.jn_round(n, self.precision, mode)
            .expect("self.precision >= 1 by invariant")
    }

    /// `Jₙ(self)` with explicit result precision.
    pub fn jn_round(
        &self,
        n: i32,
        target_precision: u32,
        mode: RoundingMode,
    ) -> Result<(Self, Status), BuildError> {
        if target_precision == 0 {
            return Err(BuildError::PrecisionZero);
        }
        Ok(bessel_j_kernel(n, self, target_precision, mode))
    }
}

#[cfg(feature = "fixed")]
impl<const PREC: u32> FixedFloat<PREC>
where
    [(); limbs_for(PREC)]:,
{
    /// `J₀(self)` for `FixedFloat`. Delegates to [`BigFloat::j0`].
    #[must_use]
    pub fn j0(&self, mode: RoundingMode) -> (Self, Status) {
        let (big, status) = self.to_big().j0(mode);
        (
            Self::try_from_big_exact(big).expect("precision matches"),
            status,
        )
    }

    /// `J₁(self)` for `FixedFloat`. Delegates to [`BigFloat::j1`].
    #[must_use]
    pub fn j1(&self, mode: RoundingMode) -> (Self, Status) {
        let (big, status) = self.to_big().j1(mode);
        (
            Self::try_from_big_exact(big).expect("precision matches"),
            status,
        )
    }

    /// `Jₙ(self)` for `FixedFloat`. Delegates to [`BigFloat::jn`].
    #[must_use]
    pub fn jn(&self, n: i32, mode: RoundingMode) -> (Self, Status) {
        let (big, status) = self.to_big().jn(n, mode);
        (
            Self::try_from_big_exact(big).expect("precision matches"),
            status,
        )
    }
}

/// `Jₙ(x)` for integer order `n` and real `x`.
///
/// Special values are handled directly; for a normal argument the
/// order and sign are reduced to `J_m(|x|)`, `m = |n|`, with one
/// parity sign, then the regime evaluator runs.
fn bessel_j_kernel(
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
            return (nan, Status::OK);
        }
        Class::Zero { .. } => {
            // J₀(±0) = 1, Jₙ(±0) = 0 for n ≠ 0 (DLMF 10.2.2); exact.
            let value = if m == 0 {
                BigFloat::try_from_i64_exact(1, target_precision).expect("precision >= 1")
            } else {
                BigFloat::try_new_zero(Sign::Positive, target_precision).expect("precision >= 1")
            };
            return (value, Status::OK);
        }
        Class::Infinity { .. } => {
            // Decaying-envelope convention (ADR-0021): bounded
            // oscillation with no limit → +0, Status::OK.
            let zero =
                BigFloat::try_new_zero(Sign::Positive, target_precision).expect("precision >= 1");
            return (zero, Status::OK);
        }
        Class::Normal { .. } => {}
    }

    // Normal argument: reduce to J_m(|x|) with one parity sign.
    // Jₙ(−x) = (−1)ⁿ Jₙ(x); J₋ₙ(x) = (−1)ⁿ Jₙ(x). Each negative
    // contributes (−1)^m, so the result is negated exactly when m is
    // odd and exactly one of {n<0, x<0} holds.
    let negate = (m % 2 == 1) && ((n < 0) ^ x.is_sign_negative());
    let ax = x.abs();

    ziv_round(
        |working_prec| {
            let v = bessel_j_eval_normal(m, &ax, working_prec);
            if negate {
                v.negated()
            } else {
                v
            }
        },
        target_precision,
        mode,
        BESSEL_J_ERROR_GUARD,
    )
}

/// Binary exponent of `v`, or `i64::MIN`/`i64::MAX` for zero /
/// non-finite (the [`super::si`] / [`super::airy`] `magnitude`
/// idiom).
fn magnitude(v: &BigFloat) -> i64 {
    match &v.class {
        Class::Normal { exponent, .. } => *exponent,
        Class::Zero { .. } => i64::MIN,
        _ => i64::MAX,
    }
}

/// `true` once `term` has fallen `working + 8` bits below the running
/// `sum`, so further terms cannot perturb the rounded result (the
/// [`super::si`] `negligible` idiom).
fn negligible(term: &BigFloat, sum: &BigFloat, working: u32) -> bool {
    match &term.class {
        Class::Zero { .. } => true,
        Class::Normal { exponent, .. } => *exponent < magnitude(sum) - i64::from(working) - 8,
        _ => false,
    }
}

/// Binary-exponent boundary below which the convergent Maclaurin
/// series ([`bessel_j_tiny`]) is used instead of Miller recurrence.
///
/// `e_x ≤ −1` ⇔ `|x| < 1`. The tiny regime exists only to keep the
/// `2k/x` recurrence away from `x → 0`; it is not a tuned crossover
/// (CLAUDE.md: no perf machinery without a bench). Miller is the
/// designed moderate-`|x|` path and carries everything `|x| ≥ 1`
/// until slice 6o.3 adds the large-`|x|` asymptotic upper cut.
/// Continuity across the boundary is pinned by a unit test.
fn bessel_j_tiny_threshold() -> i64 {
    -1
}

/// `J_m(ax)` for `m ≥ 0`, `0 < ax`, via the DLMF 10.2.2 convergent
/// Maclaurin series
///
/// ```text
/// J_m(x) = (x/2)^m · Σ_{k≥0} (−1)^k (x/2)^{2k} / (k!·(m+k)!)
/// ```
///
/// (integer order ⇒ `Γ(m+k+1) = (m+k)!`; entire, converges for all
/// `x`). Carried as a term recurrence (the [`super::si`] `si_ci_f`
/// idiom): with `t_0 = (x/2)^m / m!`,
/// `t_k = t_{k−1} · (−(x/2)²) / (k·(m+k))`. Hand-check:
/// `J_0` → `t_0 = 1, t_1 = −(x/2)², t_2 = (x/2)⁴/4 = x⁴/64`;
/// `J_1` → `t_0 = x/2, t_1 = −(x/2)³/2 = −x³/16`. Both match the
/// standard expansions. Returns the unrounded working-precision
/// value; the kernel performs the single final round.
fn bessel_j_tiny(m: u32, ax: &BigFloat, target_precision: u32) -> BigFloat {
    let working = target_precision.saturating_add(64);
    let x = ax
        .round_to_precision(working, RoundingMode::NearestEven)
        .expect("precision >= 1")
        .0;
    let two = BigFloat::try_from_i64_exact(2, working).expect("precision >= 1");
    let (half, _) = x.div(&two, RoundingMode::NearestEven);
    let (half_sq, _) = half.mul(&half, RoundingMode::NearestEven);

    // t_0 = (x/2)^m / m!  (m bounded in the tiny regime; in-module
    // recurrence-from-a-base-term, not pow.rs::pow_int).
    let mut term = BigFloat::try_from_i64_exact(1, working).expect("precision >= 1");
    for _ in 0..m {
        let (t, _) = term.mul(&half, RoundingMode::NearestEven);
        term = t;
    }
    for j in 1..=i64::from(m) {
        let d = BigFloat::try_from_i64_exact(j, working).expect("precision >= 1");
        let (t, _) = term.div(&d, RoundingMode::NearestEven);
        term = t;
    }

    let mut sum = term.clone();
    let max_iter: i64 = 1 << 22;
    for k in 1..=max_iter {
        let (t1, _) = term.mul(&half_sq, RoundingMode::NearestEven);
        let dk = BigFloat::try_from_i64_exact(k, working).expect("precision >= 1");
        let dmk = BigFloat::try_from_i64_exact(k + i64::from(m), working).expect("precision >= 1");
        let (t2, _) = t1.div(&dk, RoundingMode::NearestEven);
        let (t3, _) = t2.div(&dmk, RoundingMode::NearestEven);
        term = t3.negated();
        let (s, _) = sum.add(&term, RoundingMode::NearestEven);
        sum = s;
        if negligible(&term, &sum, working) {
            break;
        }
    }
    sum
}

/// `J_m(ax)` for `m ≥ 0`, `ax ≥ 1`, via Miller backward recurrence
/// (DLMF 10.6.1) with sum-rule normalization (DLMF 10.12.4).
///
/// The three-term recurrence `𝒞_{ν−1}+𝒞_{ν+1} = (2ν/z)𝒞_ν`
/// rearranged downward is `f_{k−1} = (2k/x)·f_k − f_{k+1}`. Started
/// at a high seed index `M` with `f_{M+1}=0`, `f_M=1`, the descent
/// converges to a fixed multiple `c·J_k(x)` of the recessive
/// solution. The DLMF 10.12.4 identity
/// `1 = J_0(x) + 2J_2(x) + 2J_4(x) + ⋯` (re-derived by setting
/// `t = 1` in the 10.12.1 generating function and using
/// `J_{−n}=(−1)^n J_n`) gives `S = f_0 + 2(f_2+f_4+⋯) = c·1 = c`, so
/// `J_m(x) = f_m / S`, every order from one descent.
///
/// `M` is derived, not guessed: DLMF 10.19.1
/// `J_M(x) ∼ (1/√(2πM))·(eX/(2M))^M` (the prefactor only shrinks the
/// bound, so dropping it is conservative). Requiring
/// `(eX/(2M))^M < 2^{−P}`, `P = target+64`, i.e. in natural logs
/// `M·(1 + ln(x/(2M))) < −P·ln2`, solved by an exponential search
/// (overshoot ≤ 2× is a deliberate robustness/cost trade; retune
/// only with a bench, CLAUDE.md) plus a small fixed step guard.
/// Working precision is boosted `≈ |x|·log₂e` bits for the recurrence
/// and sum-rule cancellation (the [`super::si`] guard idiom). Returns
/// the unrounded working-precision value.
fn bessel_j_miller(m: u32, ax: &BigFloat, target_precision: u32) -> BigFloat {
    let e_x = match &ax.class {
        Class::Normal { exponent, .. } => *exponent,
        _ => 0,
    };
    // Cancellation guard, the src/math/si.rs:163-178 idiom.
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
    let x = ax
        .round_to_precision(working, RoundingMode::NearestEven)
        .expect("precision >= 1")
        .0;

    // --- seed index M from DLMF 10.19.1 ---------------------------
    let p_bits = i64::from(target_precision) + 64;
    let lp = 64u32; // cheap precision for the M-selection test
    let x_lp = x
        .round_to_precision(lp, RoundingMode::NearestEven)
        .expect("precision >= 1")
        .0;
    let ln2 = BigFloat::try_from_i64_exact(2, lp)
        .expect("precision >= 1")
        .ln(RoundingMode::NearestEven)
        .0;
    let neg_p_ln2 = {
        let p = BigFloat::try_from_i64_exact(p_bits, lp).expect("precision >= 1");
        let (v, _) = p.mul(&ln2, RoundingMode::NearestEven);
        v.negated()
    };
    let satisfies = |big_m: i64| -> bool {
        // lhs = M·(1 + ln(x/(2M))) ; satisfied when lhs ≤ −P·ln2.
        let mm = BigFloat::try_from_i64_exact(big_m, lp).expect("precision >= 1");
        let two_m = BigFloat::try_from_i64_exact(2 * big_m, lp).expect("precision >= 1");
        let (y, _) = x_lp.div(&two_m, RoundingMode::NearestEven);
        let (lny, _) = y.ln(RoundingMode::NearestEven);
        let one = BigFloat::try_from_i64_exact(1, lp).expect("precision >= 1");
        let (s, _) = one.add(&lny, RoundingMode::NearestEven);
        let (lhs, _) = mm.mul(&s, RoundingMode::NearestEven);
        matches!(
            lhs.partial_cmp(&neg_p_ln2).0,
            Some(core::cmp::Ordering::Less | core::cmp::Ordering::Equal)
        )
    };
    let m_floor = i64::from(m) + 2;
    let start = m_floor.max(1i64 << (e_x.max(0) + 2).min(60));
    let cap: i64 = 1 << 24;
    let mut big_m = start;
    while big_m < cap && !satisfies(big_m) {
        big_m = (big_m * 2).min(cap);
    }
    big_m = (big_m + 8).min(cap).max(m_floor);

    // --- backward recurrence (DLMF 10.6.1) + sum rule (10.12.4) ---
    let two = BigFloat::try_from_i64_exact(2, working).expect("precision >= 1");
    let (inv_ax, _) = BigFloat::try_from_i64_exact(1, working)
        .expect("precision >= 1")
        .div(&x, RoundingMode::NearestEven);
    let zero = BigFloat::try_new_zero(Sign::Positive, working).expect("precision >= 1");

    let mut f_hi = zero.clone(); // f_{idx+1}
    let mut f_cur = BigFloat::try_from_i64_exact(1, working).expect("precision >= 1"); // f_M
    let mut s = zero.clone();
    let mut result = zero;
    let mut idx = big_m;
    loop {
        if idx == i64::from(m) {
            result = f_cur.clone();
        }
        if idx % 2 == 0 {
            if idx == 0 {
                let (ns, _) = s.add(&f_cur, RoundingMode::NearestEven);
                s = ns;
            } else {
                let (two_f, _) = f_cur.mul(&two, RoundingMode::NearestEven);
                let (ns, _) = s.add(&two_f, RoundingMode::NearestEven);
                s = ns;
            }
        }
        if idx == 0 {
            break;
        }
        // f_{idx−1} = (2·idx/x)·f_idx − f_{idx+1}.
        let two_idx = BigFloat::try_from_i64_exact(2 * idx, working).expect("precision >= 1");
        let (c1, _) = two_idx.mul(&inv_ax, RoundingMode::NearestEven);
        let (c2, _) = c1.mul(&f_cur, RoundingMode::NearestEven);
        let (f_lo, _) = c2.sub(&f_hi, RoundingMode::NearestEven);
        f_hi = f_cur;
        f_cur = f_lo;
        idx -= 1;
    }
    let (j_m, _) = result.div(&s, RoundingMode::NearestEven);
    j_m
}

/// Binary-exponent boundary at/above which the DLMF 10.17.3 Hankel
/// asymptotic is used instead of Miller recurrence.
///
/// The optimally-truncated Bessel-`J` asymptotic has error of order
/// `e^{−2|x|}`, so reaching `target+64` bits needs
/// `2|x|·log₂e ≥ target+64`, i.e. `|x| ≳ (target+64)·ln2/2 ≈
/// 0.347·(target+64)`. Requiring `2^{e_x} ≥ target+64` is strictly
/// more than enough: a deliberately conservative cut (Miller, always
/// correct if slower, carries everything below it; the crossover is
/// not perf-tuned without a bench, CLAUDE.md). Returns the smallest
/// such `e_x` (the [`super::erf::asymptotic_threshold_exponent`]
/// integer-loop idiom).
pub(super) fn bessel_j_threshold(target_precision: u32) -> i64 {
    let need: u64 = u64::from(target_precision) + 64;
    let mut e: i64 = 0;
    let mut pow_2: u64 = 1;
    while pow_2 < need && e < 90 {
        e += 1;
        pow_2 = pow_2.saturating_mul(2);
    }
    e
}

/// `J_m(ax)` for `m ≥ 0`, large `ax > 0`, via the DLMF 10.17.3
/// Hankel-form asymptotic
///
/// ```text
/// J_m(x) ∼ √(2/(πx))·[cos ω · Σ_{k≥0} (−1)^k a_{2k}(m)/x^{2k}
///                    − sin ω · Σ_{k≥0} (−1)^k a_{2k+1}(m)/x^{2k+1}]
/// ω = x − m·π/2 − π/4
/// ```
///
/// summed to its smallest term (the [`super::si`] / [`super::airy`]
/// `if mag > prev_mag break` optimal-truncation idiom). The
/// coefficients (DLMF 10.17.1) are **derived from the spec, not
/// recalled**: `a_0(m)=1`,
/// `a_k(m) = a_{k−1}(m)·(4m²−(2k−1)²)/(8k)`. The `8k` divisor was
/// cross-checked two independent ways against the primary source
/// (user-authorized `WebFetch` of DLMF §10.17): the closed-form ratio
/// `[(4m²−1²)…(4m²−(2k−1)²)]/(k!·8^k)` and the Pochhammer form
/// `(½−m)_k(½+m)_k/((−2)^k k!)`, agreeing at `k=1,2`
/// (`a_1=(4m²−1)/8`, `a_2=(4m²−1)(4m²−9)/128`) — the 6n Airy
/// `(2k−1)`-divisor defect is the precedent. Folding the explicit
/// `(−1)^k` into the trig assignment, the per-`j` factor on
/// `a_j(m)/x^j` is the period-4 cycle `[+cosω, −sinω, −cosω, +sinω]`
/// for `j ≡ 0,1,2,3 (mod 4)`. Returns the unrounded
/// working-precision value.
fn bessel_j_asymptotic(m: u32, ax: &BigFloat, target_precision: u32) -> BigFloat {
    let working = target_precision
        .saturating_add(64)
        .min(target_precision.saturating_add(512));
    let x = ax
        .round_to_precision(working, RoundingMode::NearestEven)
        .expect("precision >= 1")
        .0;

    // ω = x − m·(π/2) − π/4.
    let half_pi = pi_over_2_at(working);
    let two = BigFloat::try_from_i64_exact(2, working).expect("precision >= 1");
    let (quarter_pi, _) = half_pi.div(&two, RoundingMode::NearestEven);
    let m_big = BigFloat::try_from_i64_exact(i64::from(m), working).expect("precision >= 1");
    let (m_half_pi, _) = m_big.mul(&half_pi, RoundingMode::NearestEven);
    let (w0, _) = x.sub(&m_half_pi, RoundingMode::NearestEven);
    let (omega, _) = w0.sub(&quarter_pi, RoundingMode::NearestEven);
    let (cw, _) = omega.cos(RoundingMode::NearestEven);
    let (sw, _) = omega.sin(RoundingMode::NearestEven);

    // prefactor √(2/(πx)).
    let pi = pi_at(working);
    let (pi_x, _) = pi.mul(&x, RoundingMode::NearestEven);
    let (ratio, _) = two.div(&pi_x, RoundingMode::NearestEven);
    let (prefac, _) = ratio.sqrt(RoundingMode::NearestEven);

    let (inv_x, _) = BigFloat::try_from_i64_exact(1, working)
        .expect("precision >= 1")
        .div(&x, RoundingMode::NearestEven);
    let four_m2: i64 = 4 * i64::from(m) * i64::from(m);

    // j = 0: g_0 = a_0/x^0 = 1, factor +cos ω.
    let mut g = BigFloat::try_from_i64_exact(1, working).expect("precision >= 1");
    let (mut bracket, _) = cw.mul(&g, RoundingMode::NearestEven);
    let mut prev_mag = magnitude(&g);
    let max_iter: i64 = 1 << 22;
    for j in 1..=max_iter {
        // a_j/a_{j−1} = (4m²−(2j−1)²)/(8j); g_j = g_{j−1}·that·(1/x).
        let odd = 2 * j - 1;
        let num = four_m2 - odd * odd;
        let num_b = BigFloat::try_from_i64_exact(num, working).expect("precision >= 1");
        let den = BigFloat::try_from_i64_exact(8 * j, working).expect("precision >= 1");
        let (t1, _) = g.mul(&num_b, RoundingMode::NearestEven);
        let (t2, _) = t1.div(&den, RoundingMode::NearestEven);
        let (cand, _) = t2.mul(&inv_x, RoundingMode::NearestEven);
        let mag = magnitude(&cand);
        if mag > prev_mag {
            break; // smallest term passed: optimal truncation.
        }
        prev_mag = mag;
        g = cand;
        // Period-4 trig/sign cycle on a_j/x^j.
        let contribution = match j % 4 {
            0 => cw.mul(&g, RoundingMode::NearestEven).0,
            1 => sw.mul(&g, RoundingMode::NearestEven).0.negated(),
            2 => cw.mul(&g, RoundingMode::NearestEven).0.negated(),
            _ => sw.mul(&g, RoundingMode::NearestEven).0,
        };
        let (b, _) = bracket.add(&contribution, RoundingMode::NearestEven);
        bracket = b;
    }
    let (result, _) = prefac.mul(&bracket, RoundingMode::NearestEven);
    result
}

/// `J_m(ax)` for `m ≥ 0`, `ax > 0`, normal: the three-regime
/// dispatch (tiny Maclaurin / Miller recurrence / Hankel asymptotic)
/// on the binary exponent of `|x|`. Returns the unrounded
/// working-precision value; [`bessel_j_kernel`] does the final round.
fn bessel_j_eval_normal(m: u32, ax: &BigFloat, target_precision: u32) -> BigFloat {
    let e_x = match &ax.class {
        Class::Normal { exponent, .. } => *exponent,
        _ => 0,
    };
    if e_x <= bessel_j_tiny_threshold() {
        bessel_j_tiny(m, ax, target_precision)
    } else if e_x >= bessel_j_threshold(target_precision) {
        bessel_j_asymptotic(m, ax, target_precision)
    } else {
        bessel_j_miller(m, ax, target_precision)
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use core::cmp::Ordering;

    /// `|v − expected| ≤ 2^(−bits)·|expected|` (the erf/airy test
    /// helper). Reference decimals: `mpmath` `besselj(n, x)` at
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

    fn pj(s: &str, p: u32) -> BigFloat {
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

    #[test]
    fn j0_zero_is_one() {
        let z = BigFloat::try_new_zero(Sign::Positive, 53).unwrap();
        let (r, s) = z.j0(RoundingMode::NearestEven);
        let one = BigFloat::try_from_i64_exact(1, 53).unwrap();
        assert_eq!(r.partial_cmp(&one).0, Some(core::cmp::Ordering::Equal));
        assert!(!s.invalid());
    }

    #[test]
    fn jn_zero_is_zero_for_nonzero_order() {
        for n in [1i32, 2, 3, -1, -4] {
            let z = BigFloat::try_new_zero(Sign::Negative, 53).unwrap();
            let (r, _) = z.jn(n, RoundingMode::NearestEven);
            assert!(r.is_zero(), "J_{n}(0) should be 0");
        }
        // J₀(−0) is still 1.
        let z = BigFloat::try_new_zero(Sign::Negative, 53).unwrap();
        let (r, _) = z.j1(RoundingMode::NearestEven);
        assert!(r.is_zero());
    }

    #[test]
    fn jn_infinity_is_zero() {
        for n in [0i32, 1, 5, -3] {
            for s in [Sign::Positive, Sign::Negative] {
                let inf = BigFloat::try_new_infinity(s, 53).unwrap();
                let (r, st) = inf.jn(n, RoundingMode::NearestEven);
                assert!(r.is_zero() && r.is_sign_positive(), "J_{n}(±∞) = +0");
                assert!(!st.invalid());
            }
        }
    }

    #[test]
    fn jn_quiet_nan_propagates() {
        let q = BigFloat::try_new_quiet_nan(Sign::Positive, 53, &[]).unwrap();
        let (r, _) = q.jn(2, RoundingMode::NearestEven);
        assert!(r.is_quiet_nan());
    }

    #[test]
    fn jn_signaling_nan_raises_invalid() {
        let sn = BigFloat::try_new_signaling_nan(Sign::Positive, 53, &[]).unwrap();
        let (r, st) = sn.j0(RoundingMode::NearestEven);
        assert!(r.is_quiet_nan());
        assert!(st.invalid());
    }

    #[test]
    fn precision_zero_is_rejected() {
        let x = BigFloat::try_from_i64_exact(2, 53).unwrap();
        assert!(x.j0_round(0, RoundingMode::NearestEven).is_err());
        assert!(x.jn_round(3, 0, RoundingMode::NearestEven).is_err());
    }

    // Reference decimals (mpmath besselj, dps = 330).
    const J0_HALF: &str = "0.938469807240812904228404673599712625568926797096821576554705168024483425860925007342101429019355360285361507602817734654932689231019762942429695138281141842049301257720673460185045674354399576294534544926312681150467851674896190709001239711549002378317667252504582765586375476216035484008016120893020593735529906119366282097321853";
    const J1_HALF: &str = "0.242268457674873886383954576141531640800628654437959753506925305893359846884415001326989593873969791657803344005414756202472017227888834970518645480057411645603793281935366884298376839022387228059460058924704874668489145552572966351331879115284416082873175983380917044195789046356490158266097053450699953777033157360461456886474215";
    const J2_HALF: &str = "0.0306040234586826413074136309664139376335878206550174374729960555489559616767349979658569464765238063458518684188412901549553796805355769396448867819485047403658718700207940770084616817351493359433056907725068175234887305353956746963262767495886619531750366810190854111967807092099251490563720929097792213726027233224795454485750060";
    const J0_M3: &str = "0.999999750000015624999565972229003906182183160193172499526204805540539901648033580474850340233905629674281074734141489690396310470296800038290748732800624979880391262764942659652417423213757900654311294042208824834666045649269569698007629258063892397763346928063965449544480320182713897366227139286836254606553867162085459410693160";
    const J1_M3: &str = "0.000499999937500002604166612413195122612841570818899380046283907421092390895847524632454242128659746332061682958182226096175031610079246697447029450880956707613307031934293353893837691652814314828982536713111742084740510031367449642347965650467783818079301813995143619148251711124851716163069098019394037748440805918210813408480938499";
    const J0_1: &str = "0.765197686557966551449717526102663220909274289755325241861547549119278912215272440167180600098915633974929259827603576204084876855208878622462801194635043000621684026314734672618256151707142493145589701261231902194162816319453039821557896720677305893243788545259984194414064569924052713471270496393921490949691240492940714722262068";
    const J1_1: &str = "0.440050585744933515959682203718914913127372301992765251136758171780138222478015547930796592381198254162606413647919983706048911708467231602807674502243279818340473533573708213131284260839112113233616413446407818538441528048475326748299368317348297717016972242008397788984720842587370575708439842293077074095523362060185360415772992";
    const J0_52: &str = "-0.0483837764681979963272877788512034336318110200697737609317815207149902056671166513110063468770413021781191265486249083025067211750888564606148718877471724313832280470263714470734069287730240610886688681917125276807784959537460535981622749526373600359114335843856033351749743077533928906233714292727874379345955862469643895430060651";
    const J1_52: &str = "0.497094102464274038010816276264422242521234969519006818879872428918724175767054746101681708088208184949950367681287344431663438601115729542152407343804500481859007115413763080738014162947311460238019367834974845181187174079407739407691695409734588642748931012495822943117909718024270682469060599197406159844700391230218789016513871";
    const J2_52: &str = "0.446059058439617226735940799862741227648798995684979216035679463849969546280760448192351713347607850138079420693654783847837472055981440094336797762790772816870433739357381911663818259130873229279084362459692403825728235217272245124315631280425030950110578394382261689669302082172809436598619908630712365810355899231139420756217162";
    const J3_52: &str = "0.216600391039113524766689003515963721716843423576959926777214713241227098282161971006081033267964375270976705428560309724876516688454574608786469076660736025133686867558047977924095051662085706608515612100533000939978002268227852791213314638945460877427994418515795760352973613452224416088731254611733625451869047539604284193433588";
    const J5_52: &str = "0.0195016251345032198864719839258657325923572833021588195576200001314944684264085321050958195176569063682702460430875606490952208997359500848158605911632436339076403093441383752087365158839381864518143289605447540367227024564908722474949216907955796208651518536622782722160970749078689023571410731881010474253463600249104719922414505";
    const J0_7: &str = "0.300079270519555596650275377611597144343967812834484670880620264056953898212801481328259344887541899596937883592821867625581557381330208175774486639740372246563582607594400357030566044646090336861616837116925858055928471396810706458871178990455000118268644791345748528853861907433691295627134029255436445616646361113203007773835985";
    const J1_7: &str = "-0.00468282348234583269911380619631299955341731165839360755344159156363220616687699984019030366810691481929325458800252520139678589411743542482690373653588108367146798099094772108200158217760233909994526937213008943088534427443816225196380074619730281227693785250655177104649075782328853870677877698857222322168464233584105059028702811";
    const J2_7: &str = "-0.301417220085940120278593607953400858502087044736882844467317861646563099974766338425456574507001018116735956332251160540266353351078046868582173421607766841898287745020385420196852210982548148033029771223248740750467141189507324245146550632225658064633484177776191892010002123954630878114785108395028509394270544637729022228203708";

    /// Tiny regime (`|x| < 1`, DLMF 10.2.2 path): `J_{0,1,2}` vs
    /// mpmath at `p = 113`.
    #[test]
    fn tiny_regime_matches_mpmath() {
        let p = 113;
        let x = at(1, 2, p); // 0.5
        for (n, want) in [(0i32, J0_HALF), (1, J1_HALF), (2, J2_HALF)] {
            let (r, _) = x.jn(n, RoundingMode::NearestEven);
            assert!(close_at(&r, &pj(want, p), p - 8), "J_{n}(0.5)");
        }
        let x = at(1, 1000, p); // 1e-3
        for (n, want) in [(0i32, J0_M3), (1, J1_M3)] {
            let (r, _) = x.jn(n, RoundingMode::NearestEven);
            assert!(close_at(&r, &pj(want, p), p - 8), "J_{n}(1e-3)");
        }
    }

    /// Miller regime (`|x| ≥ 1`): `J_n` at `x = 2.5` and `x = 7`,
    /// orders `0,1,2,3,5`, vs mpmath at `p = 160`.
    #[test]
    fn miller_regime_matches_mpmath() {
        let p = 160;
        let x = at(5, 2, p); // 2.5
        for (n, want) in [
            (0i32, J0_52),
            (1, J1_52),
            (2, J2_52),
            (3, J3_52),
            (5, J5_52),
        ] {
            let (r, _) = x.jn(n, RoundingMode::NearestEven);
            assert!(close_at(&r, &pj(want, p), p - 12), "J_{n}(2.5)");
        }
        let x = at(7, 1, p);
        for (n, want) in [(0i32, J0_7), (1, J1_7), (2, J2_7)] {
            let (r, _) = x.jn(n, RoundingMode::NearestEven);
            assert!(close_at(&r, &pj(want, p), p - 12), "J_{n}(7)");
        }
    }

    /// Second-term-matters pin (the `derive, don't recall` reflex):
    /// large `M` at `p = 256` and `p = 1024`, validated to `p − 2`
    /// against the 330-digit references.
    #[test]
    fn high_precision_pin() {
        let x = at(5, 2, 256);
        let (r, _) = x.j0(RoundingMode::NearestEven);
        assert!(close_at(&r, &pj(J0_52, 256), 254), "J_0(2.5) p=256");

        let x = at(7, 1, 1024);
        let (r, _) = x.jn(2, RoundingMode::NearestEven);
        assert!(close_at(&r, &pj(J2_7, 1024), 1022), "J_2(7) p=1024");
    }

    /// `J_n(−x) = (−1)^n J_n(x)` and `J_{−n}(x) = (−1)^n J_n(x)`.
    #[test]
    fn parity_in_argument_and_order() {
        let p = 160;
        let pos = at(5, 2, p);
        let neg = at(-5, 2, p);
        // Even order: symmetric.
        let (a, _) = pos.jn(2, RoundingMode::NearestEven);
        let (b, _) = neg.jn(2, RoundingMode::NearestEven);
        assert_eq!(a.partial_cmp(&b).0, Some(Ordering::Equal), "J_2(−x)=J_2(x)");
        // Odd order in argument: antisymmetric.
        let (c, _) = pos.jn(3, RoundingMode::NearestEven);
        let (d, _) = neg.jn(3, RoundingMode::NearestEven);
        assert_eq!(
            c.partial_cmp(&d.negated()).0,
            Some(Ordering::Equal),
            "J_3(−x)=−J_3(x)"
        );
        // Negative odd order: antisymmetric vs positive order.
        let (e, _) = pos.jn(-3, RoundingMode::NearestEven);
        assert!(close_at(&e, &c.negated(), p - 12), "J_(-3) = -J_3");
    }

    /// Recurrence cross-tie `J_{n−1}(x)+J_{n+1}(x) = (2n/x)·J_n(x)`
    /// (DLMF 10.6.1), binding three independently descended orders.
    #[test]
    fn recurrence_spot_check() {
        let p = 160;
        let x = at(5, 2, p); // 2.5
        let (j2, _) = x.jn(2, RoundingMode::NearestEven);
        let (j3, _) = x.jn(3, RoundingMode::NearestEven);
        let (j4, _) = x.jn(4, RoundingMode::NearestEven);
        let (lhs, _) = j2.add(&j4, RoundingMode::NearestEven);
        let six = BigFloat::try_from_i64_exact(6, p).unwrap();
        let (r1, _) = six.mul(&j3, RoundingMode::NearestEven);
        let (rhs, _) = r1.div(&x, RoundingMode::NearestEven);
        assert!(close_at(&lhs, &rhs, p - 8), "J_2+J_4 = (6/x)J_3");
    }

    /// Boundary continuity: at `x = 1` the two regime evaluators
    /// (called directly) agree, and both match mpmath. Pins the
    /// `bessel_j_tiny_threshold` crossover.
    #[test]
    fn tiny_miller_continuity_at_boundary() {
        let p = 160;
        let x = at(1, 1, p);
        let t0 = bessel_j_tiny(0, &x, p);
        let m0 = bessel_j_miller(0, &x, p);
        assert!(close_at(&t0, &m0, p - 12), "tiny vs Miller J_0(1)");
        assert!(close_at(&t0, &pj(J0_1, p), p - 12), "tiny J_0(1)");
        assert!(close_at(&m0, &pj(J0_1, p), p - 12), "Miller J_0(1)");
        let t1 = bessel_j_tiny(1, &x, p);
        let m1 = bessel_j_miller(1, &x, p);
        assert!(close_at(&t1, &m1, p - 12), "tiny vs Miller J_1(1)");
        assert!(close_at(&m1, &pj(J1_1, p), p - 12), "Miller J_1(1)");
    }

    // mpmath besselj at large |x| (dps = 80).
    const J0_200: &str =
        "-0.015437439930565091591922847231344148600368768593123568900377309762820573076602530";
    const J1_200: &str =
        "-0.054304538182378222710670201774468780511705017574432993608494780629064889112381269";
    const J2_200: &str =
        "0.014894394548741309364816145213599460795251718417379238964292361956529924185478717";
    const J3_200: &str =
        "0.054602426073353048897966524678740769727610051942780578387780627868195487596090844";
    const J5_200: &str =
        "-0.055132678944014677613881610657670259235746988617144411252286985593014849978394683";
    const J0_1000: &str =
        "0.024786686152420174561330731115693708786166447133246548414806950013548862444082810";
    const J1_1000: &str =
        "0.0047283119070895239175760719012169162854180242020596368687197215361298685307353971";
    const J2_1000: &str =
        "-0.024777229528605995513495578971891274953595611084842429141069510570476602707021339";
    const J0_256: &str =
        "-0.036653498061713559638146695174670452074468458220974407801753415821915716847265129";
    const J1_256: &str =
        "-0.033884554799704389311128643638197027591868524189466937113016174597544122979420909";

    /// Asymptotic regime (DLMF 10.17.3 path): `J_n` at `x = 200`
    /// (`p = 53`) and `x = 1000` (`p = 113`) vs mpmath. Large enough
    /// that the second and later `a_k(m)` terms materially
    /// contribute, so a recalled coefficient would fail here.
    #[test]
    fn asymptotic_regime_matches_mpmath() {
        let p = 53;
        let x = at(200, 1, p);
        for (n, want) in [
            (0i32, J0_200),
            (1, J1_200),
            (2, J2_200),
            (3, J3_200),
            (5, J5_200),
        ] {
            let (r, _) = x.jn(n, RoundingMode::NearestEven);
            assert!(close_at(&r, &pj(want, p), p - 8), "J_{n}(200)");
        }
        let p = 113;
        let x = at(1000, 1, p);
        for (n, want) in [(0i32, J0_1000), (1, J1_1000), (2, J2_1000)] {
            let (r, _) = x.jn(n, RoundingMode::NearestEven);
            assert!(close_at(&r, &pj(want, p), p - 12), "J_{n}(1000)");
        }
    }

    /// Cross-regime continuity: at `x = 256` the Hankel asymptotic
    /// and Miller recurrence (called directly) agree and both match
    /// mpmath. Pins the `bessel_j_threshold` crossover and the
    /// derived `a_k(m)` recurrence against an independent path.
    #[test]
    fn asymptotic_miller_continuity() {
        let p = 113;
        let x = at(256, 1, p);
        let a0 = bessel_j_asymptotic(0, &x, p);
        let mi0 = bessel_j_miller(0, &x, p);
        assert!(close_at(&a0, &mi0, p - 12), "asymp vs Miller J_0(256)");
        assert!(close_at(&a0, &pj(J0_256, p), p - 12), "asymp J_0(256)");
        assert!(close_at(&mi0, &pj(J0_256, p), p - 12), "Miller J_0(256)");
        let a1 = bessel_j_asymptotic(1, &x, p);
        let mi1 = bessel_j_miller(1, &x, p);
        assert!(close_at(&a1, &mi1, p - 12), "asymp vs Miller J_1(256)");
        assert!(close_at(&a1, &pj(J1_256, p), p - 12), "asymp J_1(256)");
    }
}
