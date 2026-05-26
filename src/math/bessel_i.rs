//! Modified Bessel functions of the first kind `I0`, `I1`, `In`
//! (DLMF Chapter 10): integer order, real argument. Like
//! [`super::bessel_j`], `I` is entire on the real line, with no poles
//! and no domain restriction; unlike the oscillatory ordinary Bessel
//! functions it grows monotonically and diverges at infinity.
//!
//! Order and sign reduce before evaluation. The order parity is
//! **even with no sign**: `I₋ₙ(x) = Iₙ(x)` (DLMF 10.27.1), in
//! contrast to `J`/`Y` where `𝒞₋ₙ = (−1)ⁿ 𝒞ₙ`. The argument parity
//! matches `J`: `Iₙ(−x) = (−1)ⁿ Iₙ(x)` (from the `(x/2)ⁿ` prefactor
//! of the DLMF 10.25.2 series, whose remaining sum is even in `x`).
//! So the kernel evaluates `I_m(|x|)` for `m = |n| ≥ 0` and negates
//! exactly when `m` is odd and `x < 0`; the order sign never
//! contributes.
//!
//! Special cases:
//!
//! - `I₀(±0) = 1`, `Iₙ(±0) = 0` for `n ≠ 0` (exact, DLMF 10.30.1
//!   `Iν(z) ∼ (½z)ν/Γ(ν+1)`; entire, both zero signs alike).
//! - `Iₙ(+∞) = +∞`; `Iₙ(−∞) = (−1)ⁿ·∞` (so `−∞` for odd `n`). This
//!   is a **genuine infinite limit** (`I` grows like
//!   `eˣ/√(2πx) → ∞`, DLMF 10.30.4), `Status::OK`, the
//!   `exp(+∞) = +∞` precedent — explicitly **not** the
//!   decaying-envelope convention of `J`/`Y`/Airy (which covers a
//!   bounded non-converging oscillation, a different situation).
//! - `Iₙ(NaN) = NaN`; `sNaN` raises `INVALID`.
//!
//! Three regimes will dispatch on the binary exponent of `|x|` (the
//! [`super::bessel_j`] template, since `I` is the *recessive*
//! solution in order just as `J` is): tiny all-positive Maclaurin
//! (DLMF 10.25.2, slice 6q.2), Miller backward recurrence normalised
//! by the DLMF 10.35.5 sum rule `eˣ = I₀ + 2Σ_{k≥1} Iₖ` (slice 6q.3),
//! and the DLMF 10.40.1 asymptotic reusing the ADR-0023 `aₖ(ν)`
//! coefficients (slice 6q.4). ADR-0025 records the design and the
//! DLMF provenance.

use super::pi_at;
use super::ziv::ziv_round;
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
    /// `I₀(self)` rounded under `mode` to `self.precision`.
    #[must_use]
    pub fn i0(&self, mode: RoundingMode) -> (Self, Status) {
        self.i0_round(self.precision, mode)
            .expect("self.precision >= 1 by invariant")
    }

    /// `I₀(self)` with explicit result precision.
    ///
    /// Correctly rounded under every IEEE 754-2019 rounding mode
    /// via the shared `crate::math::ziv::ziv_round` driver (slice
    /// p1.33, ADR-0038).
    pub fn i0_round(
        &self,
        target_precision: u32,
        mode: RoundingMode,
    ) -> Result<(Self, Status), BuildError> {
        if target_precision == 0 {
            return Err(BuildError::PrecisionZero);
        }
        Ok(bessel_i_kernel(0, self, target_precision, mode))
    }

    /// `I₁(self)` rounded under `mode` to `self.precision`.
    #[must_use]
    pub fn i1(&self, mode: RoundingMode) -> (Self, Status) {
        self.i1_round(self.precision, mode)
            .expect("self.precision >= 1 by invariant")
    }

    /// `I₁(self)` with explicit result precision.
    ///
    /// Correctly rounded under every IEEE 754-2019 rounding mode
    /// via the shared `crate::math::ziv::ziv_round` driver (slice
    /// p1.33, ADR-0038).
    pub fn i1_round(
        &self,
        target_precision: u32,
        mode: RoundingMode,
    ) -> Result<(Self, Status), BuildError> {
        if target_precision == 0 {
            return Err(BuildError::PrecisionZero);
        }
        Ok(bessel_i_kernel(1, self, target_precision, mode))
    }

    /// `Iₙ(self)` for integer order `n` (any `i32`, including
    /// negative: `I₋ₙ = Iₙ`), rounded under `mode` to
    /// `self.precision`.
    #[must_use]
    pub fn in_(&self, n: i32, mode: RoundingMode) -> (Self, Status) {
        self.in_round(n, self.precision, mode)
            .expect("self.precision >= 1 by invariant")
    }

    /// `Iₙ(self)` with explicit result precision.
    ///
    /// Correctly rounded under every IEEE 754-2019 rounding mode
    /// via the shared `crate::math::ziv::ziv_round` driver (slice
    /// p1.33, ADR-0038).
    pub fn in_round(
        &self,
        n: i32,
        target_precision: u32,
        mode: RoundingMode,
    ) -> Result<(Self, Status), BuildError> {
        if target_precision == 0 {
            return Err(BuildError::PrecisionZero);
        }
        Ok(bessel_i_kernel(n, self, target_precision, mode))
    }
}

#[cfg(feature = "fixed")]
impl<const PREC: u32> FixedFloat<PREC>
where
    [(); limbs_for(PREC)]:,
{
    /// `I₀(self)` for `FixedFloat`. Delegates to [`BigFloat::i0`].
    #[must_use]
    pub fn i0(&self, mode: RoundingMode) -> (Self, Status) {
        let (big, status) = self.to_big().i0(mode);
        (
            Self::try_from_big_exact(big).expect("precision matches"),
            status,
        )
    }

    /// `I₁(self)` for `FixedFloat`. Delegates to [`BigFloat::i1`].
    #[must_use]
    pub fn i1(&self, mode: RoundingMode) -> (Self, Status) {
        let (big, status) = self.to_big().i1(mode);
        (
            Self::try_from_big_exact(big).expect("precision matches"),
            status,
        )
    }

    /// `Iₙ(self)` for `FixedFloat`. Delegates to [`BigFloat::in_`].
    #[must_use]
    pub fn in_(&self, n: i32, mode: RoundingMode) -> (Self, Status) {
        let (big, status) = self.to_big().in_(n, mode);
        (
            Self::try_from_big_exact(big).expect("precision matches"),
            status,
        )
    }
}

/// `Iₙ(x)` for integer order `n` and real `x`.
///
/// Special values are handled directly per the module-level domain
/// table; for a normal argument the order is reduced to `I_m(|x|)`,
/// `m = |n|`, with one argument-parity sign (`Iₙ(−x) = (−1)ⁿ Iₙ(x)`;
/// `I₋ₙ = Iₙ` adds no sign), then the regime evaluator runs.
fn bessel_i_kernel(
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
            // I₀(±0) = 1, Iₙ(±0) = 0 for n ≠ 0 (DLMF 10.30.1); exact,
            // both zero signs alike (I is entire).
            let value = if m == 0 {
                BigFloat::try_from_i64_exact(1, target_precision).expect("precision >= 1")
            } else {
                BigFloat::try_new_zero(Sign::Positive, target_precision).expect("precision >= 1")
            };
            (value, Status::OK)
        }
        Class::Infinity { sign } => {
            // Iₙ(+∞) = +∞; Iₙ(−∞) = (−1)ⁿ·∞ (the argument parity, so
            // −∞ for odd m). A genuine infinite limit (DLMF 10.30.4
            // eˣ/√(2πx) → ∞), Status::OK, the exp(+∞) precedent; not
            // the decaying-envelope convention.
            let neg = matches!(sign, Sign::Negative) && (m % 2 == 1);
            let result_sign = if neg { Sign::Negative } else { Sign::Positive };
            let inf =
                BigFloat::try_new_infinity(result_sign, target_precision).expect("precision >= 1");
            (inf, Status::OK)
        }
        Class::Normal { .. } => {
            // I is entire: no domain restriction. Reduce to I_m(|x|)
            // with one argument-parity sign. Iₙ(−x) = (−1)ⁿ Iₙ(x);
            // I₋ₙ(x) = Iₙ(x) (no order sign). Negate exactly when m is
            // odd and x < 0. Binary parity sign pinned outside Ziv.
            let negate = (m % 2 == 1) && x.is_sign_negative();
            let ax = x.abs();

            // Regime decision pinned from target_precision so it does
            // not flip across Ziv retries (slice p1.4 erf precedent).
            // Three regimes: tiny (|x| < 1, fixed cut), miller (default),
            // asymptotic (|x| >= bessel_j_threshold(target_precision)).
            // The tiny vs miller boundary is precision-independent (|x|
            // < 1); the miller vs asymptotic boundary depends on
            // target_precision via bessel_j_threshold.
            let e_x = match &ax.class {
                Class::Normal { exponent, .. } => *exponent,
                _ => 0,
            };
            let use_tiny = e_x <= bessel_i_tiny_threshold();
            let use_asymptotic =
                !use_tiny && e_x >= super::bessel_j::bessel_j_threshold(target_precision);

            let (result, status) = ziv_round(
                |w| {
                    let value = bessel_i_eval_at_w(m, &ax, w, use_tiny, use_asymptotic);
                    if negate {
                        value.negated()
                    } else {
                        value
                    }
                },
                target_precision,
                mode,
            );
            auto_raise(status);
            (result, status)
        }
    }
}

/// Same shape as [`bessel_i_eval_normal`] but takes the pinned
/// three-regime flags (`use_tiny` / `use_asymptotic`) so the regime
/// does not flip across Ziv retries. Exactly one of `use_tiny` and
/// `use_asymptotic` may be true; both false means the Miller default.
fn bessel_i_eval_at_w(
    m: u32,
    ax: &BigFloat,
    target_precision: u32,
    use_tiny: bool,
    use_asymptotic: bool,
) -> BigFloat {
    if use_tiny {
        bessel_i_tiny(m, ax, target_precision)
    } else if use_asymptotic {
        bessel_i_asymptotic(m, ax, target_precision)
    } else {
        bessel_i_miller(m, ax, target_precision)
    }
}

/// Binary-exponent boundary below which the convergent Maclaurin
/// series ([`bessel_i_tiny`]) is used instead of Miller recurrence.
///
/// `e_x ≤ −1` ⇔ `|x| < 1`. The tiny regime exists only to keep the
/// `2k/x` recurrence away from `x → 0`; it is not a tuned crossover
/// (CLAUDE.md: no perf machinery without a bench). It mirrors
/// [`super::bessel_j::bessel_j_tiny_threshold`] exactly: `I` shares
/// `J`'s recurrence-near-zero hazard. Miller carries everything
/// `|x| ≥ 1` until slice 6q.4 adds the large-`|x|` asymptotic upper
/// cut; continuity across the boundary is pinned by a unit test.
fn bessel_i_tiny_threshold() -> i64 {
    -1
}

/// `I_m(ax)` for `m ≥ 0`, normal `ax > 0`: the three-regime dispatch
/// on the binary exponent of `|x|` (tiny convergent Maclaurin /
/// Miller backward recurrence / DLMF 10.40.1 asymptotic). Returns
/// the unrounded working-precision value; [`bessel_i_kernel`] does
/// the single final round.
///
/// The asymptotic upper cut **reuses** [`super::bessel_j`]'s
/// `bessel_j_threshold`, and the reuse is *derived*, not reflexive
/// (the CLAUDE.md "derive the cut" reflex; the 6n precedent makes it
/// load-bearing): the quantity that controls accuracy is the
/// optimal-truncation **relative** error of the shared `a_k(ν)`
/// divergent series, which is `O(e^{−2x})` — *identical* to the
/// ordinary-Bessel 10.17.3/10.17.4 series, since `I`/`K` reuse the
/// same `a_k(ν)` (DLMF 10.40 ≡ §10.17(i), ADR-0023). The `eˣ/√(2πx)`
/// prefactor is computed exactly and does not enter the relative
/// error, so `bessel_j_threshold`'s conservative
/// `2^{e_x} ≥ target+64` cut is strictly more than enough for `I`
/// too. Miller (always correct, slower) carries everything below it.
/// Binary exponent of `v`, or `i64::MIN`/`i64::MAX` for zero /
/// non-finite (the [`super::bessel_y`] `magnitude` idiom; the
/// `bessel_j` copy is private, so `I` carries its own).
fn magnitude(v: &BigFloat) -> i64 {
    match &v.class {
        Class::Normal { exponent, .. } => *exponent,
        Class::Zero { .. } => i64::MIN,
        _ => i64::MAX,
    }
}

/// `true` once `term` has fallen `working + 8` bits below the running
/// `sum`, so further terms cannot perturb the rounded result (the
/// [`super::bessel_y`] `negligible` idiom). `I`'s 10.25.2 series is
/// all-positive (monotone partial sums), so this is a pure
/// tail-length cutoff with no cancellation concern.
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

/// `I_m(x)` for `m ≥ 0`, `x > 0`, via the DLMF 10.25.2 convergent
/// Maclaurin series
///
/// ```text
/// I_m(x) = (x/2)^m · Σ_{k≥0} (x/2)^{2k} / (k!·(m+k)!)
/// ```
///
/// (integer order ⇒ `Γ(m+k+1) = (m+k)!`; entire, converges for all
/// `x`). Unlike the ordinary-Bessel `J` Maclaurin (DLMF 10.2.2,
/// alternating `(−1)^k`), **every term here is positive**: the
/// modified series has no sign alternation, so the partial sums are
/// monotone and there is **no cancellation**. No `≈ x·log₂e`
/// working-precision boost is needed (contrast
/// [`super::bessel_j::bessel_j_miller`]'s cancellation guard), only
/// the `+64` base. Carried as a term recurrence (the
/// [`super::bessel_j`] `bessel_j_tiny` idiom with the negation
/// removed): with `t_0 = (x/2)^m / m!`,
/// `t_k = t_{k−1} · (x/2)² / (k·(m+k))`. Hand-check (derive, don't
/// recall): `I_0` → `t_0 = 1, t_1 = (x/2)², t_2 = (x/2)⁴/4`, i.e.
/// `Σ (x/2)^{2k}/(k!)²`; `I_1` → `t_0 = x/2, t_1 = (x/2)³/2`. Both
/// match the standard modified-Bessel expansions. Returns the
/// unrounded working-precision value.
fn bessel_i_tiny(m: u32, ax: &BigFloat, target_precision: u32) -> BigFloat {
    let working = target_precision.saturating_add(64);
    let x = ax
        .round_to_precision(working, RoundingMode::NearestEven)
        .expect("precision >= 1")
        .0;
    let two = ci(2, working);
    let (half, _) = x.div(&two, RoundingMode::NearestEven);
    let (half_sq, _) = half.mul(&half, RoundingMode::NearestEven);

    // t_0 = (x/2)^m / m!  (in-module recurrence-from-a-base-term, not
    // pow.rs::pow_int; m is small in the tiny regime).
    let mut term = ci(1, working);
    for _ in 0..m {
        let (t, _) = term.mul(&half, RoundingMode::NearestEven);
        term = t;
    }
    for j in 1..=i64::from(m) {
        let (t, _) = term.div(&ci(j, working), RoundingMode::NearestEven);
        term = t;
    }

    let mut sum = term.clone();
    let max_iter: i64 = 1 << 22;
    for k in 1..=max_iter {
        let (t1, _) = term.mul(&half_sq, RoundingMode::NearestEven);
        let (t2, _) = t1.div(&ci(k, working), RoundingMode::NearestEven);
        let (t3, _) = t2.div(&ci(k + i64::from(m), working), RoundingMode::NearestEven);
        // No negation: the modified series (DLMF 10.25.2) is
        // all-positive, vs J's alternating DLMF 10.2.2.
        term = t3;
        let (s, _) = sum.add(&term, RoundingMode::NearestEven);
        sum = s;
        if negligible(&term, &sum, working) {
            break;
        }
    }
    sum
}

/// `I_m(ax)` for `m ≥ 0`, `ax ≥ 1`, via Miller backward recurrence
/// (DLMF 10.29.1) with sum-rule normalisation (DLMF 10.35.5).
///
/// `I` is the **recessive** (minimal) solution of the
/// modified-Bessel three-term recurrence *in order* (DLMF 10.30.1
/// `Iν(z) ∼ (½z)ν/Γ(ν+1) → 0` as `ν → ∞` at fixed `z`), exactly the
/// role `J` plays for the ordinary recurrence, so the same Miller
/// backward descent applies. The DLMF 10.29.1 relation is
/// `𝒵_{ν−1}(z) − 𝒵_{ν+1}(z) = (2ν/z)·𝒵_ν(z)` — note the **minus**,
/// where the ordinary-Bessel 10.6.1 has a plus. Solving for the
/// lower order gives the descent
///
/// ```text
/// f_{k−1} = (2k/x)·f_k + f_{k+1}
/// ```
///
/// with a **plus** `f_{k+1}`, the single sign change from
/// [`super::bessel_j::bessel_j_miller`]'s `− f_{k+1}`. Started at a
/// high seed index `M` with `f_{M+1}=0`, `f_M=1`, the descent
/// converges to a fixed multiple `c·I_k(x)` of the recessive
/// solution. Because every term of the descent is **added** (the
/// `f` values are all positive and grow as the index falls), the
/// recurrence has **no subtractive cancellation** — unlike `J`'s,
/// whose `− f_{k+1}` does. The normalising constant is the DLMF
/// 10.35.5 identity `eˣ = I_0(x) + 2·Σ_{k≥1} I_k(x)` (every order
/// `k ≥ 1`, *not* the even-only `J` sum rule DLMF 10.12.4; it is the
/// modified generating function DLMF 10.35.1
/// `exp(½x(t+1/t)) = Σ I_n(x) tⁿ` at `t = 1`, using `I_{−n}=I_n`).
/// So `S = f_0 + 2·Σ_{k≥1} f_k = c·eˣ`, giving
/// `I_m(x) = f_m / c = f_m·eˣ / S`, every order from one descent and
/// one `exp(x)`.
///
/// `M` is derived, not guessed: the recessive modified-Bessel
/// solution has the same leading large-order decay
/// `I_M(x) ∼ (1/√(2πM))·(ex/(2M))^M` as the ordinary recessive
/// solution (the `J` form DLMF 10.19.1, pinned in ADR-0023). The
/// prefactor only shrinks the bound, so requiring
/// `(ex/(2M))^M < 2^{−P}`, `P = target+64`, i.e. in natural logs
/// `M·(1 + ln(x/(2M))) < −P·ln2`, is the **same conservative seed
/// criterion** as `bessel_j_miller`; it is reused verbatim (an
/// exponential search at low working precision plus a fixed step
/// guard; overshoot ≤ 2× is the deliberate robustness/cost trade,
/// CLAUDE.md). Working precision takes the same `≈ |x|·log₂e`
/// boost as `bessel_j_miller` plus the `+64` base: **not** for
/// cancellation in the recurrence (`I` has none) but to carry the
/// wide `f`-magnitude dynamic range and the `eˣ` normalisation
/// composition conservatively (CLAUDE.md: conservative, not
/// perf-tuned without a bench). Returns the unrounded
/// working-precision value.
fn bessel_i_miller(m: u32, ax: &BigFloat, target_precision: u32) -> BigFloat {
    let e_x = match &ax.class {
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
    let x = ax
        .round_to_precision(working, RoundingMode::NearestEven)
        .expect("precision >= 1")
        .0;

    // --- seed index M (same criterion as bessel_j_miller) ---------
    let p_bits = i64::from(target_precision) + 64;
    let lp = 64u32;
    let x_lp = x
        .round_to_precision(lp, RoundingMode::NearestEven)
        .expect("precision >= 1")
        .0;
    let ln2 = ci(2, lp).ln(RoundingMode::NearestEven).0;
    let neg_p_ln2 = {
        let (v, _) = ci(p_bits, lp).mul(&ln2, RoundingMode::NearestEven);
        v.negated()
    };
    let satisfies = |big_m: i64| -> bool {
        // lhs = M·(1 + ln(x/(2M))) ; satisfied when lhs ≤ −P·ln2.
        let mm = ci(big_m, lp);
        let two_m = ci(2 * big_m, lp);
        let (y, _) = x_lp.div(&two_m, RoundingMode::NearestEven);
        let (lny, _) = y.ln(RoundingMode::NearestEven);
        let (s, _) = ci(1, lp).add(&lny, RoundingMode::NearestEven);
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

    // --- backward recurrence (DLMF 10.29.1) + sum rule (10.35.5) --
    let (inv_ax, _) = ci(1, working).div(&x, RoundingMode::NearestEven);
    let zero = BigFloat::try_new_zero(Sign::Positive, working).expect("precision >= 1");

    let mut f_hi = zero.clone(); // f_{idx+1}
    let mut f_cur = ci(1, working); // f_M
    let mut s = zero.clone();
    let mut result = zero;
    let mut idx = big_m;
    loop {
        if idx == i64::from(m) {
            result = f_cur.clone();
        }
        // DLMF 10.35.5: every order contributes (f_0 once, f_{k≥1}
        // doubled), unlike J's even-only DLMF 10.12.4 sum rule.
        if idx == 0 {
            let (ns, _) = s.add(&f_cur, RoundingMode::NearestEven);
            s = ns;
            break;
        }
        let two_f = f_cur.mul(&ci(2, working), RoundingMode::NearestEven).0;
        let (ns, _) = s.add(&two_f, RoundingMode::NearestEven);
        s = ns;
        // f_{idx−1} = (2·idx/x)·f_idx + f_{idx+1}  (PLUS: DLMF
        // 10.29.1, vs bessel_j_miller's − f_{idx+1}).
        let two_idx = ci(2 * idx, working);
        let (c1, _) = two_idx.mul(&inv_ax, RoundingMode::NearestEven);
        let (c2, _) = c1.mul(&f_cur, RoundingMode::NearestEven);
        let (f_lo, _) = c2.add(&f_hi, RoundingMode::NearestEven);
        f_hi = f_cur;
        f_cur = f_lo;
        idx -= 1;
    }
    // I_m = f_m·eˣ / S  (S = c·eˣ; result = f_m = c·I_m).
    let (ex, _) = x.exp(RoundingMode::NearestEven);
    let (num, _) = result.mul(&ex, RoundingMode::NearestEven);
    let (i_m, _) = num.div(&s, RoundingMode::NearestEven);
    i_m
}

/// `I_m(ax)` for `m ≥ 0`, large `ax > 0`, via the DLMF 10.40.1
/// large-argument asymptotic
///
/// ```text
/// I_m(x) ∼ eˣ/√(2πx) · Σ_{k≥0} (−1)^k a_k(m)/x^k
/// ```
///
/// summed to its smallest term (the [`super::bessel_j`] /
/// [`super::si`] optimal-truncation idiom). There is **no trig**
/// (no `ω`, no `sin`/`cos`): `I` grows monotonically, so the
/// asymptotic is a single real exponential `eˣ/√(2πx)` times a
/// `1/x` power series — markedly simpler than `J`/`Y`'s
/// `√(2/πx)·[sinω·ΣP + cosω·ΣQ]`. The coefficients `a_k(m)` are the
/// **same** as the ordinary-Bessel 10.17.1 sequence
/// (`a_0=1`, `a_k = a_{k−1}(4m²−(2k−1)²)/(8k)`), confirmed by
/// DLMF 10.40 ≡ §10.17(i); they were already derived and
/// Pochhammer-cross-pinned at `k=1,2` in ADR-0023, so 6q reuses that
/// pin rather than re-deriving (the 6n `(2k−1)`-divisor defect is
/// the precedent). The running `g_j` carries `a_j/x^j` with the
/// recurrence's own sign; DLMF 10.40.1's explicit `(−1)^j` is
/// applied on top (the [`super::bessel_y`] sign-folding idiom, here
/// trivial since there is no trig). DLMF 10.40.5 notes a companion
/// `e^{−x}` term for `I` off the real axis; on the positive real
/// axis it is `O(e^{−2x})` relative to the `eˣ` lead, below the
/// optimal-truncation error, so the single series suffices. Returns
/// the unrounded working-precision value.
fn bessel_i_asymptotic(n: u32, ax: &BigFloat, target_precision: u32) -> BigFloat {
    let working = target_precision
        .saturating_add(64)
        .min(target_precision.saturating_add(512));
    let x = ax
        .round_to_precision(working, RoundingMode::NearestEven)
        .expect("precision >= 1")
        .0;

    // prefactor eˣ / √(2πx).
    let (ex, _) = x.exp(RoundingMode::NearestEven);
    let pi = pi_at(working);
    let two = ci(2, working);
    let (two_pi, _) = two.mul(&pi, RoundingMode::NearestEven);
    let (two_pi_x, _) = two_pi.mul(&x, RoundingMode::NearestEven);
    let (sqrt_tpx, _) = two_pi_x.sqrt(RoundingMode::NearestEven);
    let (prefac, _) = ex.div(&sqrt_tpx, RoundingMode::NearestEven);

    let (inv_x, _) = ci(1, working).div(&x, RoundingMode::NearestEven);
    let four_n2: i64 = 4 * i64::from(n) * i64::from(n);

    // j = 0: g_0 = a_0/x^0 = 1, explicit (−1)^0 = +1.
    let mut g = ci(1, working);
    let mut bracket = g.clone();
    let mut prev_mag = magnitude(&g);
    let max_iter: i64 = 1 << 22;
    for j in 1..=max_iter {
        // a_j/a_{j−1} = (4n²−(2j−1)²)/(8j); g_j = g_{j−1}·that·(1/x).
        let odd = 2 * j - 1;
        let num = four_n2 - odd * odd;
        let (t1, _) = g.mul(&ci(num, working), RoundingMode::NearestEven);
        let (t2, _) = t1.div(&ci(8 * j, working), RoundingMode::NearestEven);
        let (cand, _) = t2.mul(&inv_x, RoundingMode::NearestEven);
        let mag = magnitude(&cand);
        if mag > prev_mag {
            break; // smallest term passed: optimal truncation.
        }
        prev_mag = mag;
        g = cand;
        // DLMF 10.40.1 explicit (−1)^j on a_j/x^j (no trig cycle).
        let contribution = if j % 2 == 0 { g.clone() } else { g.negated() };
        let (b, _) = bracket.add(&contribution, RoundingMode::NearestEven);
        bracket = b;
    }
    let (result, _) = prefac.mul(&bracket, RoundingMode::NearestEven);
    result
}

#[cfg(test)]
mod tests {
    use super::*;
    use core::cmp::Ordering;

    /// `|v − expected| ≤ 2^(−bits)·|expected|` (the `bessel_j`/
    /// `bessel_y` test helper). Reference decimals: `mpmath`
    /// `besseli(n, x)` at `mp.dps = 330`
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

    // mpmath besseli(n, x), dps = 330 (truncated to fit p = 113).
    const I0_HALF: &str = "1.0634833707413235192631844154453565293295231748211049891695720746879267185056918544345638318413338632140262062685854965842115745278400339764545357831903650048756518532239906976239160110977076739869409855388779460377122111653287321663731483128995723917359295616246216272807698482061839385267638373181295164392770720209276";
    const I1_HALF: &str = "0.25789430539089631636247965952320963418774314964079457273094519087056586338943968672536228301582848151057763899364628714853943592929796482473494201469626935091697781469896391107400830696094154538084267866284323038046698724066499433700734845701239890709248465780426768018991327156274423210473243761778181862075267797104367";
    const I2_HALF: &str = "0.031906149177738253813265777352517992578550576257926698245791311205663264947933107533114699778019937171715650294000347990053830810648174677514767724405287601207740594428135053327882783253941492463570270887505024515844262202668754818343754484849976763365990930407550906521116761955207010107834086847002241956266360136752864";
    const I0_M3: &str = "1.0000002500000156250004340277845594618733723963042836154429310157005298814807852214627858798978125575203137027779731355520389368525875727930175622870915975274808774949175131459213950006210729609766529579761255440783007377345436645584398399793814634997573562954218419951217288510387549623390874312651559458644797876109271";
    const I1_M3: &str = "0.00050000006250000260416672092013956705729731807005678745399166577219035323182580872495967316870712609911518156781046862548209247659471322803652228085818111311007268387292938807391088260084818729361434529510384540765354317903680136378764049646696799553657051450026313904110998344356862769536348026706961849807071345642209947";
    const I3_09: &str = "0.015972113178804609529461790187571543905141741067083162419700313547801708401962198772146567393232180854094617058760119899142348716140103805064328880191029242076916568624283819842047920532285312768835653662378464614201343534994513683166896581076317023213369599754384556341742335404928764895724006073215211346429966867932024";

    /// Tiny / convergent-Maclaurin regime (DLMF 10.25.2 path):
    /// `I_{0,1,2}(0.5)`, `I_{0,1}(1e-3)` (the
    /// monotone-near-the-origin shape), and `I_3(0.9)` where the
    /// second and later `t_k` terms materially contribute, vs mpmath
    /// at `p = 113`. Because the series is all-positive a coefficient
    /// error cannot hide in cancellation (the `derive, don't recall`
    /// reflex: a leading-term-only check would not catch a wrong
    /// `k·(m+k)` divisor).
    #[test]
    fn tiny_regime_matches_mpmath() {
        let p = 113;
        let x = at(1, 2, p); // 0.5
        for (n, want) in [(0i32, I0_HALF), (1, I1_HALF), (2, I2_HALF)] {
            let (r, _) = x.in_(n, RoundingMode::NearestEven);
            assert!(close_at(&r, &py(want, p), p - 8), "I_{n}(0.5)");
        }
        let x = at(1, 1000, p); // 1e-3
        for (n, want) in [(0i32, I0_M3), (1, I1_M3)] {
            let (r, _) = x.in_(n, RoundingMode::NearestEven);
            assert!(close_at(&r, &py(want, p), p - 8), "I_{n}(1e-3)");
        }
        let x = at(9, 10, p); // 0.9
        let (r, _) = x.in_(3, RoundingMode::NearestEven);
        assert!(close_at(&r, &py(I3_09, p), p - 8), "I_3(0.9)");
    }

    // mpmath besseli(n, x), dps = 330.
    const I0_52: &str = "3.2898391440501230357059082299060560261118015753483941612552870534405381281069497332407051341783812221934084659616546791811212606383909101136377412978747376249803972150614381637175830861753964515785207689723598032584993035238239289703790819954960397879528423835134676611617142099196745844735391351097153851303617722988494";
    const I1_52: &str = "2.5167162452886984415281917481223776723889473033969492189945395536862123621998489805478611328623654879854309633948532751443897799743868030305513678468672073092043537834766850107277309570952811649702000233303932592586843992778535656938613510512689669963141673599204696295648188684564962922311601663177882163200383719996789";
    const I2_52: &str = "1.2764661478191642824833548314081538882006437326308347860596554104915682383470705488024162278884888318050636952457720590656094366588814676891966470203809717776169141882800901551353983204991715196023607503080451958515517841015410764152900011544808661909015084955770919575098591151544775506886110020554848120743310746991063";
    const I3_52: &str = "0.47437040877803558955482401786933145126791733118761356129909089689970318084453610246399516824078335709732905100161798063941468132017645472783673261425765246501729108222854076251109364429660673360642282283752094589620154471538784342939734920409958109087175376699712249754904428420933221112938256302901251700110865248110878";
    const I5_52: &str = "0.032843475172023213389137014599704554763462490289814396685211671516405247019947613219746131086835277828612337907573482940169236550710930432193568626537313707975961191247445322163018206897197584976195701243937183654063699004635036438240987397245592057881995511688329014694155218442675430399485440514277248931763664498884024";
    const I0_7: &str = "168.59390851028969885732662718750084037652267923453171419319405566855416412467567826523169181672021486929560740073878851136031864873853514579691924920051719620176847781535776336671740853295939586010914436537083026003191102696772561838030422796757453778910033451199739338748330017676012004587611339859666944366147702231554";
    const I1_7: &str = "156.03909286995545346239058066071115563003105204154494317062487461487033616940013352003661452334856107081952187307407510752031015044268041681675441168719517175559714279703872634857336249277754166744896266292847517491713297098270599426677882319174023951343080424877167156982090527896085153837147505053744589452455698324444";
    const I2_7: &str = "124.01131054744528358235788985586908162508523579409030185872980577859121093341849725950694481004919742049002972271762419492594432004062645527784656014703286141445500844477527012426787639216581252655229789024840878148415874954409533430408170705564875507097724758377691579610589866848559103491283481272882775951160359853142";
    const I0_1: &str = "1.2660658777520083355982446252147175376076703113549622068081353312135750161227754703948183571472801018710361346890561387866044362393033515852308035320408747903825409142441362953734395725600995545926438383227198691560588492206548377304558882075871377934066816603165868599396338832481107108607514361134215851892896832534291";
    const I1_1: &str = "0.56515910399248502720769602760986330732889962162109200948029448947925564096437113409266499776681441006467788605552630267685763768491717981204113120812102680026939328297640534978443200181007060344916618001109074264750052022176527680673939532145235807641855528642638633569454900860852508478096649245186461937460814440257591";

    /// Miller backward-recurrence regime (`|x| ≥ 1`, DLMF 10.29.1
    /// descent normalised by the DLMF 10.35.5 `eˣ` sum rule):
    /// `I_n(2.5)` orders `0,1,2,3,5` and `I_n(7)` orders `0,1,2`,
    /// where the recurrence depth and the `eˣ` normalisation both
    /// materially matter, vs mpmath at `p = 160`. A wrong recurrence
    /// sign (the `+ f_{k+1}` vs `− f_{k+1}` fork) or an even-only
    /// sum rule fails here.
    #[test]
    fn miller_regime_matches_mpmath() {
        let p = 160;
        let x = at(5, 2, p); // 2.5
        for (n, want) in [
            (0i32, I0_52),
            (1, I1_52),
            (2, I2_52),
            (3, I3_52),
            (5, I5_52),
        ] {
            let (r, _) = x.in_(n, RoundingMode::NearestEven);
            assert!(close_at(&r, &py(want, p), p - 12), "I_{n}(2.5)");
        }
        let x = at(7, 1, p);
        for (n, want) in [(0i32, I0_7), (1, I1_7), (2, I2_7)] {
            let (r, _) = x.in_(n, RoundingMode::NearestEven);
            assert!(close_at(&r, &py(want, p), p - 12), "I_{n}(7)");
        }
    }

    /// Boundary continuity: at `x = 1` the two regime evaluators
    /// (called directly) agree, and both match mpmath. Pins the
    /// `bessel_i_tiny_threshold` crossover and the all-positive
    /// 10.25.2 series against the independent Miller path.
    #[test]
    fn tiny_miller_continuity_at_boundary() {
        let p = 160;
        let x = at(1, 1, p);
        let t0 = bessel_i_tiny(0, &x, p);
        let m0 = bessel_i_miller(0, &x, p);
        assert!(close_at(&t0, &m0, p - 12), "tiny vs Miller I_0(1)");
        assert!(close_at(&t0, &py(I0_1, p), p - 12), "tiny I_0(1)");
        assert!(close_at(&m0, &py(I0_1, p), p - 12), "Miller I_0(1)");
        let t1 = bessel_i_tiny(1, &x, p);
        let m1 = bessel_i_miller(1, &x, p);
        assert!(close_at(&t1, &m1, p - 12), "tiny vs Miller I_1(1)");
        assert!(close_at(&m1, &py(I1_1, p), p - 12), "Miller I_1(1)");
    }

    // mpmath besseli at large x (dps = 120).
    const I0_200: &str = "20396871734097246195416731267794596223326757361483433789432837835530189904402084547243.5871";
    const I1_200: &str = "20345815493320627034274279771390695038966116168112296415921960606931666871258425531892.7136";
    const I2_200: &str =
        "20193413579164039925073988470080689272937096199802310825273618229460873235689500291924.66";
    const I3_200: &str = "19941947221737346235772800001989081253507374244116250199416488242342449406544635526054.2204";
    const I5_200: &str = "19158141015236869454252767823188240580094099245217097266644843299054825416404909077008.499";
    const I0_1000: &str = "2.48568609607586417456277148414567563132911034374842167789492204933350605326104812721297733851453673308429390778241907277e+432";
    const I1_1000: &str = "2.4844429420058669729947094283340842192656945142670877097696307799870865072981949999169252725636201630993234673884229269e+432";
    const I2_1000: &str = "2.48071721019185244061678206528900746289057895471988750247538278777353188024645173721314348796940949275809526084764222692e+432";

    /// DLMF 10.40.1 asymptotic path, called directly so the reused
    /// ADR-0023 `a_k(n)` recurrence and the explicit `(−1)^j` sign
    /// (no trig) are pinned independently of the regime dispatch:
    /// `I_n(200)` (`p = 53`) and `I_n(1000)` (`p = 113`) vs mpmath.
    /// Large enough that the second and later `a_k` terms materially
    /// contribute, so a recalled coefficient or a wrong sign fails
    /// here. The public path at `x = 200`, `p = 53` must route
    /// through the dispatch into the asymptotic (pins the *derived*
    /// `bessel_j_threshold` reuse).
    #[test]
    fn asymptotic_regime_matches_mpmath() {
        let p = 53;
        let x = at(200, 1, p);
        for (n, want) in [
            (0u32, I0_200),
            (1, I1_200),
            (2, I2_200),
            (3, I3_200),
            (5, I5_200),
        ] {
            let r = bessel_i_asymptotic(n, &x, p);
            assert!(close_at(&r, &py(want, p), p - 8), "I_{n}(200)");
        }
        let p = 113;
        let x = at(1000, 1, p);
        for (n, want) in [(0u32, I0_1000), (1, I1_1000), (2, I2_1000)] {
            let r = bessel_i_asymptotic(n, &x, p);
            assert!(close_at(&r, &py(want, p), p - 12), "I_{n}(1000)");
        }

        let x = at(200, 1, 53);
        let (r0, _) = x.i0(RoundingMode::NearestEven);
        assert!(close_at(&r0, &py(I0_200, 53), 45), "i0(200) via dispatch");
        let (r1, _) = x.i1(RoundingMode::NearestEven);
        assert!(close_at(&r1, &py(I1_200, 53), 45), "i1(200) via dispatch");
    }

    // mpmath besseli (dps = 120) for the regime-boundary check.
    const I0_256: &str = "37704216813925512248187965940740515799724074745159949233554531401110625525761751932277313607778930620432331127.5288431961";
    const I1_256: &str = "37630503567727454079825096501224087285338948500695164037686115528501114205565327708769760904440709640999968454.0273413245";

    /// Cross-regime continuity: at `x = 256` the DLMF 10.40.1
    /// asymptotic and the DLMF 10.29.1 Miller recurrence (called
    /// directly) agree and both match mpmath. Pins the *derived*
    /// `bessel_j_threshold` crossover and the reused `a_k(n)`
    /// recurrence against the independent Miller path (the
    /// `bessel_j` `asymptotic_miller_continuity` analog).
    #[test]
    fn asymptotic_miller_continuity() {
        let p = 113;
        let x = at(256, 1, p);
        let a0 = bessel_i_asymptotic(0, &x, p);
        let mi0 = bessel_i_miller(0, &x, p);
        assert!(close_at(&a0, &mi0, p - 12), "asymp vs Miller I_0(256)");
        assert!(close_at(&a0, &py(I0_256, p), p - 12), "asymp I_0(256)");
        assert!(close_at(&mi0, &py(I0_256, p), p - 12), "Miller I_0(256)");
        let a1 = bessel_i_asymptotic(1, &x, p);
        let mi1 = bessel_i_miller(1, &x, p);
        assert!(close_at(&a1, &mi1, p - 12), "asymp vs Miller I_1(256)");
        assert!(close_at(&a1, &py(I1_256, p), p - 12), "asymp I_1(256)");
    }

    // mpmath besseli (dps = 340) for the p = 1024 second-term pin.
    const I0_52_BIG: &str = "3.28983914405012303570590822990605602611180157534839416125528705344053812810694973324070513417838122219340846596165467918112126063839091011363774129787473762498039721506143816371758308617539645157852076897235980325849930352382392897037908199549603978795284238351346766116171420991967458447353913510971538513036177229884938048029043";
    const I2_256_BIG: &str = "37410228504802641513189332374324702617807364209998268264510108623544210571030772809552549850712987576362018873.9817545920191533882152982129453659415258064928448971088079928843708892549110307533090426032502310242193371456222198148772771860159740743468648991188382406203461068619103845253344226863032223166060549000013756896507611957";

    /// Second-term-matters pin at `p = 1024` (the `derive, don't
    /// recall` reflex): the Miller path `I_0(2.5)` and the deeper
    /// Miller path `I_2(256)` (at `p = 1024` the
    /// precision-scaled `bessel_j_threshold` keeps `x = 256` in
    /// Miller, not the asymptotic), validated to `p − 2` against the
    /// 330-digit references. A coefficient, recurrence-sign, or
    /// sum-rule error invisible at low precision fails here.
    #[test]
    fn high_precision_pin() {
        let x = at(5, 2, 1024);
        let (r, _) = x.i0(RoundingMode::NearestEven);
        assert!(close_at(&r, &py(I0_52_BIG, 1024), 1022), "I_0(2.5) p=1024");

        let x = at(256, 1, 1024);
        let (r, _) = x.in_(2, RoundingMode::NearestEven);
        assert!(close_at(&r, &py(I2_256_BIG, 1024), 1022), "I_2(256) p=1024");
    }

    #[test]
    fn i0_zero_is_one() {
        for s in [Sign::Positive, Sign::Negative] {
            let z = BigFloat::try_new_zero(s, 53).unwrap();
            let (r, st) = z.i0(RoundingMode::NearestEven);
            let one = BigFloat::try_from_i64_exact(1, 53).unwrap();
            assert_eq!(r.partial_cmp(&one).0, Some(Ordering::Equal), "I0(±0) = 1");
            assert!(!st.invalid());
        }
    }

    #[test]
    fn in_zero_is_zero_for_nonzero_order() {
        for n in [1i32, 2, 3, -1, -4] {
            for s in [Sign::Positive, Sign::Negative] {
                let z = BigFloat::try_new_zero(s, 53).unwrap();
                let (r, _) = z.in_(n, RoundingMode::NearestEven);
                assert!(r.is_zero(), "I_{n}(±0) = 0");
            }
        }
    }

    #[test]
    fn i_positive_infinity_is_positive_infinity() {
        let inf = BigFloat::try_new_infinity(Sign::Positive, 53).unwrap();
        for n in [0i32, 1, 2, 5, -3] {
            let (r, st) = inf.in_(n, RoundingMode::NearestEven);
            assert!(r.is_infinite() && r.is_sign_positive(), "I_{n}(+∞) = +∞");
            assert!(!st.invalid());
        }
    }

    #[test]
    fn i_negative_infinity_is_signed_by_parity() {
        let ninf = BigFloat::try_new_infinity(Sign::Negative, 53).unwrap();
        // Even order: I_n(−∞) = +∞. Odd order: −∞.
        for n in [0i32, 2, 4, -2] {
            let (r, st) = ninf.in_(n, RoundingMode::NearestEven);
            assert!(r.is_infinite() && r.is_sign_positive(), "I_{n}(−∞) = +∞");
            assert!(!st.invalid());
        }
        for n in [1i32, 3, -1, -3] {
            let (r, st) = ninf.in_(n, RoundingMode::NearestEven);
            assert!(r.is_infinite() && r.is_sign_negative(), "I_{n}(−∞) = −∞");
            assert!(!st.invalid());
        }
    }

    #[test]
    fn i_quiet_nan_propagates() {
        let q = BigFloat::try_new_quiet_nan(Sign::Positive, 53, &[]).unwrap();
        let (r, _) = q.in_(2, RoundingMode::NearestEven);
        assert!(r.is_quiet_nan());
    }

    #[test]
    fn i_signaling_nan_raises_invalid() {
        let sn = BigFloat::try_new_signaling_nan(Sign::Positive, 53, &[]).unwrap();
        let (r, st) = sn.i0(RoundingMode::NearestEven);
        assert!(r.is_quiet_nan());
        assert!(st.invalid());
    }

    #[test]
    fn precision_zero_is_rejected() {
        let x = BigFloat::try_from_i64_exact(2, 53).unwrap();
        assert!(x.i0_round(0, RoundingMode::NearestEven).is_err());
        assert!(x.in_round(3, 0, RoundingMode::NearestEven).is_err());
    }
}
