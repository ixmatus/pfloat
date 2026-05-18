//! Modified Bessel functions of the second kind `K0`, `K1`, `Kn`
//! (DLMF Chapter 10): integer order, real argument. Like
//! [`super::bessel_y`], `K` is real-valued only for `x > 0`: it has a
//! logarithmic branch point at the origin and is complex for
//! `x < 0`. The domain convention follows the [`super::ci`] /
//! [`super::li`] / [`super::bessel_y`] precedent (real-only, complex
//! off the positive axis):
//!
//! - `Kₙ(+0) = +∞`, raising `DIV_BY_ZERO` (a pole: DLMF 10.30.2
//!   `Kν(z) ∼ ½Γ(ν)(½z)⁻ν` for `ν > 0` and DLMF 10.30.3
//!   `K₀(z) ∼ −ln z` both diverge to **+∞** as `x → 0⁺`). Note the
//!   sign is the opposite of `Yₙ(+0) = −∞`.
//! - `x < 0` (and `−0`, `−∞`) ⇒ `NaN` + `INVALID` (`K` is complex in
//!   the reals there).
//! - `Kₙ(+∞) = +0` for every order, `Status::OK`. This is a
//!   **genuine exponential-decay limit** (DLMF 10.40.2
//!   `√(π/2x)·e⁻ˣ → 0`), **not** the decaying-envelope *convention*
//!   used by `J`/`Y`/Airy. There the function oscillates with a
//!   shrinking but non-converging envelope and `+0` is a
//!   conservative choice that keeps the function total; here `K`
//!   actually converges to `0`, so `+0` is the true mathematical
//!   limit. ADR-0025 records the distinction.
//! - `Kₙ(NaN) = NaN`; `sNaN` raises `INVALID`.
//!
//! The order parity is **even with no sign**: `K₋ₙ(x) = Kₙ(x)`
//! (DLMF 10.27.3), in contrast to `Y₋ₙ = (−1)ⁿ Yₙ`. There is no
//! argument-parity reduction (the domain is `x > 0`; a negative
//! argument is `INVALID`, not folded to `|x|`). So the kernel
//! evaluates `K_m(x)` for `m = |n| ≥ 0` with no sign applied at all.
//!
//! `K` is the *dominant* solution of the modified-Bessel three-term
//! recurrence in order (DLMF 10.30.2 `Kν → ∞` as `ν → ∞` at fixed
//! `z`), so the kernel computes `K₀` and `K₁` directly and climbs to
//! `Kₙ` by **upward** recurrence (DLMF 10.29.1; slice 6q.6), the
//! [`super::bessel_y`] template. The recurrence is
//! `K_{k+1}(x) = (2k/x)·K_k(x) + K_{k−1}(x)`: DLMF 10.29.1 is the
//! unified `𝒵_{ν−1} − 𝒵_{ν+1} = (2ν/z)𝒵_ν` with the §10.25(ii)
//! convention `𝒵_ν = I_ν` **or `e^{νπi} K_ν`**, and the `e^{νπi}`
//! factor flips `K`'s sign relative to a naive reading (and
//! relative to `I`); see [`bessel_k_eval_normal`] for the
//! derivation. The base pair
//! `K₀`/`K₁` uses two regimes on the binary exponent of `x`: the
//! DLMF 10.31.1 logarithmic series below the cut (slice 6q.5), the
//! DLMF 10.40.2 asymptotic at/above it (slice 6q.7). ADR-0025
//! records the design and the DLMF provenance.

use super::euler_gamma_at;
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
    /// `K₀(self)` rounded under `mode` to `self.precision`.
    #[must_use]
    pub fn k0(&self, mode: RoundingMode) -> (Self, Status) {
        self.k0_round(self.precision, mode)
            .expect("self.precision >= 1 by invariant")
    }

    /// `K₀(self)` with explicit result precision.
    pub fn k0_round(
        &self,
        target_precision: u32,
        mode: RoundingMode,
    ) -> Result<(Self, Status), BuildError> {
        if target_precision == 0 {
            return Err(BuildError::PrecisionZero);
        }
        Ok(bessel_k_kernel(0, self, target_precision, mode))
    }

    /// `K₁(self)` rounded under `mode` to `self.precision`.
    #[must_use]
    pub fn k1(&self, mode: RoundingMode) -> (Self, Status) {
        self.k1_round(self.precision, mode)
            .expect("self.precision >= 1 by invariant")
    }

    /// `K₁(self)` with explicit result precision.
    pub fn k1_round(
        &self,
        target_precision: u32,
        mode: RoundingMode,
    ) -> Result<(Self, Status), BuildError> {
        if target_precision == 0 {
            return Err(BuildError::PrecisionZero);
        }
        Ok(bessel_k_kernel(1, self, target_precision, mode))
    }

    /// `Kₙ(self)` for integer order `n` (any `i32`, including
    /// negative: `K₋ₙ = Kₙ`), rounded under `mode` to
    /// `self.precision`.
    #[must_use]
    pub fn kn(&self, n: i32, mode: RoundingMode) -> (Self, Status) {
        self.kn_round(n, self.precision, mode)
            .expect("self.precision >= 1 by invariant")
    }

    /// `Kₙ(self)` with explicit result precision.
    pub fn kn_round(
        &self,
        n: i32,
        target_precision: u32,
        mode: RoundingMode,
    ) -> Result<(Self, Status), BuildError> {
        if target_precision == 0 {
            return Err(BuildError::PrecisionZero);
        }
        Ok(bessel_k_kernel(n, self, target_precision, mode))
    }
}

#[cfg(feature = "fixed")]
impl<const PREC: u32> FixedFloat<PREC>
where
    [(); limbs_for(PREC)]:,
{
    /// `K₀(self)` for `FixedFloat`. Delegates to [`BigFloat::k0`].
    #[must_use]
    pub fn k0(&self, mode: RoundingMode) -> (Self, Status) {
        let (big, status) = self.to_big().k0(mode);
        (
            Self::try_from_big_exact(big).expect("precision matches"),
            status,
        )
    }

    /// `K₁(self)` for `FixedFloat`. Delegates to [`BigFloat::k1`].
    #[must_use]
    pub fn k1(&self, mode: RoundingMode) -> (Self, Status) {
        let (big, status) = self.to_big().k1(mode);
        (
            Self::try_from_big_exact(big).expect("precision matches"),
            status,
        )
    }

    /// `Kₙ(self)` for `FixedFloat`. Delegates to [`BigFloat::kn`].
    #[must_use]
    pub fn kn(&self, n: i32, mode: RoundingMode) -> (Self, Status) {
        let (big, status) = self.to_big().kn(n, mode);
        (
            Self::try_from_big_exact(big).expect("precision matches"),
            status,
        )
    }
}

/// `Kₙ(x)` for integer order `n` and real `x`.
///
/// Special values are handled directly per the module-level domain
/// table; for a normal positive argument the order is reduced to
/// `K_m(x)`, `m = |n|`, with no sign applied (`K₋ₙ = Kₙ`; the domain
/// is `x > 0`, so there is no argument parity), then the regime
/// evaluator runs.
fn bessel_k_kernel(
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
                // −0: K is complex off the positive axis (the Ci/li
                // convention; −0 groups with x < 0).
                let nan = BigFloat::try_new_quiet_nan(Sign::Positive, target_precision, &[])
                    .expect("precision >= 1");
                auto_raise(Status::INVALID);
                return (nan, Status::INVALID);
            }
            // Kₙ(+0) = +∞ + DIV_BY_ZERO (a pole, DLMF 10.30.2/10.30.3;
            // +∞, the opposite sign of Yₙ(+0) = −∞).
            let pinf = BigFloat::try_new_infinity(Sign::Positive, target_precision)
                .expect("precision >= 1");
            auto_raise(Status::DIV_BY_ZERO);
            (pinf, Status::DIV_BY_ZERO)
        }
        Class::Infinity { sign } => {
            if matches!(sign, Sign::Negative) {
                // −∞: complex (groups with x < 0).
                let nan = BigFloat::try_new_quiet_nan(Sign::Positive, target_precision, &[])
                    .expect("precision >= 1");
                auto_raise(Status::INVALID);
                return (nan, Status::INVALID);
            }
            // Kₙ(+∞) = +0, the genuine exponential-decay limit
            // (DLMF 10.40.2 √(π/2x)·e⁻ˣ → 0), Status::OK. This is a
            // true limit, not the decaying-envelope convention.
            let zero =
                BigFloat::try_new_zero(Sign::Positive, target_precision).expect("precision >= 1");
            (zero, Status::OK)
        }
        Class::Normal { sign, .. } => {
            if matches!(sign, Sign::Negative) {
                // x < 0: Kₙ(−x) is complex in the reals.
                let nan = BigFloat::try_new_quiet_nan(Sign::Positive, target_precision, &[])
                    .expect("precision >= 1");
                auto_raise(Status::INVALID);
                return (nan, Status::INVALID);
            }

            // K₋ₙ(x) = Kₙ(x) (DLMF 10.27.3): even in order, no sign.
            // No argument parity (x > 0 only). Evaluate K_m(x).
            let value = bessel_k_eval_normal(m, x, target_precision);

            let (rounded, status) = value
                .round_to_precision(target_precision, mode)
                .expect("precision >= 1");
            auto_raise(status);
            (rounded, status)
        }
    }
}

/// `K_m(x)` for `m ≥ 0`, normal `x > 0`: the base pair plus upward
/// recurrence. Returns the unrounded working-precision value;
/// [`bessel_k_kernel`] does the single final round.
///
/// The base pair `K₀`/`K₁` goes through [`bessel_k_series`] (slice
/// 6q.5; slice 6q.7 adds the two-regime series/asymptotic base
/// dispatch). `Kₙ (n ≥ 2)` climbs from that pair by the **upward**
/// recurrence (DLMF 10.29.1)
///
/// ```text
/// K_{k+1}(x) = (2k/x)·K_k(x) + K_{k−1}(x)
/// ```
///
/// **derived, not recalled** (the 6n `(2k−1)` defect is the
/// precedent): DLMF 10.29.1 is the unified
/// `𝒵_{ν−1} − 𝒵_{ν+1} = (2ν/z)·𝒵_ν` where, by the §10.25(ii)
/// standard-solution convention, `𝒵_ν = I_ν(z)` **or
/// `e^{νπi} K_ν(z)`**. Substituting the latter and dividing by
/// `e^{νπi}` turns `e^{−πi}K_{ν−1} − e^{πi}K_{ν+1} = (2ν/z)K_ν`
/// into `−K_{ν−1} + K_{ν+1} = (2ν/z)K_ν`, i.e. the **`+ K_{k−1}`**
/// form above — the opposite sign to a naive reading of the
/// unified relation, and the opposite of `I`'s
/// `f_{k−1} = (2k/x)f_k + f_{k+1}` rearrangement. Numerically
/// pinned by `K₂ − K₀ = (2/x)K₁` at `x = 1` (`recurrence_spot_check`).
/// `K` is the *dominant* solution in order (DLMF 10.30.2
/// `Kν ∼ ½Γ(ν)(2/z)ν → ∞` as `ν → ∞`), so this forward climb is
/// numerically stable (it grows the dominant solution; the
/// all-positive `(2k/x)·K_k + K_{k−1}` has no cancellation). The
/// climb is a rolling pair, `O(m)` time and `O(1)` space, no `Vec`,
/// the [`super::bessel_y`] template (with the sign change and the
/// `+` recurrence rather than `Y`'s `Y_{k+1}=(2k/x)Y_k−Y_{k−1}`).
fn bessel_k_eval_normal(m: u32, x: &BigFloat, target_precision: u32) -> BigFloat {
    if m <= 1 {
        return bessel_k_series(m, x, target_precision);
    }

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

    let xw = x
        .round_to_precision(working, RoundingMode::NearestEven)
        .expect("precision >= 1")
        .0;
    let (inv_x, _) = ci(1, working).div(&xw, RoundingMode::NearestEven);

    // Seeds K₀, K₁ at the recurrence working precision.
    let mut k_prev = bessel_k_series(0, x, working)
        .round_to_precision(working, RoundingMode::NearestEven)
        .expect("precision >= 1")
        .0;
    let mut k_cur = bessel_k_series(1, x, working)
        .round_to_precision(working, RoundingMode::NearestEven)
        .expect("precision >= 1")
        .0;

    // K_{k+1} = (2k/x)·K_k + K_{k−1}, k = 1 … m−1.
    for k in 1..i64::from(m) {
        let two_k = ci(2 * k, working);
        let (c1, _) = two_k.mul(&inv_x, RoundingMode::NearestEven);
        let (c2, _) = c1.mul(&k_cur, RoundingMode::NearestEven);
        let (k_next, _) = c2.add(&k_prev, RoundingMode::NearestEven);
        k_prev = k_cur;
        k_cur = k_next;
    }
    k_cur
}

/// Binary exponent of `v`, or `i64::MIN`/`i64::MAX` for zero /
/// non-finite (the [`super::bessel_y`] `magnitude` idiom; `K`
/// carries its own, as `bessel_j`/`bessel_y`/`bessel_i` do).
fn magnitude(v: &BigFloat) -> i64 {
    match &v.class {
        Class::Normal { exponent, .. } => *exponent,
        Class::Zero { .. } => i64::MIN,
        _ => i64::MAX,
    }
}

/// `true` once `term` has fallen `working + 8` bits below the running
/// `sum`, so further terms cannot perturb the rounded result (the
/// [`super::bessel_y`] `negligible` idiom).
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

/// `K_n(x)` for integer `n ≥ 0`, `x > 0`, via the DLMF 10.31.1
/// logarithmic series
///
/// ```text
/// K_n(x) = ½(x/2)^{−n} Σ_{k=0}^{n−1} (n−k−1)!/k! (−x²/4)^k
///        + (−1)^{n+1} ln(x/2) · I_n(x)
///        + (−1)^n ½(x/2)^n Σ_{k≥0} (ψ(k+1)+ψ(n+k+1)) (x²/4)^k/(k!(n+k)!)
/// ```
///
/// the modified-Bessel analog of the `Y` log series (DLMF 10.8.1),
/// with three sign differences derived from DLMF 10.31.1 (not
/// recalled; the 6n `(2k−1)` defect is the precedent): the finite
/// head carries `+½` and an **alternating** `(−x²/4)^k` (where `Y`'s
/// head is `−1/π` and a plain `(x²/4)^k`); the log term is
/// `(−1)^{n+1}` with **no `2/π`** (where `Y`'s is `+2/π`); the tail
/// carries `(−1)^n ½` and a **positive** `(x²/4)^k` (where `Y`'s
/// tail is `−1/π` and an alternating `(−x²/4)^k`). The finite head
/// sum is empty for `n = 0` (the `n ≥ 1` pole source). The digamma
/// terms reduce to running harmonic sums exactly as for `Y`:
/// `ψ(k+1) = −γ + H_k`, `ψ(n+k+1) = −γ + H_{n+k}`, so
/// `ψ(k+1)+ψ(n+k+1) = −2γ + H_k + H_{n+k}` (`H_0 = 0`) — **no
/// digamma kernel**, only harmonic partial sums and the in-tree
/// Euler–Mascheroni `γ` ([`super::euler_gamma_at`], slice 6m0). The
/// `I_n(x)` piece is the 6q `bessel_i` kernel
/// ([`BigFloat::in_round`]).
///
/// Cross-check (the `derive, don't recall` reflex, the DLMF
/// 10.31.1 ↔ 10.31.2 worked identity): for `n = 0` the finite head
/// is empty, `(−1)^{0+1} ln(x/2) I_0 = −ln(x/2) I_0`, and the tail
/// is `½ Σ 2ψ(k+1)(x²/4)^k/(k!)² = Σ ψ(k+1)(x²/4)^k/(k!)²`.
/// Substituting `ψ(k+1) = −γ + H_k`, using
/// `Σ (x²/4)^k/(k!)² = I_0(x)` and `H_0 = 0`, the total collapses to
/// `K_0(x) = −(ln(x/2)+γ) I_0(x) + Σ_{k≥1} H_k (x²/4)^k/(k!)²`,
/// exactly DLMF 10.31.2. The two DLMF forms agreeing is the worked
/// check (pinned by `k0_dlmf_10_31_2_crosscheck`).
///
/// The head's `(−x²/4)^k` alternation and the log/head/tail
/// composition cancel, so working precision is boosted `≈ x·log₂e`
/// (the [`super::bessel_y`] `ci.rs`/`si.rs` capped guard). Returns
/// the unrounded working-precision value.
fn bessel_k_series(n: u32, x: &BigFloat, target_precision: u32) -> BigFloat {
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
    let gamma = euler_gamma_at(working);
    let one = ci(1, working);
    let half_ci = one.div(&two, RoundingMode::NearestEven).0; // ½

    // (x/2)^n by repeated multiply (n small; in-module recurrence
    // from a base term, not pow.rs::pow_int).
    let mut half_pow_n = one.clone();
    for _ in 0..n {
        half_pow_n = half_pow_n.mul(&half, RoundingMode::NearestEven).0;
    }

    // Log term: (−1)^{n+1} ln(x/2) · I_n(x). (−1)^{n+1} = +1 for odd
    // n, −1 for even n. No 2/π factor (contrast Y's 10.8.1).
    let (ln_half_x, _) = half.ln(RoundingMode::NearestEven);
    let (i_n, _) = xw
        .in_round(n as i32, working, RoundingMode::NearestEven)
        .expect("precision >= 1");
    let (lx_in, _) = ln_half_x.mul(&i_n, RoundingMode::NearestEven);
    let log_term = if n % 2 == 1 { lx_in } else { lx_in.negated() };

    // Finite head +½(x/2)^{−n} Σ_{k=0}^{n−1} (n−k−1)!/k! (−x²/4)^k.
    // h_0 = (n−1)!; h_k = h_{k−1}·(−x²/4)/(k·(n−k)). Empty for n = 0.
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
            term = term.mul(&neg_qxs, RoundingMode::NearestEven).0;
            term = term
                .div(&ci(k * (nn - k), working), RoundingMode::NearestEven)
                .0;
            hd = hd.add(&term, RoundingMode::NearestEven).0;
        }
        let (inv_hp, _) = one.div(&half_pow_n, RoundingMode::NearestEven); // (x/2)^{−n}
        let (a, _) = hd.mul(&inv_hp, RoundingMode::NearestEven);
        a.mul(&half_ci, RoundingMode::NearestEven).0 // +½ · …
    };

    // Tail (−1)^n ½(x/2)^n Σ_{k≥0} (−2γ + H_k + H_{n+k}) c_k,
    // c_0 = 1/n!, c_k = c_{k−1}·(x²/4)/(k·(n+k))  (POSITIVE x²/4,
    // contrast Y's −x²/4).
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
        c = c.mul(&qxs, RoundingMode::NearestEven).0;
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
    let (tp1, _) = tp0.mul(&half_ci, RoundingMode::NearestEven); // ½ · …
    let tail_term = if n % 2 == 1 { tp1.negated() } else { tp1 }; // (−1)^n

    let (s1, _) = head_term.add(&log_term, RoundingMode::NearestEven);
    let (k_val, _) = s1.add(&tail_term, RoundingMode::NearestEven);
    k_val
}

#[cfg(test)]
mod tests {
    use super::*;
    use core::cmp::Ordering;

    /// `|v − expected| ≤ 2^(−bits)·|expected|` (the `bessel_j`/
    /// `bessel_y`/`bessel_i` test helper). Reference decimals:
    /// `mpmath` `besselk(n, x)` at `mp.dps = 330`
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

    // mpmath besselk(n, x), dps = 330 (truncated to fit p = 113).
    const K0_HALF: &str = "0.92441907122766586178192416753021698953876831195352968481501974063291996009501604867818076098235912000427368840729058351991770641488764243451672765636339533291397998639934907236026386948677230468237592737790977883002994546654825607094160547875637743429736784627856231110527836389079862970859130931388929445638628327371615";
    const K1_HALF: &str = "1.6564411200033008936964454031740915115341007594640774460554278145261965895145190494662987924191674125535071795407350533285497439222917759436984298078978886519109188566482265764992625931087886254516334162416038604462247270788447321549733571254845780489233396611700563796431383011337659885036898595269150267772017950278808";
    const K2_HALF: &str = "7.5501835512408694365677057802265830356751713498098394690367309987377063181530922465433759306590287702183024065702307968341166821040547462093104468879549499405576554129922553783573142419219268064889095923443252206149288537819271846908350339806946896299907264909587878296778315684258625837233507474215494015651934633852393";
    const K0_52: &str = "0.062347553200366186029169529476013925996005578743445303838599172078372612679992899458805237455015772383234650436222307193733385924874664762854707475585477522914665853136001287110210662495170045365059918105322606104140571847060926644897405882307279238484634954137887406874971337052133924153978853479760979621356943415886755";
    const K1_52: &str = "0.073890816347747063648993540591217582101975744210962005306723010118640815423264789079998669903769446192709979004231962186722904993570366173945941847023125465297962409843223129226161716463143121265289812018988432110004162870359315151908540027633085782112292566248103032238256452998185913523743132471215697189809655576115526";
    const K2_52: &str = "0.12146020627856383694836436194898799167758617411221490808397758017328526501860473072280417337803132933740263363960787694311170991973095770201146095320397789515303578101057979049114003566568454237729176772051335179214390214334837876642423790441374786417446900713636983266557649945068265497297335945673353737320466787677918";
    const K0_7: &str = "0.00042479574186923180685159865280657229397931752463243151357157944598225095002003175131102192637869339162928957505052080071428655942756877188480724097957734454681984001374809296247474819146725015388359624973830793466403113827256838132763217970871380153407717390543264830296218677328044508642002661301005696955299043061126624";
    const K1_7: &str = "0.00045418248688489697123995940271024363126799741919564224934062310511312502497023217346402141962579812124144905779053397536579617082739569706271086103049567037209874117639561728778405315524542892357946837306989706465927472359774959741339129189364362862918895555361032059573782353608243367967495720792939623359103473402586464";
    const K2_7: &str = "0.00055456216669348808434872991072378476005588821583118644195461461887171524286866951515788518912892142626970359156210193653308546537825325675986748698829039322456233749271826933041304909296594413204915864204399281028096677358621112344574397739261198114241687549217845418745870778358971185204144295813274160772185749747579899";
    const K0_M3: &str = "7.0236888005623813436120800630140551147904217215113641495598227153519412581158941161072858437366484115301350390760361520149619794407508448324440924964988768279653955526630381497792016604102238019358857312218785573692672185734096932143775900974500473783684717352347088010659927078707528628821387875091928113944349324310431";
    const K1_M3: &str = "999.99623815608557427795340401620491430757739517990119191272555082855272043700974597141908872781027251247264327790341305248295819815723540543922349946607690893732796182434680861134870411968711710070520015398619451671937650195301472462110240158215861035367554706577716352330844503098990602073076965497981667540059270153659";
    const K2_M3: &str = "1999999.5000009717109372504201124728426702695807815238951896006614798207928152776078369542847414642816733568166908459021411179313582939115617232794430246503167514838890442462802608471874410346444252123361937036109119961222711246028589354191807544146707547294626032895617554179560546876827943244214487471425436125798380056";

    /// DLMF 10.31.1 log-series path vs mpmath at `p = 113`: orders
    /// `0,1,2` over small (`0.5`), moderate (`2.5`, `7`) `x`, plus
    /// the pole approach (`1e-3`, where the `(x/2)^{−n}` head
    /// dominates for `n ≥ 1`). Large enough orders/arguments that the
    /// second tail term and the alternating head both materially
    /// contribute (the `derive, don't recall` reflex: a low-`x`-only
    /// check would hide the head/tail sign forks).
    #[test]
    fn series_matches_mpmath() {
        let p = 113;
        let x = at(1, 2, p); // 0.5
        for (n, want) in [(0i32, K0_HALF), (1, K1_HALF), (2, K2_HALF)] {
            let (r, _) = x.kn(n, RoundingMode::NearestEven);
            assert!(close_at(&r, &py(want, p), p - 8), "K_{n}(0.5)");
        }
        let x = at(5, 2, p); // 2.5
        for (n, want) in [(0i32, K0_52), (1, K1_52), (2, K2_52)] {
            let (r, _) = x.kn(n, RoundingMode::NearestEven);
            assert!(close_at(&r, &py(want, p), p - 10), "K_{n}(2.5)");
        }
        let x = at(7, 1, p);
        for (n, want) in [(0i32, K0_7), (1, K1_7), (2, K2_7)] {
            let (r, _) = x.kn(n, RoundingMode::NearestEven);
            assert!(close_at(&r, &py(want, p), p - 12), "K_{n}(7)");
        }
        let x = at(1, 1000, p); // 1e-3, pole approach
        for (n, want) in [(0i32, K0_M3), (1, K1_M3), (2, K2_M3)] {
            let (r, _) = x.kn(n, RoundingMode::NearestEven);
            assert!(close_at(&r, &py(want, p), p - 12), "K_{n}(1e-3)");
        }
    }

    /// The DLMF 10.31.1 ↔ 10.31.2 worked cross-check (the
    /// `derive, don't recall` reflex, the doc's reduction made
    /// runnable): `bessel_k_series(0, x)` (the general 10.31.1 path,
    /// `n = 0`) must equal the independent DLMF 10.31.2 closed form
    /// `K_0(x) = −(ln(x/2)+γ) I_0(x) + Σ_{k≥1} H_k (x²/4)^k/(k!)²`,
    /// evaluated here by a separate series. Agreement pins the
    /// harmonic / `−2γ` reduction and the head/log/tail signs at
    /// `n = 0` against a second derivation, not against mpmath.
    #[test]
    fn k0_dlmf_10_31_2_crosscheck() {
        let p = 240;
        // The independent 10.31.2 form has a `−(ln(x/2)+γ)I_0 + Σ…`
        // catastrophic cancellation (≈ hundreds → ~1e−4 at x = 7);
        // give the hand-rolled side ample headroom so the comparison
        // reflects formula agreement, not the test's own rounding
        // (`bessel_k_series` boosts internally; this must too).
        let pp = p + 192;
        for &(num, den) in &[(1i64, 2i64), (5, 2), (7, 1)] {
            let x = at(num, den, p);
            let via_10_31_1 = bessel_k_series(0, &x, p);

            // Independent DLMF 10.31.2 evaluation at boosted pp:
            // K_0 = −(ln(x/2)+γ) I_0 + Σ_{k≥1} H_k (x²/4)^k/(k!)².
            let xw = x
                .round_to_precision(pp, RoundingMode::NearestEven)
                .unwrap()
                .0;
            let two = ci(2, pp);
            let (half, _) = xw.div(&two, RoundingMode::NearestEven);
            let (x_sq, _) = xw.mul(&xw, RoundingMode::NearestEven);
            let (qxs, _) = x_sq.div(&ci(4, pp), RoundingMode::NearestEven);
            let (ln_half, _) = half.ln(RoundingMode::NearestEven);
            let gamma = euler_gamma_at(pp);
            let (lg, _) = ln_half.add(&gamma, RoundingMode::NearestEven);
            let (i0, _) = xw.i0(RoundingMode::NearestEven);
            let (lead, _) = lg.mul(&i0, RoundingMode::NearestEven);
            let mut acc = lead.negated();
            let one = ci(1, pp);
            let mut t = one.clone();
            let mut h = BigFloat::try_new_zero(Sign::Positive, pp).unwrap();
            for k in 1..=(1i64 << 20) {
                t = t.mul(&qxs, RoundingMode::NearestEven).0;
                t = t.div(&ci(k * k, pp), RoundingMode::NearestEven).0;
                let (rk, _) = one.div(&ci(k, pp), RoundingMode::NearestEven);
                h = h.add(&rk, RoundingMode::NearestEven).0;
                let term = t.mul(&h, RoundingMode::NearestEven).0;
                acc = acc.add(&term, RoundingMode::NearestEven).0;
                if negligible(&term, &acc, pp) {
                    break;
                }
            }
            assert!(
                close_at(&via_10_31_1, &acc, p - 12),
                "10.31.1 vs 10.31.2 at x={num}/{den}"
            );
        }
    }

    /// Upward-recurrence cross-tie `K_{n+1}(x) − K_{n−1}(x) =
    /// (2n/x)·K_n(x)` (the verified DLMF 10.29.1 form for `K`, the
    /// `+ K_{k−1}` rearrangement), binding three orders that all
    /// climb from the same `(K₀, K₁)` base pair. This is the numeric
    /// pin for the recurrence **sign** (a wrong sign — the plan's
    /// original `− K_{k−1}` — fails here; the derive-don't-recall
    /// catch).
    #[test]
    fn recurrence_spot_check() {
        let p = 160;
        let x = at(5, 2, p); // 2.5
        let (k2, _) = x.kn(2, RoundingMode::NearestEven);
        let (k3, _) = x.kn(3, RoundingMode::NearestEven);
        let (k4, _) = x.kn(4, RoundingMode::NearestEven);
        let (lhs, _) = k4.sub(&k2, RoundingMode::NearestEven);
        let six = ci(6, p);
        let (r1, _) = six.mul(&k3, RoundingMode::NearestEven);
        let (rhs, _) = r1.div(&x, RoundingMode::NearestEven);
        assert!(close_at(&lhs, &rhs, p - 10), "K_4−K_2 = (6/x)K_3");
    }

    /// Recurrence agrees with the direct DLMF 10.31.1 series for
    /// `n ≥ 2` (independent paths: `bessel_k_eval_normal` climbs from
    /// `(K₀,K₁)`, `bessel_k_series` computes `K_n` directly), at
    /// moderate `x`.
    #[test]
    fn recurrence_matches_series() {
        let p = 160;
        let x = at(7, 2, p); // 3.5
        for n in [2u32, 3, 5] {
            let recur = bessel_k_eval_normal(n, &x, p);
            let series = bessel_k_series(n, &x, p);
            assert!(close_at(&recur, &series, p - 12), "K_{n}(3.5) recur=series");
        }
    }

    #[test]
    fn k_positive_zero_is_pole() {
        let z = BigFloat::try_new_zero(Sign::Positive, 53).unwrap();
        let (r, s) = z.k0(RoundingMode::NearestEven);
        assert!(r.is_infinite() && r.is_sign_positive(), "K0(+0) = +∞");
        assert!(s.div_by_zero());
        let (r1, s1) = z.k1(RoundingMode::NearestEven);
        assert!(r1.is_infinite() && r1.is_sign_positive(), "K1(+0) = +∞");
        assert!(s1.div_by_zero());
        let (rn, sn) = z.kn(3, RoundingMode::NearestEven);
        assert!(rn.is_infinite() && rn.is_sign_positive(), "K3(+0) = +∞");
        assert!(sn.div_by_zero());
    }

    #[test]
    fn k_negative_zero_is_invalid() {
        let z = BigFloat::try_new_zero(Sign::Negative, 53).unwrap();
        let (r, s) = z.k0(RoundingMode::NearestEven);
        assert!(r.is_quiet_nan(), "K0(−0) = NaN");
        assert!(s.invalid());
    }

    #[test]
    fn k_negative_argument_is_invalid() {
        let x = BigFloat::try_from_i64_exact(-3, 53).unwrap();
        for n in [0i32, 1, 2, -2] {
            let (r, s) = x.kn(n, RoundingMode::NearestEven);
            assert!(r.is_quiet_nan(), "K_{n}(−3) = NaN (complex)");
            assert!(s.invalid());
        }
    }

    #[test]
    fn k_positive_infinity_is_zero() {
        let inf = BigFloat::try_new_infinity(Sign::Positive, 53).unwrap();
        for n in [0i32, 1, 5, -3] {
            let (r, st) = inf.kn(n, RoundingMode::NearestEven);
            assert!(r.is_zero() && r.is_sign_positive(), "K_{n}(+∞) = +0");
            assert!(!st.invalid());
        }
    }

    #[test]
    fn k_negative_infinity_is_invalid() {
        let ninf = BigFloat::try_new_infinity(Sign::Negative, 53).unwrap();
        let (r, s) = ninf.k1(RoundingMode::NearestEven);
        assert!(r.is_quiet_nan(), "K1(−∞) = NaN (complex)");
        assert!(s.invalid());
    }

    #[test]
    fn k_quiet_nan_propagates() {
        let q = BigFloat::try_new_quiet_nan(Sign::Positive, 53, &[]).unwrap();
        let (r, _) = q.kn(2, RoundingMode::NearestEven);
        assert!(r.is_quiet_nan());
    }

    #[test]
    fn k_signaling_nan_raises_invalid() {
        let sn = BigFloat::try_new_signaling_nan(Sign::Positive, 53, &[]).unwrap();
        let (r, st) = sn.k0(RoundingMode::NearestEven);
        assert!(r.is_quiet_nan());
        assert!(st.invalid());
    }

    #[test]
    fn precision_zero_is_rejected() {
        let x = BigFloat::try_from_i64_exact(2, 53).unwrap();
        assert!(x.k0_round(0, RoundingMode::NearestEven).is_err());
        assert!(x.kn_round(3, 0, RoundingMode::NearestEven).is_err());
    }
}
