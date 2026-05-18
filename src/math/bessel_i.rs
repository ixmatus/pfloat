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
            // odd and x < 0.
            let negate = (m % 2 == 1) && x.is_sign_negative();
            let ax = x.abs();

            let value = bessel_i_eval_normal(m, &ax, target_precision);
            let value = if negate { value.negated() } else { value };

            let (rounded, status) = value
                .round_to_precision(target_precision, mode)
                .expect("precision >= 1");
            auto_raise(status);
            (rounded, status)
        }
    }
}

/// `I_m(ax)` for `m ≥ 0`, normal `ax > 0`: the regime evaluator.
/// Returns the unrounded working-precision value; [`bessel_i_kernel`]
/// does the single final round.
///
/// Slice 6q.2 wires the DLMF 10.25.2 convergent Maclaurin series
/// ([`bessel_i_tiny`]), which is entire and converges for every `x`,
/// so it carries the whole real line on its own. Slices 6q.3 (Miller
/// backward recurrence normalised by the DLMF 10.35.5 sum rule, for
/// moderate `|x|`) and 6q.4 (DLMF 10.40.1 asymptotic plus the
/// binary-exponent regime dispatch, for large `|x|`) layer the
/// faster regimes on top; until then the series is the only path and
/// there is no dead code.
fn bessel_i_eval_normal(m: u32, ax: &BigFloat, target_precision: u32) -> BigFloat {
    bessel_i_tiny(m, ax, target_precision)
}

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
