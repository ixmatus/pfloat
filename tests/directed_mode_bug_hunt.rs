//! Slice 1 (pf-3rtr.2): adversarial directed-mode bug-hunt.
//!
//! pfloat-ball's soundness rests on directed-mode rounding
//! (`TowardNegative` / `TowardPositive`, and the `TowardZero` /
//! `NearestAway` siblings). This lane attacks the directed paths that
//! the existing coverage leaves thin, with `rug`/MPFR as the
//! independent oracle (the `NearestAway` arm synthesized via
//! `mpfr_oracle_for_mode`, since MPFR has no roundTiesToAway).
//!
//! Targets, in the priority the phase plan sets:
//!
//! 1. `log2` / `log10`. These are the only transcendentals on the
//!    surface still rounding by a fixed `target + 64` guard then a
//!    single directed round, with no Ziv interval-test certification
//!    (`src/math/log2.rs`, `log10.rs`). The differential harness's own
//!    `NEAREST_EVEN_ROUNDING_MODES` doc records that this pattern can
//!    diverge up to one ULP from the correctly-rounded result under a
//!    non-`NearestEven` mode at a tie. Both a general sweep and an
//!    adversarial near-power-of-two sweep (outputs land next to an
//!    integer grid point, the directed hard-to-round zone) probe it.
//!
//! 2. `tanh` general range. Expected clean: `tanh` is already
//!    `ziv_round`-wrapped. A green result here confirms the lane
//!    discriminates an uncertified kernel from a certified one rather
//!    than flagging everything.
//!
//! 3. The negate-without-mirror / odd-function class. The 1.2.0 fix
//!    routed the irrational-constant special cases through
//!    `signed_constant_at_round`; the odd-function kernels negate
//!    inside their `ziv_round` eval closure so the interval test
//!    re-certifies the signed value. This lane pins both: the
//!    special-case constants under all five modes, and a negative
//!    -domain sweep of the odd functions.
//!
//! Any mismatch panics with the reproducing `(input, mode, precision,
//! got, want)`. A confirmed finding is re-derived against a second
//! oracle before it drives a kernel edit (verify-the-verdict); this
//! lane is the first oracle.

#![cfg(all(unix, feature = "differential-mpfr"))]

mod differential;

use differential::{
    bigfloat_to_rug, mpfr_oracle_for_mode, mpfr_round_of, next_u64, round_ties_to_away,
    BIT_EXACT_ROUNDING_MODES,
};
use pfloat::{BigFloat, RoundingMode, Sign};
use rug::float::{Constant, Round};
use rug::Float;

/// Precisions probed. p=24 is the f32 grid (coarse, ties common),
/// p=53 the f64 grid, p=113 a wider transcendental target.
const PRECISIONS: &[u32] = &[24, 53, 113];

/// Per-(precision, mode) sample count. Kept modest so the lane runs in
/// the per-push budget; `PFLOAT_DEEP=1` widens it for a local soak.
fn sample() -> u32 {
    if std::env::var("PFLOAT_DEEP").is_ok() {
        200_000
    } else {
        4_000
    }
}

/// Exact in-place `x * 2^e` (power-of-two scaling only shifts the
/// binary exponent, so no rounding occurs at any precision).
fn scale2(x: &BigFloat, e: i64) -> BigFloat {
    let two = BigFloat::try_from_i64_exact(2, x.precision()).expect("precision >= 1");
    let mut v = x.clone();
    for _ in 0..e.unsigned_abs() {
        v = if e < 0 {
            v.div(&two, RoundingMode::NearestEven).0
        } else {
            v.mul(&two, RoundingMode::NearestEven).0
        };
    }
    v
}

/// A generic positive off-grid `BigFloat` at precision `p`: the ratio
/// of two random integers (so the value carries a full `p`-bit
/// significand after the division) scaled by a random power of two.
fn rand_pos(state: &mut u64, p: u32) -> BigFloat {
    let bits = p.min(53).saturating_sub(1).max(1);
    let modulus = 1u64 << bits;
    let a = (next_u64(state) % modulus).max(1) as i64;
    let b = (next_u64(state) % modulus).max(1) as i64;
    let af = BigFloat::try_from_i64_exact(a, p).expect("a < 2^(p-1)");
    let bf = BigFloat::try_from_i64_exact(b, p).expect("b < 2^(p-1)");
    let (x, _) = af.div(&bf, RoundingMode::NearestEven);
    let e = (next_u64(state) % 81) as i64 - 40;
    scale2(&x, e)
}

/// `2^k * (1 + small)`: a value just above the power of two `2^k`, so
/// `log2(x)` sits just above the integer grid point `k`. This is the
/// directed-mode hard-to-round zone for `log2` (the rounding boundary
/// is the grid point itself, not the midpoint). Not a power of two, so
/// the exact-dispatch does not fire.
fn near_pow2(state: &mut u64, p: u32) -> BigFloat {
    let k = (next_u64(state) % 60) as i64 - 30;
    let base = scale2(
        &BigFloat::try_from_i64_exact(1, p).expect("precision >= 1"),
        k,
    );
    // One ULP of `base` is `2^(k - p + 1)`.
    let ulp = scale2(
        &BigFloat::try_from_i64_exact(1, p).expect("precision >= 1"),
        k - i64::from(p) + 1,
    );
    let m = (next_u64(state) % ((1u64 << (p.min(20))) - 1) + 1) as i64;
    let mf = BigFloat::try_from_i64_exact(m, p).expect("m < 2^20 <= 2^(p) for p>=20; p>=24 here");
    let (delta, _) = ulp.mul(&mf, RoundingMode::NearestEven);
    base.add(&delta, RoundingMode::NearestEven).0
}

/// A positive `BigFloat` at precision `p` with `|x|` bounded to roughly
/// `[2^-max_abs_exp, 2^(max_abs_exp+1))`. Used by the asymptote-function
/// controls (`tanh`, `erf`, …) to stay in the Ziv-converging regime,
/// well short of the saturation tail that
/// `saturation_directed_sweep` and `saturation_limit_reproducers` cover. A mantissa in `[1, 2)` times a bounded power of two.
fn rand_bounded(state: &mut u64, p: u32, max_abs_exp: i64) -> BigFloat {
    let bits = p.min(53).saturating_sub(1).max(1);
    let a = (next_u64(state) % (1u64 << bits)) as i64;
    let one = BigFloat::try_from_i64_exact(1, p).expect("precision >= 1");
    let af = BigFloat::try_from_i64_exact(a, p).expect("a < 2^(p-1)");
    let denom = scale2(&one, i64::from(bits)); // 2^bits
    let (frac, _) = af.div(&denom, RoundingMode::NearestEven); // [0, 1)
    let mantissa = one.add(&frac, RoundingMode::NearestEven).0; // [1, 2)
    let span = 2 * max_abs_exp as u64 + 1;
    let e = (next_u64(state) % span) as i64 - max_abs_exp;
    scale2(&mantissa, e)
}

/// Round a high-precision constant `hp` to precision `p` under a pfloat
/// rounding mode, synthesizing `NearestAway` (MPFR has no
/// roundTiesToAway). `hp` must carry far more than `p` bits.
fn round_const(hp: &Float, mode: RoundingMode, p: u32) -> Float {
    match mpfr_round_of(mode) {
        Some(r) => Float::with_val_round(p, hp, r).0,
        None => round_ties_to_away(hp, p),
    }
}

/// Assert pfloat's unary kernel agrees with the MPFR oracle bit-for-bit
/// across all five modes on `x`. `op` names the rug reference; `kernel`
/// the pfloat method. Panics with the reproducer on the first mismatch.
fn check_unary(
    label: &str,
    x: &BigFloat,
    kernel: impl Fn(&BigFloat, RoundingMode) -> BigFloat,
    op: impl Fn(&Float, u32, Round) -> Float,
) {
    let p = x.precision();
    let rug_x = bigfloat_to_rug(x);
    for &mode in BIT_EXACT_ROUNDING_MODES {
        let got = bigfloat_to_rug(&kernel(x, mode));
        let want = mpfr_oracle_for_mode(|prec, round| op(&rug_x, prec, round), mode, p);
        // Bit-for-bit on finite values; NaN handled by callers (none
        // of the probed inputs map to NaN).
        assert!(
            got == want,
            "{label} mismatch at p={p}, mode={mode:?}: input={rug_x:?} \
             pfloat={got:?} mpfr={want:?}"
        );
    }
}

#[test]
fn log2_directed_modes_match_mpfr() {
    let mut state = u64::from_le_bytes(*b"log2hunt");
    let n = sample();
    for &p in PRECISIONS {
        for _ in 0..n {
            let x = rand_pos(&mut state, p);
            check_unary(
                "log2",
                &x,
                |bf, m| bf.log2(m).0,
                |rx, prec, round| Float::with_val_round(prec, rx.log2_ref(), round).0,
            );
        }
    }
}

#[test]
fn log10_directed_modes_match_mpfr() {
    let mut state = u64::from_le_bytes(*b"log10hnt");
    let n = sample();
    for &p in PRECISIONS {
        for _ in 0..n {
            let x = rand_pos(&mut state, p);
            check_unary(
                "log10",
                &x,
                |bf, m| bf.log10(m).0,
                |rx, prec, round| Float::with_val_round(prec, rx.log10_ref(), round).0,
            );
        }
    }
}

#[test]
fn log2_near_power_of_two_boundary_directed() {
    let mut state = u64::from_le_bytes(*b"log2near");
    let n = sample();
    for &p in PRECISIONS {
        for _ in 0..n {
            let x = near_pow2(&mut state, p);
            check_unary(
                "log2-near-pow2",
                &x,
                |bf, m| bf.log2(m).0,
                |rx, prec, round| Float::with_val_round(prec, rx.log2_ref(), round).0,
            );
            check_unary(
                "log10-near-pow2",
                &x,
                |bf, m| bf.log10(m).0,
                |rx, prec, round| Float::with_val_round(prec, rx.log10_ref(), round).0,
            );
        }
    }
}

#[test]
fn tanh_directed_modes_match_mpfr_control() {
    // Expected clean: tanh is ziv_round-wrapped. Confirms the lane
    // discriminates a certified kernel from the uncertified log family.
    let mut state = u64::from_le_bytes(*b"tanhctrl");
    let n = sample();
    for &p in PRECISIONS {
        for _ in 0..n {
            // tanh over the whole real line, bounded to |x| < 16 (the
            // Ziv-converging regime; the saturation tail is in the
            // dedicated saturation tests below).
            let mut x = rand_bounded(&mut state, p, 3);
            if next_u64(&mut state) & 1 == 1 {
                x = x.negated();
            }
            check_unary(
                "tanh",
                &x,
                |bf, m| bf.tanh(m).0,
                |rx, prec, round| Float::with_val_round(prec, rx.tanh_ref(), round).0,
            );
        }
    }
}

#[test]
fn odd_functions_negative_domain_directed() {
    // The odd-function kernels negate inside their ziv_round eval
    // closure; this confirms the signed value is correctly rounded
    // under every directed mode (a negate-after-final-round bug would
    // surface here as a 1-ULP directed-mode divergence).
    let mut state = u64::from_le_bytes(*b"oddnegfn");
    let n = sample() / 2;
    for &p in PRECISIONS {
        for _ in 0..n {
            // Bounded to |x| < 16 so erf stays in the converging regime
            // (its saturation tail is covered by the saturation tests below).
            let neg = rand_bounded(&mut state, p, 3).negated();
            // sin / tan / sinh / cbrt / erf accept any negative input.
            check_unary(
                "sin",
                &neg,
                |b, m| b.sin(m).0,
                |rx, pr, r| Float::with_val_round(pr, rx.sin_ref(), r).0,
            );
            check_unary(
                "sinh",
                &neg,
                |b, m| b.sinh(m).0,
                |rx, pr, r| Float::with_val_round(pr, rx.sinh_ref(), r).0,
            );
            check_unary(
                "cbrt",
                &neg,
                |b, m| b.cbrt(m).0,
                |rx, pr, r| Float::with_val_round(pr, rx.cbrt_ref(), r).0,
            );
            check_unary(
                "erf",
                &neg,
                |b, m| b.erf(m).0,
                |rx, pr, r| Float::with_val_round(pr, rx.erf_ref(), r).0,
            );
            check_unary(
                "asinh",
                &neg,
                |b, m| b.asinh(m).0,
                |rx, pr, r| Float::with_val_round(pr, rx.asinh_ref(), r).0,
            );

            // atanh needs |x| < 1: build a negative value in (-1, 0).
            let small = scale2(&rand_pos(&mut state, p), -64).negated();
            check_unary(
                "atanh",
                &small,
                |b, m| b.atanh(m).0,
                |rx, pr, r| Float::with_val_round(pr, rx.atanh_ref(), r).0,
            );
        }
    }
}

#[test]
fn special_case_constants_directed_modes() {
    // The irrational-constant special cases the 1.2.0 negate-without
    // -mirror fix touched, re-pinned under all five modes vs an
    // independent constant oracle. A regression here would mean a
    // directed constant landed on the wrong side of its grid.
    for &p in PRECISIONS {
        let pi = Float::with_val(512, Constant::Pi);
        let pi_2 = Float::with_val(512, &pi) / 2.0f64;
        let pi_4 = Float::with_val(512, &pi) / 4.0f64;
        let three_pi_4 = Float::with_val(512, &pi) * 3.0f64 / 4.0f64;

        let pos_one = BigFloat::from_f64(1.0)
            .round_to_precision(p, RoundingMode::NearestEven)
            .unwrap()
            .0;
        let neg_one = pos_one.negated();
        let zero = BigFloat::try_new_zero(Sign::Positive, p).unwrap();
        let pos_inf = BigFloat::try_new_infinity(Sign::Positive, p).unwrap();
        let neg_inf = BigFloat::try_new_infinity(Sign::Negative, p).unwrap();

        for &mode in BIT_EXACT_ROUNDING_MODES {
            // asin(±1) = ±π/2
            assert_eq!(
                bigfloat_to_rug(&pos_one.asin(mode).0),
                round_const(&pi_2, mode, p),
                "asin(1) p={p} mode={mode:?}"
            );
            assert_eq!(
                bigfloat_to_rug(&neg_one.asin(mode).0),
                round_const(&Float::with_val(512, -&pi_2), mode, p),
                "asin(-1) p={p} mode={mode:?}"
            );
            // acos(0) = π/2, acos(-1) = π
            assert_eq!(
                bigfloat_to_rug(&zero.acos(mode).0),
                round_const(&pi_2, mode, p),
                "acos(0) p={p} mode={mode:?}"
            );
            assert_eq!(
                bigfloat_to_rug(&neg_one.acos(mode).0),
                round_const(&pi, mode, p),
                "acos(-1) p={p} mode={mode:?}"
            );
            // atan(±∞) = ±π/2
            assert_eq!(
                bigfloat_to_rug(&pos_inf.atan(mode).0),
                round_const(&pi_2, mode, p),
                "atan(+inf) p={p} mode={mode:?}"
            );
            assert_eq!(
                bigfloat_to_rug(&neg_inf.atan(mode).0),
                round_const(&Float::with_val(512, -&pi_2), mode, p),
                "atan(-inf) p={p} mode={mode:?}"
            );
            // Si(±∞) = ±π/2
            assert_eq!(
                bigfloat_to_rug(&pos_inf.si(mode).0),
                round_const(&pi_2, mode, p),
                "Si(+inf) p={p} mode={mode:?}"
            );
            assert_eq!(
                bigfloat_to_rug(&neg_inf.si(mode).0),
                round_const(&Float::with_val(512, -&pi_2), mode, p),
                "Si(-inf) p={p} mode={mode:?}"
            );
            // atan2(+∞, +∞) = π/4, atan2(+∞, -∞) = 3π/4
            assert_eq!(
                bigfloat_to_rug(&pos_inf.atan2(&pos_inf, mode).0),
                round_const(&pi_4, mode, p),
                "atan2(+inf,+inf) p={p} mode={mode:?}"
            );
            assert_eq!(
                bigfloat_to_rug(&pos_inf.atan2(&neg_inf, mode).0),
                round_const(&three_pi_4, mode, p),
                "atan2(+inf,-inf) p={p} mode={mode:?}"
            );
        }
    }
}

#[test]
fn airy_at_zero_directed_modes() {
    // Ai(0), Bi(0), Ai′(0), Bi′(0). Ai′(0) is a NEGATIVE irrational
    // constant: the textbook negate-without-mirror location. Oracle is
    // an independent rug computation from gamma + 3^a (Ai(0) =
    // 1/(3^{2/3} Γ(2/3)), Ai′(0) = -1/(3^{1/3} Γ(1/3)), Bi(0) =
    // 1/(3^{1/6} Γ(2/3)), Bi′(0) = 3^{1/6}/Γ(1/3)).
    let w = 512u32;
    let ln3 = Float::with_val(w, Float::with_val(w, 3.0f64).ln_ref());
    let pow3 = |a: &Float| -> Float {
        let t = Float::with_val(w, a * &ln3);
        Float::with_val(w, t.exp_ref())
    };
    let two_thirds = Float::with_val(w, 2.0f64) / 3.0f64;
    let one_third = Float::with_val(w, 1.0f64) / 3.0f64;
    let one_sixth = Float::with_val(w, 1.0f64) / 6.0f64;
    let g23 = two_thirds.clone().gamma();
    let g13 = one_third.clone().gamma();

    let ai0 = Float::with_val(w, 1.0f64) / (pow3(&two_thirds) * &g23);
    let aip0 = -(Float::with_val(w, 1.0f64) / (pow3(&one_third) * &g13));
    let bi0 = Float::with_val(w, 1.0f64) / (pow3(&one_sixth) * &g23);
    let bip0 = pow3(&one_sixth) / &g13;

    for &p in PRECISIONS {
        let zero = BigFloat::try_new_zero(Sign::Positive, p).unwrap();
        for &mode in BIT_EXACT_ROUNDING_MODES {
            assert_eq!(
                bigfloat_to_rug(&zero.ai(mode).0),
                round_const(&ai0, mode, p),
                "Ai(0) p={p} mode={mode:?}"
            );
            assert_eq!(
                bigfloat_to_rug(&zero.ai_prime(mode).0),
                round_const(&aip0, mode, p),
                "Ai'(0) p={p} mode={mode:?}"
            );
            assert_eq!(
                bigfloat_to_rug(&zero.bi(mode).0),
                round_const(&bi0, mode, p),
                "Bi(0) p={p} mode={mode:?}"
            );
            assert_eq!(
                bigfloat_to_rug(&zero.bi_prime(mode).0),
                round_const(&bip0, mode, p),
                "Bi'(0) p={p} mode={mode:?}"
            );
        }
    }
}

/// A `p = 24` `BigFloat` from an `f64`. The reproducer inputs are exact
/// integers / round f64s, so the lift is lossless.
fn bf24(v: f64) -> BigFloat {
    BigFloat::from_f64(v)
        .round_to_precision(24, RoundingMode::NearestEven)
        .expect("precision >= 1")
        .0
}

/// Regression pin for the directed-mode saturation bug (Slice 1 finding,
/// pf-3rtr.2; fixed in pf-3rtr.11).
///
/// Five kernels approach a nonzero, on-grid asymptotic limit (`tanh`,
/// `erf` → ±1; `erfc` → 2 as x→−∞; `expm1` → −1 as x→−∞; `zeta` → 1
/// from above). For `|x|` past where the residual underflows the maximum
/// Ziv working precision (`target + ZIV_GUARD_CAP`), the kernel used to
/// evaluate to the exact limit at every iteration, the interval test never
/// converged, the loop capped at `ZIV_MAX_ITERS`, and the cap-fallback
/// returned the limit — correct under `NearestEven` / `NearestAway` and
/// the outward directed mode, but wrong under `TowardZero` and the inward
/// directed mode, where the correctly-rounded result is the interior
/// neighbour, over an unbounded input tail. The fix short-circuits each
/// kernel to the mode-aware rounding of the limit with the residual's
/// infinitesimal. (`atan` → π/2 and `log2`/`log10` are NOT affected: their
/// limits are off-grid, so the final directed round moves off them
/// correctly.)
#[test]
fn saturation_limit_reproducers() {
    // tanh → ±1.
    check_unary(
        "tanh(8066)",
        &bf24(8066.0),
        |b, m| b.tanh(m).0,
        |rx, p, r| Float::with_val_round(p, rx.tanh_ref(), r).0,
    );
    check_unary(
        "tanh(-8066)",
        &bf24(-8066.0),
        |b, m| b.tanh(m).0,
        |rx, p, r| Float::with_val_round(p, rx.tanh_ref(), r).0,
    );
    // erf → ±1.
    check_unary(
        "erf(-4.5e8)",
        &bf24(-4.5e8),
        |b, m| b.erf(m).0,
        |rx, p, r| Float::with_val_round(p, rx.erf_ref(), r).0,
    );
    // erfc → 2 as x → −∞.
    check_unary(
        "erfc(-30)",
        &bf24(-30.0),
        |b, m| b.erfc(m).0,
        |rx, p, r| Float::with_val_round(p, rx.erfc_ref(), r).0,
    );
    // expm1 → −1 as x → −∞.
    check_unary(
        "expm1(-1500)",
        &bf24(-1500.0),
        |b, m| b.expm1(m).0,
        |rx, p, r| Float::with_val_round(p, rx.exp_m1_ref(), r).0,
    );

    // zeta → 1 from above. No rug ζ; the limit behaviour fixes the
    // correctly-rounded values by hand: the true value is 1 + tiny, so
    // TowardPositive must give the next grid point above 1 (1 + ulp) and
    // TowardZero / NearestEven must give 1.0.
    let z = bf24(5000.0);
    let one24 = Float::with_val(24, 1.0f64);
    let one_plus_ulp = Float::with_val(
        24,
        Float::with_val(53, 1.0f64) + (Float::with_val(53, 1.0f64) >> 23u32),
    );
    assert_eq!(
        bigfloat_to_rug(&z.zeta(RoundingMode::TowardPositive).0),
        one_plus_ulp,
        "zeta(5000) TowardPositive should round up to 1 + ulp"
    );
    assert_eq!(
        bigfloat_to_rug(&z.zeta(RoundingMode::TowardZero).0),
        one24,
        "zeta(5000) TowardZero should be 1.0"
    );
}

/// Sweep the saturation boundary and tail for the rug-orable kernels
/// (`tanh`, `erf`, `erfc`, `expm1`), all five modes, vs MPFR. Guards both
/// failure shapes of the fix: a *missed* saturation (the short-circuit
/// fires too late, so the Ziv loop caps and returns the on-grid limit) and
/// an *over-eager* one (it fires before the residual is within half a ULP,
/// so the infinitesimal model is invalid). The exponent sweep crosses from
/// the Ziv-converging regime, through the short-circuit threshold, into the
/// deep tail. `zeta` has no rug oracle: its tail is pinned in
/// `saturation_limit_reproducers` and its interior covered by
/// `differential_zeta` and the kernel unit tests.
#[test]
fn saturation_directed_sweep() {
    let mut state = u64::from_le_bytes(*b"satsweep");
    let n = if std::env::var("PFLOAT_DEEP").is_ok() {
        200
    } else {
        8
    };
    for &p in PRECISIONS {
        for e in 2..=28i64 {
            for _ in 0..n {
                // A magnitude in [2^e, 2^(e+1)): mantissa in [1, 2) scaled.
                let mag = scale2(&rand_bounded(&mut state, p, 0), e);
                // tanh and erf are odd: exercise both signs (each side has
                // its own inward directed mode).
                for x in [mag.clone(), mag.negated()] {
                    check_unary(
                        "tanh-tail",
                        &x,
                        |b, m| b.tanh(m).0,
                        |rx, pr, r| Float::with_val_round(pr, rx.tanh_ref(), r).0,
                    );
                    check_unary(
                        "erf-tail",
                        &x,
                        |b, m| b.erf(m).0,
                        |rx, pr, r| Float::with_val_round(pr, rx.erf_ref(), r).0,
                    );
                }
                // erfc saturates to 2 only on the negative side; expm1 to
                // −1 only on the negative side.
                check_unary(
                    "erfc-tail",
                    &mag.negated(),
                    |b, m| b.erfc(m).0,
                    |rx, pr, r| Float::with_val_round(pr, rx.erfc_ref(), r).0,
                );
                check_unary(
                    "expm1-tail",
                    &mag.negated(),
                    |b, m| b.expm1(m).0,
                    |rx, pr, r| Float::with_val_round(pr, rx.exp_m1_ref(), r).0,
                );
            }
        }
    }
}
