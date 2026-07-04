//! Regression guards for the 2026-07-04 remediation arc R5.1
//! (epic pf-8iji, pf-1vzg, ADR-0125): near a function zero the large-x
//! oscillatory asymptotic `f·sin ∓ g·cos` (Ci) / `cos φ·P ± sin φ·Q`
//! (Airy) / `[±cosω, ±sinω, …]` bracket (Bessel J/Y) cancels
//! catastrophically. The asymptotic is a DIVERGENT series with an
//! irreducible truncation floor `≈ e^{−c·x}`; below it the value cannot be
//! computed at any working precision, so the pre-fix kernel CERTIFIED THE
//! FLOOR VALUE — the certified-wrong-at-a-zero family (RC2/ADR-0097). The
//! fix detects the below-floor condition and hands off to a convergent
//! series (Ci→`ci_series`, J→`bessel_j_tiny` Maclaurin, Y→`bessel_y_series`
//! log series), which has no floor and resolves any depth; Airy's floor is
//! deeper than the guard cap so its `cancellation_boosted` wrapper already
//! resolves the reachable gap.
//!
//! Each kernel gets two rows sharing one function zero:
//! - `*_deep_certifies` (release-gated): input parsed ~1300 bits from the
//!   zero, so `|f(input)| ~ 2^-1300` — past every kernel's floor. The
//!   pre-fix asymptotic returned `≈ floor` (Ci `2^-148`, Bessel `2^-547`),
//!   wrong by hundreds of orders / often the wrong sign; the convergent
//!   fallback recovers the true value. (Airy's floor `2^-1593` is past
//!   this D, so its wrapper alone suffices and this row exercises it.)
//! - `*_shallow_control` (always runs): input ~245 bits from the zero,
//!   within every floor, so the asymptotic (or its recovery) is correct.
//!   The control proves the fix leaves the common near-zero path right.
//!
//! Oracles: `mpmath` 1.4.1 at `/usr/bin/python3`, 5000-bit working, on
//! the BIT-IDENTICAL dyadic input (integer mantissa × 2^exp — the same
//! value pfloat reconstructs via `parse_str` + `scale_by_pow2`, so no
//! decimal-rounding step sits between pfloat and the oracle). Depths and
//! constants: `scratchpad/gen_r51_oracles_v2.py`.
//!
//! Run (release, where the deep rows are active): `cargo test --release
//! --test regression_review_2026_07_04_r5 --features
//! std,big,trig,integrals,airy,bessel`. The CI MPFR full-union job
//! (`--release --features=differential-mpfr`) exercises it.

#![cfg(all(
    feature = "big",
    feature = "integrals",
    feature = "airy",
    feature = "bessel"
))]

use core::cmp::Ordering;

use pfloat::{BigFloat, RoundingMode};

const NE: RoundingMode = RoundingMode::NearestEven;

/// Reconstruct the exact dyadic `mantissa · 2^exp2` bit-for-bit. The
/// mantissa is an integer that fits in `prec` bits, so `parse_str` holds
/// it exactly and `scale_by_pow2` only shifts the exponent; the result
/// is the identical value the mpmath oracle was evaluated at.
fn dyadic(mantissa: &str, exp2: i64, prec: u32) -> BigFloat {
    let (m, pst) = BigFloat::parse_str(mantissa, prec, NE).unwrap();
    assert!(
        pst.is_ok(),
        "mantissa integer must parse exactly at {prec} bits"
    );
    let (x, sst) = m.scale_by_pow2(exp2);
    assert!(sst.is_ok(), "scale_by_pow2 is exact (power-of-two shift)");
    x
}

/// `got` bit-identical to `reference` rounded to `p` bits under `mode`
/// (the R1 `assert_bit_exact` idiom).
fn assert_bit_exact(label: &str, got: &BigFloat, reference: &str, p: u32, mode: RoundingMode) {
    let expected = BigFloat::parse_str(reference, p, mode).unwrap().0;
    assert_eq!(
        got.total_cmp(&expected),
        Ordering::Equal,
        "{label}: got {got}, want {expected}"
    );
}

// ---------------------------------------------------------------------
// Ci: the named site (ci.rs, DLMF 6.12.4). Zero ~ 100.5409068603591 in
// the asymptotic regime (e_x = 6 ≥ threshold 6 at target 53).
// ---------------------------------------------------------------------

#[cfg(not(debug_assertions))]
const CI_ZERO_DEEP_MAN: &str = "17144593471777766391103144555246797526259976398238155235390221541021307194687664997994397453248013338662552787016457590818281236214048969693537926852239617168749471889176115188895418739622029363192039477176321304074672905221440189218896574432311742500363245927339624317078934384918725311015045966005081038910731903869432991761634191637707532584937621534793994147176070018572784911560425034579";

/// Ci near a large-x zero, ~1302-bit proximity: D exceeds the 1024 cap,
/// so the un-boosted asymptotic path certified a wrong value.
/// Release-gated (~seconds).
#[cfg(not(debug_assertions))]
#[test]
fn ci_asymptotic_near_zero_deep_certifies() {
    let x = dyadic(CI_ZERO_DEEP_MAN, -1293, 1400);
    let (r, st) = x.ci_round(53, NE).unwrap();
    assert!(
        !r.is_sign_negative() && !r.is_zero() && !r.is_infinite(),
        "Ci(deep zero) must be a tiny positive normal, got {r}"
    );
    assert_bit_exact(
        "Ci(deep zero)",
        &r,
        "1.0953495180322616669467626098655598310304937857209e-392",
        53,
        NE,
    );
    assert!(st.inexact());
}

/// Ci near the same zero at ~252-bit proximity (recovered by the legacy
/// guard doubling): the boost must leave it correct.
#[test]
fn ci_asymptotic_near_zero_shallow_control() {
    let x = dyadic(
        "1421123249411058809752764872096889576279965292226880261464773671365819311615",
        -243,
        400,
    );
    let (r, st) = x.ci_round(53, NE).unwrap();
    assert_bit_exact(
        "Ci(shallow zero)",
        &r,
        "1.8262011909521957157559334590568133377979463611235e-76",
        53,
        NE,
    );
    assert!(st.inexact());
}

// ---------------------------------------------------------------------
// Airy Ai on the negative axis (airy_asymptotic_neg, DLMF 9.7.9). Zero
// ~ -139.5294620819492 (|x| ≥ threshold 128 at target 53). Input
// negative, so the positive mantissa is negated.
// ---------------------------------------------------------------------

#[cfg(not(debug_assertions))]
const AI_ZERO_DEEP_MAN: &str = "11896530374712727142169882224157771188164063490421272054322116055745298106734385871016869933328398455296661633645510980692772399875871964071284482893194876842085765450579993795755893420277638368241304318476948088149386714838986928655645412198283760081632493183232286503963775358134932584560971681802118284711304639549719053525455588055433539909009016247760913684239686556912794997843359529893";

/// Ai(−t) near a large negative-axis zero, ~1294-bit proximity: the
/// oscillatory `cos φ·Pu + sin φ·Qu` cancels past the cap.
/// Release-gated.
#[cfg(not(debug_assertions))]
#[test]
fn ai_asymptotic_near_zero_deep_certifies() {
    let x = dyadic(AI_ZERO_DEEP_MAN, -1292, 1400).negated();
    let (r, st) = x.ai_round(53, NE).unwrap();
    assert!(
        !r.is_sign_negative() && !r.is_zero() && !r.is_infinite(),
        "Ai(deep neg zero) must be a tiny positive normal, got {r}"
    );
    assert_bit_exact(
        "Ai(deep neg zero)",
        &r,
        "3.7667058372709013283719563680797866317930732444075e-390",
        53,
        NE,
    );
    assert!(st.inexact());
}

/// Ai(−t) near the same zero at ~242-bit proximity (recovered).
#[test]
fn ai_asymptotic_near_zero_shallow_control() {
    let x = dyadic(
        "123263610264799418318823510914390753649868494355360361809939630052190167057",
        -239,
        400,
    )
    .negated();
    let (r, st) = x.ai_round(53, NE).unwrap();
    assert_bit_exact(
        "Ai(shallow neg zero)",
        &r,
        "-1.2523692857549240985092580923659998559362294859798e-73",
        53,
        NE,
    );
    assert!(st.inexact());
}

// ---------------------------------------------------------------------
// Bessel J0 (bessel_j_asymptotic, DLMF 10.17.3). Zero ~ 150.0118824569
// (e_x = 7 ≥ threshold 7 at target 53).
// ---------------------------------------------------------------------

#[cfg(not(debug_assertions))]
const J0_ZERO_DEEP_MAN: &str = "6395140100120375282372916601820345333476318323350462096075038286535816279798966625549124879204499344556690492978652865349155151747536911776941015760525700203534591783666725584197796140467066421660161994887351611900475944502774338430734340955479974217437184943686152133653990333421484490843696928134289520605231232719349203264033506841327435724995495800935052882062500974895164229057574458919";

/// J0 near a large zero, ~1298-bit proximity: the `[+cosω, −sinω, …]`
/// bracket cancels past the cap. Release-gated.
#[cfg(not(debug_assertions))]
#[test]
fn j0_asymptotic_near_zero_deep_certifies() {
    let x = dyadic(J0_ZERO_DEEP_MAN, -1291, 1400);
    let (r, st) = x.j0_round(53, NE).unwrap();
    assert!(
        r.is_sign_negative() && !r.is_zero() && !r.is_infinite(),
        "J0(deep zero) must be a tiny negative normal, got {r}"
    );
    assert_bit_exact(
        "J0(deep zero)",
        &r,
        "-1.8087591060810815773762375553715845116080867661123e-391",
        53,
        NE,
    );
    assert!(st.inexact());
}

/// J0 near the same zero at ~248-bit proximity (recovered).
#[test]
fn j0_asymptotic_near_zero_shallow_control() {
    let x = dyadic(
        "530096108401901036418758653721902090518058637028735652397231676648780552125",
        -241,
        400,
    );
    let (r, st) = x.j0_round(53, NE).unwrap();
    assert_bit_exact(
        "J0(shallow zero)",
        &r,
        "2.9262788138201282224970637035086754673183550342854e-75",
        53,
        NE,
    );
    assert!(st.inexact());
}

// ---------------------------------------------------------------------
// Bessel Y0 (bessel_y_asymptotic, DLMF 10.17.4). Zero ~ 148.4410949470
// (e_x = 7 ≥ threshold 6 at target 53).
// ---------------------------------------------------------------------

#[cfg(not(debug_assertions))]
const Y0_ZERO_DEEP_MAN: &str = "6328176030150676961761385507768860958729574439643510697303734825645390398939554645674187003552279567776140360861936596959149385443314723965569199299678923834095732967395924909787531433967594515565816492339811174496577035470059431806151041577997384867718802305819911206221435563551678231579751552607776955673035837291251370354545106841946688286103997560626703652131445554487785401755341293085";

/// Y0 near a large zero, ~1298-bit proximity: the DLMF 10.17.4 bracket
/// cancels past the cap. Release-gated.
#[cfg(not(debug_assertions))]
#[test]
fn y0_asymptotic_near_zero_deep_certifies() {
    let x = dyadic(Y0_ZERO_DEEP_MAN, -1291, 1400);
    let (r, st) = x.y0_round(53, NE).unwrap();
    assert!(
        !r.is_sign_negative() && !r.is_zero() && !r.is_infinite(),
        "Y0(deep zero) must be a tiny positive normal, got {r}"
    );
    assert_bit_exact(
        "Y0(deep zero)",
        &r,
        "1.7546974349779011811214319165273389947164176333217e-391",
        53,
        NE,
    );
    assert!(st.inexact());
}

/// Y0 near the same zero at ~248-bit proximity (recovered).
#[test]
fn y0_asymptotic_near_zero_shallow_control() {
    let x = dyadic(
        "1049090851598989273872152425314966008395572028100259158758535461598266362973",
        -242,
        400,
    );
    let (r, st) = x.y0_round(53, NE).unwrap();
    assert_bit_exact(
        "Y0(shallow zero)",
        &r,
        "3.0664562385457375180025091160954898486746787787157e-75",
        53,
        NE,
    );
    assert!(st.inexact());
}
