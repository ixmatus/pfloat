//! Shared helpers for the MPFR differential test lane.
//!
//! Each `tests/differential_<op>.rs` integration-test file does
//! `mod differential;` to pull this module in. Cargo treats the
//! file as a shared submodule of every test crate that imports it
//! (not as a standalone test binary).
//!
//! ADR-0014 records the comparison strategy. Binary radix means
//! one canonical normalized form per finite value, so bit-for-bit
//! equality of `rug::Float` values is the right test.
//!
//! Slice 7a (ADR-0016) replaces the original `Display + parse`
//! converter with a bit-exact one built on the public
//! [`pfloat::Parts`] accessor. The new converter is rounding-mode
//! independent: it reads pfloat's raw mantissa limbs and exponent
//! and constructs the corresponding [`rug::Float`] via
//! [`rug::Integer::from_digits`] and `mul_2si`. The differential
//! sweep accordingly exercises all five IEEE rounding modes.

#![cfg(all(unix, feature = "differential-mpfr"))]
// Each `tests/differential_*.rs` test crate uses a different
// subset of the helpers below; the unused subset would otherwise
// generate `dead_code` warnings under that crate's compilation.
#![allow(dead_code)]

use pfloat::{BigFloat, Parts, RoundingMode, Sign};
use rug::float::{Round, Special};
use rug::integer::Order;
use rug::{Float, Integer};

/// Map pfloat's [`RoundingMode`] to MPFR's [`Round`], or `None` when
/// MPFR has no equivalent.
///
/// `None` is returned only for [`RoundingMode::NearestAway`]: MPFR
/// has no roundTiesToAway mode. `MPFR_RNDA` (rug [`Round::AwayZero`])
/// is *directed* round-away-from-zero, which takes the farther
/// neighbour of every inexact value, not only ties; it differs from
/// IEEE 754 roundTiesToAway on every inexact non-tie value (e.g.
/// `63^-3` at p=53: roundTiesToAway ...263, `MPFR_RNDA` ...271). A
/// lane that sweeps `NearestAway` must synthesize the oracle with
/// [`round_ties_to_away`]. Returning `None` makes the absent
/// equivalent impossible to use by accident (the previous total
/// mapping silently aliased `NearestAway` to the wrong directed
/// mode; it was masked only because the four lanes that reach it
/// sweep exact arithmetic — pf-suo).
pub fn mpfr_round_of(mode: RoundingMode) -> Option<Round> {
    match mode {
        RoundingMode::NearestEven => Some(Round::Nearest),
        RoundingMode::NearestAway => None,
        RoundingMode::TowardZero => Some(Round::Zero),
        RoundingMode::TowardPositive => Some(Round::Up),
        RoundingMode::TowardNegative => Some(Round::Down),
    }
}

/// Compute an MPFR oracle value for `op` at precision `p` under the
/// requested pfloat `mode`.
///
/// For the four IEEE modes MPFR has a direct equivalent for
/// ([`Round::Nearest`], [`Round::Zero`], [`Round::Up`],
/// [`Round::Down`]), calls `op(p, round)` once at the target
/// precision. For [`RoundingMode::NearestAway`] — which MPFR has no
/// equivalent for, per [`mpfr_round_of`] returning `None` — calls
/// `op(p + 128, Round::Nearest)` to land a high-precision result
/// then routes through [`round_ties_to_away`] to synthesize the
/// IEEE roundTiesToAway value at `p`. The `+128` guard captures
/// every realistic value's two `p`-bit neighbours exactly and
/// resolves every non-tie unambiguously.
///
/// Lifted into the shared module at Phase 1f slice p1.23
/// (ADR-0038); the canonical caller pattern is `differential_pow.rs`
/// before this helper landed. Every differential lane widened to
/// [`BIT_EXACT_ROUNDING_MODES`] should route its `NearestAway` arm
/// through this helper rather than open-coding the synthesis.
pub fn mpfr_oracle_for_mode<F: Fn(u32, Round) -> Float>(
    op: F,
    mode: RoundingMode,
    p: u32,
) -> Float {
    match mpfr_round_of(mode) {
        Some(round) => op(p, round),
        None => {
            let hp = op(p + 128, Round::Nearest);
            round_ties_to_away(&hp, p)
        }
    }
}

/// IEEE 754 roundTiesToAway of a high-precision value `hp` to
/// precision `p`, synthesized because MPFR offers no such mode (see
/// [`mpfr_round_of`]).
///
/// `hp` must carry enough precision that both `p`-bit neighbours and
/// their distances to the true value are represented without further
/// rounding (the caller's responsibility — typically the operation
/// evaluated at `p + guard`). The two neighbours are obtained with
/// MPFR's directed modes; on an exact tie the away-from-zero
/// neighbour wins, and [`Round::AwayZero`] yields the correctly
/// signed larger-magnitude value so this is sign-correct for
/// negatives too.
pub fn round_ties_to_away(hp: &Float, p: u32) -> Float {
    let (lo, _) = Float::with_val_round(p, hp, Round::Zero);
    let (hi, _) = Float::with_val_round(p, hp, Round::AwayZero);
    if lo == hi {
        return lo; // exactly representable at p
    }
    let g = hp.prec();
    let d_lo = Float::with_val(g, hp - &lo).abs();
    let d_hi = Float::with_val(g, &hi - hp).abs();
    if d_hi < d_lo {
        hi
    } else if d_lo < d_hi {
        lo
    } else {
        hi // exact tie → away from zero
    }
}

/// Convert a [`BigFloat`] to a [`rug::Float`] at the same precision,
/// bit-exact regardless of which rounding mode produced the value.
///
/// Reads pfloat's raw representation via [`BigFloat::parts`] and
/// builds the corresponding [`Float`] directly: a [`rug::Integer`]
/// from the little-endian mantissa limbs, then `mul_2si` to apply
/// the signed exponent shift. The construction is exact because the
/// pfloat mantissa-as-integer has exactly `precision` significant
/// bits (top-bit set) and the destination Float is built at
/// precision `precision`.
///
/// Specials map to MPFR's special constructors. NaN payload is not
/// preserved (MPFR does not expose payload bits via the public
/// API); differential tests do not compare NaN values for bit
/// equality (NaN != NaN under IEEE) so the loss is intentional.
pub fn bigfloat_to_rug(value: &BigFloat) -> Float {
    let p = value.precision();
    match value.parts() {
        Parts::Zero { sign } => signed(Float::with_val(p, 0u32), sign),
        Parts::Infinity { sign } => signed(Float::with_val(p, Special::Infinity), sign),
        Parts::Nan { .. } => Float::with_val(p, Special::Nan),
        Parts::Normal {
            sign,
            exponent,
            mantissa,
            precision: _precision,
        } => {
            // pfloat's mantissa Vec<u64> is left-aligned within the
            // most significant limb (top bit of the highest limb is
            // 1, bits below the precision are zero). The integer
            // read out of the limbs is therefore the precision-bit
            // mantissa value shifted left by
            // `limbs * 64 - precision` bits. Setting a `p`-precision
            // Float to that integer truncates the storage padding
            // exactly (the trailing zeros below the precision do
            // not contribute), and the corrected shift below
            // accounts for the storage alignment so the result
            // recovers the original value bit-for-bit.
            let int = Integer::from_digits(mantissa, Order::Lsf);
            let mut f = Float::with_val(p, &int);
            let stored_bits = (mantissa.len() as i64) * 64;
            let shift: i64 = exponent + 1 - stored_bits;
            mul_2si_chunked(&mut f, shift);
            signed(f, sign)
        }
    }
}

fn signed(f: Float, sign: Sign) -> Float {
    if matches!(sign, Sign::Negative) {
        -f
    } else {
        f
    }
}

/// In-place `f *= 2^shift`, splitting the exponent into i32-sized
/// chunks so the helper is correct for any `i64` shift even though
/// MPFR's `mpfr_mul_2si` takes a `long` (i32 on rug's interface).
/// For the precisions and operand magnitudes the differential lane
/// exercises, this loop runs once.
fn mul_2si_chunked(f: &mut Float, shift: i64) {
    let mut remaining = shift;
    while remaining != 0 {
        let step = if remaining >= 0 {
            remaining.min(i64::from(i32::MAX)) as i32
        } else {
            remaining.max(i64::from(i32::MIN)) as i32
        };
        // Float << i32 is exact mul_2si.
        *f <<= step;
        remaining -= i64::from(step);
    }
}

/// Construct a [`BigFloat`] at the given precision from an [`i64`].
///
/// Exact for any `n` whose magnitude fits in `p` bits, which is the
/// only case slice 6a exercises.
pub fn bigfloat_from_i64(n: i64, p: u32) -> BigFloat {
    BigFloat::try_from_i64_exact(n, p).expect("i64 fits in precision")
}

/// Construct a [`rug::Float`] at the given precision from an [`i64`].
pub fn rug_from_i64(n: i64, p: u32) -> Float {
    Float::with_val(p, n)
}

/// Construct a [`BigFloat`] equal to `2^exp` at precision `p`.
///
/// A power of two is a single set bit, so repeated multiply/divide by
/// two only shifts the stored exponent and never rounds; the value is
/// exact at any precision. The integer-input lanes cannot reach the
/// tiny-x activation band of the small-argument fast-paths (ADR-0059);
/// this and [`rug_pow2`] feed the dedicated band sweep.
pub fn bigfloat_pow2(exp: i64, p: u32) -> BigFloat {
    let two = BigFloat::try_from_i64_exact(2, p).expect("precision >= 1");
    let mut v = BigFloat::try_from_i64_exact(1, p).expect("precision >= 1");
    for _ in 0..exp.unsigned_abs() {
        v = if exp < 0 {
            v.div(&two, RoundingMode::NearestEven).0
        } else {
            v.mul(&two, RoundingMode::NearestEven).0
        };
    }
    v
}

/// Construct a [`rug::Float`] equal to `2^exp` at precision `p`, exact
/// via `mul_2si` (split into i32 chunks like [`mul_2si_chunked`]).
/// Companion to [`bigfloat_pow2`] for the tiny-x band sweep.
pub fn rug_pow2(exp: i64, p: u32) -> Float {
    let mut f = Float::with_val(p, 1u32);
    let mut remaining = exp;
    while remaining != 0 {
        let step = if remaining >= 0 {
            remaining.min(i64::from(i32::MAX)) as i32
        } else {
            remaining.max(i64::from(i32::MIN)) as i32
        };
        f <<= step;
        remaining -= i64::from(step);
    }
    f
}

/// Splitmix64 step. Used by each `differential_*` test for
/// deterministic input generation; consolidated here so the
/// helper isn't duplicated 23 times and so the i64 range math is
/// fixed in one place.
pub fn next_u64(state: &mut u64) -> u64 {
    *state = state.wrapping_add(0x9E37_79B9_7F4A_7C15);
    let mut z = *state;
    z = (z ^ (z >> 30)).wrapping_mul(0xBF58_476D_1CE4_E5B9);
    z = (z ^ (z >> 27)).wrapping_mul(0x94D0_49BB_1331_11EB);
    z ^ (z >> 31)
}

/// Draw an i64 uniformly from `[lo, hi]` using the splitmix64
/// state. Uses i128 arithmetic so the span can cover the full i64
/// range without overflow.
pub fn next_i64_in(state: &mut u64, lo: i64, hi: i64) -> i64 {
    debug_assert!(lo <= hi);
    let span = (i128::from(hi) - i128::from(lo) + 1) as u64;
    let offset = next_u64(state) % span;
    (i128::from(lo) + i128::from(offset)) as i64
}

/// Number of random pairs exercised per `(precision, mode)` cell in
/// each cargo-test run. The deep sweep (10⁶) runs locally under
/// `PFLOAT_DEEP=1` per ADR-0014.
pub fn sweep_size() -> u32 {
    if std::env::var("PFLOAT_DEEP").is_ok() {
        1_000_000
    } else {
        10_000
    }
}

/// The four precisions exercised by the CI sweep for arithmetic and
/// parsing, where pfloat is expected to be bit-exact against MPFR at
/// every tested precision.
pub const SWEEP_PRECISIONS: &[u32] = &[53, 113, 256, 1024];

/// Precisions used by transcendental and tier-1 special function
/// differential tests.
///
/// Slice 6h originally capped this at 256 bits because pfloat's
/// transcendentals used hardcoded 1024-bit constants (`ln(2)`,
/// `π`, `2/π`, etc.) for argument reduction and ADR-0014 attributed
/// p>256 divergence to the 64-bit guard exhausting the constant's
/// reach. Slice 7b (ADR-0017) discovered the underlying problem was
/// a faulty 1024-bit `LN2_LIMBS_1024` encoding (correct only to
/// ~450 bits) and replaced it with AGM-based on-the-fly computation
/// (Brent–Salamin for π, atanh series for ln(2)). The kernels can
/// now correctly compute transcendentals at any target precision.
///
/// Slice 7b deferred the p=1024 lane because AGM-based constants were
/// recomputed on every call, making the sweep hour-scale. That
/// follow-up has landed: `math::agm_constants` now memoizes `π`,
/// `ln(2)`, etc. per `(kind, precision)` in a thread-local table, so
/// each distinct working precision is computed once per sweep rather
/// than per iteration. The lane is restored to include 1024, matching
/// `SWEEP_PRECISIONS`.
pub const TRANSCENDENTAL_PRECISIONS: &[u32] = &[53, 113, 256, 1024];

/// All five IEEE 754-2019 rounding modes. Used by the differential
/// tests for operations whose pfloat kernel is bit-exact correctly
/// rounded under any mode: add, sub, mul, div, sqrt, fma, and
/// decimal parsing. Slice 7a (ADR-0016) unblocks this list by
/// landing the bit-exact `BigFloat`→`rug::Float` converter; the
/// old `Display + parse` route lost up to 1 ULP under non-NE rounding
/// and obscured the bit-exact kernel agreement.
pub const BIT_EXACT_ROUNDING_MODES: &[RoundingMode] = &[
    RoundingMode::NearestEven,
    RoundingMode::NearestAway,
    RoundingMode::TowardZero,
    RoundingMode::TowardPositive,
    RoundingMode::TowardNegative,
];

/// Rounding modes exercised by differential tests for operations
/// whose pfloat kernel is correctly rounded only under `NearestEven`:
/// all transcendentals (exp, ln, sin, cos, tan, atan2, sinh, cosh,
/// asinh), all tier-1 specials (erf, erfc, gamma, lgamma, digamma,
/// beta), and AGM.
///
/// These kernels compute at working precision `target + 64` under
/// `NearestEven` and apply the user's mode only at the final round
/// to target precision. The fixed 64-bit guard is not Ziv-strategy
/// retry: at tie cases the final round under non-NearestEven modes
/// can diverge from MPFR's correctly-rounded result by up to 1 ULP.
/// Slice 7c (ADR-0022) shipped the Ziv recompute-and-compare driver
/// for `pow`, which is therefore the first kernel off this tier and
/// uses [`BIT_EXACT_ROUNDING_MODES`]; the rest stay under
/// `NearestEven` only until later slices extend the driver, the same
/// correctness floor that ships today.
pub const NEAREST_EVEN_ROUNDING_MODES: &[RoundingMode] = &[RoundingMode::NearestEven];

/// Back-compat alias preserved during the slice 7a rollout. New
/// tests should pick [`BIT_EXACT_ROUNDING_MODES`] or
/// [`NEAREST_EVEN_ROUNDING_MODES`] explicitly based on what the
/// pfloat kernel under test guarantees.
pub const ALL_ROUNDING_MODES: &[RoundingMode] = NEAREST_EVEN_ROUNDING_MODES;

#[cfg(test)]
mod converter_tests {
    use super::*;

    // Sanity test: converting a BigFloat built by arithmetic round
    // trips through rug bit-exactly when both sides perform the same
    // operation under NearestEven. If this passes but
    // `differential_div` fails under non-NearestEven, the converter
    // is sound and the divergence is a kernel correctness gap.
    #[test]
    fn converter_round_trips_under_nearest_even() {
        let a = BigFloat::try_from_i64_exact(3, 53).unwrap();
        let b = BigFloat::try_from_i64_exact(7, 53).unwrap();
        let (q_bf, _) = a.div(&b, RoundingMode::NearestEven);
        let q_rug = bigfloat_to_rug(&q_bf);
        let a_r = rug_from_i64(3, 53);
        let b_r = rug_from_i64(7, 53);
        let (q_direct, _) = Float::with_val_round(53, &a_r / &b_r, Round::Nearest);
        assert_eq!(
            q_rug, q_direct,
            "3/7 NE: bf->rug={q_rug}, direct={q_direct}"
        );
    }

    // The slice-6h problem case under NearestEven: must round-trip.
    #[test]
    fn converter_round_trips_slice_6h_case_ne() {
        let a = BigFloat::try_from_i64_exact(-966132233652331, 53).unwrap();
        let b = BigFloat::try_from_i64_exact(1233101814760529, 53).unwrap();
        let (q_bf, _) = a.div(&b, RoundingMode::NearestEven);
        let q_rug = bigfloat_to_rug(&q_bf);
        let a_r = rug_from_i64(-966132233652331, 53);
        let b_r = rug_from_i64(1233101814760529, 53);
        let (q_direct, _) = Float::with_val_round(53, &a_r / &b_r, Round::Nearest);
        assert_eq!(
            q_rug, q_direct,
            "slice-6h NE: bf->rug={q_rug}, direct={q_direct}"
        );
    }

    // Integer round-trips at several precisions.
    #[test]
    fn converter_round_trips_integers() {
        for n in [1i64, -1, 7, -42, 1_000_000, -(1 << 50)] {
            for p in [53u32, 113, 256, 1024] {
                let bf = bigfloat_from_i64(n, p);
                let bf_as_rug = bigfloat_to_rug(&bf);
                let direct = rug_from_i64(n, p);
                assert_eq!(
                    bf_as_rug, direct,
                    "i64 {n} at p={p}: bf->rug={bf_as_rug}, direct={direct}"
                );
            }
        }
    }
}
