//! MPFR differential for `BigFloat::zeta` (Riemann zeta, real
//! argument; the `differential_yn` bit-exact precedent).
//!
//! rug 1.30 exposes MPFR `mpfr_zeta` (`Float::zeta_ref()`), a
//! genuine external oracle, so this is a **bit-exact** lane
//! (`assert_eq!`, the `differential_ei` / `differential_jn` idiom)
//! under `NearestEven` across the full `TRANSCENDENTAL_PRECISIONS`
//! including `p = 1024` (the user-confirmed posture; the oracle
//! fork resolved to a native MPFR primitive, unlike 6q's tiered
//! `differential_ik`).
//!
//! Domain: every sweep argument is an exact dyadic rational
//! (power-of-two denominator) so the constructed `s` is
//! bit-identical in pfloat and rug at every precision. The points
//! span the three code paths — `s > 1` and `0 < s < 1` (the
//! Borwein eta-acceleration core, DLMF 25.2.3) and `s < 0` (the
//! functional equation DLMF 25.4.2 reflecting into the Borwein
//! region) — and stay clear of the pole at `s = 1` and of the
//! negative even integers (the exactly-zero trivial zeros, where
//! both sides are `0` but which carry no rounding information).
//!
//! Cost note (`feedback_differential_lane_cost`): the negative
//! axis composes `Γ`+`sin`+`pow`+Borwein, so each `s < 0` point at
//! `p = 1024` is the cost driver and this lane is a deliberately
//! slow CI tier (the `differential_yn` / `differential_si`
//! posture). The sweep is a small bounded dyadic set, not a wide
//! random sweep. Validate a fast subset locally (one region at
//! `p = 53`); the full lane is the CI slow tier.

#![cfg(all(unix, feature = "differential-mpfr"))]

mod differential;

use differential::{
    bigfloat_from_i64, bigfloat_to_rug, mpfr_oracle_for_mode, rug_from_i64,
    BIT_EXACT_ROUNDING_MODES, NEAREST_EVEN_ROUNDING_MODES, TRANSCENDENTAL_PRECISIONS,
};
use pfloat::RoundingMode;

/// Exact dyadic-rational inputs `(numerator, denominator)`, the
/// denominator a power of two so the argument is bit-identical in
/// pfloat and rug at every precision. Clear of `s = 1` (pole) and
/// the negative even integers (trivial zeros). Spans `s > 1`
/// (`5/2, 9/4, 17/8, 3, 7`), the critical strip `0 < s < 1`
/// (`1/2, 1/4, 3/4, 7/8`), and `s < 0` exercising the functional
/// equation including the negative odd integers (`-1, -3`) and
/// non-integers (`-1/2, -5/4, -9/4, -7/2`).
const DYADIC: &[(i64, i64)] = &[
    (5, 2),
    (9, 4),
    (17, 8),
    (3, 1),
    (7, 1),
    (1, 2),
    (1, 4),
    (3, 4),
    (7, 8),
    (-1, 2),
    (-5, 4),
    (-9, 4),
    (-7, 2),
    (-1, 1),
    (-3, 1),
];

#[test]
fn zeta_matches_mpfr() {
    for &p in TRANSCENDENTAL_PRECISIONS {
        // The `s < 0` points at p = 1024 compose Γ + sin + pow + Borwein and
        // are the differential suite's wall-clock floor: ~53 min of a ~76 min
        // release CI job, all in this one non-parallelizable test (ADR-0014
        // amendment, slice ci-green-zeta). The directed modes at p = 1024 were
        // the Phase 1f widening (ADR-0038) over zeta's original, user-confirmed
        // NearestEven-only-at-p=1024 posture; they mostly re-exercise the
        // final-round converter already covered by the five-mode tier at
        // p <= 256. Keep all five modes where ties are common (p <= 256);
        // NearestEven-only at p = 1024 (one of the three deliberate
        // NearestEven-only differential lanes named in ADR-0079, with
        // beta and parse).
        let modes: &[RoundingMode] = if p >= 1024 {
            NEAREST_EVEN_ROUNDING_MODES
        } else {
            BIT_EXACT_ROUNDING_MODES
        };
        for &mode in modes {
            for &(num, den) in DYADIC {
                let x_bf = if den == 1 {
                    bigfloat_from_i64(num, p)
                } else {
                    let (q, _) = bigfloat_from_i64(num, p)
                        .div(&bigfloat_from_i64(den, p), RoundingMode::NearestEven);
                    q
                };

                let bf_r = {
                    let (r, _s) = x_bf.zeta(mode);
                    bigfloat_to_rug(&r)
                };
                let rug_r = mpfr_oracle_for_mode(
                    |prec, round| {
                        let xr = if den == 1 {
                            rug_from_i64(num, prec)
                        } else {
                            rug_from_i64(num, prec) / rug_from_i64(den, prec)
                        };
                        rug::Float::with_val_round(prec, xr.zeta_ref(), round).0
                    },
                    mode,
                    p,
                );
                assert_eq!(bf_r, rug_r, "ζ({num}/{den}) at p={p}, mode={mode:?}");
            }
        }
    }
}
