//! MPFR differential: `BigFloat::parse_str` agrees with
//! `rug::Float::parse` for canonical decimal strings.
//!
//! Both pfloat and rug round to the requested precision under
//! `NearestEven`; the bit-for-bit comparison is the standard test.
//!
//! This lane is `NearestEven`-only by design: it pins the canonical
//! decimal round-trip bit-for-bit. It is one of the three deliberate
//! `NearestEven`-only differential lanes named in ADR-0079 (with
//! `beta`, a loose two-ULP oracle, and `zeta` at p = 1024 for cost);
//! directed-mode decimal parsing is a separate concern from the
//! correct-rounding sweep this phase verifies.

#![cfg(all(unix, feature = "differential-mpfr"))]

mod differential;

use differential::{bigfloat_to_rug, mpfr_round_of, NEAREST_EVEN_ROUNDING_MODES, SWEEP_PRECISIONS};

use pfloat::{BigFloat, RoundingMode};

/// A hand-curated battery of decimal strings that exercise the
/// parser's edge cases (sign, exponent forms, leading zeros,
/// trailing zeros, fractional only, integer only).
const STRINGS: &[&str] = &[
    "0",
    "-0",
    "1",
    "-1",
    "0.5",
    "-0.5",
    "1.5",
    "0.0",
    "10",
    "100",
    "1e10",
    "1e-10",
    "1.23456789",
    "3.14159265358979323846",
    "2.71828182845904523536",
    "0.000001",
    "1000000.000001",
    "1.5e20",
    "-1.5e-20",
];

/// Large in-range decimal exponents: well inside the
/// `MAX_DECIMAL_EXPONENT = 10^6` cost-budget cap (ADR-0031, amended
/// for the parse-oom slice), so they parse to correct finite values
/// rather than saturating. The property the test pins is that pfloat
/// and rug/MPFR agree bit-exact at a large exponent that exercises
/// the multi-limb `pow5` and the Algorithm-D divmod, not just the
/// small canonical strings.
///
/// `10^5` is deliberate: MPFR's `strtod` cost climbs steeply with the
/// decimal exponent (a parse near the old `~3 * 10^6` band ran into
/// the tens of minutes locally), and `10^5` keeps both sides fast
/// while staying an order of magnitude into the large-`pow5` regime.
/// Tested separately from [`STRINGS`] so they hit only the bit-exact
/// comparison and not the Display round-trip.
const BOUNDARY_STRINGS: &[&str] = &["1e100000", "1e-100000"];

#[test]
fn parse_matches_mpfr_on_canonical_strings() {
    for &p in SWEEP_PRECISIONS {
        for &s in STRINGS {
            for &mode in NEAREST_EVEN_ROUNDING_MODES {
                let bf_r = {
                    let (parsed, _status) = BigFloat::parse_str(s, p, mode)
                        .unwrap_or_else(|e| panic!("pfloat parse failed for {s:?}: {e:?}"));
                    bigfloat_to_rug(&parsed)
                };
                let rug_r = {
                    let parsed = rug::Float::parse(s)
                        .unwrap_or_else(|e| panic!("rug parse failed for {s:?}: {e:?}"));
                    let (r, _ord) = rug::Float::with_val_round(
                        p,
                        parsed,
                        mpfr_round_of(mode)
                            .expect("NE-only lane: NearestEven has an MPFR equivalent (pf-suo)"),
                    );
                    r
                };
                assert_eq!(
                    bf_r, rug_r,
                    "parse({s:?}) at p={p}, mode={mode:?}: pfloat={bf_r}, rug={rug_r}"
                );
            }
        }
    }
}

/// Bit-exact parse comparison at a large in-range decimal exponent
/// (`10^5`, well inside the `MAX_DECIMAL_EXPONENT = 10^6` cap). Single
/// precision `p = 113`: the property the test pins is that pfloat and
/// rug/MPFR produce the same correctly rounded value when the
/// conversion goes through a multi-limb `pow5` and the Algorithm-D
/// divmod; matching across the full sweep precision ladder would
/// multiply the per-string cost without adding coverage.
#[test]
fn parse_matches_mpfr_at_large_in_range_exponent() {
    let p: u32 = 113;
    let mode = RoundingMode::NearestEven;
    for &s in BOUNDARY_STRINGS {
        let bf_r = {
            let (parsed, _status) = BigFloat::parse_str(s, p, mode)
                .unwrap_or_else(|e| panic!("pfloat parse failed for {s:?}: {e:?}"));
            bigfloat_to_rug(&parsed)
        };
        let rug_r = {
            let parsed = rug::Float::parse(s)
                .unwrap_or_else(|e| panic!("rug parse failed for {s:?}: {e:?}"));
            let (r, _ord) = rug::Float::with_val_round(
                p,
                parsed,
                mpfr_round_of(mode).expect("NE has an MPFR equivalent"),
            );
            r
        };
        assert_eq!(
            bf_r, rug_r,
            "parse({s:?}) at p={p}: pfloat={bf_r}, rug={rug_r}"
        );
    }
}

/// Parse-format round-trip for `BigFloat::parse_str(s, p, NE)` →
/// Display → re-parse at the same precision: should reproduce
/// numerically equal values.
#[test]
fn parse_display_roundtrip_preserves_value() {
    let prec: u32 = 113;
    for &s in STRINGS {
        let (v, _status) = BigFloat::parse_str(s, prec, RoundingMode::NearestEven)
            .unwrap_or_else(|e| panic!("pfloat parse failed for {s:?}: {e:?}"));
        let rendered = v.to_string();
        let (reparsed, _) = BigFloat::parse_str(&rendered, prec, RoundingMode::NearestEven)
            .unwrap_or_else(|e| panic!("pfloat re-parse failed for {rendered:?}: {e:?}"));
        let (cmp, _) = v.partial_cmp(&reparsed);
        assert_eq!(
            cmp,
            Some(std::cmp::Ordering::Equal),
            "round-trip mismatch for {s:?} → {rendered:?} → {reparsed}"
        );
    }
}
