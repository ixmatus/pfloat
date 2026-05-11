//! MPFR differential: `BigFloat::parse_str` agrees with
//! `rug::Float::parse` for canonical decimal strings.
//!
//! Both pfloat and rug round to the requested precision under
//! NearestEven; the bit-for-bit comparison is the standard test.

#![cfg(all(unix, feature = "differential-mpfr"))]

mod differential;

use differential::{bigfloat_to_rug, mpfr_round_of, ALL_ROUNDING_MODES, SWEEP_PRECISIONS};

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

#[test]
fn parse_matches_mpfr_on_canonical_strings() {
    for &p in SWEEP_PRECISIONS {
        for &s in STRINGS {
            for &mode in ALL_ROUNDING_MODES {
                let bf_r = {
                    let (parsed, _status) = BigFloat::parse_str(s, p, mode)
                        .unwrap_or_else(|e| panic!("pfloat parse failed for {s:?}: {e:?}"));
                    bigfloat_to_rug(&parsed)
                };
                let rug_r = {
                    let parsed = rug::Float::parse(s)
                        .unwrap_or_else(|e| panic!("rug parse failed for {s:?}: {e:?}"));
                    let (r, _ord) = rug::Float::with_val_round(p, parsed, mpfr_round_of(mode));
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
