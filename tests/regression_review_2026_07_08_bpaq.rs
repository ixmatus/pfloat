//! Regression guard for pf-bpaq (ADR-0128): the parse/Display round-trip
//! near the `±MAX_DECIMAL_EXPONENT` (10^6) cap.
//!
//! `parse` used to cap on the point-baked exponent `e_part - frac_digits`,
//! so a value's own `round_trip_digit_count`-digit `Display` output (whose
//! baked re-parse exponent sits `frac_digits` below the value magnitude)
//! could exceed the cap and be saturated to `0`/`inf` — a finite value that
//! did NOT round-trip (the libFuzzer `parse` round-trip assertion caught it).
//! An independent `Fraction` oracle confirmed those values DO round-trip
//! through 36 digits, so parse was wrong.
//!
//! The fix caps representability on the value MAGNITUDE and widens the pow5
//! cost budget (parse) / render cap (fmt) by `round_trip_digit_count`, so
//! every value parse produces renders exactly and round-trips.
//!
//! Run: `cargo test --release --features std,big,fmt
//! --test regression_review_2026_07_08_bpaq`.

#![cfg(all(feature = "big", feature = "fmt"))]

use core::cmp::Ordering;

use pfloat::{BigFloat, RoundingMode};

const NE: RoundingMode = RoundingMode::NearestEven;

/// Assert `s` parses to a finite value that round-trips through Display.
fn assert_roundtrips(s: &str, p: u32) {
    let (parsed, _) = BigFloat::parse_str(s, p, NE).unwrap();
    assert!(
        !parsed.is_infinite() && !parsed.is_zero() && !parsed.is_nan(),
        "{s}: expected a finite nonzero value, got {parsed}"
    );
    let rendered = format!("{parsed}");
    let (reparsed, _) = BigFloat::parse_str(&rendered, p, NE).unwrap();
    assert_eq!(
        parsed.partial_cmp(&reparsed).0,
        Some(Ordering::Equal),
        "{s} -> Display {rendered} -> re-parse must equal the original"
    );
}

#[test]
fn tiny_near_cap_roundtrips() {
    // The exact fuzz failure family: tiny values whose 36-digit render has a
    // baked re-parse exponent past the old cap. All must round-trip now.
    for e in [-999_990i64, -999_999, -1_000_000] {
        for m in ["1", "3", "7", "9"] {
            assert_roundtrips(&format!("{m}e{e}"), 113);
        }
    }
    // The `1e-1000000` case: the log10 estimate off-by-one wrongly routed it
    // to the approximate token pre-fix.
    assert_roundtrips("1e-1000000", 113);
}

#[test]
fn value_magnitude_past_cap_saturates_consistently() {
    // 99e1000000 = 9.9e1000001: value magnitude 1000001 > cap, so parse now
    // overflows to inf (consistent with fmt, which cannot render it exactly).
    // Display("inf") re-parses to inf, so there is no finite mismatch.
    for s in ["99e1000000", "100e1000000", "999e999999", "12345e999998"] {
        let (v, st) = BigFloat::parse_str(s, 113, NE).unwrap();
        assert!(
            v.is_infinite() && st.overflow(),
            "{s}: magnitude past cap must overflow to inf, got {v}"
        );
    }
}

#[test]
fn cap_boundary_is_on_magnitude_not_baked_exponent() {
    // A value at exactly magnitude ±cap stays finite regardless of how many
    // significant digits the input carries (precision-independent cap).
    // Same magnitude 10^6 as a 1-digit and a 7-digit integer input:
    assert_roundtrips("1e1000000", 113);
    assert_roundtrips("1000000e999994", 113);
    // The same magnitude at 40 significant digits (its baked exponent lands
    // far past the old cap) still parses finite and round-trips.
    let mut long = String::from("1.");
    for _ in 0..39 {
        long.push('0');
    }
    long.push_str("e1000000");
    assert_roundtrips(&long, 113);
}

#[test]
fn clearly_past_cap_still_saturates_mode_aware() {
    // Existing pf-mw6u behaviour is preserved (2x past the cap).
    let (a, _) = BigFloat::parse_str("1e2000000", 113, NE).unwrap();
    assert!(a.is_infinite(), "1e2000000 NE -> inf");
    let (b, _) = BigFloat::parse_str("1e2000000", 113, RoundingMode::TowardZero).unwrap();
    assert!(!b.is_infinite(), "1e2000000 TZ -> largest finite, not inf");
    let (c, _) = BigFloat::parse_str("1e-2000000", 113, NE).unwrap();
    assert!(c.is_zero(), "1e-2000000 NE -> 0");
    let (d, _) = BigFloat::parse_str("1e-2000000", 113, RoundingMode::TowardPositive).unwrap();
    assert!(
        !d.is_zero() && d.is_sign_positive(),
        "1e-2000000 TP -> min positive"
    );
}

#[test]
fn roundtrip_sweep_near_cap() {
    // Sweep both cap boundaries with 2-digit mantissas; a finite parse must
    // round-trip through its own Display. This mirrors the libFuzzer `parse`
    // invariant, which rounds to nearest-even for BOTH the parse and the
    // re-parse — the only mode Display targets (its digit count is
    // round-trip-safe under NearestEven, not the directed modes).
    for e in [-1_000_000i64, -999_999, -999_998, 999_997, 999_998, 999_999] {
        for m in ["11", "37", "95"] {
            let s = format!("{m}e{e}");
            let (parsed, _) = BigFloat::parse_str(&s, 113, NE).unwrap();
            if parsed.is_infinite() || parsed.is_zero() || parsed.is_nan() {
                continue; // magnitude past the cap saturates; not a mismatch
            }
            let rendered = format!("{parsed}");
            let (reparsed, _) = BigFloat::parse_str(&rendered, 113, NE).unwrap();
            assert_eq!(
                parsed.partial_cmp(&reparsed).0,
                Some(Ordering::Equal),
                "{s}: {rendered} must re-parse equal"
            );
        }
    }
}
