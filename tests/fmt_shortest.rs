//! Verification for the Dragon4 shortest-output formatter (ADR-0071).
//!
//! Oracle: for an f64 value, pfloat's shortest output must (a) parse
//! back to exactly that f64 and (b) have the same number of significant
//! digits as Rust's own shortest f64 formatting (which is minimal).
//! Exact dyadic-halfway ties (two equal-length outputs equidistant from
//! the value) are accepted whichever neighbor each formatter picks:
//! pfloat rounds half to even on the decimal digit, Rust does not, and
//! both are valid shortest round-trips. A wrong digit on a non-tie value
//! fails the parse-back check, so real bugs are still caught.

#![cfg(all(feature = "big", feature = "std"))]

use pfloat::{BigFloat, RoundingMode};

const NE: RoundingMode = RoundingMode::NearestEven;

/// Significant digit count: digits with sign, exponent, point, and
/// leading/trailing zeros removed (zero counts as one).
fn sig_digit_count(s: &str) -> usize {
    sig_digits(s).len()
}

/// Significant digit sequence of a value.
fn sig_digits(s: &str) -> String {
    let s = s.trim_start_matches('-');
    let mantissa = s.split(['e', 'E']).next().unwrap_or(s);
    let digits: String = mantissa.chars().filter(char::is_ascii_digit).collect();
    let trimmed = digits.trim_start_matches('0').trim_end_matches('0');
    if trimmed.is_empty() {
        "0".to_string()
    } else {
        trimmed.to_string()
    }
}

/// pfloat's shortest output `s` for f64 `x` is correct and minimal.
fn assert_shortest_ok(s: &str, x: f64) {
    assert_eq!(
        s.parse::<f64>().expect("reparse"),
        x,
        "round-trip {s} -> {x}"
    );
    let rust = format!("{x}");
    assert_eq!(
        sig_digit_count(s),
        sig_digit_count(&rust),
        "minimality x={x}: pfloat={s} rust={rust}"
    );
}

#[test]
fn matches_rust_f64_shortest() {
    let cases = [
        0.5f64,
        1.0,
        -1.0,
        3.5,
        0.1,
        0.2,
        0.3,
        1.0 / 3.0,
        2.0 / 3.0,
        123.456,
        1e20,
        1e-20,
        2.5,
        1_234_567_890.0,
        9.999_999,
        100.0,
        0.001,
        core::f64::consts::PI,
        core::f64::consts::E,
        1.5e300,
        -7.25e-100,
        f64::MIN_POSITIVE, // smallest normal f64
        f64::MAX,
        123_456.789,
        // Exact dyadic halfway tie (-201562347225087.625); via bits to
        // avoid a clippy excessive-precision rewrite that would change it.
        f64::from_bits(0xc2e6_ea3c_8367_fff4),
    ];
    // f64 subnormals are excluded: they carry fewer than 53 bits of
    // precision, but pfloat has no subnormals, so a precision-53 value
    // always carries 53 bits — its shortest round-trip is correct but
    // differs from f64's subnormal shortest (different precision models).
    for &x in &cases {
        let bf = BigFloat::from_f64(x).round_to_precision(53, NE).unwrap().0;
        let s = bf.to_shortest_decimal_string();
        // Round-trips through pfloat's own parser to the exact value.
        let back = BigFloat::parse_str(&s, 53, NE).expect("parse").0;
        assert_eq!(back, bf, "pfloat round-trip x={x}: {s}");
        assert_shortest_ok(&s, x);
    }
}

#[test]
fn random_f64_sweep() {
    // Sweep random f64 bit patterns; each finite normal value must round
    // back to itself and be minimal.
    let mut state = 0x1234_5678_9abc_def0u64;
    let mut checked = 0usize;
    for _ in 0..200_000 {
        state ^= state << 13;
        state ^= state >> 7;
        state ^= state << 17;
        let x = f64::from_bits(state);
        if !x.is_finite() || x == 0.0 || x.is_subnormal() {
            continue;
        }
        let bf = BigFloat::from_f64(x).round_to_precision(53, NE).unwrap().0;
        let s = bf.to_shortest_decimal_string();
        assert_shortest_ok(&s, x);
        let back = BigFloat::parse_str(&s, 53, NE).expect("parse").0;
        assert_eq!(back, bf, "pfloat round-trip x={x}: {s}");
        checked += 1;
    }
    assert!(checked > 1_000, "only {checked} normals checked");
}

#[test]
fn round_trips_at_arbitrary_precision() {
    for &p in &[24u32, 53, 64, 113, 200, 256] {
        let one = BigFloat::try_from_i64_exact(1, p).unwrap();
        let three = BigFloat::try_from_i64_exact(3, p).unwrap();
        let seven = BigFloat::try_from_i64_exact(7, p).unwrap();
        for v in [one.div(&three, NE).0, one.div(&seven, NE).0] {
            let s = v.to_shortest_decimal_string();
            let back = BigFloat::parse_str(&s, p, NE).expect("parse").0;
            assert_eq!(back, v, "round-trip at p={p}: {s}");
            // Shortest must not exceed the round-trip-safe digit count.
            assert!(
                sig_digit_count(&s) <= BigFloat::round_trip_digit_count(p) as usize,
                "p={p}: {s} has {} sig digits, cap {}",
                sig_digit_count(&s),
                BigFloat::round_trip_digit_count(p)
            );
        }
    }
}

#[test]
fn exact_integers_and_powers_of_two() {
    for &n in &[0i64, 1, -1, 2, 8, 10, 100, 1024, -4096, 1_000_000] {
        for &p in &[53u32, 113] {
            let v = BigFloat::try_from_i64_exact(n, p).unwrap();
            let s = v.to_shortest_decimal_string();
            let back = BigFloat::parse_str(&s, p, NE).expect("parse").0;
            assert_eq!(back, v, "integer {n} at p={p}: {s}");
            if n != 0 {
                assert_eq!(
                    sig_digits(&s),
                    sig_digits(&format!("{n}")),
                    "integer {n}: {s}"
                );
            }
        }
    }
}

#[test]
fn special_values() {
    let p = 53;
    let pos_zero = BigFloat::try_from_i64_exact(0, p).unwrap();
    assert_eq!(pos_zero.to_shortest_decimal_string(), "0");
    assert_eq!(pos_zero.negated().to_shortest_decimal_string(), "-0");
}
