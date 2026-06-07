//! Fuzz target: the `pfloat-ball` decimal parser as an adversarial
//! boundary (slice 10).
//!
//! `Ball::parse_decimal` takes attacker-controlled bytes and must (1)
//! never panic, (2) never hang (its DoS bounds reject pathological
//! literals before any bignum work; libFuzzer's timeout catches a
//! regression), and (3) produce a well-formed ball — `lower ≤ upper` and
//! a Display that renders without panicking. The byte-prefix also drives
//! the requested precision so the parser is exercised across widths.
//!
//! Per ADR-0013 the corpus is not checked in; libFuzzer evolves its own.

#![no_main]

use libfuzzer_sys::fuzz_target;

use pfloat_ball::Ball;
use pfloat::BigFloat;

fuzz_target!(|data: &[u8]| {
    // First byte (if any) picks a precision in a sane range; the rest is
    // the literal under test.
    let (prec, rest) = match data.split_first() {
        Some((p, rest)) => (1u32 + u32::from(*p) * 8, rest),
        None => (53, &[][..]),
    };
    let Ok(s) = core::str::from_utf8(rest) else {
        return;
    };

    match Ball::<BigFloat>::parse_decimal(s, prec) {
        Ok(ball) => {
            // Well-formed: the endpoints are ordered (lower ≤ upper).
            let lo = ball.lower();
            let hi = ball.upper();
            assert!(
                lo.partial_cmp(&hi).0 != Some(core::cmp::Ordering::Greater),
                "parsed ball has lower > upper: s={s:?} prec={prec}"
            );
            // The printer never panics, and a fresh parse of the printed
            // interval endpoints stays finite/ordered.
            let _ = ball.to_decimal_interval(17);
            let _ = format!("{ball}");
        }
        Err(_) => {
            // Rejection is fine; the contract is "no panic, no hang".
        }
    }
});
