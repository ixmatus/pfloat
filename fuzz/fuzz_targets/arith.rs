//! Fuzz target: arithmetic-core dispatch and identity invariants.
//!
//! Feeds an [`Op`] tag plus an `(a, b, c)` i64 triple plus a
//! rounding-mode byte through each of `add`, `sub`, `mul`, `div`,
//! `sqrt`, and `fma`. The body asserts panic-freedom and the cheap
//! identity invariants that hold for finite operands:
//!
//! - `a + 0 ≡ a`
//! - `a × 1 ≡ a`
//! - `a − a` is `+0` (or `−0` under `TowardNegative`)
//! - `sqrt(a × a) ≥ |a|` is left to MPFR differential; this fuzz
//!   target tracks the cheaper structural invariants.
//!
//! Per ADR-0013 the corpus is not checked into the repo; libFuzzer
//! evolves its own. Counterexamples that survive get promoted to
//! `.proptest-regressions` against the relevant property test.

#![no_main]

use arbitrary::Arbitrary;
use libfuzzer_sys::fuzz_target;

use pfloat::{BigFloat, RoundingMode};

#[derive(Arbitrary, Debug)]
struct Input {
    op: u8,
    a: i64,
    b: i64,
    c: i64,
    mode: u8,
}

fn rounding_mode(byte: u8) -> RoundingMode {
    match byte % 5 {
        0 => RoundingMode::NearestEven,
        1 => RoundingMode::NearestAway,
        2 => RoundingMode::TowardZero,
        3 => RoundingMode::TowardPositive,
        4 => RoundingMode::TowardNegative,
        _ => unreachable!(),
    }
}

fn try_bigfloat_from_i64(n: i64, prec: u32) -> Option<BigFloat> {
    BigFloat::try_from_i64_exact(n, prec).ok()
}

fuzz_target!(|input: Input| {
    let prec: u32 = 113;
    let mode = rounding_mode(input.mode);

    // Construct operands. If a value does not fit at this precision
    // (mostly i64::MIN at p < 64), skip the iteration.
    let Some(a) = try_bigfloat_from_i64(input.a, prec) else {
        return;
    };
    let Some(b) = try_bigfloat_from_i64(input.b, prec) else {
        return;
    };
    let Some(c) = try_bigfloat_from_i64(input.c, prec) else {
        return;
    };

    // Primary dispatch — exercise each op for panic-freedom.
    match input.op % 6 {
        0 => {
            let _ = a.add(&b, mode);
        }
        1 => {
            let _ = a.sub(&b, mode);
        }
        2 => {
            let _ = a.mul(&b, mode);
        }
        3 => {
            let _ = a.div(&b, mode);
        }
        4 => {
            let _ = a.sqrt(mode);
        }
        5 => {
            let _ = a.fma(&b, &c, mode);
        }
        _ => unreachable!(),
    }

    // Identity invariants under round-to-nearest-even (the IEEE
    // default; identity holds for representable operands).
    let rm = RoundingMode::NearestEven;
    let zero = BigFloat::try_new_zero(pfloat::Sign::Positive, prec).expect("precision >= 1");
    let one = BigFloat::try_from_i64_exact(1, prec).expect("1 fits");

    // a + 0 ≡ a (for finite a; trivially holds for non-finite via
    // dispatch elsewhere — restrict to finite to keep the assertion
    // tight against `partial_cmp`).
    if a.is_finite() {
        let (sum, _) = a.add(&zero, rm);
        let (cmp, _) = sum.partial_cmp(&a);
        assert_eq!(cmp, Some(core::cmp::Ordering::Equal), "a + 0 ≡ a for a={a}");
    }

    // a × 1 ≡ a
    if a.is_finite() {
        let (prod, _) = a.mul(&one, rm);
        let (cmp, _) = prod.partial_cmp(&a);
        assert_eq!(cmp, Some(core::cmp::Ordering::Equal), "a × 1 ≡ a for a={a}");
    }

    // a − a is ±0 for finite a. The sign depends on rounding mode;
    // here we assert it is zero (sign is the special-case-dispatch
    // job, covered by Kani harnesses).
    if a.is_finite() {
        let (diff, _) = a.sub(&a, rm);
        assert!(diff.is_zero(), "a − a is zero for a={a}; got {diff}");
    }
});
