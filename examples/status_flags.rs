//! Tutorial 2: read and act on the status flags.
//!
//! Run with: `cargo run --example status_flags`
//!
//! Companion to `docs/guides/02-status-flags.md`. Every line here is the code
//! the guide walks through.

use pfloat::{flags, BigFloat, RoundingMode};

const NE: RoundingMode = RoundingMode::NearestEven;

fn main() {
    // --- An inexact result -------------------------------------------------
    // sqrt(2) is irrational, so no finite precision holds it exactly. The
    // kernel rounds and reports INEXACT for this one call.
    let two = BigFloat::try_from_i64_exact(2, 200).unwrap();
    let (_root2, s_inexact) = two.sqrt(NE);
    assert!(s_inexact.inexact());
    assert!(!s_inexact.overflow());
    assert!(!s_inexact.invalid());

    // --- An exact result ---------------------------------------------------
    // sqrt(4) is 2 exactly, representable at 200 bits. No flag is raised, so
    // is_ok() is true.
    let four = BigFloat::try_from_i64_exact(4, 200).unwrap();
    let (_root4, s_exact) = four.sqrt(NE);
    assert!(!s_exact.inexact());
    assert!(s_exact.is_ok());

    // --- An overflow -------------------------------------------------------
    // pfloat carries the exponent in an i64, so OVERFLOW is reached only by
    // pushing a value past i64::MAX. scale_by_pow2(i64::MAX) on a value whose
    // exponent is already positive saturates to the largest finite value and
    // raises OVERFLOW. The result is finite, never infinity: pfloat has no
    // emax to clamp against.
    let three = BigFloat::try_from_i64_exact(3, 53).unwrap();
    let (saturated, s_over) = three.scale_by_pow2(i64::MAX);
    assert!(s_over.overflow());
    assert!(!s_over.underflow());
    assert!(!saturated.is_infinite());

    // --- The thread-local accumulator --------------------------------------
    // Under std, every flag-producing op also OR-accumulates its Status into a
    // per-thread sticky cell, the IEEE 754 "global flags" model. Clear it,
    // run two ops, then read what accumulated across both.
    //
    // flags::clear() returns the previous value and resets the cell to OK;
    // flags::test() reads without clearing. (Names confirmed by reading
    // src/status.rs: the module exposes test, clear, set, raise.)
    let _ = flags::clear();

    // First op is inexact, second is exact. The accumulator unions them, so
    // the sticky set after both ops still carries INEXACT.
    let _ = two.sqrt(NE); // raises INEXACT
    let _ = four.sqrt(NE); // raises nothing
    let accumulated = flags::test();
    assert!(accumulated.inexact());

    // clear() hands back the accumulated set and resets the cell, so a fresh
    // exact-only computation starts clean.
    let cleared = flags::clear();
    assert!(cleared.inexact());
    let after = flags::test();
    assert!(after.is_ok());

    println!("inexact call flags  : inexact={}", s_inexact.inexact());
    println!("exact call flags    : is_ok={}", s_exact.is_ok());
    println!("overflow call flags : overflow={}", s_over.overflow());
    println!("accumulated flags   : inexact={}", accumulated.inexact());
}
