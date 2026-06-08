//! How-to: run pfloat in a `no_std` or embedded build with `FixedFloat<PREC>`.
//!
//! Run with: `cargo run --example fixed_precision --features fixed`
//!
//! Companion to `docs/guides/08-no-std-embedded.md`. This example itself is a
//! `std` binary (it prints), but every numeric value it builds lives in
//! stack-allocated `FixedFloat<113>` storage: binary128 width, no heap. The
//! same `FixedFloat<113>` type compiles in a `default-features = false` build
//! with only the `fixed` feature, the configuration an embedded target uses.

// FixedFloat<PREC> carries the const-generic bound `[(); limbs_for(PREC)]:`,
// so any crate that names the type (including this example) declares the same
// nightly feature the library declares. The pinned toolchain provides it.
#![feature(generic_const_exprs)]
#![allow(incomplete_features)]

use pfloat::{FixedFloat, RoundingMode};

// binary128 width: the IEEE 754 quad mantissa is 113 bits. The precision is a
// type parameter, fixed at compile time, so the storage size is known to the
// compiler and the value lives on the stack.
type Quad = FixedFloat<113>;

const NE: RoundingMode = RoundingMode::NearestEven;

fn main() {
    // The precision rides on the type, so the constructor takes no precision
    // argument. Compare BigFloat::try_from_i64_exact(2, 200), where 200 is a
    // runtime bit count, against Quad::try_from_i64_exact(2), where 113 is
    // baked into the type.
    let two = Quad::try_from_i64_exact(2).unwrap();
    assert_eq!(Quad::PRECISION, 113);

    // sqrt takes a rounding mode and returns the value paired with a Status,
    // exactly like BigFloat. sqrt(2) is irrational, so it is inexact at any
    // finite width.
    let (root, status) = two.sqrt(NE);
    assert!(status.inexact());

    // A perfect square is exact at this width: sqrt(4) = 2, no rounding, so the
    // INEXACT flag stays clear. This is the per-call exactness witness, the
    // same contract the BigFloat surface offers.
    let four = Quad::try_from_i64_exact(4).unwrap();
    let (two_again, exact_status) = four.sqrt(NE);
    assert!(!exact_status.inexact());

    // To print, or to reach the decimal formatter and the elementary
    // functions, widen to BigFloat with to_big(). FixedFloat carries
    // arithmetic and roots without alloc; the decimal string formatter lives
    // on BigFloat.
    let root_big = root.to_big();
    assert_eq!(root_big.precision(), 113);
    let digits = root_big.to_decimal_string(30, NE);
    println!("sqrt(2) at binary128 width = {digits}");
    // The leading digits of sqrt(2) are known; pin a prefix so the run checks
    // itself. 113 bits pins about 34 decimal digits, so all 30 here are
    // correct. (sqrt(2) = 1.41421356237309504880168872420969807...)
    assert!(digits.starts_with("1.4142135623730950488016887242"));

    // Round-trip the exact case through BigFloat and confirm it reads back as
    // two. From<FixedFloat<PREC>> for BigFloat is the same widening as
    // to_big().
    let two_back = two_again.to_big();
    let two_target = Quad::try_from_i64_exact(2).unwrap().to_big();
    let (ord, cmp_status) = two_back.partial_cmp(&two_target);
    assert!(cmp_status.is_ok());
    assert_eq!(ord, Some(core::cmp::Ordering::Equal));

    println!(
        "sqrt(4) at binary128 width = {}",
        two_back.to_decimal_string(2, NE)
    );
}
