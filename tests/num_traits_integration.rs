//! Integration tests for the `num-traits` feature: generic numeric
//! code over `FixedFloat<PREC>`, plus direct trait-method checks.
//! ADR-0070.

// `FixedFloat<PREC>`'s `[(); limbs_for(PREC)]` bound needs the same
// nightly feature the library declares.
#![feature(generic_const_exprs)]
#![allow(incomplete_features)]
#![cfg(feature = "num-traits")]

use core::cmp::Ordering;

use num_traits::{FromPrimitive, Inv, Num, NumCast, One, Signed, ToPrimitive, Zero};
use pfloat::{FixedFloat, RadixParseError, RoundingMode};

type F = FixedFloat<53>;

const NE: RoundingMode = RoundingMode::NearestEven;

/// Value equality via pfloat's `partial_cmp`; expands where `PREC` is
/// concrete so no const-generic bound is needed.
macro_rules! val_eq {
    ($a:expr, $b:expr) => {
        matches!(($a).partial_cmp(&($b)).0, Some(Ordering::Equal))
    };
}

fn fi(n: i64) -> F {
    F::from_i64(n).unwrap()
}

/// Generic accumulation over `T: Zero + Add`, the canonical num-traits
/// use case.
fn sum<T>(xs: &[T]) -> T
where
    T: Zero + Copy + core::ops::Add<Output = T>,
{
    xs.iter().fold(T::zero(), |acc, &x| acc + x)
}

/// Generic polynomial over `T: Num`, exercising the `NumOps` bound.
fn poly<T: Num + Copy>(x: T) -> T {
    x * x + x + T::one()
}

#[test]
fn generic_sum_and_poly() {
    assert!(val_eq!(sum(&[fi(1), fi(2), fi(3)]), fi(6)));
    // 2*2 + 2 + 1 = 7
    assert!(val_eq!(poly(fi(2)), fi(7)));
}

#[test]
fn zero_and_one() {
    assert!(F::zero().is_zero());
    // Negative zero is still zero: the is_zero override (not `== zero()`)
    // catches it.
    assert!(F::zero().negated().is_zero());
    assert!(F::one().is_one());
    assert!(!fi(2).is_one());
    // Identities round-trip exactly.
    assert!(val_eq!(fi(7) + F::zero(), fi(7)));
    assert!(val_eq!(fi(7) * F::one(), fi(7)));
}

#[test]
fn signed() {
    assert!(fi(-5).is_negative() && !fi(-5).is_positive());
    assert!(fi(5).is_positive() && !fi(5).is_negative());
    assert!(!F::zero().is_positive() && !F::zero().is_negative());
    assert!(val_eq!(Signed::abs(&fi(-5)), fi(5)));
    assert!(val_eq!(fi(-5).signum(), fi(-1)));
    // abs_sub: max(a - b, 0).
    assert!(val_eq!(fi(7).abs_sub(&fi(3)), fi(4)));
    assert!(val_eq!(fi(3).abs_sub(&fi(7)), F::zero()));
}

#[test]
fn primitives_round_trip() {
    assert_eq!(fi(42).to_i64(), Some(42));
    assert_eq!(fi(42).to_f64(), Some(42.0));
    assert_eq!(F::from_f64(0.5).unwrap().to_f64(), Some(0.5));
    assert_eq!(F::from_f32(0.25).unwrap().to_f64(), Some(0.25));
    // u64 above i64::MAX uses the split path and stays positive.
    let big = F::from_u64(u64::MAX).expect("from_u64");
    assert!(big.is_positive());
    // to_i64 out of range returns None.
    assert_eq!(F::from_f64(1e40).unwrap().to_i64(), None);
}

#[test]
fn num_cast_and_inv() {
    let c: F = NumCast::from(3.5_f64).expect("num_cast");
    assert_eq!(c.to_f64(), Some(3.5));
    assert_eq!(fi(2).inv().to_f64(), Some(0.5));
}

#[test]
fn from_str_radix_decimal_and_radix_error() {
    let v = F::from_str_radix("3.5", 10).expect("decimal parse");
    assert_eq!(v.to_f64(), Some(3.5));
    assert_eq!(
        F::from_str_radix("11", 2),
        Err(RadixParseError::UnsupportedRadix(2))
    );
}

#[test]
fn from_u64_rounds_once() {
    // Regression (adversarial review): the `n >= 2^63` path must round
    // the full magnitude once. The prior split rounded the low part to
    // PREC and then added 2^63, a double rounding that lands 1 ULP off
    // (16.7% of top-bit-set u64 at PREC=53). The independent reference
    // is the exact decimal parsed once, which the parser rounds correctly.
    let cases = [
        0xfe80_0d65_69fa_1bff_u64, // the reported 1-ULP trigger at PREC=53
        u64::MAX,
        0x8000_0000_0000_0001,
        0xffff_ffff_ffff_fc00,
        0xc000_0000_0000_07ff,
        0x9249_2492_4924_9249,
    ];
    for &n in &cases {
        let got = F::from_u64(n).expect("from_u64");
        let want = F::parse_str(&n.to_string(), NE).expect("parse").0;
        assert!(
            val_eq!(got, want),
            "from_u64({n}): got {:?}, want {:?}",
            got.to_f64(),
            want.to_f64()
        );
    }
    // Deterministic xorshift sweep over top-bit-set values, where the
    // double rounding used to bite.
    let mut s: u64 = 0x1234_5678_9abc_def0;
    for _ in 0..20_000 {
        s ^= s << 13;
        s ^= s >> 7;
        s ^= s << 17;
        let n = s | (1 << 63);
        let got = F::from_u64(n).expect("from_u64");
        let want = F::parse_str(&n.to_string(), NE).expect("parse").0;
        assert!(val_eq!(got, want), "from_u64({n}) mismatch");
    }
    // PREC=1 made the double rounding a 2x error (3*2^62 rounds to 2^64
    // instead of 2^63). One round gives 2^63.
    type H = FixedFloat<1>;
    let n = 0xbfff_ffff_ffff_ffff_u64;
    let got = H::from_u64(n).expect("from_u64");
    let want = H::parse_str(&n.to_string(), NE).expect("parse").0;
    assert!(
        matches!(got.partial_cmp(&want).0, Some(Ordering::Equal)),
        "PREC=1 from_u64 double rounding: got {:?}, want {:?}",
        got.to_f64(),
        want.to_f64()
    );
}

#[test]
fn abs_sub_propagates_nan() {
    // Regression (adversarial review): a NaN operand must propagate, not
    // collapse to +0. This matches the num-traits f32/f64 reference and
    // pfloat's own NaN discipline.
    let nan = F::from_f64(f64::NAN).expect("nan");
    assert!(nan.abs_sub(&fi(3)).is_nan());
    assert!(fi(3).abs_sub(&nan).is_nan());
    // The ordered branches are unchanged.
    assert!(val_eq!(fi(7).abs_sub(&fi(3)), fi(4)));
    assert!(fi(3).abs_sub(&fi(7)).is_zero());
}

#[test]
fn works_at_another_precision() {
    type G = FixedFloat<128>;
    assert!(G::one().is_one());
    let s = sum(&[G::from_i64(1).unwrap(), G::from_i64(2).unwrap()]);
    assert!(val_eq!(s, G::from_i64(3).unwrap()));
}
