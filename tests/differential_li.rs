//! MPFR differential for `BigFloat::li` via the identity
//! `li(x) = Ei(ln x)`.
//!
//! MPFR has no `li`; the oracle composes `mpfr_eint(mpfr_log(x))`,
//! which is valid only for `x > 1` (there `ln x > 0`, the domain
//! where `eint == Ei`). That lane compares to within a few ULP
//! because both sides double-round the composition. The `0 < x < 1`
//! range (no oracle) is covered by a high-precision self-consistency
//! lane.

#![cfg(all(unix, feature = "differential-mpfr"))]

mod differential;

use differential::{bigfloat_to_rug, next_i64_in, sweep_size, TRANSCENDENTAL_PRECISIONS};
use pfloat::{BigFloat, RoundingMode};

fn close_within(a: &BigFloat, b: &BigFloat, bits: u32) -> bool {
    let (diff, _) = a.sub(b, RoundingMode::NearestEven);
    let abs = diff.abs();
    if abs.is_zero() {
        return true;
    }
    let p = a.precision().max(b.precision());
    let two = BigFloat::try_from_i64_exact(2, p).unwrap();
    let mut bound = b.abs();
    if bound.is_zero() {
        bound = BigFloat::try_from_i64_exact(1, p).unwrap();
    }
    for _ in 0..bits {
        bound = bound.div(&two, RoundingMode::NearestEven).0;
    }
    matches!(
        abs.partial_cmp(&bound).0,
        Some(core::cmp::Ordering::Less | core::cmp::Ordering::Equal)
    )
}

#[test]
fn li_matches_mpfr_eint_of_ln_for_x_gt_1() {
    let mut state: u64 = u64::from_le_bytes(*b"pfl6mliA");
    let cases = sweep_size().min(300);

    for &p in TRANSCENDENTAL_PRECISIONS {
        for _ in 0..cases {
            let a = next_i64_in(&mut state, 2, 50);
            let bf = {
                let (r, _) = BigFloat::try_from_i64_exact(a, p)
                    .unwrap()
                    .li(RoundingMode::NearestEven);
                r
            };
            // Oracle: eint(ln x) at a generous working precision.
            // li is a composition, so both sides double-round;
            // compare in rug to within ≤ ~2 ULP at precision p.
            let work = p + 256;
            let t = rug::Float::with_val(work, a).ln();
            let want = rug::Float::with_val(work, t.eint_ref());
            let got = bigfloat_to_rug(&bf);
            let err = rug::Float::with_val(work, &got - &want).abs();
            let mut tol = want.clone().abs();
            if tol == 0 {
                tol = rug::Float::with_val(work, 1);
            }
            tol >>= p - 2; // ≤ 2 ULP at precision p
            assert!(
                err <= tol,
                "li({a}) at p={p}: got {bf}, |err|={err:e}, tol={tol:e}"
            );
        }
    }
}

#[test]
fn li_self_consistent_on_unit_interval() {
    let mut state: u64 = u64::from_le_bytes(*b"pfl6mliB");
    let cases = sweep_size().min(150);

    for &p in TRANSCENDENTAL_PRECISIONS {
        for _ in 0..cases {
            // x = 1/d ∈ (0, 1): no MPFR oracle (ln x < 0).
            let d = next_i64_in(&mut state, 2, 50);
            let x_lo = BigFloat::try_from_i64_exact(1, p)
                .unwrap()
                .div(
                    &BigFloat::try_from_i64_exact(d, p).unwrap(),
                    RoundingMode::NearestEven,
                )
                .0;
            let x_hi = BigFloat::try_from_i64_exact(1, p + 96)
                .unwrap()
                .div(
                    &BigFloat::try_from_i64_exact(d, p + 96).unwrap(),
                    RoundingMode::NearestEven,
                )
                .0;
            let (lo, _) = x_lo.li(RoundingMode::NearestEven);
            let hi = {
                let (r, _) = x_hi.li(RoundingMode::NearestEven);
                r.round_to_precision(p, RoundingMode::NearestEven)
                    .expect("p >= 1")
                    .0
            };
            assert!(
                close_within(&lo, &hi, p - 8),
                "li(1/{d}) self-consistency at p={p}: lo={lo}, hi={hi}"
            );
        }
    }
}
