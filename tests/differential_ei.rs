//! MPFR differential for `BigFloat::ei`.
//!
//! MPFR's `mpfr_eint` equals `Ei(x)` only for `x > 0` (it returns
//! NaN for `x < 0`), so the bit-exact MPFR lane sweeps positive
//! inputs. Negative `x` (a genuine real value for `Ei`, with no MPFR
//! oracle) is covered by a high-precision self-consistency lane:
//! `Ei` at `p` must agree with `Ei` at `p + 96` to within a few ULP,
//! which catches regime-dispatch and series bugs even though it does
//! not pin correct rounding.

#![cfg(all(unix, feature = "differential-mpfr"))]

mod differential;

use differential::{
    bigfloat_from_i64, bigfloat_to_rug, mpfr_round_of, next_i64_in, rug_from_i64, sweep_size,
    ALL_ROUNDING_MODES, TRANSCENDENTAL_PRECISIONS,
};
use pfloat::RoundingMode;

#[test]
fn ei_matches_mpfr_on_positive_inputs() {
    let mut state: u64 = u64::from_le_bytes(*b"pfloat6m");
    let cases = sweep_size().min(400);

    for &p in TRANSCENDENTAL_PRECISIONS {
        for _ in 0..cases {
            // x > 0: MPFR eint == Ei. Cap |x| so eˣ in the
            // asymptotic regime stays a sane size.
            let a = next_i64_in(&mut state, 1, 40);
            for &mode in ALL_ROUNDING_MODES {
                let bf_r = {
                    let a_bf = bigfloat_from_i64(a, p);
                    let (r, _status) = a_bf.ei(mode);
                    bigfloat_to_rug(&r)
                };
                let rug_r = {
                    let a_rg = rug_from_i64(a, p);
                    let (r, _ord) =
                        rug::Float::with_val_round(p, a_rg.eint_ref(), mpfr_round_of(mode));
                    r
                };
                assert_eq!(bf_r, rug_r, "Ei({a}) at p={p}, mode={mode:?}");
            }
        }
    }
}

#[test]
fn ei_self_consistent_on_negative_inputs() {
    let mut state: u64 = u64::from_le_bytes(*b"pfl6mneg");
    let cases = sweep_size().min(200);

    for &p in TRANSCENDENTAL_PRECISIONS {
        for _ in 0..cases {
            let a = next_i64_in(&mut state, -40, -1);
            let lo = {
                let (r, _) = bigfloat_from_i64(a, p).ei(RoundingMode::NearestEven);
                r
            };
            let hi = {
                let (r, _) = bigfloat_from_i64(a, p + 96).ei(RoundingMode::NearestEven);
                r.round_to_precision(p, RoundingMode::NearestEven)
                    .expect("p >= 1")
                    .0
            };
            // Agreement to p-8 bits: catches algorithm/regime bugs
            // without demanding correct rounding (no oracle here).
            let (diff, _) = lo.sub(&hi, RoundingMode::NearestEven);
            let abs = diff.abs();
            if abs.is_zero() {
                continue;
            }
            let two = pfloat::BigFloat::try_from_i64_exact(2, p).unwrap();
            let mut bound = hi.abs();
            for _ in 0..(p - 8) {
                bound = bound.div(&two, RoundingMode::NearestEven).0;
            }
            assert!(
                matches!(
                    abs.partial_cmp(&bound).0,
                    Some(core::cmp::Ordering::Less | core::cmp::Ordering::Equal)
                ),
                "Ei({a}) self-consistency at p={p}: lo={lo}, hi={hi}"
            );
        }
    }
}
