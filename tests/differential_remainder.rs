//! MPFR differential: `BigFloat::remainder` matches `rug::Float`'s
//! IEEE 754-2019 remainder (`mpfr_remainder`) bit-for-bit.
//!
//! `remainder` is exact (it never rounds), so there is no rounding-mode
//! sweep; the result must equal MPFR's at a precision of
//! `max(px, py)`. The sweep covers random integer pairs, power-of-two
//! exponent gaps (exercising the modular-exponentiation reduction and
//! the `2|x| < |y|` early exit), explicit round-to-nearest-even ties,
//! and fractional operands. ADR-0069, slice C3c prerequisite (pf-2138).

#![cfg(all(unix, feature = "differential-mpfr"))]

mod differential;

use differential::{bigfloat_from_i64, bigfloat_pow2, bigfloat_to_rug, next_i64_in, sweep_size};
use pfloat::BigFloat;
use rug::Float;

/// Compare pfloat's remainder against MPFR for one operand pair, both
/// converted exactly to `rug::Float`.
fn check(x: &BigFloat, y: &BigFloat, label: &str) {
    let target = x.precision().max(y.precision());
    let (pf, _status) = x.remainder(y);
    let pf_rug = bigfloat_to_rug(&pf);

    let x_r = bigfloat_to_rug(x);
    let y_r = bigfloat_to_rug(y);
    let oracle = Float::with_val(target, x_r.remainder_ref(&y_r));

    assert_eq!(pf_rug, oracle, "remainder mismatch: {label}");
}

#[test]
fn remainder_matches_mpfr_on_integers() {
    let mut state: u64 = u64::from_le_bytes(*b"pf-rem01");
    let cases = sweep_size().min(2_000);
    for &p in &[53u32, 113, 256] {
        let cap = if p >= 62 {
            i64::MAX
        } else {
            (1_i64 << (p - 1)) - 1
        };
        for _ in 0..cases {
            let a = next_i64_in(&mut state, -cap, cap);
            let mut b = next_i64_in(&mut state, -cap, cap);
            if b == 0 {
                b = 1;
            }
            check(
                &bigfloat_from_i64(a, p),
                &bigfloat_from_i64(b, p),
                &format!("rem({a}, {b}) p={p}"),
            );
        }
    }
}

#[test]
fn remainder_matches_mpfr_on_pow2_gaps() {
    // x = ±2^e dominates y: exercises the modexp reduction (ax > ay)
    // for a range of exponent gaps. The reduction is O(log k), so these
    // moderate exponents drive the same code path a 10^9 gap would.
    let p = 128u32;
    for &e in &[60i64, 64, 100, 200, 500, 1000, 3000] {
        for &b in &[3i64, -3, 5, 7, -11, 1_000_003] {
            let x = bigfloat_pow2(e, p);
            let y = bigfloat_from_i64(b, p);
            check(&x, &y, &format!("rem(2^{e}, {b})"));
            // Negative dividend: remainder is odd in x.
            let neg_x = bigfloat_from_i64(-1, p)
                .mul(&x, pfloat::RoundingMode::NearestEven)
                .0;
            check(&neg_x, &y, &format!("rem(-2^{e}, {b})"));
        }
    }
}

#[test]
fn remainder_matches_mpfr_on_y_dominant() {
    // |y| >> |x|: exercises the early exit (2|x| < |y| → x) and the
    // y-dominant divmod branch near the boundary.
    let p = 128u32;
    for &e in &[60i64, 100, 500, 2000] {
        for &a in &[1i64, 3, -7, 12345] {
            let x = bigfloat_from_i64(a, p);
            let y = bigfloat_pow2(e, p);
            check(&x, &y, &format!("rem({a}, 2^{e})"));
        }
    }
}

#[test]
fn remainder_matches_mpfr_on_ties() {
    // Half-integer quotients drive the round-to-nearest-even tie path.
    // y a power of two, x an odd multiple of y/2 → x/y = k + 0.5.
    let p = 64u32;
    for &y in &[2i64, 4, 8, 16] {
        for k in -6i64..=6 {
            // x = (2k + 1) * y / 2  → x / y = k + 0.5
            let x = (2 * k + 1) * y / 2;
            check(
                &bigfloat_from_i64(x, p),
                &bigfloat_from_i64(y, p),
                &format!("tie rem({x}, {y})"),
            );
        }
    }
}

#[test]
fn remainder_matches_mpfr_on_fractionals() {
    // Non-integer operands at mixed precisions via exact construction:
    // x = na / 2^sa, y = nb / 2^sb.
    let mut state: u64 = u64::from_le_bytes(*b"pf-rem02");
    for &(px, py) in &[(53u32, 53u32), (53, 113), (200, 64)] {
        for _ in 0..500 {
            let na = next_i64_in(&mut state, -(1 << 40), 1 << 40);
            let mut nb = next_i64_in(&mut state, -(1 << 40), 1 << 40);
            if nb == 0 {
                nb = 1;
            }
            let sa = next_i64_in(&mut state, 0, 30);
            let sb = next_i64_in(&mut state, 0, 30);
            // x = na * 2^-sa at px, y = nb * 2^-sb at py.
            let x = bigfloat_from_i64(na, px)
                .mul(&bigfloat_pow2(-sa, px), pfloat::RoundingMode::NearestEven)
                .0;
            let y = bigfloat_from_i64(nb, py)
                .mul(&bigfloat_pow2(-sb, py), pfloat::RoundingMode::NearestEven)
                .0;
            check(&x, &y, &format!("frac rem({na}/2^{sa}, {nb}/2^{sb})"));
        }
    }
}
