//! MPFR differential: `BigFloat::add` matches `rug::Float`
//! addition bit-for-bit at every tested precision and rounding mode.
//!
//! Slice 6a ships the canonical example: integer-valued operands
//! produced by `try_from_i64_exact`. Slice 6b extends to
//! finite-finite operands generated from random mantissa + exponent
//! pairs across the full normal range. ADR-0014 records the lane's
//! gating and sweep-size policy.

#![cfg(all(unix, feature = "differential-mpfr"))]

mod differential;

use differential::{
    bigfloat_from_i64, bigfloat_to_rug, mpfr_round_of, rug_from_i64, sweep_size,
    ALL_ROUNDING_MODES, SWEEP_PRECISIONS,
};

/// Splitmix64. Seeded deterministically; no `rand` dependency.
fn next_u64(state: &mut u64) -> u64 {
    *state = state.wrapping_add(0x9E37_79B9_7F4A_7C15);
    let mut z = *state;
    z = (z ^ (z >> 30)).wrapping_mul(0xBF58_476D_1CE4_E5B9);
    z = (z ^ (z >> 27)).wrapping_mul(0x94D0_49BB_1331_11EB);
    z ^ (z >> 31)
}

fn next_i64_in(state: &mut u64, lo: i64, hi: i64) -> i64 {
    debug_assert!(lo <= hi);
    let span = (hi - lo) as u64 + 1;
    lo + (next_u64(state) % span) as i64
}

#[test]
fn add_matches_mpfr_on_i64_pairs() {
    let mut state: u64 = u64::from_le_bytes(*b"pfloat6a");
    let cases = sweep_size();

    for &p in SWEEP_PRECISIONS {
        // Cap operand magnitude so each operand and the sum fit in
        // p bits exactly. The cap is 2^(p-2) so sum stays inside
        // 2^(p-1). For p >= 64 the natural i64 range covers it.
        let cap = if p >= 64 {
            i64::MAX
        } else {
            (1_i64 << (p as i64 - 2)) - 1
        };
        for _ in 0..cases {
            let a = next_i64_in(&mut state, -cap, cap);
            let b = next_i64_in(&mut state, -cap, cap);
            for &mode in ALL_ROUNDING_MODES {
                let bf_sum = {
                    let a_bf = bigfloat_from_i64(a, p);
                    let b_bf = bigfloat_from_i64(b, p);
                    let (sum, _status) = a_bf.add(&b_bf, mode);
                    bigfloat_to_rug(&sum)
                };
                let rug_sum = {
                    let a_rg = rug_from_i64(a, p);
                    let b_rg = rug_from_i64(b, p);
                    let (sum, _ord) =
                        rug::Float::with_val_round(p, &a_rg + &b_rg, mpfr_round_of(mode));
                    sum
                };
                assert_eq!(
                    bf_sum, rug_sum,
                    "add({a}, {b}) at p={p}, mode={mode:?}: pfloat={bf_sum} mpfr={rug_sum}"
                );
            }
        }
    }
}

#[test]
fn add_matches_mpfr_at_boundary_cases() {
    // A handful of hand-curated pairs that touch sign-of-zero,
    // exponent boundaries, and cancellation. These complement the
    // randomized sweep.
    // Boundary values must fit in the smallest tested precision
    // (53 bits) under try_from_i64_exact; |v| < 2^52.
    let pairs: &[(i64, i64)] = &[
        (0, 0),
        (0, 1),
        (1, -1),
        (1_000_000_000, -1_000_000_000),
        (1_i64 << 50, 1),
        (-(1_i64 << 50), -1),
        (1, 1),
        (-1, -1),
        (1_000_000, 999_999),
        (-999_999, 1_000_000),
    ];
    for &p in SWEEP_PRECISIONS {
        for &(a, b) in pairs {
            for &mode in ALL_ROUNDING_MODES {
                let a_bf = bigfloat_from_i64(a, p);
                let b_bf = bigfloat_from_i64(b, p);
                let (sum_bf, _status) = a_bf.add(&b_bf, mode);
                let bf_as_rug = bigfloat_to_rug(&sum_bf);
                let a_rg = rug_from_i64(a, p);
                let b_rg = rug_from_i64(b, p);
                let (rug_sum, _ord) =
                    rug::Float::with_val_round(p, &a_rg + &b_rg, mpfr_round_of(mode));
                assert_eq!(
                    bf_as_rug, rug_sum,
                    "boundary add({a}, {b}) at p={p}, mode={mode:?}: pfloat={bf_as_rug} mpfr={rug_sum}"
                );
            }
        }
    }
}
