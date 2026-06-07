//! Kani proof harnesses over the `Vec`-free [`Mag`](crate::Mag)
//! invariants.
//!
//! Compiled only under `cfg(kani)`. `Mag` is `Copy`, single-`u64`-limb,
//! and `Vec`-free, so its round-up invariants discharge cleanly under
//! CBMC — exactly the property the 64-bit single-limb choice was made
//! for (ADR-0074). Anything `BigFloat`-backed (`Ball`, the conversions)
//! is unverifiable at the heap-allocation level (CBMC is `Vec`-hostile,
//! ADR-0062), so those soundness claims rest on the property and
//! differential lanes instead. ADR-0078.
//!
//! Run locally with `cargo kani --features=kani` (or via the manual
//! Kani CI lane). The lane is advisory.

use crate::mag::Mag;

/// A nondeterministic finite `Mag` with the canonical top-bit-set
/// mantissa and a bounded exponent (bounding keeps CBMC tractable; the
/// invariants are exponent-translation-invariant).
fn any_finite() -> Mag {
    let mantissa: u64 = kani::any();
    kani::assume(mantissa & (1u64 << 63) != 0);
    let exponent: i64 = kani::any();
    kani::assume(exponent.abs() < 1_000);
    Mag::Finite { exponent, mantissa }
}

fn is_canonical(m: Mag) -> bool {
    match m {
        Mag::Zero | Mag::Infinity => true,
        Mag::Finite { mantissa, .. } => mantissa & (1u64 << 63) != 0,
    }
}

#[kani::proof]
fn mag_add_is_monotone_and_canonical() {
    let a = any_finite();
    let b = any_finite();
    let s = a.add(b);
    // Round-up: the sum is at least each operand (Ord is value order).
    assert!(s >= a);
    assert!(s >= b);
    assert!(is_canonical(s));
}

#[kani::proof]
fn mag_add_is_commutative() {
    let a = any_finite();
    let b = any_finite();
    assert!(a.add(b) == b.add(a));
}

#[kani::proof]
fn mag_zero_is_additive_identity() {
    let a = any_finite();
    assert!(Mag::ZERO.add(a) == a);
    assert!(a.add(Mag::ZERO) == a);
    // Exact-in-exact-out at the Mag level: 0 + 0 = 0.
    assert!(Mag::ZERO.add(Mag::ZERO) == Mag::ZERO);
}

#[kani::proof]
fn mag_zero_annihilates_under_mul() {
    let a = any_finite();
    assert!(Mag::ZERO.mul(a) == Mag::ZERO);
    assert!(a.mul(Mag::ZERO) == Mag::ZERO);
}

#[kani::proof]
fn mag_mul_is_canonical_and_commutative() {
    let a = any_finite();
    let b = any_finite();
    let p = a.mul(b);
    assert!(is_canonical(p));
    assert!(p == b.mul(a));
}

#[kani::proof]
fn mag_infinity_absorbs() {
    let a = any_finite();
    assert!(Mag::INFINITY.add(a) == Mag::INFINITY);
    assert!(a.add(Mag::INFINITY) == Mag::INFINITY);
    assert!(Mag::INFINITY.mul(a) == Mag::INFINITY);
    assert!(a.mul(Mag::INFINITY) == Mag::INFINITY);
    // 0 · ∞ saturates to ∞ (never a false exact zero).
    assert!(Mag::ZERO.mul(Mag::INFINITY) == Mag::INFINITY);
}

#[kani::proof]
fn mag_from_pow2_is_canonical_power() {
    let k: i64 = kani::any();
    kani::assume(k.abs() < 1_000);
    let m = Mag::from_pow2(k);
    assert!(is_canonical(m));
    assert!(
        m == Mag::Finite {
            exponent: k,
            mantissa: 1u64 << 63
        }
    );
}

#[kani::proof]
fn mag_ordering_is_total_across_variants() {
    let a = any_finite();
    // Zero < every finite < Infinity.
    assert!(Mag::ZERO < a);
    assert!(a < Mag::INFINITY);
    assert!(Mag::ZERO < Mag::INFINITY);
}
