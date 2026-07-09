//! Regression guard for the pf-nt21 review-tail digamma pole-sign finding
//! (`specials/type/8`), adjudicated CONFIRMED and fixed.
//!
//! `ψ(x) ~ −1/x` is a DIRECTIONAL pole, so `digamma(-0)` must be `+∞`
//! (the one-sided limit from `x → 0⁻`), not `−∞`; the pre-fix zero arm
//! discarded the sign and always returned `−∞`.
//!
//! The sibling `specials/type/11` overflow fix (unchecked
//! `target_precision + 60` → `saturating_add`) is guarded by in-module
//! unit tests on `z_min_for_target` in `digamma.rs` / `lgamma.rs`: the
//! full kernel cannot be driven at a target near `u32::MAX` (the honest
//! result would need billions of mantissa bits), so the private sizing
//! helper is tested directly.
//!
//! Run: `cargo test --release --features std,big,specials
//! --test regression_nt21_specials`.

#![cfg(all(feature = "big", feature = "specials"))]

use pfloat::{BigFloat, RoundingMode, Sign};

const NE: RoundingMode = RoundingMode::NearestEven;

#[test]
fn digamma_pole_is_sign_directional() {
    // ψ(+0) = −∞ (x → 0⁺), ψ(−0) = +∞ (x → 0⁻), both DIV_BY_ZERO.
    let (dp, sp) = BigFloat::try_new_zero(Sign::Positive, 53)
        .unwrap()
        .digamma(NE);
    assert!(
        dp.is_infinite() && dp.is_sign_negative(),
        "ψ(+0) = −∞, got {dp}"
    );
    assert!(sp.div_by_zero());
    let (dn, sn) = BigFloat::try_new_zero(Sign::Negative, 53)
        .unwrap()
        .digamma(NE);
    assert!(
        dn.is_infinite() && dn.is_sign_positive(),
        "ψ(−0) = +∞, got {dn}"
    );
    assert!(sn.div_by_zero());
}
