//! Adjacent-representable values and unit-in-last-place for
//! [`BigFloat`] (IEEE 754-2019 §5.3.1 `nextUp`/`nextDown`, plus
//! `ulp`). ADR-0073.
//!
//! # pfloat's representable set, and why the boundaries differ from
//! IEEE binary formats
//!
//! At precision `p` pfloat represents `±0`, `±∞`, NaN, and finite
//! normals `±m·2^(e−p+1)` where `m` is a `p`-bit integer with the top
//! bit set (`m ∈ [2^(p−1), 2^p−1]`) and `e` is an `i64`. There are **no
//! subnormals and no `emin`/`emax`**: the exponent `e` saturates within
//! the `i64` range instead of an IEEE format's bounded one. Two derived
//! values name the edges of the finite range at precision `p`:
//!
//! - **`MinPos`** — the smallest positive value: `m = 2^(p−1)`,
//!   `e = i64::MIN` (the top-bit-only mantissa at the lowest exponent).
//!   Its value is `2^(i64::MIN)`.
//! - **`MaxFinite`** — the largest finite magnitude: `m = 2^p − 1`
//!   (all-ones mantissa), `e = i64::MAX`.
//!
//! The exponent floor and ceiling are where `nextUp`/`nextDown` meet
//! `±0` and `±∞`:
//!
//! - `next_up(MaxFinite)` is `+∞` (the exponent cannot increment past
//!   `i64::MAX`, so the successor is infinity, exactly as IEEE
//!   `nextUp(maxFinite) = +∞`).
//! - `next_up(−MinPos)` is `−0` (nothing lies between the
//!   smallest-magnitude negative and zero; IEEE 754-2019 §5.3.1 fixes
//!   the sign of this zero as negative).
//! - `next_up(±0)` is `+MinPos` (the gap from zero to the nearest
//!   representable; pfloat has no subnormals, so this gap is the full
//!   `2^(i64::MIN)`, not a subnormal step).
//!
//! `next_down(x)` is defined by the IEEE identity `−next_up(−x)`, which
//! makes the negative-side boundaries fall out of the positive ones.
//!
//! # `ulp`
//!
//! `ulp(x)` is the positive weight of the last mantissa bit,
//! `2^(e−p+1)` — the distance from a finite non-zero `x` to its
//! larger-magnitude neighbour. `ulp(±0) = MinPos`, `ulp(±∞) = +∞`,
//! `ulp(NaN) = NaN`. Because pfloat has no `emin`, the true `ulp` of a
//! value near the exponent floor (`e − p + 1 < i64::MIN`) is itself
//! below `MinPos` and not representable; `ulp` then **saturates upward**
//! to `MinPos` and raises [`Status::UNDERFLOW`]. An upward-saturated
//! `ulp` is an over-estimate, which is the sound direction for the
//! radius and `from_interval` consumers in `pfloat-ball`.
//!
//! # Signaling NaN
//!
//! `nextUp`/`nextDown` are general-computational operations and signal
//! [`Status::INVALID`] on a signaling NaN, returning a quiet NaN
//! (§5.3.1). `ulp` follows the same rule. A quiet NaN passes through
//! unchanged with [`Status::OK`]. The result NaN carries the input's
//! sign; IEEE 754-2019 leaves a NaN result's sign unspecified, and
//! preserving it keeps the three operations mutually consistent under
//! the `next_down(x) = −next_up(−x)` identity (a caller must not rely
//! on a NaN's sign regardless).

use alloc::vec;

use crate::big::BigFloat;
use crate::class::Class;
use crate::mantissa::limbs_for;
use crate::sign::Sign;
use crate::status::{auto_raise, Status};

impl BigFloat {
    /// IEEE 754-2019 §5.3.1 `nextUp(self)`: the least representable
    /// value (at `self`'s precision) that compares greater than
    /// `self`.
    ///
    /// `next_up(+∞) = +∞`, `next_up(−∞) = −MaxFinite`,
    /// `next_up(±0) = +MinPos`, and `next_up(−MinPos) = −0` (see the
    /// module docs for `MinPos`/`MaxFinite` and the saturating-exponent
    /// boundaries). A quiet NaN passes through; a signaling NaN raises
    /// [`Status::INVALID`] and is quieted.
    #[must_use]
    pub fn next_up(&self) -> (Self, Status) {
        next_up_kernel(self)
    }

    /// [`next_up`](Self::next_up) accumulating into a caller-supplied
    /// flag bag.
    #[must_use]
    pub fn next_up_with_flags(&self, flags: &mut Status) -> Self {
        let (value, status) = self.next_up();
        *flags |= status;
        value
    }

    /// IEEE 754-2019 §5.3.1 `nextDown(self)`: the greatest
    /// representable value (at `self`'s precision) that compares less
    /// than `self`. Defined by the identity `nextDown(x) = −nextUp(−x)`.
    ///
    /// `next_down(−∞) = −∞`, `next_down(+∞) = +MaxFinite`,
    /// `next_down(±0) = −MinPos`, and `next_down(+MinPos) = +0`. A quiet
    /// NaN passes through; a signaling NaN raises [`Status::INVALID`]
    /// and is quieted.
    #[must_use]
    pub fn next_down(&self) -> (Self, Status) {
        // nextDown(x) = -nextUp(-x). Negation is exact and preserves
        // NaN payload/signaling, so the single INVALID (if any) is
        // raised exactly once inside next_up.
        let (up, status) = self.negated().next_up();
        (up.negated(), status)
    }

    /// [`next_down`](Self::next_down) accumulating into a
    /// caller-supplied flag bag.
    #[must_use]
    pub fn next_down_with_flags(&self, flags: &mut Status) -> Self {
        let (value, status) = self.next_down();
        *flags |= status;
        value
    }

    /// Unit in the last place: the positive distance from `self` to its
    /// larger-magnitude representable neighbour, `2^(e−p+1)` for a
    /// finite non-zero value.
    ///
    /// `ulp(±0) = MinPos`, `ulp(±∞) = +∞`, `ulp(NaN) = NaN`. The result
    /// always has the same precision as `self` and a positive sign. For
    /// a value near the exponent floor the true `ulp` can fall below
    /// `MinPos`; `ulp` then saturates upward to `MinPos` and raises
    /// [`Status::UNDERFLOW`] (a sound over-estimate). A signaling NaN
    /// raises [`Status::INVALID`].
    #[must_use]
    pub fn ulp(&self) -> (Self, Status) {
        ulp_kernel(self)
    }

    /// [`ulp`](Self::ulp) accumulating into a caller-supplied flag bag.
    #[must_use]
    pub fn ulp_with_flags(&self, flags: &mut Status) -> Self {
        let (value, status) = self.ulp();
        *flags |= status;
        value
    }
}

fn next_up_kernel(x: &BigFloat) -> (BigFloat, Status) {
    let p = x.precision;
    match &x.class {
        Class::Nan {
            quiet,
            sign,
            payload,
        } => {
            if *quiet {
                let nan = BigFloat::try_new_quiet_nan(*sign, p, payload)
                    .expect("BigFloat invariant: precision >= 1");
                (nan, Status::OK)
            } else {
                // Preserve the input sign when quieting (IEEE leaves a
                // NaN result's sign unspecified; preserving it keeps
                // next_up/next_down/ulp mutually consistent under the
                // next_down = -next_up(-x) identity).
                let nan = BigFloat::try_new_quiet_nan(*sign, p, &[])
                    .expect("BigFloat invariant: precision >= 1");
                auto_raise(Status::INVALID);
                (nan, Status::INVALID)
            }
        }
        Class::Infinity { sign } => match sign {
            // next_up(+inf) = +inf (already the maximum).
            Sign::Positive => (
                BigFloat::try_new_infinity(Sign::Positive, p).unwrap(),
                Status::OK,
            ),
            // next_up(-inf) = the most-negative finite value.
            Sign::Negative => (max_finite(Sign::Negative, p), Status::OK),
        },
        // next_up(±0) = the smallest positive value.
        Class::Zero { .. } => (min_positive(p), Status::OK),
        Class::Normal {
            sign,
            exponent,
            mantissa,
        } => match sign {
            Sign::Positive => next_up_positive(*exponent, mantissa, p),
            Sign::Negative => next_up_negative(*exponent, mantissa, p),
        },
    }
}

/// `next_up` of a positive finite value: increase the magnitude by one
/// ulp, carrying into the exponent (and into `+∞` at the ceiling).
fn next_up_positive(exponent: i64, mantissa: &[u64], p: u32) -> (BigFloat, Status) {
    let limbs = limbs_for(p);
    let lz = low_zero(p, limbs);
    let mut storage = mantissa.to_vec();
    let carried_out = add_one_ulp(&mut storage, lz);
    if carried_out {
        // The mantissa was all-ones; the magnitude crosses a power of
        // two upward, so the new mantissa is the top bit only at e+1.
        if exponent == i64::MAX {
            // No exponent above the ceiling: the successor is +∞.
            return (
                BigFloat::try_new_infinity(Sign::Positive, p).unwrap(),
                Status::OK,
            );
        }
        let mut top = vec![0u64; limbs];
        top[limbs - 1] = 1u64 << 63;
        return (normal(Sign::Positive, exponent + 1, top, p), Status::OK);
    }
    (normal(Sign::Positive, exponent, storage, p), Status::OK)
}

/// `next_up` of a negative finite value: decrease the magnitude by one
/// ulp toward zero, borrowing across a power of two (and reaching `−0`
/// at the floor).
fn next_up_negative(exponent: i64, mantissa: &[u64], p: u32) -> (BigFloat, Status) {
    let limbs = limbs_for(p);
    let lz = low_zero(p, limbs);
    if is_power_of_two(mantissa, limbs) {
        // Value is −2^e. Decreasing the magnitude crosses a power of
        // two downward: the predecessor magnitude is the all-ones
        // mantissa at e−1.
        if exponent == i64::MIN {
            // −MinPos: nothing lies between it and zero. IEEE fixes the
            // sign of this zero as negative.
            return (
                BigFloat::try_new_zero(Sign::Negative, p).unwrap(),
                Status::OK,
            );
        }
        return (max_mantissa(Sign::Negative, exponent - 1, p), Status::OK);
    }
    // Interior of the binade: m − 1 stays top-bit-set, exponent fixed.
    let mut storage = mantissa.to_vec();
    sub_one_ulp(&mut storage, lz);
    (normal(Sign::Negative, exponent, storage, p), Status::OK)
}

fn ulp_kernel(x: &BigFloat) -> (BigFloat, Status) {
    let p = x.precision;
    match &x.class {
        Class::Nan {
            quiet,
            sign,
            payload,
        } => {
            if *quiet {
                let nan = BigFloat::try_new_quiet_nan(*sign, p, payload)
                    .expect("BigFloat invariant: precision >= 1");
                (nan, Status::OK)
            } else {
                // Preserve the input sign when quieting (IEEE leaves a
                // NaN result's sign unspecified; preserving it keeps
                // next_up/next_down/ulp mutually consistent under the
                // next_down = -next_up(-x) identity).
                let nan = BigFloat::try_new_quiet_nan(*sign, p, &[])
                    .expect("BigFloat invariant: precision >= 1");
                auto_raise(Status::INVALID);
                (nan, Status::INVALID)
            }
        }
        // ulp(±∞) = +∞ (the gap above the largest finite is unbounded).
        Class::Infinity { .. } => (
            BigFloat::try_new_infinity(Sign::Positive, p).unwrap(),
            Status::OK,
        ),
        // ulp(±0) = the smallest positive value (the gap above zero).
        Class::Zero { .. } => (min_positive(p), Status::OK),
        Class::Normal { exponent, .. } => {
            // ulp = 2^(e − p + 1): a power of two, top-bit-only mantissa,
            // positive sign. e − p + 1 ≤ e ≤ i64::MAX cannot overflow
            // upward; it can only underflow below i64::MIN, in which
            // case the true ulp is below MinPos and we saturate upward.
            let wide = i128::from(*exponent) - i128::from(p) + 1;
            let limbs = limbs_for(p);
            let mut top = vec![0u64; limbs];
            top[limbs - 1] = 1u64 << 63;
            if wide < i128::from(i64::MIN) {
                let saturated = normal(Sign::Positive, i64::MIN, top, p);
                auto_raise(Status::UNDERFLOW);
                (saturated, Status::UNDERFLOW)
            } else {
                (normal(Sign::Positive, wide as i64, top, p), Status::OK)
            }
        }
    }
}

// ---------- internal builders ----------

/// Assemble a `Normal` directly from validated components. Internal to
/// the crate (no raw-parts constructor is exposed); callers here supply
/// a canonical top-bit-set mantissa of `limbs_for(precision)` limbs.
fn normal(sign: Sign, exponent: i64, mantissa: alloc::vec::Vec<u64>, precision: u32) -> BigFloat {
    debug_assert_eq!(mantissa.len(), limbs_for(precision));
    debug_assert!(
        mantissa[mantissa.len() - 1] & (1u64 << 63) != 0,
        "top bit must be set"
    );
    BigFloat {
        class: Class::Normal {
            sign,
            exponent,
            mantissa,
        },
        precision,
    }
}

/// The smallest positive value at `precision`: top-bit-only mantissa at
/// the exponent floor `i64::MIN` (value `2^(i64::MIN)`).
fn min_positive(precision: u32) -> BigFloat {
    let limbs = limbs_for(precision);
    let mut mantissa = vec![0u64; limbs];
    mantissa[limbs - 1] = 1u64 << 63;
    normal(Sign::Positive, i64::MIN, mantissa, precision)
}

/// An all-ones (`2^precision − 1`) mantissa at the given sign and
/// exponent. With `exponent = i64::MAX` this is `±MaxFinite`.
fn max_mantissa(sign: Sign, exponent: i64, precision: u32) -> BigFloat {
    let limbs = limbs_for(precision);
    let mut mantissa = vec![u64::MAX; limbs];
    clear_low(&mut mantissa, low_zero(precision, limbs));
    normal(sign, exponent, mantissa, precision)
}

/// `±MaxFinite`: the largest finite magnitude at `precision`
/// (all-ones mantissa, exponent ceiling `i64::MAX`).
fn max_finite(sign: Sign, precision: u32) -> BigFloat {
    max_mantissa(sign, i64::MAX, precision)
}

// ---------- internal bit helpers ----------

/// Absolute storage bit position of the mantissa's least-significant
/// (value) bit: the `limbs*64 − precision` low bits are always zero in
/// canonical form, so the value's LSB sits at this position.
#[inline]
fn low_zero(precision: u32, limbs: usize) -> usize {
    ((limbs as u64) * 64 - u64::from(precision)) as usize
}

/// `true` when `mantissa` is a power of two in canonical form: the top
/// bit of the most-significant limb set and every other bit zero.
fn is_power_of_two(mantissa: &[u64], limbs: usize) -> bool {
    if mantissa[limbs - 1] != 1u64 << 63 {
        return false;
    }
    mantissa[..limbs - 1].iter().all(|&l| l == 0)
}

/// Add one ulp (a `1` at storage bit `lz`), propagating the carry
/// upward. Returns `true` when the carry leaves the top of the storage
/// (the mantissa was all-ones), signalling a power-of-two renormalize.
fn add_one_ulp(storage: &mut [u64], lz: usize) -> bool {
    let limb_idx = lz / 64;
    let bit = lz % 64;
    let mut carry = 1u64 << bit;
    for limb in &mut storage[limb_idx..] {
        let (sum, overflowed) = limb.overflowing_add(carry);
        *limb = sum;
        if overflowed {
            carry = 1;
        } else {
            return false;
        }
    }
    carry != 0
}

/// Subtract one ulp (a `1` at storage bit `lz`), propagating the borrow
/// upward. Precondition: the result stays top-bit-set (the caller only
/// invokes this on a non-power-of-two mantissa), so the borrow never
/// leaves the top of the storage.
fn sub_one_ulp(storage: &mut [u64], lz: usize) {
    let limb_idx = lz / 64;
    let bit = lz % 64;
    let mut borrow = 1u64 << bit;
    for limb in &mut storage[limb_idx..] {
        let (diff, underflowed) = limb.overflowing_sub(borrow);
        *limb = diff;
        if underflowed {
            borrow = 1;
        } else {
            return;
        }
    }
    debug_assert!(
        false,
        "sub_one_ulp borrowed past the top of a non-power-of-two mantissa"
    );
}

/// Zero the `count` lowest bits of `storage`.
fn clear_low(storage: &mut [u64], count: usize) {
    let len = storage.len();
    let full = count / 64;
    let partial = count % 64;
    for limb in &mut storage[..full.min(len)] {
        *limb = 0;
    }
    if partial > 0 && full < len {
        storage[full] &= !((1u64 << partial) - 1);
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::rounding::RoundingMode;
    use core::cmp::Ordering;

    fn from_i64(n: i64, p: u32) -> BigFloat {
        BigFloat::try_from_i64_exact(n, p).unwrap()
    }

    fn eq(a: &BigFloat, b: &BigFloat) -> bool {
        a.partial_cmp(b).0 == Some(Ordering::Equal)
    }

    fn less(a: &BigFloat, b: &BigFloat) -> bool {
        a.partial_cmp(b).0 == Some(Ordering::Less)
    }

    // ---------- adjacency / strict monotonicity ----------

    #[test]
    fn next_up_is_strictly_greater() {
        for &(n, p) in &[(1i64, 53), (3, 53), (-3, 53), (42, 113), (-42, 256), (7, 8)] {
            let x = from_i64(n, p);
            let (up, s) = x.next_up();
            assert!(s.is_ok());
            assert!(less(&x, &up), "next_up({n}) must exceed {n} at p={p}");
            let (down, _) = x.next_down();
            assert!(less(&down, &x), "next_down({n}) must be below {n} at p={p}");
            assert_eq!(up.precision(), p);
            assert_eq!(down.precision(), p);
        }
    }

    #[test]
    fn nothing_strictly_between_x_and_next_up() {
        // next_down(next_up(x)) == x for an interior value.
        for &(n, p) in &[(5i64, 53), (-5, 53), (1, 64), (-1, 8)] {
            let x = from_i64(n, p);
            let (up, _) = x.next_up();
            let (back, _) = up.next_down();
            assert!(eq(&back, &x), "round trip up/down failed for {n} at p={p}");
            let (down, _) = x.next_down();
            let (back2, _) = down.next_up();
            assert!(eq(&back2, &x), "round trip down/up failed for {n} at p={p}");
        }
    }

    #[test]
    fn next_up_equals_add_one_ulp_for_interior() {
        // For an interior value, next_up(x) == x + ulp(x) under any
        // rounding mode (the neighbour is exactly representable).
        for &(n, p) in &[(5i64, 53), (-5, 53), (3, 113), (-7, 8), (1, 200)] {
            let x = from_i64(n, p);
            let (u, _) = x.ulp();
            let (sum, ss) = x.add(&u, RoundingMode::NearestEven);
            assert!(ss.is_ok(), "x + ulp must be exact for {n} at p={p}");
            let (up, _) = x.next_up();
            assert!(eq(&sum, &up), "x + ulp != next_up for {n} at p={p}");
        }
    }

    // ---------- power-of-two binade crossings ----------

    #[test]
    fn next_up_crosses_power_of_two_upward() {
        // 2 at precision 2: mantissa 0b10 (top-bit-only), exponent 1.
        // next_up should be 3 = 0b11 at exponent 1.
        let two = from_i64(2, 2);
        let (up, _) = two.next_up();
        assert!(eq(&up, &from_i64(3, 2)));
        // next_up(3) crosses upward: 3 is all-ones (0b11) at e=1, so the
        // successor is 4 = top-bit-only at e=2.
        let (up2, _) = up.next_up();
        assert!(eq(&up2, &from_i64(4, 2)));
    }

    #[test]
    fn next_up_negative_crosses_power_of_two_downward() {
        // -2 at precision 2 (top-bit-only). next_up toward zero is the
        // largest magnitude below 2, i.e. -1.5 = -(0b11 at e=0).
        let neg_two = from_i64(-2, 2);
        let (up, _) = neg_two.next_up();
        // -1.5 = -3 * 2^-1; build it as next_down of -1 ... instead
        // check numerically: up should be -1.5, strictly between -2 and
        // -1.
        assert!(less(&neg_two, &up));
        assert!(less(&up, &from_i64(-1, 2)));
        // And exactly -1.5: -3 scaled by 2^-1.
        let (neg_one_point_five, _) = from_i64(-3, 2).scale_by_pow2(-1);
        assert!(eq(&up, &neg_one_point_five));
    }

    // ---------- exponent-floor / ceiling boundaries ----------

    #[test]
    fn next_up_of_zero_is_min_positive() {
        for sign in [Sign::Positive, Sign::Negative] {
            let z = BigFloat::try_new_zero(sign, 53).unwrap();
            let (up, s) = z.next_up();
            assert!(s.is_ok());
            assert!(up.is_normal() && up.is_sign_positive());
            assert!(eq(&up, &min_positive(53)));
        }
    }

    #[test]
    fn next_down_of_zero_is_neg_min_positive() {
        for sign in [Sign::Positive, Sign::Negative] {
            let z = BigFloat::try_new_zero(sign, 53).unwrap();
            let (down, s) = z.next_down();
            assert!(s.is_ok());
            assert!(down.is_normal() && down.is_sign_negative());
            assert!(eq(&down, &min_positive(53).negated()));
        }
    }

    #[test]
    fn next_up_of_neg_min_positive_is_neg_zero() {
        let neg_min = min_positive(53).negated();
        let (up, s) = neg_min.next_up();
        assert!(s.is_ok());
        assert!(up.is_zero() && up.is_sign_negative(), "must be -0");
    }

    #[test]
    fn next_up_of_max_finite_is_pos_infinity() {
        let max = max_finite(Sign::Positive, 53);
        let (up, s) = max.next_up();
        assert!(s.is_ok());
        assert!(up.is_infinite() && up.is_sign_positive());
    }

    #[test]
    fn next_up_of_neg_infinity_is_neg_max_finite() {
        let ni = BigFloat::try_new_infinity(Sign::Negative, 53).unwrap();
        let (up, s) = ni.next_up();
        assert!(s.is_ok());
        assert!(eq(&up, &max_finite(Sign::Negative, 53)));
        assert!(up.is_sign_negative());
    }

    #[test]
    fn next_up_of_pos_infinity_is_pos_infinity() {
        let pi = BigFloat::try_new_infinity(Sign::Positive, 53).unwrap();
        let (up, s) = pi.next_up();
        assert!(s.is_ok());
        assert!(up.is_infinite() && up.is_sign_positive());
    }

    #[test]
    fn next_down_of_pos_infinity_is_max_finite() {
        let pi = BigFloat::try_new_infinity(Sign::Positive, 53).unwrap();
        let (down, s) = pi.next_down();
        assert!(s.is_ok());
        assert!(eq(&down, &max_finite(Sign::Positive, 53)));
    }

    #[test]
    fn boundary_round_trips() {
        // next_down(next_up(MaxFinite)) == MaxFinite (via +∞).
        let max = max_finite(Sign::Positive, 53);
        let (up, _) = max.next_up();
        let (back, _) = up.next_down();
        assert!(eq(&back, &max));
    }

    // ---------- ulp ----------

    #[test]
    fn ulp_of_one_is_two_pow_minus_p_minus_1() {
        // 1.0 at precision p: exponent 0, ulp = 2^(0 - p + 1) = 2^(1-p).
        let one = from_i64(1, 53);
        let (u, s) = one.ulp();
        assert!(s.is_ok());
        let (expected, _) = from_i64(1, 53).scale_by_pow2(1 - 53);
        assert!(eq(&u, &expected));
        assert_eq!(u.precision(), 53);
        assert!(u.is_sign_positive());
    }

    #[test]
    fn ulp_is_sign_invariant() {
        let pos = from_i64(7, 53);
        let neg = from_i64(-7, 53);
        let (up, _) = pos.ulp();
        let (un, _) = neg.ulp();
        assert!(eq(&up, &un));
        assert!(up.is_sign_positive() && un.is_sign_positive());
    }

    #[test]
    fn ulp_of_zero_is_min_positive() {
        let z = BigFloat::try_new_zero(Sign::Positive, 53).unwrap();
        let (u, s) = z.ulp();
        assert!(s.is_ok());
        assert!(eq(&u, &min_positive(53)));
    }

    #[test]
    fn ulp_of_infinity_is_infinity() {
        let pi = BigFloat::try_new_infinity(Sign::Positive, 53).unwrap();
        let (u, s) = pi.ulp();
        assert!(s.is_ok());
        assert!(u.is_infinite() && u.is_sign_positive());
        let ni = BigFloat::try_new_infinity(Sign::Negative, 53).unwrap();
        let (u2, _) = ni.ulp();
        assert!(u2.is_infinite() && u2.is_sign_positive());
    }

    #[test]
    fn ulp_saturates_upward_near_exponent_floor() {
        // A value whose exponent − p + 1 falls below i64::MIN: its true
        // ulp is below MinPos and must saturate upward with UNDERFLOW.
        let (near_floor, _) = from_i64(1, 53).scale_by_pow2(i64::MIN + 10);
        let (u, s) = near_floor.ulp();
        assert!(s.underflow());
        assert!(eq(&u, &min_positive(53)));
    }

    #[test]
    fn ulp_at_floor_is_exact_when_representable() {
        // MinPos itself: exponent i64::MIN, p=1, ulp = 2^(i64::MIN). At
        // p=1 the value is the mantissa weight, so ulp == MinPos exactly
        // and no saturation occurs.
        let mp = min_positive(1);
        let (u, s) = mp.ulp();
        assert!(s.is_ok(), "p=1 ulp at floor is exactly representable");
        assert!(eq(&u, &mp));
    }

    // ---------- NaN handling ----------

    #[test]
    fn quiet_nan_passes_through() {
        let q = BigFloat::try_new_quiet_nan(Sign::Negative, 53, &[9]).unwrap();
        for (v, s) in [q.next_up(), q.next_down(), q.ulp()] {
            assert!(s.is_ok());
            assert!(v.is_quiet_nan());
        }
    }

    #[test]
    fn signaling_nan_raises_invalid() {
        let sn = BigFloat::try_new_signaling_nan(Sign::Positive, 53, &[]).unwrap();
        for (v, s) in [sn.next_up(), sn.next_down(), sn.ulp()] {
            assert!(s.invalid());
            assert!(v.is_quiet_nan());
        }
    }

    #[test]
    fn nan_sign_preserved_consistently() {
        // The three operations agree on the quieted NaN's sign for a
        // signaling input (the cross-op asymmetry the adversarial review
        // flagged). IEEE leaves it unspecified; pfloat preserves the
        // input sign so next_up/next_down/ulp do not disagree.
        for sign in [Sign::Positive, Sign::Negative] {
            let sn = BigFloat::try_new_signaling_nan(sign, 53, &[]).unwrap();
            let want_neg = matches!(sign, Sign::Negative);
            for (v, _) in [sn.next_up(), sn.next_down(), sn.ulp()] {
                assert!(v.is_quiet_nan());
                assert_eq!(v.is_sign_negative(), want_neg);
            }
            // Quiet NaN likewise keeps its sign across all three.
            let qn = BigFloat::try_new_quiet_nan(sign, 53, &[]).unwrap();
            for (v, _) in [qn.next_up(), qn.next_down(), qn.ulp()] {
                assert_eq!(v.is_sign_negative(), want_neg);
            }
        }
    }

    // ---------- multi-limb precision ----------

    #[test]
    fn multi_limb_increment_and_borrow() {
        // precision 200 spans 4 limbs. Exercise an interior increment
        // and the round trip at a precision where the mantissa is not a
        // single limb.
        let x = from_i64(123_456_789, 200);
        let (up, _) = x.next_up();
        assert!(less(&x, &up));
        let (back, _) = up.next_down();
        assert!(eq(&back, &x));
        // ulp consistency at multi-limb precision.
        let (u, _) = x.ulp();
        let (sum, ss) = x.add(&u, RoundingMode::NearestEven);
        assert!(ss.is_ok());
        assert!(eq(&sum, &up));
    }

    #[test]
    fn with_flags_siblings_accumulate() {
        let sn = BigFloat::try_new_signaling_nan(Sign::Positive, 53, &[]).unwrap();
        let mut flags = Status::OK;
        let _ = sn.next_up_with_flags(&mut flags);
        assert!(flags.invalid());
        let _ = sn.next_down_with_flags(&mut flags);
        let (near_floor, _) = from_i64(1, 53).scale_by_pow2(i64::MIN + 10);
        let _ = near_floor.ulp_with_flags(&mut flags);
        assert!(flags.invalid() && flags.underflow());
    }
}
