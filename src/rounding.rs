//! Rounding mode and the universal rounding pipeline.
//!
//! [`RoundingMode`] enumerates the five IEEE 754-2019 §4.3 rounding
//! attributes. [`round_finite_to_precision`] is the funnel every
//! arithmetic kernel routes through to discharge IEEE 754-2019 §6.5
//! rounding behavior. The pipeline matches ferrodec's
//! `src/ops/round.rs::round_and_pack_finite` adapted for
//! arbitrary-precision binary mantissas.
//!
//! Slice 1b (this slice) ships the funnel and exercises it through
//! [`BigFloat::try_from_i64_round`](crate::big::BigFloat::try_from_i64_round)
//! and [`BigFloat::round_to_precision`](crate::big::BigFloat::round_to_precision).
//! Slices 1c–1f connect the arithmetic kernels.

#[cfg(feature = "alloc")]
extern crate alloc;

#[cfg(feature = "big")]
use alloc::vec;
#[cfg(feature = "big")]
use alloc::vec::Vec;

#[cfg(feature = "big")]
use crate::big::BigFloat;
#[cfg(feature = "big")]
use crate::class::Class;
use crate::mantissa::limbs_for;
use crate::sign::Sign;
use crate::status::Status;

/// IEEE 754-2019 §4.3 rounding-direction attribute.
///
/// Five modes, each implementing one of §4.3.1's rounding rules.
/// The default is [`NearestEven`](Self::NearestEven), matching IEEE
/// 754's default rounding for arithmetic operations.
#[derive(Clone, Copy, Debug, Default, Eq, PartialEq, Hash)]
#[cfg_attr(feature = "serde", derive(serde::Serialize, serde::Deserialize))]
pub enum RoundingMode {
    /// Round to nearest, ties to even (IEEE 754 default).
    /// Also called "banker's rounding."
    #[default]
    NearestEven,

    /// Round to nearest, ties away from zero.
    NearestAway,

    /// Round toward zero (truncate).
    TowardZero,

    /// Round toward `+∞` (ceiling for positive values, floor for
    /// negative).
    TowardPositive,

    /// Round toward `-∞` (floor for positive values, ceiling for
    /// negative).
    TowardNegative,
}

/// Round an intermediate finite mantissa to a target precision.
///
/// This is the funnel every arithmetic kernel routes through. The
/// caller provides:
///
/// - `sign`, `exponent`: the value's sign and the unbiased binary
///   exponent of the intermediate mantissa.
/// - `intermediate_mantissa`: little-endian limbs holding the
///   intermediate's `intermediate_precision`-bit mantissa, top-bit-
///   set per ADR-0001's normalization rule. Bits below the
///   precision must be zero in canonical form.
/// - `intermediate_precision`: the precision of the intermediate
///   (the number of meaningful bits). Must satisfy
///   `intermediate_precision >= 1`.
/// - `pre_sticky`: any rounding-bit information already accumulated
///   upstream (e.g., during exponent alignment in add/sub).
/// - `target_precision`: the precision to round to. Must satisfy
///   `target_precision >= 1`.
/// - `mode`: the [`RoundingMode`].
///
/// Returns the rounded [`BigFloat`] at `target_precision` and a
/// [`Status`] carrying [`Status::INEXACT`] if rounding discarded
/// any non-zero bit, plus [`Status::OVERFLOW`] if a round-up caused
/// the exponent to saturate `i64::MAX`.
///
/// Behavior is undefined (debug-asserts) when:
/// - `intermediate_precision == 0` or `target_precision == 0`.
/// - `intermediate_mantissa` is not top-bit-set normalized.
/// - `intermediate_mantissa.len() != limbs_for(intermediate_precision)`.
///
/// # Direction of rounding
///
/// The pipeline computes a guard bit (the bit immediately below the
/// kept mantissa) and a sticky bit (the OR of all bits below the
/// guard, plus `pre_sticky`). Together with the lowest kept bit and
/// the rounding mode, these determine whether to round up.
///
/// IEEE 754-2019 §4.3.1 rules:
/// - `NearestEven`: round up iff `guard && (sticky || lowest_kept)`.
/// - `NearestAway`: round up iff `guard`.
/// - `TowardZero`: never round up.
/// - `TowardPositive`: round up iff inexact and `sign == Positive`.
/// - `TowardNegative`: round up iff inexact and `sign == Negative`.
///
/// A round-up that overflows the mantissa width (e.g., `0b11 + 1 =
/// 0b100` at 2-bit precision) renormalizes by shifting right one bit
/// and incrementing the exponent.
#[cfg(feature = "big")]
#[allow(clippy::too_many_arguments)] // matches ferrodec's round_and_pack_finite
pub(crate) fn round_finite_to_precision(
    sign: Sign,
    exponent: i64,
    intermediate_mantissa: &[u64],
    intermediate_precision: u32,
    pre_sticky: bool,
    target_precision: u32,
    mode: RoundingMode,
) -> (BigFloat, Status) {
    debug_assert!(intermediate_precision >= 1);
    debug_assert!(target_precision >= 1);
    debug_assert_eq!(
        intermediate_mantissa.len(),
        limbs_for(intermediate_precision),
        "intermediate mantissa storage must match its precision"
    );

    if target_precision >= intermediate_precision {
        // Exact extension: copy the intermediate into a wider
        // storage with low bits zero. INEXACT depends only on
        // pre_sticky.
        return extend_exact(
            sign,
            exponent,
            intermediate_mantissa,
            intermediate_precision,
            pre_sticky,
            target_precision,
        );
    }

    let drop_bits = intermediate_precision - target_precision;
    let intermediate_low_zero = (intermediate_mantissa.len() as u32) * 64 - intermediate_precision;

    // Absolute storage bit positions:
    //   - intermediate's mantissa LSB:    intermediate_low_zero
    //   - intermediate's mantissa MSB:    intermediate_low_zero + intermediate_precision - 1
    // Within the kept high portion:
    //   - kept LSB position (within intermediate storage):
    //     intermediate_low_zero + drop_bits
    //   - guard bit position: kept LSB - 1 = intermediate_low_zero + drop_bits - 1
    //   - sticky range: bits at storage positions [intermediate_low_zero, guard - 1]
    //     plus any bits below intermediate_low_zero (which are zero by invariant)
    let kept_lsb_storage_pos = intermediate_low_zero + drop_bits;
    let guard_storage_pos = kept_lsb_storage_pos - 1;

    let guard = bit_at(intermediate_mantissa, guard_storage_pos as usize);
    let sticky_below = any_bit_below(intermediate_mantissa, guard_storage_pos as usize);
    let sticky = pre_sticky || sticky_below;

    let lowest_kept = bit_at(intermediate_mantissa, kept_lsb_storage_pos as usize);

    let inexact = guard || sticky;
    let round_up = if inexact {
        match mode {
            RoundingMode::NearestEven => guard && (sticky || lowest_kept),
            RoundingMode::NearestAway => guard,
            RoundingMode::TowardZero => false,
            RoundingMode::TowardPositive => matches!(sign, Sign::Positive),
            RoundingMode::TowardNegative => matches!(sign, Sign::Negative),
        }
    } else {
        false
    };

    // Build the truncated mantissa in target storage by copying the
    // top target_limbs limbs from the intermediate, then masking
    // off the now-redundant low bits (which were guard / sticky in
    // the intermediate frame but should be zero in the target frame).
    let target_limbs = limbs_for(target_precision);
    let mut storage: Vec<u64> = vec![0u64; target_limbs];

    // Source range in intermediate storage: top target_limbs limbs.
    let source_start = intermediate_mantissa.len() - target_limbs;
    storage[..].copy_from_slice(&intermediate_mantissa[source_start..]);

    // Mask the low (target_limbs * 64 - target_precision) bits to
    // zero in the target storage.
    let target_low_zero = (target_limbs as u32) * 64 - target_precision;
    clear_low_bits(&mut storage, target_low_zero as usize);

    let mut new_exponent = exponent;
    let mut status = Status::OK;

    if round_up {
        let lsb_storage_pos = target_low_zero as usize;
        let overflowed = increment_at(&mut storage, lsb_storage_pos);
        if overflowed {
            // Mantissa was all-ones at target_precision; +1 wraps
            // to 0 and we need to renormalize. The new mantissa is
            // 2^(target_precision - 1) (top bit only), exponent +1.
            storage.fill(0);
            storage[target_limbs - 1] = 1u64 << 63;
            match new_exponent.checked_add(1) {
                Some(e) => new_exponent = e,
                None => {
                    // Astronomically unreachable but documented in
                    // the OVERFLOW flag's contract.
                    status |= Status::OVERFLOW;
                    new_exponent = i64::MAX;
                }
            }
        }
    }

    if inexact {
        status |= Status::INEXACT;
    }

    (
        BigFloat {
            class: Class::Normal {
                sign,
                exponent: new_exponent,
                mantissa: storage,
            },
            precision: target_precision,
        },
        status,
    )
}

/// Pad a finite mantissa with trailing zeros to reach a wider
/// target precision. No rounding required.
#[cfg(feature = "big")]
fn extend_exact(
    sign: Sign,
    exponent: i64,
    intermediate_mantissa: &[u64],
    intermediate_precision: u32,
    pre_sticky: bool,
    target_precision: u32,
) -> (BigFloat, Status) {
    debug_assert!(target_precision >= intermediate_precision);

    let target_limbs = limbs_for(target_precision);
    let mut storage: Vec<u64> = vec![0u64; target_limbs];

    // Place the intermediate into the top of the target storage,
    // top-bit-aligned: source's MSL bit 63 maps to target's MSL bit
    // 63. Both storages have their high data at the top of their
    // respective MSLs (top-bit-set normalization), so the shift
    // between source and target is the difference in limb count
    // times 64 (always a whole-limb shift).
    let source = intermediate_mantissa;
    let limb_diff = (target_limbs - source.len()) as i64;
    let shift_bits = limb_diff * 64;
    place_top_aligned(&mut storage, source, shift_bits);

    let status = if pre_sticky {
        Status::INEXACT
    } else {
        Status::OK
    };

    (
        BigFloat {
            class: Class::Normal {
                sign,
                exponent,
                mantissa: storage,
            },
            precision: target_precision,
        },
        status,
    )
}

/// Read a single bit at an absolute storage position.
#[cfg(feature = "big")]
fn bit_at(storage: &[u64], position: usize) -> bool {
    let limb_idx = position / 64;
    let bit_in_limb = position % 64;
    if limb_idx >= storage.len() {
        return false;
    }
    (storage[limb_idx] >> bit_in_limb) & 1 == 1
}

/// Returns `true` if any bit at storage positions `[0, position)`
/// is set.
#[cfg(feature = "big")]
fn any_bit_below(storage: &[u64], position: usize) -> bool {
    if position == 0 {
        return false;
    }
    let last_full_limb = position / 64;
    let bits_in_last_limb = position % 64;

    for &limb in &storage[..last_full_limb.min(storage.len())] {
        if limb != 0 {
            return true;
        }
    }
    if bits_in_last_limb > 0 && last_full_limb < storage.len() {
        let mask = (1u64 << bits_in_last_limb) - 1;
        if storage[last_full_limb] & mask != 0 {
            return true;
        }
    }
    false
}

/// Zero the bits in `storage` at absolute positions `[0, count)`.
#[cfg(feature = "big")]
fn clear_low_bits(storage: &mut [u64], count: usize) {
    if count == 0 {
        return;
    }
    let full_zero_limbs = count / 64;
    let partial_bits = count % 64;
    let zero_end = full_zero_limbs.min(storage.len());
    for limb in &mut storage[..zero_end] {
        *limb = 0;
    }
    if partial_bits > 0 && full_zero_limbs < storage.len() {
        let mask = !((1u64 << partial_bits) - 1);
        storage[full_zero_limbs] &= mask;
    }
}

/// Add `1` at absolute bit position `position`, propagating carries
/// upward. Returns `true` if the carry propagated past the top of
/// the storage (i.e., the increment overflowed the available bits).
#[cfg(feature = "big")]
fn increment_at(storage: &mut [u64], position: usize) -> bool {
    let limb_idx = position / 64;
    let bit_in_limb = position % 64;
    debug_assert!(limb_idx < storage.len());

    let mut carry: u64 = 1u64 << bit_in_limb;
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

/// Copy `source` into `target` shifted by `shift_bits` bits. A
/// positive `shift_bits` shifts the source bits up (toward the MSL);
/// a negative shift moves them down. The source is treated as
/// occupying its own storage layout; this function realigns it into
/// the target's frame.
///
/// For slice 1b's [`extend_exact`] the shift is always non-negative
/// (target storage is at least as wide as source storage); we
/// implement only that direction here.
#[cfg(feature = "big")]
fn place_top_aligned(target: &mut [u64], source: &[u64], shift_bits: i64) {
    debug_assert!(shift_bits >= 0, "extend_exact uses non-negative shift only");
    let shift = shift_bits as u64;
    let limb_shift = (shift / 64) as usize;
    let bit_shift = (shift % 64) as u32;

    if bit_shift == 0 {
        // Whole-limb shift only.
        for (i, &s) in source.iter().enumerate() {
            if let Some(slot) = target.get_mut(i + limb_shift) {
                *slot = s;
            }
        }
    } else {
        // Mixed shift: each source limb spans two target limbs.
        let inv_shift = 64 - bit_shift;
        for (i, &s) in source.iter().enumerate() {
            let lo = s << bit_shift;
            let hi = s >> inv_shift;
            if let Some(slot) = target.get_mut(i + limb_shift) {
                *slot |= lo;
            }
            if let Some(slot) = target.get_mut(i + limb_shift + 1) {
                *slot |= hi;
            }
        }
    }
}

/// Correctly round `base` perturbed by an infinitesimal residue to
/// `target_precision` under `mode`.
///
/// The residue is a single bit placed strictly below both `base`'s and
/// the target's least-significant bit, so it carries only direction and
/// stickiness for the rounding and never reaches the guard position.
/// The signed result is `result_sign · (|base| ± ε)`; the result sign
/// is supplied explicitly because callers (the add/sub huge-gap path)
/// track the effective sign separately from the operand's stored sign.
/// `subtracts_magnitude` selects whether the magnitude shrinks (an
/// opposite-sign add whose smaller operand lies far past the rounding
/// boundary, or `log1p` of a tiny positive argument) or grows (a
/// same-sign add, or `expm1` of a tiny positive argument). The result
/// is always inexact.
///
/// Implemented as an exact wide add of the residue bit followed by a
/// single mode-aware round, so it inherits the rounding pipeline's
/// correctness for every mode, including the borrow renormalisation
/// when `base` is a power of two. `base` must be a finite non-zero
/// `Normal`. Review 2026-05-29 (add/sub huge-gap directed rounding;
/// expm1/log1p tiny-x collapse).
#[cfg(feature = "big")]
pub(crate) fn round_with_infinitesimal(
    base: &BigFloat,
    result_sign: Sign,
    subtracts_magnitude: bool,
    target_precision: u32,
    mode: RoundingMode,
) -> (BigFloat, Status) {
    let e = match &base.class {
        Class::Normal { exponent, .. } => *exponent,
        _ => unreachable!("round_with_infinitesimal requires a finite non-zero base"),
    };

    // Operate on `result_sign · |base|` so the directed-mode rounding
    // sees the correct sign even when it differs from `base`'s stored
    // sign (the huge-gap subtraction case).
    let mut signed_base = base.clone();
    if let Class::Normal { sign, .. } = &mut signed_base.class {
        *sign = result_sign;
    }

    // The residue bit sits three places below max(precision, target):
    // strictly under the target rounding boundary, with room for a
    // sticky bit beneath the guard after any borrow renormalisation.
    let wide = base.precision.max(target_precision).saturating_add(3);
    let delta_sign = if subtracts_magnitude {
        result_sign.flip()
    } else {
        result_sign
    };
    let mut delta_mant = vec![0u64; limbs_for(1)];
    let top = delta_mant.len() - 1;
    delta_mant[top] = 1u64 << 63;
    let delta = BigFloat {
        class: Class::Normal {
            sign: delta_sign,
            exponent: e.saturating_sub(i64::from(wide)).saturating_add(1),
            mantissa: delta_mant,
        },
        precision: 1,
    };

    // The sum spans at most `wide` significant bits, so the add is
    // exact; the single subsequent round applies the caller's mode and
    // auto-raises the resulting INEXACT.
    let (sum, _) = signed_base
        .add_round(&delta, wide, RoundingMode::NearestEven)
        .expect("wide >= 1");
    sum.round_to_precision(target_precision, mode)
        .expect("target_precision >= 1")
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn rounding_mode_default_is_nearest_even() {
        assert_eq!(RoundingMode::default(), RoundingMode::NearestEven);
    }

    #[test]
    fn bit_at_helper() {
        let storage = [0xFFu64, 0x0u64];
        assert!(bit_at(&storage, 0));
        assert!(bit_at(&storage, 7));
        assert!(!bit_at(&storage, 8));
        assert!(!bit_at(&storage, 64));
        assert!(!bit_at(&storage, 200)); // out of range, returns false
    }

    #[test]
    fn any_bit_below_helper() {
        let storage = [0b1010u64, 0u64];
        assert!(!any_bit_below(&storage, 0));
        assert!(!any_bit_below(&storage, 1));
        assert!(any_bit_below(&storage, 2)); // bit 1 set
        assert!(any_bit_below(&storage, 64));
        assert!(!any_bit_below(&[0u64, 0u64], 128));
    }

    #[test]
    fn clear_low_bits_helper() {
        // `clear_low_bits(&mut s, n)` zeroes bit positions [0, n).
        let mut s = [0xFFFF_FFFF_FFFF_FFFFu64, 0xFFFF_FFFF_FFFF_FFFFu64];
        clear_low_bits(&mut s, 65);
        // Positions 0..63 = limb 0 (all cleared).
        assert_eq!(s[0], 0);
        // Position 64 = bit 0 of limb 1 (cleared); positions
        // 65..127 untouched.
        assert_eq!(s[1] & 0b1, 0);
        assert_eq!(s[1] >> 1, 0xFFFF_FFFF_FFFF_FFFFu64 >> 1);
    }

    #[test]
    fn increment_at_no_overflow() {
        let mut s = [0u64, 0u64];
        let overflow = increment_at(&mut s, 0);
        assert!(!overflow);
        assert_eq!(s[0], 1);

        let mut s = [0u64, 0u64];
        let overflow = increment_at(&mut s, 64);
        assert!(!overflow);
        assert_eq!(s[1], 1);
    }

    #[test]
    fn increment_at_with_carry() {
        let mut s = [0xFFFF_FFFF_FFFF_FFFFu64, 0u64];
        let overflow = increment_at(&mut s, 0);
        assert!(!overflow);
        assert_eq!(s[0], 0);
        assert_eq!(s[1], 1);
    }

    #[test]
    fn increment_at_overflow_past_top() {
        let mut s = [0xFFFF_FFFF_FFFF_FFFFu64];
        let overflow = increment_at(&mut s, 0);
        assert!(overflow);
        assert_eq!(s[0], 0);
    }
}
