//! Multi-precision integer primitives shared by the arithmetic
//! kernels.
//!
//! The buffer layout these helpers operate on is the "integer view":
//! limb 0 holds bits 0..63, limb 1 holds bits 64..127, and so on
//! (little-endian limb order). This is distinct from pfloat's
//! storage layout for mantissas, which is top-bit-set (top of an
//! integer-as-mantissa lives at the MSL bit 63). Conversion between
//! the two views is the job of [`extract_value_limbs`] and
//! [`or_left_shifted_into`].
//!
//! All helpers here are `pub(crate)`; they support
//! [`crate::ops::addsub`], [`crate::ops::mul`] (slice 1d), and
//! later kernels. Slice 1d split them out of `addsub.rs` so mul
//! could share them without circular imports.

use alloc::vec;
use alloc::vec::Vec;

/// Extract the `bits` consecutive bits at storage position
/// `low_bit_pos` from `src`, packing them as a little-endian limb
/// array (bit 0 of the extracted value at `out[0]` bit 0).
///
/// Used to convert a top-aligned mantissa into the integer view.
pub(crate) fn extract_value_limbs(src: &[u64], low_bit_pos: u32, bits: u32, out: &mut [u64]) {
    if bits == 0 {
        return;
    }
    let limb_offset = (low_bit_pos / 64) as usize;
    let bit_offset = low_bit_pos % 64;

    if bit_offset == 0 {
        let copy_count = bits.div_ceil(64) as usize;
        for i in 0..copy_count {
            if let Some(&v) = src.get(limb_offset + i) {
                if let Some(slot) = out.get_mut(i) {
                    *slot = v;
                }
            }
        }
        let last_idx = ((bits - 1) / 64) as usize;
        let last_bit = bits % 64;
        if last_bit != 0 {
            if let Some(slot) = out.get_mut(last_idx) {
                let mask = (1u64 << last_bit) - 1;
                *slot &= mask;
            }
        }
        return;
    }

    let inv_offset = 64 - bit_offset;
    let mut i = 0usize;
    let mut bits_remaining = bits;
    while bits_remaining > 0 {
        let lo = src.get(limb_offset + i).copied().unwrap_or(0);
        let hi = src.get(limb_offset + i + 1).copied().unwrap_or(0);
        let combined = (lo >> bit_offset) | (hi << inv_offset);
        let take = bits_remaining.min(64);
        let mask = if take == 64 {
            !0u64
        } else {
            (1u64 << take) - 1
        };
        if let Some(slot) = out.get_mut(i) {
            *slot = combined & mask;
        }
        bits_remaining = bits_remaining.saturating_sub(64);
        i += 1;
    }
}

/// Convenience: extract a top-aligned mantissa as a bottom-aligned
/// integer of `precision` bits. The result Vec is
/// `limbs_for(precision)` long.
pub(crate) fn extract_as_integer(mantissa: &[u64], precision: u32) -> Vec<u64> {
    let storage_bits = (mantissa.len() as u32) * 64;
    let low_bit_pos = storage_bits - precision;
    let limb_count = crate::mantissa::limbs_for(precision);
    let mut out = vec![0u64; limb_count];
    extract_value_limbs(mantissa, low_bit_pos, precision, &mut out);
    out
}

/// OR `value` (a little-endian integer-view limb array of
/// `value_bits` bits) into `dst` shifted left by `left_shift`
/// storage-bit positions.
pub(crate) fn or_left_shifted_into(
    dst: &mut [u64],
    value: &[u64],
    value_bits: u32,
    left_shift: u32,
) {
    if value_bits == 0 {
        return;
    }
    let limb_shift = (left_shift / 64) as usize;
    let bit_shift = left_shift % 64;

    if bit_shift == 0 {
        for (i, &v) in value.iter().enumerate() {
            if let Some(slot) = dst.get_mut(i + limb_shift) {
                *slot |= v;
            }
        }
        return;
    }

    let inv_shift = 64 - bit_shift;
    for (i, &v) in value.iter().enumerate() {
        let lo = v << bit_shift;
        let hi = v >> inv_shift;
        if let Some(slot) = dst.get_mut(i + limb_shift) {
            *slot |= lo;
        }
        if let Some(slot) = dst.get_mut(i + limb_shift + 1) {
            *slot |= hi;
        }
    }
}

/// `dst += src`, returning the final carry (true if the addition
/// produced a carry out of the top of `dst`).
pub(crate) fn limbs_add_assign(dst: &mut [u64], src: &[u64]) -> bool {
    let mut carry: u64 = 0;
    for (i, &s) in src.iter().enumerate() {
        let Some(slot) = dst.get_mut(i) else {
            if s != 0 {
                return true;
            }
            continue;
        };
        let (sum1, c1) = slot.overflowing_add(s);
        let (sum2, c2) = sum1.overflowing_add(carry);
        *slot = sum2;
        carry = u64::from(c1) + u64::from(c2);
    }
    if carry == 0 {
        return false;
    }
    for slot in dst.iter_mut().skip(src.len()) {
        let (sum, c) = slot.overflowing_add(carry);
        *slot = sum;
        if c {
            carry = 1;
        } else {
            return false;
        }
    }
    carry != 0
}

/// `dst -= src`. Caller-supplied invariant: `dst >= src`.
pub(crate) fn limbs_sub_assign(dst: &mut [u64], src: &[u64]) {
    let mut borrow: u64 = 0;
    for (i, &s) in src.iter().enumerate() {
        let Some(slot) = dst.get_mut(i) else {
            debug_assert_eq!(s, 0, "sub overflow: src has bits past dst");
            continue;
        };
        let (diff1, b1) = slot.overflowing_sub(s);
        let (diff2, b2) = diff1.overflowing_sub(borrow);
        *slot = diff2;
        borrow = u64::from(b1) + u64::from(b2);
    }
    for slot in dst.iter_mut().skip(src.len()) {
        if borrow == 0 {
            return;
        }
        let (diff, b) = slot.overflowing_sub(borrow);
        *slot = diff;
        borrow = u64::from(b);
    }
    debug_assert_eq!(borrow, 0, "sub overflow: dst < src");
}

/// Returns the (zero-indexed) position of the most-significant set
/// bit of `buffer`, or `None` if all limbs are zero.
pub(crate) fn top_set_bit(buffer: &[u64]) -> Option<usize> {
    for (i, &limb) in buffer.iter().enumerate().rev() {
        if limb != 0 {
            let leading = limb.leading_zeros();
            let bit_in_limb = 63 - (leading as usize);
            return Some(i * 64 + bit_in_limb);
        }
    }
    None
}

/// Schoolbook multiplication of two little-endian limb arrays.
///
/// Returns the product as a `Vec<u64>` of length `a.len() + b.len()`.
/// O(`a.len()` × `b.len()`).
pub(crate) fn multiply_limbs_schoolbook(a: &[u64], b: &[u64]) -> Vec<u64> {
    let mut result = vec![0u64; a.len() + b.len()];
    for i in 0..a.len() {
        let ai = u128::from(a[i]);
        if ai == 0 {
            continue;
        }
        let mut carry: u64 = 0;
        for j in 0..b.len() {
            let bj = u128::from(b[j]);
            let prod = ai * bj + u128::from(result[i + j]) + u128::from(carry);
            result[i + j] = prod as u64;
            carry = (prod >> 64) as u64;
        }
        // Add final carry into the next limb.
        let mut k = i + b.len();
        let mut c = carry;
        while c != 0 && k < result.len() {
            let (sum, overflow) = result[k].overflowing_add(c);
            result[k] = sum;
            c = u64::from(overflow);
            k += 1;
        }
    }
    result
}

/// Karatsuba threshold (operand limb count below which we fall
/// back to schoolbook). Phase 7 will tune empirically; for now we
/// match MPFR's default ballpark.
pub(crate) const KARATSUBA_THRESHOLD: usize = 30;

/// Karatsuba-with-schoolbook-fallback multiplication.
///
/// Recurses with the identity
/// `(a_hi · B + a_lo)(b_hi · B + b_lo) = a_hi·b_hi · B² + ((a_hi+a_lo)(b_hi+b_lo) − a_hi·b_hi − a_lo·b_lo) · B + a_lo·b_lo`,
/// reducing four multiplications to three. Falls back to
/// [`multiply_limbs_schoolbook`] when either operand has at most
/// [`KARATSUBA_THRESHOLD`] limbs.
pub(crate) fn multiply_limbs_karatsuba(a: &[u64], b: &[u64]) -> Vec<u64> {
    if a.len() <= KARATSUBA_THRESHOLD || b.len() <= KARATSUBA_THRESHOLD {
        return multiply_limbs_schoolbook(a, b);
    }
    // Choose split point based on the longer operand.
    let split = a.len().max(b.len()) / 2;

    let (a_lo, a_hi) = split_at_or_zero(a, split);
    let (b_lo, b_hi) = split_at_or_zero(b, split);

    let z0 = multiply_limbs_karatsuba(&a_lo, &b_lo);
    let z2 = multiply_limbs_karatsuba(&a_hi, &b_hi);

    // (a_lo + a_hi) and (b_lo + b_hi).
    let a_sum = add_owned(&a_lo, &a_hi);
    let b_sum = add_owned(&b_lo, &b_hi);
    let z1_full = multiply_limbs_karatsuba(&a_sum, &b_sum);

    // z1 = z1_full - z0 - z2.
    let mut z1 = z1_full;
    sub_owned_in_place(&mut z1, &z0);
    sub_owned_in_place(&mut z1, &z2);

    // Combine: result = z0 + (z1 << (split*64)) + (z2 << (2*split*64)).
    let result_len = a.len() + b.len();
    let mut result = vec![0u64; result_len];
    // z0 at offset 0.
    add_into_at_offset(&mut result, &z0, 0);
    // z1 at offset `split` limbs.
    add_into_at_offset(&mut result, &z1, split);
    // z2 at offset `2 * split` limbs.
    add_into_at_offset(&mut result, &z2, 2 * split);

    result
}

/// Multi-precision multiplication dispatcher: schoolbook for small
/// inputs, Karatsuba for larger ones.
pub(crate) fn multiply_limbs(a: &[u64], b: &[u64]) -> Vec<u64> {
    if a.len().min(b.len()) <= KARATSUBA_THRESHOLD {
        multiply_limbs_schoolbook(a, b)
    } else {
        multiply_limbs_karatsuba(a, b)
    }
}

// --- Karatsuba support helpers ---

/// Split `a` into two halves of `split` limbs (`lo`) and the rest
/// (`hi`). If `a.len() <= split`, the high half is empty (returns
/// an empty Vec).
fn split_at_or_zero(a: &[u64], split: usize) -> (Vec<u64>, Vec<u64>) {
    if a.len() <= split {
        (a.to_vec(), Vec::new())
    } else {
        (a[..split].to_vec(), a[split..].to_vec())
    }
}

/// Owned-output addition: returns `a + b` as a Vec.
fn add_owned(a: &[u64], b: &[u64]) -> Vec<u64> {
    let mut result = vec![0u64; a.len().max(b.len()) + 1];
    for (i, &v) in a.iter().enumerate() {
        result[i] = v;
    }
    // Now add b into result.
    let _carry = limbs_add_assign(&mut result, b);
    // Trim trailing zero (result is up to 1 limb wider than max).
    while result.len() > 1 && *result.last().unwrap() == 0 {
        result.pop();
    }
    result
}

/// In-place subtraction: `dst -= src`. Caller-supplied invariant:
/// `dst >= src`. Trims trailing zeros.
fn sub_owned_in_place(dst: &mut Vec<u64>, src: &[u64]) {
    limbs_sub_assign(dst, src);
    while dst.len() > 1 && *dst.last().unwrap() == 0 {
        dst.pop();
    }
}

/// `dst[offset..] += src`. `dst` is expected to be wide enough to
/// hold the result.
fn add_into_at_offset(dst: &mut [u64], src: &[u64], offset: usize) {
    if offset >= dst.len() {
        return;
    }
    let _carry = limbs_add_assign(&mut dst[offset..], src);
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn top_set_bit_works() {
        assert_eq!(top_set_bit(&[0u64, 0u64]), None);
        assert_eq!(top_set_bit(&[1u64, 0u64]), Some(0));
        assert_eq!(top_set_bit(&[0u64, 1u64]), Some(64));
        assert_eq!(top_set_bit(&[1u64 << 63, 0u64]), Some(63));
        assert_eq!(top_set_bit(&[0u64, 1u64 << 63]), Some(127));
    }

    #[test]
    fn schoolbook_small() {
        // 5 * 7 = 35
        let p = multiply_limbs_schoolbook(&[5], &[7]);
        assert_eq!(p, vec![35, 0]);
    }

    #[test]
    fn schoolbook_carries() {
        // 0xFFFF_FFFF_FFFF_FFFF * 0xFFFF_FFFF_FFFF_FFFF = 0xFFFF_FFFF_FFFF_FFFE_0000_0000_0000_0001
        let p = multiply_limbs_schoolbook(&[u64::MAX], &[u64::MAX]);
        assert_eq!(p, vec![1, 0xFFFF_FFFF_FFFF_FFFEu64]);
    }

    #[test]
    fn schoolbook_multi_limb() {
        // (2^64) * 3 = 3 * 2^64
        let p = multiply_limbs_schoolbook(&[0, 1], &[3]);
        assert_eq!(p, vec![0, 3, 0]);
    }

    #[test]
    fn karatsuba_matches_schoolbook_small() {
        // Should fall through to schoolbook.
        let a = vec![1u64, 2u64];
        let b = vec![3u64, 4u64];
        let p_k = multiply_limbs_karatsuba(&a, &b);
        let p_s = multiply_limbs_schoolbook(&a, &b);
        assert_eq!(p_k, p_s);
    }

    #[test]
    fn karatsuba_matches_schoolbook_large() {
        // Build operands above the threshold to exercise the
        // recursive path. Use deterministic-looking values.
        let n = KARATSUBA_THRESHOLD + 4;
        let mut a = vec![0u64; n];
        let mut b = vec![0u64; n];
        for i in 0..n {
            a[i] = (i as u64).wrapping_mul(0x9E37_79B9_7F4A_7C15) ^ 0xCAFE_BABE;
            b[i] = (i as u64).wrapping_mul(0xC6BC_2796_31E1_4F69) ^ 0xDEAD_BEEF;
        }
        let p_k = multiply_limbs_karatsuba(&a, &b);
        let p_s = multiply_limbs_schoolbook(&a, &b);
        assert_eq!(p_k, p_s, "Karatsuba must agree with schoolbook bit-for-bit");
    }
}
