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

/// Read the bit at absolute storage position `pos` from `buf`.
pub(crate) fn bit_at(buf: &[u64], pos: usize) -> bool {
    let limb_idx = pos / 64;
    let bit_in_limb = pos % 64;
    if limb_idx >= buf.len() {
        return false;
    }
    (buf[limb_idx] >> bit_in_limb) & 1 == 1
}

/// Shift `buf` left by one bit in place. Bits shifted off the top
/// are lost (caller is responsible for sizing `buf` to fit any
/// expected growth).
pub(crate) fn shift_left_one_bit(buf: &mut [u64]) {
    let mut carry: u64 = 0;
    for limb in buf.iter_mut() {
        let new_carry = *limb >> 63;
        *limb = (*limb << 1) | carry;
        carry = new_carry;
    }
}

/// Compare two little-endian limb arrays as unsigned integers. The
/// arrays may differ in length; the shorter one is treated as
/// zero-extended at the top.
pub(crate) fn cmp_limbs(a: &[u64], b: &[u64]) -> core::cmp::Ordering {
    use core::cmp::Ordering;
    let len = a.len().max(b.len());
    for i in (0..len).rev() {
        let av = a.get(i).copied().unwrap_or(0);
        let bv = b.get(i).copied().unwrap_or(0);
        let ord = av.cmp(&bv);
        if ord != Ordering::Equal {
            return ord;
        }
    }
    Ordering::Equal
}

/// Integer square root: returns `(s, r)` such that
/// `n = s² + r` and `0 <= r <= 2s`, equivalently
/// `s = floor(sqrt(n))`.
///
/// `s` is sized to roughly `n.bit_length() / 2` bits;
/// `r` is sized to at most `n`'s limb count.
///
/// Algorithm: classic bit-by-bit / pair-at-a-time digit-recurrence.
/// Process the bits of `n` two at a time from MSB toward LSB; at
/// each step, slide two bits into the running remainder, test
/// against `(s << 2) | 1`, and either subtract-and-set-bit or
/// shift-left-zero. Complexity O(`n.bit_length()` × `n.len()`).
/// Phase 7 may replace with Karatsuba sqrt or a Newton-iteration
/// path; correctness ships in slice 1f.
pub(crate) fn isqrt_limbs(n: &[u64]) -> (Vec<u64>, Vec<u64>) {
    let Some(n_top) = top_set_bit(n) else {
        return (vec![0u64], vec![0u64; n.len().max(1)]);
    };

    let s_bits = n_top / 2 + 1;
    let s_limbs = s_bits.div_ceil(64).max(1);
    let r_limbs = n.len() + 2; // headroom for the shift-by-2 inside the loop

    let mut s = vec![0u64; s_limbs];
    let mut r = vec![0u64; r_limbs];

    let max_pair = (n_top + 1).div_ceil(2);

    let mut test = vec![0u64; r_limbs];

    for pair_idx in (0..max_pair).rev() {
        // r <<= 2; OR in the next two bits of `n`.
        shift_left_one_bit(&mut r);
        shift_left_one_bit(&mut r);
        let high = bit_at(n, 2 * pair_idx + 1);
        let low = bit_at(n, 2 * pair_idx);
        r[0] |= (u64::from(high) << 1) | u64::from(low);

        // test = (s << 2) | 1, computed in a scratch buffer the
        // same size as r so the shift cannot overflow.
        test.fill(0);
        for (i, &sv) in s.iter().enumerate() {
            if let Some(slot) = test.get_mut(i) {
                *slot = sv;
            }
        }
        shift_left_one_bit(&mut test);
        shift_left_one_bit(&mut test);
        test[0] |= 1;

        if matches!(
            cmp_limbs(&r, &test),
            core::cmp::Ordering::Greater | core::cmp::Ordering::Equal
        ) {
            limbs_sub_assign(&mut r, &test);
            shift_left_one_bit(&mut s);
            s[0] |= 1;
        } else {
            shift_left_one_bit(&mut s);
        }
    }

    r.truncate(n.len().max(1));
    (s, r)
}

/// Divisor size (in significant limbs) at or above which
/// [`divmod_limbs`] switches from Knuth Algorithm D to the recursive
/// (Burnikel-Ziegler) divider. Below it Knuth's `O(qlen·dlen)` has the
/// smaller constant; at or above it the recursive divider's
/// `O(M(n)·log n)` wins (ADR-0052). The decimal base conversion in
/// `fmt` and high-precision `div` / `parse` are the callers whose
/// divisors reach this size.
const RECURSIVE_DIV_DISPATCH: usize = 64;

/// Divisor-limb count at or below which the recursive divider falls
/// back to Knuth Algorithm D for its base case. The recursive divider
/// normalizes the divisor to a power-of-two limb count, so this is
/// effectively the smallest power-of-two block it recurses down to.
const RECURSIVE_DIV_BASE: usize = 32;

/// Integer `(quotient, remainder)` of `dividend / divisor`.
///
/// `quotient` is returned with `dividend.len()` limbs and `remainder`
/// with `divisor.len()` limbs (each value zero-extended into that
/// width); `0 ≤ remainder < divisor`. `divisor` must be non-zero.
///
/// Dispatches on divisor size: a large divisor (`≥
/// RECURSIVE_DIV_DISPATCH` significant limbs, with a strictly larger
/// dividend) takes the sub-quadratic recursive path
/// ([`divmod_recursive`]); everything else takes the Knuth Algorithm D
/// core ([`divmod_knuth`]). Both honour the identical contract, and the
/// recursive path is differentially checked against the bit-at-a-time
/// oracle (`divmod_limbs_bitwise`) in the tests.
pub(crate) fn divmod_limbs(dividend: &[u64], divisor: &[u64]) -> (Vec<u64>, Vec<u64>) {
    debug_assert!(top_set_bit(divisor).is_some(), "divisor must be non-zero");
    let n = effective_len(divisor);
    let dn = effective_len(dividend);
    if n >= RECURSIVE_DIV_DISPATCH && dn >= n {
        return divmod_recursive(
            &dividend[..dn],
            &divisor[..n],
            dividend.len().max(1),
            divisor.len().max(1),
        );
    }
    divmod_knuth(dividend, divisor)
}

/// Integer `(quotient, remainder)` of `dividend / divisor`, both
/// little-endian limb slices. `quotient` is returned with
/// `dividend.len()` limbs and `remainder` with `divisor.len()` limbs
/// (each value zero-extended into that width).
///
/// Knuth Algorithm D (TAOCP vol. 2, §4.3.1), base `2^64`. The earlier
/// implementation was a bit-at-a-time long division whose cost was
/// `O(dividend_bits × divisor_limbs)`. For the decimal-conversion
/// callers (`fmt::compute_scaled`, the negative-exponent
/// `parse::finite_to_bigfloat`) the dividend and divisor are both
/// near the parse exponent budget — tens of millions of bits — while
/// the quotient is only a handful of limbs, so the bit loop was
/// quadratic in the exponent and let a 10-byte input drive
/// multi-minute, multi-gigabyte work (slice parse-oom; the libFuzzer
/// `parse` out-of-memory). Algorithm D is `O(quotient_limbs ×
/// divisor_limbs)`, linear in the operand size, which bounds that
/// work to the same magnitude the exact `pow5` itself costs. The
/// bit-at-a-time routine survives as a test-only oracle
/// (`divmod_limbs_bitwise`) that this routine is differentially
/// checked against. This is the Knuth core; [`divmod_limbs`] dispatches
/// here for small divisors and to [`divmod_recursive`] for large ones,
/// and the recursive divider uses this as its base case.
fn divmod_knuth(dividend: &[u64], divisor: &[u64]) -> (Vec<u64>, Vec<u64>) {
    debug_assert!(top_set_bit(divisor).is_some(), "divisor must be non-zero");

    let q_len = dividend.len().max(1);
    let r_len = divisor.len().max(1);
    let mut quotient = vec![0u64; q_len];

    // Zero dividend: 0 / d = 0 remainder 0.
    if top_set_bit(dividend).is_none() {
        return (quotient, vec![0u64; r_len]);
    }

    let n = effective_len(divisor); // significant divisor limbs, >= 1
    let dn = effective_len(dividend); // significant dividend limbs, >= 1

    // dividend < divisor (strictly fewer significant limbs): quotient
    // 0, remainder = dividend.
    if dn < n {
        let mut remainder = vec![0u64; r_len];
        remainder[..dividend.len().min(r_len)]
            .copy_from_slice(&dividend[..dividend.len().min(r_len)]);
        return (quotient, remainder);
    }

    // Single-limb divisor: a plain base-2^64 long division, no
    // normalization needed.
    if n == 1 {
        let d = u128::from(divisor[0]);
        let mut rem: u128 = 0;
        for i in (0..dn).rev() {
            let cur = (rem << 64) | u128::from(dividend[i]);
            quotient[i] = (cur / d) as u64;
            rem = cur % d;
        }
        let mut remainder = vec![0u64; r_len];
        remainder[0] = rem as u64;
        return (quotient, remainder);
    }

    // --- Knuth Algorithm D ---
    const B: u128 = 1u128 << 64;
    const MASK: u128 = B - 1;

    // D1. Normalize so the divisor's top limb has its high bit set.
    let shift = divisor[n - 1].leading_zeros();
    let v = shl_limbs(&divisor[..n], shift); // length n (top limb high-bit set)
    let v = &v[..n];
    let mut u = shl_limbs(&dividend[..dn], shift); // length dn + 1
    debug_assert_eq!(u.len(), dn + 1);
    let m = dn - n; // quotient has m + 1 limbs

    let vn1 = u128::from(v[n - 1]);
    let vn2 = u128::from(v[n - 2]);

    // D2-D7. Main loop, one quotient limb per iteration, top down.
    for j in (0..=m).rev() {
        // D3. Estimate qhat = floor((u[j+n]·B + u[j+n-1]) / v[n-1]).
        let num = (u128::from(u[j + n]) << 64) | u128::from(u[j + n - 1]);
        let mut qhat = num / vn1;
        let mut rhat = num % vn1;
        // The invariant u[j+n] <= v[n-1] bounds qhat <= B, so the
        // `qhat >= B` test trips at most twice; `||` short-circuits
        // before `qhat * vn2` can overflow u128.
        while qhat >= B || qhat * vn2 > (rhat << 64) | u128::from(u[j + n - 2]) {
            qhat -= 1;
            rhat += vn1;
            if rhat >= B {
                break;
            }
        }

        // D4. Multiply and subtract: u[j..=j+n] -= qhat · v.
        let mut borrow: i128 = 0;
        for i in 0..n {
            let p = qhat * u128::from(v[i]);
            let t = i128::from(u[j + i]) - borrow - (p & MASK) as i128;
            u[j + i] = t as u64;
            borrow = (p >> 64) as i128 - (t >> 64);
        }
        let t = i128::from(u[j + n]) - borrow;
        u[j + n] = t as u64;

        // D5/D6. If we subtracted too much, qhat was one too large:
        // add the divisor back and decrement qhat.
        if t < 0 {
            qhat -= 1;
            let mut carry: u128 = 0;
            for i in 0..n {
                let sum = u128::from(u[j + i]) + u128::from(v[i]) + carry;
                u[j + i] = sum as u64;
                carry = sum >> 64;
            }
            u[j + n] = (u128::from(u[j + n]) + carry) as u64; // final carry discarded
        }

        quotient[j] = qhat as u64;
    }

    // D8. Denormalize the remainder (u[0..n] >> shift).
    let r_norm = shr_limbs(&u[..n], shift);
    let mut remainder = vec![0u64; r_len];
    let copy = r_norm.len().min(r_len);
    remainder[..copy].copy_from_slice(&r_norm[..copy]);
    (quotient, remainder)
}

/// Recursive (Burnikel-Ziegler) integer division, `O(M(n)·log n)`.
///
/// `dividend` and `divisor` are significant-length slices (no high zero
/// limbs) with `dividend.len() ≥ divisor.len() ≥ 1`. The result is
/// zero-extended to `q_len` / `r_len` limbs to honour [`divmod_limbs`]'s
/// contract.
///
/// The divisor is normalized to a power-of-two limb count `m` with its
/// top bit set (a combined whole-limb + sub-limb left shift `s`); the
/// same shift applied to the dividend leaves the quotient unchanged and
/// scales the remainder by `2^s`, which the final right shift undoes.
/// The dividend is then processed in `m`-limb blocks from the top, each
/// `[remainder | block]` divided by the normalized divisor via
/// [`div_2n_1n`]. Because `m` is a power of two, every recursion halves
/// to an even size down to [`RECURSIVE_DIV_BASE`].
fn divmod_recursive(
    dividend: &[u64],
    divisor: &[u64],
    q_len: usize,
    r_len: usize,
) -> (Vec<u64>, Vec<u64>) {
    let n = divisor.len();
    let m = n.next_power_of_two();

    // Normalize: shift the divisor up to exactly `m` limbs, top bit set.
    let clz = divisor[n - 1].leading_zeros() as usize;
    let s = (m - n) * 64 + clz;
    let b_full = shift_left_bits(divisor, s);
    let mut b = vec![0u64; m];
    let copy = b_full.len().min(m);
    b[..copy].copy_from_slice(&b_full[..copy]);
    debug_assert_eq!(effective_len(&b), m, "normalized divisor must fill m limbs");

    // Apply the same shift to the dividend and pad to a multiple of `m`.
    let a = shift_left_bits(dividend, s);
    let t = a.len().div_ceil(m).max(1);
    let mut a_pad = vec![0u64; t * m];
    a_pad[..a.len()].copy_from_slice(&a);

    let mut quotient = vec![0u64; t * m];
    let mut r = vec![0u64; m]; // running remainder, always < b
    for i in (0..t).rev() {
        let zi = &a_pad[i * m..(i + 1) * m];
        let mut dd = Vec::with_capacity(2 * m);
        dd.extend_from_slice(zi); // low m limbs
        dd.extend_from_slice(&r); // high m limbs (r·β^m)
        let (qi, ri) = div_2n_1n(&dd, &b);
        let qc = qi.len().min(m);
        quotient[i * m..i * m + qc].copy_from_slice(&qi[..qc]);
        r = ri;
    }

    // Denormalize the remainder; the quotient is shift-invariant.
    let rem = shift_right_bits(&r, s);

    let mut q_out = vec![0u64; q_len];
    let qc = quotient.len().min(q_len);
    debug_assert!(
        quotient[qc..].iter().all(|&l| l == 0),
        "quotient exceeds dividend width"
    );
    q_out[..qc].copy_from_slice(&quotient[..qc]);
    let mut r_out = vec![0u64; r_len];
    let rc = rem.len().min(r_len);
    r_out[..rc].copy_from_slice(&rem[..rc]);
    (q_out, r_out)
}

/// Divide a `2n`-limb dividend by an `n`-limb normalized divisor (top
/// bit set), with the precondition `a < β^n · b` so the quotient fits
/// `n` limbs. Returns `(Q, R)` with `Q` exactly `n` limbs and `R`
/// exactly `n` limbs, `a = Q·b + R`, `0 ≤ R < b`.
fn div_2n_1n(a: &[u64], b: &[u64]) -> (Vec<u64>, Vec<u64>) {
    let n = b.len();
    debug_assert_eq!(a.len(), 2 * n);
    if n <= RECURSIVE_DIV_BASE {
        let (mut q, r) = divmod_knuth(a, b);
        q.truncate(n); // value < β^n
        return (q, r);
    }
    let half = n / 2;
    // a = [A4 | A3 | A2 | A1] in `half`-limb quarters, A1 most significant.
    let top3 = &a[half..2 * n]; // [A3 | A2 | A1], 3·half limbs
    let (q1, r1) = div_3n_2n(top3, b);
    let mut dd2 = Vec::with_capacity(3 * half);
    dd2.extend_from_slice(&a[0..half]); // A4 (low half)
    dd2.extend_from_slice(&r1); // r1 (n limbs)
    let (q0, r) = div_3n_2n(&dd2, b);
    let mut q = vec![0u64; n];
    let c0 = q0.len().min(half);
    q[..c0].copy_from_slice(&q0[..c0]);
    let c1 = q1.len().min(half);
    q[half..half + c1].copy_from_slice(&q1[..c1]);
    (q, r)
}

/// Divide a `3·half`-limb dividend by an `n = 2·half`-limb normalized
/// divisor, with the precondition `a < β^half · b` so the quotient fits
/// `half` limbs. Returns `(Q, R)` with `Q` exactly `half` limbs and `R`
/// exactly `n` limbs, `a = Q·b + R`, `0 ≤ R < b`. The Burnikel-Ziegler
/// "three halves by two" step.
fn div_3n_2n(a: &[u64], b: &[u64]) -> (Vec<u64>, Vec<u64>) {
    use core::cmp::Ordering;
    let n = b.len();
    let half = n / 2;
    debug_assert_eq!(a.len(), 3 * half);
    let a3 = &a[0..half];
    let a12 = &a[half..3 * half]; // n limbs
    let a1 = &a[2 * half..3 * half]; // top half
    let b1 = &b[half..n]; // high half of b (inherits b's top-bit-set)
    let b2 = &b[0..half]; // low half of b

    let mut q: Vec<u64>;
    let r1: Vec<u64>;
    if cmp_limbs(a1, b1) == Ordering::Less {
        // Quotient estimate fits: divide the top n limbs by b1.
        let (qq, rr) = div_2n_1n(a12, b1);
        q = qq;
        r1 = rr;
    } else {
        // Estimate would overflow; clamp to β^half − 1 and form the
        // partial remainder r1 = A12 − (β^half − 1)·b1 = A12 + b1 − b1·β^half.
        q = vec![u64::MAX; half];
        let mut tmp = vec![0u64; n + 1];
        tmp[..a12.len()].copy_from_slice(a12);
        let carry = limbs_add_assign(&mut tmp, b1);
        debug_assert!(!carry, "A12 + b1 overflow");
        let b1_shifted = shift_up_limbs(b1, half); // b1·β^half, n limbs
        limbs_sub_assign(&mut tmp, &b1_shifted); // valid: A1 ≥ b1 ⟹ A12 ≥ b1·β^half
        r1 = tmp;
    }

    // D = Q·b2; candidate C = r1·β^half + A3; then R = C − D with ≤ 2
    // corrections (decrement Q, add b back) when the estimate ran high.
    let d = multiply_limbs(&q, b2);
    let mut c = vec![0u64; 3 * half + 2];
    c[half..half + r1.len()].copy_from_slice(&r1);
    let carry = limbs_add_assign(&mut c, a3);
    debug_assert!(!carry, "C low-half add overflow");

    let mut corrections = 0;
    while cmp_limbs(&c, &d) == Ordering::Less {
        sub_one(&mut q);
        let carry = limbs_add_assign(&mut c, b);
        debug_assert!(!carry, "correction add overflow");
        corrections += 1;
        debug_assert!(corrections <= 2, "more than two BZ corrections");
    }
    limbs_sub_assign(&mut c, &d);
    c.truncate(n);
    debug_assert!(cmp_limbs(&c, b) == Ordering::Less, "remainder >= divisor");

    q.resize(half, 0); // value < β^half
    (q, c)
}

/// `x · β^k`: prepend `k` zero limbs (a whole-limb left shift).
fn shift_up_limbs(x: &[u64], k: usize) -> Vec<u64> {
    let mut out = vec![0u64; k + x.len()];
    out[k..].copy_from_slice(x);
    out
}

/// In-place `x -= 1`; `x` must be non-zero.
fn sub_one(x: &mut [u64]) {
    for limb in x.iter_mut() {
        let (v, borrow) = limb.overflowing_sub(1);
        *limb = v;
        if !borrow {
            return;
        }
    }
    debug_assert!(false, "sub_one underflow");
}

/// Logical left shift by `s` bits (any `s`), combining a whole-limb
/// shift with a sub-limb [`shl_limbs`].
fn shift_left_bits(x: &[u64], s: usize) -> Vec<u64> {
    let limb_shift = s / 64;
    let bit_shift = (s % 64) as u32;
    let body = if bit_shift == 0 {
        let mut v = x.to_vec();
        v.push(0);
        v
    } else {
        shl_limbs(x, bit_shift)
    };
    let mut out = vec![0u64; limb_shift];
    out.extend_from_slice(&body);
    out
}

/// Logical right shift by `s` bits (any `s`), combining a whole-limb
/// drop with a sub-limb [`shr_limbs`].
fn shift_right_bits(x: &[u64], s: usize) -> Vec<u64> {
    let limb_shift = s / 64;
    let bit_shift = (s % 64) as u32;
    if limb_shift >= x.len() {
        return vec![0u64; 1];
    }
    let hi = &x[limb_shift..];
    if bit_shift == 0 {
        hi.to_vec()
    } else {
        shr_limbs(hi, bit_shift)
    }
}

/// Significant little-endian limb count (high zero limbs stripped),
/// always at least 1.
fn effective_len(x: &[u64]) -> usize {
    let mut n = x.len();
    while n > 1 && x[n - 1] == 0 {
        n -= 1;
    }
    n.max(1)
}

/// Logical left shift of a little-endian limb slice by `s` bits
/// (`0 ≤ s < 64`). Returns `src.len() + 1` limbs so the shifted-out
/// high bits are never lost.
fn shl_limbs(src: &[u64], s: u32) -> Vec<u64> {
    let mut out = vec![0u64; src.len() + 1];
    if s == 0 {
        out[..src.len()].copy_from_slice(src);
        return out;
    }
    let mut carry: u64 = 0;
    for (i, &x) in src.iter().enumerate() {
        out[i] = (x << s) | carry;
        carry = x >> (64 - s);
    }
    out[src.len()] = carry;
    out
}

/// Logical right shift of a little-endian limb slice by `s` bits
/// (`0 ≤ s < 64`), returning `src.len()` limbs.
fn shr_limbs(src: &[u64], s: u32) -> Vec<u64> {
    if s == 0 {
        return src.to_vec();
    }
    let n = src.len();
    let mut out = vec![0u64; n];
    for i in 0..n {
        let lo = src[i] >> s;
        let hi = if i + 1 < n { src[i + 1] << (64 - s) } else { 0 };
        out[i] = lo | hi;
    }
    out
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

/// Karatsuba threshold: operand limb count at or below which the
/// dispatcher (and the recursion base case) uses schoolbook.
///
/// Calibrated empirically in slice 7d against
/// `benches/mul_thresholds.rs` on the equal-size sweep; ADR-0027
/// records the methodology and the curves. The MPFR-ballpark value
/// of 30 was too low: this in-tree Karatsuba allocates several
/// `Vec`s per recursion node, so its constant factor only pays off
/// past ~48 limbs. At 30, multiplications in the 32..48-limb band
/// were routed to Karatsuba and ran ~20% slower than schoolbook;
/// 48 keeps that band on schoolbook and still wins large-n via a
/// better recursion base case (n=512 ~15% faster than at 30). The
/// threshold is host- and arch-dependent (calibrated on
/// aarch64-apple-darwin); 48 is a deliberate, measured point, not
/// an asymptotic guess.
pub(crate) const KARATSUBA_THRESHOLD: usize = 48;

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

    /// The pre-Algorithm-D bit-at-a-time long division, retained as a
    /// differential oracle for [`divmod_limbs`]. Simple and obviously
    /// correct (it is the schoolbook restoring algorithm), at the cost
    /// of being quadratic in the dividend bit length — which is exactly
    /// why the production routine moved to Algorithm D.
    fn divmod_limbs_bitwise(dividend: &[u64], divisor: &[u64]) -> (Vec<u64>, Vec<u64>) {
        let mut quotient = vec![0u64; dividend.len()];
        let mut remainder = vec![0u64; divisor.len() + 1];
        let dividend_top = match top_set_bit(dividend) {
            Some(t) => t,
            None => return (quotient, vec![0u64; divisor.len()]),
        };
        for bit_idx in (0..=dividend_top).rev() {
            shift_left_one_bit(&mut remainder);
            if bit_at(dividend, bit_idx) {
                remainder[0] |= 1;
            }
            if cmp_limbs(&remainder, divisor) != core::cmp::Ordering::Less {
                limbs_sub_assign(&mut remainder, divisor);
                let q_idx = bit_idx / 64;
                if q_idx < quotient.len() {
                    quotient[q_idx] |= 1u64 << (bit_idx % 64);
                }
            }
        }
        remainder.truncate(divisor.len());
        (quotient, remainder)
    }

    // Deterministic xorshift64 for the differential sweeps.
    fn xorshift(state: &mut u64) -> u64 {
        *state ^= *state << 13;
        *state ^= *state >> 7;
        *state ^= *state << 17;
        *state
    }

    #[test]
    fn divmod_matches_bitwise_reference() {
        let mut state: u64 = 0x9E37_79B9_7F4A_7C15;
        for _ in 0..5000 {
            let dlen = 1 + (xorshift(&mut state) % 8) as usize;
            let vlen = 1 + (xorshift(&mut state) % 8) as usize;
            let dividend: Vec<u64> = (0..dlen).map(|_| xorshift(&mut state)).collect();
            let mut divisor: Vec<u64> = (0..vlen).map(|_| xorshift(&mut state)).collect();
            // Sometimes zero high divisor limbs (trailing-zero form) to
            // exercise effective-length handling and the single-limb path.
            if xorshift(&mut state) & 1 == 0 && vlen > 1 {
                let keep = 1 + (xorshift(&mut state) % vlen as u64) as usize;
                for limb in divisor.iter_mut().skip(keep) {
                    *limb = 0;
                }
            }
            if divisor.iter().all(|&x| x == 0) {
                divisor[0] = 1;
            }
            let (q1, r1) = divmod_limbs(&dividend, &divisor);
            let (q2, r2) = divmod_limbs_bitwise(&dividend, &divisor);
            assert_eq!(q1, q2, "quotient: {dividend:?} / {divisor:?}");
            assert_eq!(r1, r2, "remainder: {dividend:?} / {divisor:?}");
        }
    }

    #[test]
    fn divmod_reconstructs_large_nearly_equal() {
        // The decimal-conversion shape: huge divisor, small quotient.
        // Verify q·v + r == dividend and r < v at multi-limb scale,
        // where the bitwise oracle would be too slow to use directly.
        let mut state: u64 = 0xD1B5_4A32_D192_ED03;
        for _ in 0..300 {
            let n = 2 + (xorshift(&mut state) % 40) as usize;
            let mut divisor: Vec<u64> = (0..n).map(|_| xorshift(&mut state)).collect();
            *divisor.last_mut().unwrap() |= 1u64 << 63; // top limb non-zero
            let qfac = 1 + (xorshift(&mut state) % 1_000_000);
            let add = xorshift(&mut state) % divisor[0].max(1);
            let mut dividend = multiply_limbs(&divisor, &[qfac]);
            let _ = limbs_add_assign(&mut dividend, &[add]);
            let (q, r) = divmod_limbs(&dividend, &divisor);
            assert_eq!(
                cmp_limbs(&r, &divisor),
                core::cmp::Ordering::Less,
                "remainder >= divisor"
            );
            let mut recon = multiply_limbs(&q, &divisor);
            let _ = limbs_add_assign(&mut recon, &r);
            assert_eq!(
                cmp_limbs(&recon, &dividend),
                core::cmp::Ordering::Equal,
                "q·v + r != dividend"
            );
        }
    }

    #[test]
    fn divmod_recursive_matches_bitwise_reference() {
        // Force the recursive (Burnikel-Ziegler) path: divisors at and
        // above the dispatch threshold, spanning the power-of-two
        // normalization boundary and several recursion levels. The
        // bit-at-a-time routine is the independent ground truth.
        let mut state: u64 = 0x0123_4567_89AB_CDEF;
        for _ in 0..50 {
            let vlen = 64 + (xorshift(&mut state) % 65) as usize; // 64..=128
            let mut divisor: Vec<u64> = (0..vlen).map(|_| xorshift(&mut state)).collect();
            *divisor.last_mut().unwrap() |= 1u64 << 63;
            let extra = (xorshift(&mut state) % (vlen as u64 + 1)) as usize;
            let dlen = vlen + extra;
            let mut dividend: Vec<u64> = (0..dlen).map(|_| xorshift(&mut state)).collect();
            *dividend.last_mut().unwrap() |= 1; // dn == dlen ≥ vlen → recursive path
            let (q1, r1) = divmod_limbs(&dividend, &divisor);
            let (q2, r2) = divmod_limbs_bitwise(&dividend, &divisor);
            assert_eq!(q1, q2, "quotient mismatch (vlen={vlen}, dlen={dlen})");
            assert_eq!(r1, r2, "remainder mismatch (vlen={vlen}, dlen={dlen})");
        }
    }

    #[test]
    fn divmod_recursive_reconstructs() {
        // Cheap structural sweep at recursive sizes: construct
        // dividend = q·divisor + r with r < divisor by construction
        // (at most n−1 limbs, and divisor has its top bit set so
        // divisor ≥ β^(n−1)), then check exact recovery.
        use core::cmp::Ordering;
        let mut state: u64 = 0x5DEE_CE66_2A9E_1B07;
        for _ in 0..200 {
            let n = 64 + (xorshift(&mut state) % 201) as usize; // 64..=264
            let mut divisor: Vec<u64> = (0..n).map(|_| xorshift(&mut state)).collect();
            *divisor.last_mut().unwrap() |= 1u64 << 63;
            let rlen = (xorshift(&mut state) % n as u64) as usize; // 0..=n−1
            let r: Vec<u64> = (0..rlen).map(|_| xorshift(&mut state)).collect();
            let qlen = 1 + (xorshift(&mut state) % n as u64) as usize; // 1..=n
            let mut q: Vec<u64> = (0..qlen).map(|_| xorshift(&mut state)).collect();
            *q.last_mut().unwrap() |= 1; // nonzero quotient → dn ≥ n → recursive path
            let mut dividend = multiply_limbs(&q, &divisor);
            let _ = limbs_add_assign(&mut dividend, &r);
            let (qq, rr) = divmod_limbs(&dividend, &divisor);
            assert_eq!(
                cmp_limbs(&rr, &divisor),
                Ordering::Less,
                "remainder >= divisor (n={n})"
            );
            assert_eq!(cmp_limbs(&rr, &r), Ordering::Equal, "remainder mismatch (n={n})");
            assert_eq!(cmp_limbs(&qq, &q), Ordering::Equal, "quotient mismatch (n={n})");
            let mut recon = multiply_limbs(&qq, &divisor);
            let _ = limbs_add_assign(&mut recon, &rr);
            assert_eq!(
                cmp_limbs(&recon, &dividend),
                Ordering::Equal,
                "q·v + r != dividend (n={n})"
            );
        }
    }

    #[test]
    fn divmod_recursive_adversarial() {
        // Structured divisors that stress the Burnikel-Ziegler
        // correction branches at recursive sizes (the all-ones divisor
        // in particular drives quotient estimates high), checked against
        // the bitwise oracle.
        let n = 80usize;
        let cases: Vec<Vec<u64>> = vec![
            vec![u64::MAX; n], // all ones (max divisor)
            {
                let mut v = vec![0u64; n];
                v[n - 1] = 1 << 63;
                v
            }, // single top bit: 2^(64n−1)
            {
                let mut v = vec![0xAAAA_AAAA_AAAA_AAAAu64; n];
                v[n - 1] |= 1 << 63;
                v
            }, // alternating bits
            {
                let mut v = vec![1u64; n];
                v[n - 1] |= 1 << 63;
                v
            }, // sparse low bits
        ];
        let mut state: u64 = 0xF00D_BABE_1234_5678;
        for divisor in &cases {
            for _ in 0..6 {
                let dlen = n + (xorshift(&mut state) % (n as u64 + 1)) as usize;
                let mut dividend: Vec<u64> = (0..dlen).map(|_| xorshift(&mut state)).collect();
                *dividend.last_mut().unwrap() |= 1;
                let (q1, r1) = divmod_limbs(&dividend, divisor);
                let (q2, r2) = divmod_limbs_bitwise(&dividend, divisor);
                assert_eq!(q1, q2, "quotient mismatch (structured divisor)");
                assert_eq!(r1, r2, "remainder mismatch (structured divisor)");
            }
        }
    }

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

    #[test]
    fn bit_at_works() {
        let buf = [0b1010u64, 1u64 << 63];
        assert!(!bit_at(&buf, 0));
        assert!(bit_at(&buf, 1));
        assert!(!bit_at(&buf, 2));
        assert!(bit_at(&buf, 3));
        assert!(!bit_at(&buf, 64)); // bit 0 of limb 1
        assert!(bit_at(&buf, 127)); // top bit of limb 1
        assert!(!bit_at(&buf, 200)); // out of range
    }

    #[test]
    fn shift_left_one_bit_works() {
        let mut buf = [0b1010u64, 0u64];
        shift_left_one_bit(&mut buf);
        assert_eq!(buf, [0b10100u64, 0u64]);

        let mut buf = [1u64 << 63, 0u64];
        shift_left_one_bit(&mut buf);
        // Bit 63 of limb 0 carries into bit 0 of limb 1.
        assert_eq!(buf, [0u64, 1u64]);
    }

    #[test]
    fn cmp_limbs_works() {
        use core::cmp::Ordering;
        assert_eq!(cmp_limbs(&[0u64], &[0u64]), Ordering::Equal);
        assert_eq!(cmp_limbs(&[1u64], &[0u64]), Ordering::Greater);
        assert_eq!(cmp_limbs(&[0u64], &[1u64]), Ordering::Less);
        // High limb dominates.
        assert_eq!(
            cmp_limbs(&[0u64, 1u64], &[u64::MAX, 0u64]),
            Ordering::Greater
        );
        // Different lengths.
        assert_eq!(cmp_limbs(&[5u64], &[5u64, 0u64]), Ordering::Equal);
        assert_eq!(cmp_limbs(&[5u64], &[5u64, 1u64]), Ordering::Less);
    }

    #[test]
    fn divmod_small_exact() {
        // 35 / 5 = 7 remainder 0.
        let (q, r) = divmod_limbs(&[35u64], &[5u64]);
        assert_eq!(q[0], 7);
        assert_eq!(r[0], 0);
    }

    #[test]
    fn divmod_small_with_remainder() {
        // 23 / 5 = 4 remainder 3.
        let (q, r) = divmod_limbs(&[23u64], &[5u64]);
        assert_eq!(q[0], 4);
        assert_eq!(r[0], 3);
    }

    #[test]
    fn divmod_dividend_smaller_than_divisor() {
        // 3 / 5 = 0 remainder 3.
        let (q, r) = divmod_limbs(&[3u64], &[5u64]);
        assert_eq!(q[0], 0);
        assert_eq!(r[0], 3);
    }

    #[test]
    fn divmod_multi_limb_dividend() {
        // 2^64 / 3 = 6148914691236517205 (= 0x5555_5555_5555_5555) remainder 1.
        let (q, r) = divmod_limbs(&[0u64, 1u64], &[3u64]);
        assert_eq!(q, vec![0x5555_5555_5555_5555u64, 0u64]);
        assert_eq!(r[0], 1);
    }

    #[test]
    fn divmod_zero_dividend() {
        // 0 / anything = 0 r 0.
        let (q, r) = divmod_limbs(&[0u64, 0u64], &[5u64]);
        assert_eq!(q, vec![0u64, 0u64]);
        assert_eq!(r, vec![0u64]);
    }

    #[test]
    fn isqrt_perfect_squares() {
        for &(input, expected) in &[
            (0u64, 0u64),
            (1, 1),
            (4, 2),
            (9, 3),
            (16, 4),
            (25, 5),
            (100, 10),
            (10000, 100),
            (2147483648u64 * 2147483648u64, 2147483648), // exact 2^62
        ] {
            let (s, r) = isqrt_limbs(&[input]);
            assert_eq!(s[0], expected, "isqrt({input})");
            // Verify s² <= input < (s+1)².
            let s_squared = s[0].checked_mul(s[0]).expect("no overflow at this scale");
            assert!(s_squared <= input);
            if let Some(next) = (s[0] + 1).checked_mul(s[0] + 1) {
                assert!(input < next);
            }
            assert_eq!(r[0], input - s_squared);
        }
    }

    #[test]
    fn isqrt_non_perfect_squares() {
        for &input in &[2u64, 3, 5, 7, 10, 26, 999] {
            let (s, r) = isqrt_limbs(&[input]);
            // s² <= input < (s+1)²
            let s_val = s[0];
            assert!(s_val * s_val <= input, "isqrt({input}): s² > n");
            assert!(
                (s_val + 1) * (s_val + 1) > input,
                "isqrt({input}): (s+1)² <= n"
            );
            assert_eq!(r[0], input - s_val * s_val);
        }
    }

    #[test]
    fn isqrt_multi_limb_perfect_square() {
        // (2^64) squared = 2^128, occupying limb 2.
        // sqrt should be 2^64 = [0, 1].
        let n = [0u64, 0u64, 1u64];
        let (s, r) = isqrt_limbs(&n);
        assert_eq!(s, vec![0u64, 1u64]);
        assert!(r.iter().all(|&v| v == 0));
    }

    #[test]
    fn isqrt_multi_limb_near_square() {
        // Pick a known perfect square: 10^18 squared.
        let base: u64 = 1_000_000_000_000_000_000;
        let squared = multiply_limbs_schoolbook(&[base], &[base]);
        let (s, r) = isqrt_limbs(&squared);
        assert_eq!(s[0], base, "isqrt of (10^18)² should be 10^18");
        assert!(
            r.iter().all(|&v| v == 0),
            "remainder should be zero for perfect square"
        );
    }

    #[test]
    fn divmod_random_round_trip() {
        // For random a and non-zero b, verify a == q*b + r, r < b.
        use core::cmp::Ordering;
        let a: Vec<u64> = vec![0xDEAD_BEEF_CAFE_BABE, 0x1234_5678_9ABC_DEF0];
        let b: Vec<u64> = vec![0x0000_0000_0001_2345];
        let (q, r) = divmod_limbs(&a, &b);
        // Verify q × b + r == a.
        let qb = multiply_limbs_schoolbook(&q, &b);
        let mut reconstructed = qb;
        // Pad r to match.
        while reconstructed.len() < r.len() {
            reconstructed.push(0);
        }
        let mut r_padded = r.clone();
        while r_padded.len() < reconstructed.len() {
            r_padded.push(0);
        }
        let _ = limbs_add_assign(&mut reconstructed, &r_padded);
        // Trim leading zeros to match a.
        while reconstructed.len() > a.len() {
            assert_eq!(*reconstructed.last().unwrap(), 0);
            reconstructed.pop();
        }
        assert_eq!(reconstructed, a);
        // Verify r < b.
        assert_eq!(cmp_limbs(&r, &b), Ordering::Less);
    }
}
