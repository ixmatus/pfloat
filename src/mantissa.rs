//! Mantissa storage abstractions shared by `BigFloat` and (in slice
//! 1g) `FixedFloat<PREC>`.
//!
//! The trait is `pub(crate)`. It is the internal shape that
//! arithmetic kernels in slices 1c onward read through; users never
//! see it. The two implementations live with the types they back:
//! `BigFloat` carries a `Vec<u64>` mantissa, `FixedFloat<PREC>`
//! carries a `[u64; limbs_for(PREC)]` mantissa, both little-endian
//! limb order with the most-significant bit of the most-significant
//! limb set for normalized non-zero values per ADR-0001 and ADR-0002.

/// Number of `u64` limbs needed to hold a `precision`-bit mantissa.
///
/// `((precision + 63) / 64)`. Returns `0` if `precision == 0`, which
/// is invalid for a normalized non-zero value; the public
/// constructors validate `precision >= 1` before calling this.
#[allow(dead_code)] // unused under `--no-default-features`; consumed by `big` and (slice 1g) `fixed`
#[inline]
#[must_use]
pub const fn limbs_for(precision: u32) -> usize {
    (precision as usize).div_ceil(64)
}

/// Internal abstraction over mantissa storage.
///
/// Slice 1a defines the trait; slice 1c (Add/Sub) is the first
/// kernel to consume it. Methods land here as kernels need them. The
/// initial surface is read-only because slice 1a's classification
/// and comparison only need to walk the limbs.
#[allow(dead_code)] // consumed by slice 1c onward
pub(crate) trait Mantissa {
    /// Limb storage, little-endian (limb 0 is least significant).
    fn as_limbs(&self) -> &[u64];

    /// Number of meaningful precision bits in the mantissa. The
    /// limb storage rounds up to whole `u64`s; bits below the
    /// precision are zero in canonical form.
    fn precision_bits(&self) -> u32;

    /// `true` when the most-significant bit of the most-significant
    /// limb is `1`, i.e. the mantissa is canonically normalized for
    /// a non-zero value (ADR-0001's top-bit-set rule).
    fn is_top_bit_set(&self) -> bool;
}

/// Computes the bit-shift needed to place a `value_bits`-bit
/// magnitude as the top `value_bits` bits of a `precision`-bit
/// mantissa stored in `limbs` `u64`s.
///
/// In other words: if the magnitude is laid out in the storage so
/// that its most-significant bit lands at the very top of the
/// most-significant limb (bit 63 of `mantissa[limbs - 1]`), this
/// function returns `(limbs * 64) - value_bits`.
///
/// Used by [`crate::convert::storage_layout_from_magnitude`] in
/// slice 1a's exact integer constructors and by the rounding
/// pipeline in slice 1b.
#[allow(dead_code)] // consumed by slice 1a's `convert.rs`
#[inline]
#[must_use]
pub(crate) const fn storage_shift(limbs: usize, value_bits: u32) -> u32 {
    (limbs as u32) * 64 - value_bits
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn limbs_for_round_up_to_64() {
        assert_eq!(limbs_for(1), 1);
        assert_eq!(limbs_for(53), 1);
        assert_eq!(limbs_for(64), 1);
        assert_eq!(limbs_for(65), 2);
        assert_eq!(limbs_for(113), 2);
        assert_eq!(limbs_for(128), 2);
        assert_eq!(limbs_for(129), 3);
        assert_eq!(limbs_for(256), 4);
        assert_eq!(limbs_for(257), 5);
    }

    #[test]
    fn storage_shift_places_top_bit_at_msb() {
        // For precision=53 in a 1-limb storage, a 1-bit value
        // (e.g. magnitude 1) needs 63 bits of left-shift to reach
        // bit 63 of the MSL.
        assert_eq!(storage_shift(1, 1), 63);
        // For precision=53 in a 1-limb storage, a 53-bit value
        // (e.g. magnitude with top bit at position 52 already)
        // needs 11 bits of left-shift.
        assert_eq!(storage_shift(1, 53), 11);
        // For precision=128 in a 2-limb storage, a 1-bit value
        // needs 127 bits of left-shift.
        assert_eq!(storage_shift(2, 1), 127);
        // For precision=64 in a 1-limb storage, a 64-bit value
        // needs 0 bits.
        assert_eq!(storage_shift(1, 64), 0);
    }
}
