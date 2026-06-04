//! Sub-range sharding for the exhaustive sweep.
//!
//! The capability pfloat's pf-hcz4 runner lacked: split a single
//! function's `[0, total)` input space into `count` contiguous shards so
//! the 2^32 `f32` grid fans out across instances. All arithmetic is
//! `u64` so the `f32` case (`total = 2^32`) never wraps a `u32` at the
//! final shard.

#![cfg(all(unix, feature = "differential-mpfr"))]

/// The `[start, end)` sub-range of `[0, total)` for shard `index` of
/// `count`. The shards partition `[0, total)` with no gaps or overlaps;
/// shards whose start reaches `total` are empty (`start == end == total`).
pub fn shard_range(index: u64, count: u64, total: u64) -> (u64, u64) {
    let count = count.max(1);
    let per = total.div_ceil(count);
    let start = index.saturating_mul(per).min(total);
    let end = index.saturating_add(1).saturating_mul(per).min(total);
    (start, end)
}

#[cfg(test)]
mod tests {
    use super::shard_range;

    const T: u64 = 1 << 32;

    #[test]
    fn shards_partition_without_gaps_or_overlap() {
        for &count in &[1u64, 2, 3, 7, 16, 257, 65536] {
            let mut prev_end = 0u64;
            for idx in 0..count {
                let (s, e) = shard_range(idx, count, T);
                assert_eq!(s, prev_end, "gap or overlap at shard {idx}/{count}");
                assert!(e >= s, "empty-or-forward range at {idx}/{count}");
                prev_end = e;
            }
            assert_eq!(prev_end, T, "shards must cover [0, 2^32) for count={count}");
        }
    }

    #[test]
    fn last_shard_ends_exactly_at_two_pow_32_no_u32_wrap() {
        let (s, e) = shard_range(65535, 65536, T);
        assert_eq!(e, T, "final shard must end at 2^32");
        assert!(s < e);
        // The final input is 2^32 - 1, which fits u32 without wrapping.
        assert_eq!((e - 1) as u32, u32::MAX);
    }
}
