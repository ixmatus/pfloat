//! Sticky exception flags per IEEE 754-2019 §7.
//!
//! # Slice 1a stub
//!
//! This file is a **stub for slice 1a**. It defines the [`Status`]
//! shape, the IEEE 754-2019 §7.2 `INVALID` flag (the only flag
//! [`partial_cmp`](crate::cmp) can raise without arithmetic), and the
//! `OK` constant. The full enumeration of all five IEEE flags
//! (`INVALID`, `DIV_BY_ZERO`, `OVERFLOW`, `UNDERFLOW`, `INEXACT`),
//! the `BitOr` / `BitOrAssign` plumbing, and the
//! `std`-thread-local-vs-`no_std`-passed-context flag-storage policy
//! land in slice 1b alongside the rounding pipeline.
//!
//! Until slice 1b ships, the only flag-producing API in pfloat is
//! [`partial_cmp`](crate::cmp), which raises `INVALID` on a
//! signaling NaN comparand per IEEE 754-2019 §5.11. ADR-0007 records
//! the full design.

/// Accumulating bag of IEEE 754-2019 §7 sticky exception flags.
///
/// The bit layout matches ferrodec's `src/status.rs::Status` so the
/// differential lane in Phase 6 has a 1:1 translation. Slice 1a only
/// exposes [`INVALID`](Status::INVALID); the remaining four flags
/// (`DIV_BY_ZERO`, `OVERFLOW`, `UNDERFLOW`, `INEXACT`) are reserved
/// bit positions and surface in slice 1b.
#[derive(Clone, Copy, Debug, Default, Eq, PartialEq, Hash)]
#[repr(transparent)]
pub struct Status(u8);

impl Status {
    /// No flags set.
    pub const OK: Self = Self(0);

    /// IEEE 754-2019 §7.2 *invalid operation*: the operation has no
    /// useful real-number result and the result is a NaN.
    ///
    /// In slice 1a, only [`partial_cmp`](crate::cmp) raises this
    /// flag, and only when one of the operands is a signaling NaN
    /// per IEEE 754-2019 §5.11.
    pub const INVALID: Self = Self(0b0000_0001);

    /// Returns `true` when the [`INVALID`](Status::INVALID) flag is
    /// set in this status.
    #[inline]
    #[must_use]
    pub const fn invalid(self) -> bool {
        self.0 & Self::INVALID.0 != 0
    }

    /// Returns `true` when no flags are set.
    #[inline]
    #[must_use]
    pub const fn is_ok(self) -> bool {
        self.0 == 0
    }

    /// Returns the bitwise union of two statuses.
    ///
    /// Equivalent to `self | rhs`. Provided as a `const fn` so it
    /// can be used in const contexts (slice 1b's rounding pipeline
    /// computes flag accumulations at compile time for some
    /// const-precision paths).
    #[inline]
    #[must_use]
    pub const fn merge(self, rhs: Self) -> Self {
        Self(self.0 | rhs.0)
    }
}

impl core::ops::BitOr for Status {
    type Output = Self;

    #[inline]
    fn bitor(self, rhs: Self) -> Self {
        Self(self.0 | rhs.0)
    }
}

impl core::ops::BitOrAssign for Status {
    #[inline]
    fn bitor_assign(&mut self, rhs: Self) {
        self.0 |= rhs.0;
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn ok_has_no_flags() {
        assert!(Status::OK.is_ok());
        assert!(!Status::OK.invalid());
    }

    #[test]
    fn invalid_constructor_sets_invalid() {
        assert!(Status::INVALID.invalid());
        assert!(!Status::INVALID.is_ok());
    }

    #[test]
    fn merge_is_union() {
        let s = Status::OK | Status::INVALID;
        assert!(s.invalid());
        let s2 = Status::OK.merge(Status::INVALID);
        assert_eq!(s, s2);
    }

    #[test]
    fn merge_is_idempotent() {
        let s = Status::INVALID;
        assert_eq!(s.merge(s), s);
    }

    #[test]
    fn bitor_assign_accumulates() {
        let mut s = Status::OK;
        s |= Status::INVALID;
        assert!(s.invalid());
        s |= Status::OK;
        assert!(s.invalid()); // OK is identity
    }

    #[test]
    fn default_is_ok() {
        assert_eq!(Status::default(), Status::OK);
    }
}
