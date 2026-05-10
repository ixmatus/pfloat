//! Sticky exception flags per IEEE 754-2019 §7.
//!
//! [`Status`] packs the five IEEE 754-2019 sticky exception flags
//! (`INVALID`, `DIV_BY_ZERO`, `OVERFLOW`, `UNDERFLOW`, `INEXACT`)
//! into a single byte. The bit layout matches ferrodec's
//! `src/status.rs::Status` so the differential lane in Phase 6 has a
//! 1:1 translation.
//!
//! Two access patterns coexist:
//!
//! - **Explicit, always available**: every flag-producing operation
//!   returns `(value, Status)`. Callers thread flags through long
//!   expressions via [`Status::merge`] or `|=`. This is the only
//!   form available under `no_std`.
//! - **Thread-local, std-only**: under the `std` feature, every
//!   flag-producing operation also OR-accumulates its `Status` into
//!   a thread-local `Cell<Status>`. The accessor functions in the
//!   [`flags`] submodule (`test`, `clear`, `set`, `raise`) read and
//!   write that cell. Useful for callers who want IEEE 754's
//!   "global" sticky-flag semantics; matches `fenv.h` shape without
//!   the cross-thread footgun.
//!
//! Both patterns return identical bit patterns for the same input,
//! so callers can mix them freely. ADR-0007 records the design.

/// Accumulating bag of IEEE 754-2019 §7 sticky exception flags.
///
/// `Status` is `Copy`, `repr(transparent)` over `u8`, and OR-monoid
/// under [`merge`](Self::merge): combine two statuses by union of
/// their flags.
#[derive(Clone, Copy, Debug, Default, Eq, PartialEq, Hash)]
#[repr(transparent)]
pub struct Status(u8);

impl Status {
    /// No flags set.
    pub const OK: Self = Self(0);

    /// IEEE 754-2019 §7.2 *invalid operation*: the operation has no
    /// useful real-number result and the result is a NaN.
    ///
    /// Raised by signaling-NaN operands in comparisons and
    /// arithmetic, by `0/0`, by `∞ - ∞`, and by `√(negative)`.
    pub const INVALID: Self = Self(0b0000_0001);

    /// IEEE 754-2019 §7.3 *division by zero*: a finite non-zero
    /// dividend divided by zero. The result is the appropriately
    /// signed infinity.
    pub const DIV_BY_ZERO: Self = Self(0b0000_0010);

    /// IEEE 754-2019 §7.4 *overflow*: the rounded result's
    /// magnitude exceeds the largest representable finite value.
    /// Result depends on rounding mode (infinity for round-to-
    /// nearest and round-toward-positive/negative away from zero;
    /// the largest finite for round-toward-zero and the
    /// nearest-toward-zero of the two extremes).
    ///
    /// pfloat at arbitrary precision rarely reaches this flag
    /// because the exponent is `i64`; only catastrophic chains of
    /// scaling operations approach the limit. Slice 1b's pipeline
    /// raises `OVERFLOW` only when an exponent-add saturation
    /// produces `i64::MAX + 1` during round-up.
    pub const OVERFLOW: Self = Self(0b0000_0100);

    /// IEEE 754-2019 §7.5 *underflow*: the rounded result's
    /// magnitude is less than the smallest representable normal
    /// value AND the result is inexact.
    ///
    /// pfloat at arbitrary precision has no implicit minimum
    /// exponent (`i64` covers any realistic working range), so
    /// `UNDERFLOW` is essentially unreachable for normal use.
    /// Defined for IEEE compatibility and for the differential
    /// lane.
    pub const UNDERFLOW: Self = Self(0b0000_1000);

    /// IEEE 754-2019 §7.6 *inexact*: the rounded result differs
    /// from the unrounded result. Set whenever rounding discarded
    /// any non-zero bit.
    pub const INEXACT: Self = Self(0b0001_0000);

    /// Returns `true` when the [`INVALID`](Self::INVALID) flag is
    /// set.
    #[inline]
    #[must_use]
    pub const fn invalid(self) -> bool {
        self.0 & Self::INVALID.0 != 0
    }

    /// Returns `true` when the [`DIV_BY_ZERO`](Self::DIV_BY_ZERO)
    /// flag is set.
    #[inline]
    #[must_use]
    pub const fn div_by_zero(self) -> bool {
        self.0 & Self::DIV_BY_ZERO.0 != 0
    }

    /// Returns `true` when the [`OVERFLOW`](Self::OVERFLOW) flag is
    /// set.
    #[inline]
    #[must_use]
    pub const fn overflow(self) -> bool {
        self.0 & Self::OVERFLOW.0 != 0
    }

    /// Returns `true` when the [`UNDERFLOW`](Self::UNDERFLOW) flag
    /// is set.
    #[inline]
    #[must_use]
    pub const fn underflow(self) -> bool {
        self.0 & Self::UNDERFLOW.0 != 0
    }

    /// Returns `true` when the [`INEXACT`](Self::INEXACT) flag is
    /// set.
    #[inline]
    #[must_use]
    pub const fn inexact(self) -> bool {
        self.0 & Self::INEXACT.0 != 0
    }

    /// Returns `true` when no flags are set.
    #[inline]
    #[must_use]
    pub const fn is_ok(self) -> bool {
        self.0 == 0
    }

    /// Returns the bitwise union of two statuses. Equivalent to
    /// `self | rhs` with [`BitOr`](core::ops::BitOr); provided as a
    /// `const fn` for use in const contexts.
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

/// Thread-local sticky flag accessors (std-only).
///
/// Under `std`, every flag-producing pfloat operation also
/// OR-accumulates its [`Status`] into a thread-local cell. The
/// functions here read and write that cell.
///
/// Use the per-call `(value, Status)` return when local flag
/// accumulation is preferred (e.g., inside a function that should
/// not observe its caller's flag state). Use the thread-local for
/// the IEEE 754 "global state" mental model where flags accumulate
/// across an entire computation.
///
/// Cross-thread visibility: each thread has its own
/// `Cell<Status>`. Spawning a worker thread starts it with
/// [`Status::OK`].
#[cfg(feature = "std")]
pub mod flags {
    use super::Status;
    use std::cell::Cell;

    std::thread_local! {
        static FLAGS: Cell<Status> = const { Cell::new(Status::OK) };
    }

    /// Reads the current thread-local flag set.
    #[must_use]
    pub fn test() -> Status {
        FLAGS.with(Cell::get)
    }

    /// Clears the thread-local flag set, returning its previous
    /// value (so the caller can inspect what was cleared).
    pub fn clear() -> Status {
        FLAGS.with(|f| f.replace(Status::OK))
    }

    /// Overwrites the thread-local flag set with `s`. Use
    /// [`raise`] for OR-accumulation, which is the IEEE 754 sticky
    /// semantics.
    pub fn set(s: Status) {
        FLAGS.with(|f| f.set(s));
    }

    /// OR-accumulates `s` into the thread-local flag set. This is
    /// the operation pfloat's flag-producing ops use internally.
    pub fn raise(s: Status) {
        FLAGS.with(|f| f.set(f.get() | s));
    }
}

/// Internal helper invoked by every flag-producing pfloat
/// operation. Under `std`, OR-accumulates `s` into the thread-local
/// flag set. Under `no_std`, a no-op.
#[allow(dead_code)] // unused under `--no-default-features` (no flag-producers)
#[cfg(feature = "std")]
#[inline]
pub(crate) fn auto_raise(s: Status) {
    flags::raise(s);
}

/// Internal helper invoked by every flag-producing pfloat
/// operation. No-op under `no_std`; the explicit `(value, Status)`
/// return is the only flag transport.
#[allow(dead_code)] // unused under `--no-default-features` (no flag-producers)
#[cfg(not(feature = "std"))]
#[inline]
pub(crate) fn auto_raise(_: Status) {}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn ok_has_no_flags() {
        assert!(Status::OK.is_ok());
        assert!(!Status::OK.invalid());
        assert!(!Status::OK.div_by_zero());
        assert!(!Status::OK.overflow());
        assert!(!Status::OK.underflow());
        assert!(!Status::OK.inexact());
    }

    #[test]
    fn each_flag_is_independent() {
        for (s, name) in [
            (Status::INVALID, "INVALID"),
            (Status::DIV_BY_ZERO, "DIV_BY_ZERO"),
            (Status::OVERFLOW, "OVERFLOW"),
            (Status::UNDERFLOW, "UNDERFLOW"),
            (Status::INEXACT, "INEXACT"),
        ] {
            // Exactly one query method returns true.
            let count = u32::from(s.invalid())
                + u32::from(s.div_by_zero())
                + u32::from(s.overflow())
                + u32::from(s.underflow())
                + u32::from(s.inexact());
            assert_eq!(count, 1, "{name} should set exactly one flag");
        }
    }

    #[test]
    fn merge_is_union() {
        let s = Status::INVALID | Status::INEXACT;
        assert!(s.invalid());
        assert!(s.inexact());
        assert!(!s.overflow());
        let s2 = Status::INVALID.merge(Status::INEXACT);
        assert_eq!(s, s2);
    }

    #[test]
    fn merge_is_idempotent() {
        let s = Status::INEXACT;
        assert_eq!(s.merge(s), s);
    }

    #[test]
    fn bitor_assign_accumulates() {
        let mut s = Status::OK;
        s |= Status::INEXACT;
        assert!(s.inexact());
        s |= Status::OVERFLOW;
        assert!(s.inexact());
        assert!(s.overflow());
    }

    #[test]
    fn default_is_ok() {
        assert_eq!(Status::default(), Status::OK);
    }

    #[cfg(feature = "std")]
    #[test]
    fn thread_local_cleared_at_start() {
        // Each test gets a fresh thread; thread-local starts OK.
        assert_eq!(flags::test(), Status::OK);
    }

    #[cfg(feature = "std")]
    #[test]
    fn thread_local_raise_accumulates() {
        flags::clear();
        flags::raise(Status::INEXACT);
        assert!(flags::test().inexact());
        flags::raise(Status::OVERFLOW);
        let s = flags::test();
        assert!(s.inexact());
        assert!(s.overflow());
    }

    #[cfg(feature = "std")]
    #[test]
    fn thread_local_clear_returns_previous() {
        flags::set(Status::INVALID);
        let prev = flags::clear();
        assert!(prev.invalid());
        assert_eq!(flags::test(), Status::OK);
    }

    #[cfg(feature = "std")]
    #[test]
    fn thread_local_isolated_per_thread() {
        flags::clear();
        flags::set(Status::INVALID);
        let h = std::thread::spawn(|| {
            // Other thread sees OK initially.
            let initial = flags::test();
            flags::raise(Status::OVERFLOW);
            (initial, flags::test())
        });
        let (initial, ours_in_other) = h.join().unwrap();
        assert_eq!(initial, Status::OK);
        assert!(ours_in_other.overflow());
        // Our thread's flags are unaffected.
        assert!(flags::test().invalid());
        assert!(!flags::test().overflow());
        flags::clear();
    }
}
