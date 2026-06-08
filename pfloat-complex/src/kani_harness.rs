//! Kani proof harnesses for pfloat-complex (ADR-0092, advisory).
//!
//! Compiled only under `cfg(kani)`. The verifiable surface here is deliberately
//! narrow: every `Complex<T>` operation runs in `BigFloat`, which is
//! `Vec`-backed and unverifiable at the heap-allocation level (CBMC is
//! `Vec`-hostile, ADR-0062), so the elementary kernels' soundness rests on the
//! enumerated, identity, and acb differential lanes instead. The one `Vec`-free
//! invariant the public API depends on is the componentwise [`Status`] merge:
//! every `add`/`sub`/`mul`/`div`/`abs`/… returns `s_re | s_im`, and the crate's
//! contract is that this is exactly the union of the two component statuses.
//! These harnesses prove that union monoid for ALL flag combinations, lifting
//! it from pfloat's example tests (`status.rs`) to proof tier.
//!
//! Run with `cargo kani --no-default-features` (from this crate's directory),
//! NOT the default-feature build: the default profile pulls `pfloat/big`, which
//! activates pfloat's own `#[cfg(all(kani, feature = "big"))]` verify suite and
//! that suite needs math features (`tanh`, `bessel`, …) this crate deliberately
//! does not forward, so it fails to compile. The merge harness needs only
//! `pfloat::Status`, which is available in the bare build. The lane is advisory.

use pfloat::Status;

/// A nondeterministic `Status`: each of the five IEEE 754-2019 §7 flags is
/// independently present or absent, so the harness ranges over all 32 reachable
/// flag sets.
fn any_status() -> Status {
    let mut s = Status::OK;
    if kani::any() {
        s |= Status::INVALID;
    }
    if kani::any() {
        s |= Status::DIV_BY_ZERO;
    }
    if kani::any() {
        s |= Status::OVERFLOW;
    }
    if kani::any() {
        s |= Status::UNDERFLOW;
    }
    if kani::any() {
        s |= Status::INEXACT;
    }
    s
}

/// The componentwise merge a `Complex` operation performs on its two component
/// statuses. Defined here so the harness proves the exact expression the public
/// methods use (`s_re | s_im`).
fn complex_merge(s_re: Status, s_im: Status) -> Status {
    s_re | s_im
}

/// The load-bearing contract: the merged complex status sets a flag iff EITHER
/// component set it. A `Complex` operation must never drop a component's
/// exception flag, and must never invent one.
#[kani::proof]
fn complex_status_merge_is_exact_union() {
    let a = any_status();
    let b = any_status();
    let m = complex_merge(a, b);
    assert!(m.invalid() == (a.invalid() || b.invalid()));
    assert!(m.div_by_zero() == (a.div_by_zero() || b.div_by_zero()));
    assert!(m.overflow() == (a.overflow() || b.overflow()));
    assert!(m.underflow() == (a.underflow() || b.underflow()));
    assert!(m.inexact() == (a.inexact() || b.inexact()));
}

/// `Status::OK` is the identity of the componentwise merge, so a component that
/// raised nothing never perturbs the other component's flags.
#[kani::proof]
fn complex_status_ok_is_identity() {
    let a = any_status();
    assert!(complex_merge(a, Status::OK) == a);
    assert!(complex_merge(Status::OK, a) == a);
}

/// The merge is commutative: the result does not depend on which component is
/// the real and which the imaginary part.
#[kani::proof]
fn complex_status_merge_is_commutative() {
    let a = any_status();
    let b = any_status();
    assert!(complex_merge(a, b) == complex_merge(b, a));
}

/// The merge is associative, so chained merges (e.g. `to_polar`'s `s_r | s_t`,
/// or a multi-step expression) are unambiguous regardless of grouping.
#[kani::proof]
fn complex_status_merge_is_associative() {
    let a = any_status();
    let b = any_status();
    let c = any_status();
    assert!(complex_merge(complex_merge(a, b), c) == complex_merge(a, complex_merge(b, c)));
}

/// The merge is idempotent: merging a status with itself is a no-op, so a
/// component status folded in twice cannot accumulate spurious flags.
#[kani::proof]
fn complex_status_merge_is_idempotent() {
    let a = any_status();
    assert!(complex_merge(a, a) == a);
}
