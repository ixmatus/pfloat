//! Slice 7f.0: allocation profiling pass (zero new dependency for
//! pfloat).
//!
//! ADR-0004 deferred smallvec-style inline `BigFloat` storage until
//! "a hot path shows the allocation cost dominating". No profiling
//! was done in Phase 6, so that trigger was unmeasured. This harness
//! measures heap allocations per operation on a representative kernel
//! set so the slice 7f disposition rests on data, not a guess.
//!
//! It is a standalone crate, not a pfloat example or lib module:
//! pfloat is `#![forbid(unsafe_code)]` package-wide (and `forbid`
//! cannot be locally overridden), while a counting global allocator
//! is necessarily `unsafe impl GlobalAlloc`. Living here keeps
//! pfloat's unsafe-free invariant intact and adds no dependency to
//! pfloat. ADR-0028 records the rationale.
//!
//! Run:
//! `cargo run --release --manifest-path tools/alloc-profile/Cargo.toml`

use std::alloc::{GlobalAlloc, Layout, System};
use std::hint::black_box;
use std::sync::atomic::{AtomicUsize, Ordering};

use pfloat::{BigFloat, RoundingMode};

static ALLOC_CALLS: AtomicUsize = AtomicUsize::new(0);
static ALLOC_BYTES: AtomicUsize = AtomicUsize::new(0);

/// Forwards every call verbatim to the system allocator and, on the
/// allocating paths, bumps two relaxed counters.
struct Counting;

// SAFETY: every method delegates unchanged to `System`, which upholds
// the `GlobalAlloc` contract (valid pointers, correct layouts, no
// aliasing). The wrapper only increments two relaxed atomics, which
// introduces no new pointer, lifetime, or aliasing obligation, so the
// safety contract is exactly `System`'s.
unsafe impl GlobalAlloc for Counting {
    unsafe fn alloc(&self, layout: Layout) -> *mut u8 {
        ALLOC_CALLS.fetch_add(1, Ordering::Relaxed);
        ALLOC_BYTES.fetch_add(layout.size(), Ordering::Relaxed);
        System.alloc(layout)
    }
    unsafe fn dealloc(&self, ptr: *mut u8, layout: Layout) {
        System.dealloc(ptr, layout);
    }
    unsafe fn alloc_zeroed(&self, layout: Layout) -> *mut u8 {
        ALLOC_CALLS.fetch_add(1, Ordering::Relaxed);
        ALLOC_BYTES.fetch_add(layout.size(), Ordering::Relaxed);
        System.alloc_zeroed(layout)
    }
    unsafe fn realloc(&self, ptr: *mut u8, layout: Layout, new_size: usize) -> *mut u8 {
        ALLOC_CALLS.fetch_add(1, Ordering::Relaxed);
        ALLOC_BYTES.fetch_add(new_size, Ordering::Relaxed);
        System.realloc(ptr, layout, new_size)
    }
}

#[global_allocator]
static GLOBAL: Counting = Counting;

fn snapshot() -> (usize, usize) {
    (
        ALLOC_CALLS.load(Ordering::Relaxed),
        ALLOC_BYTES.load(Ordering::Relaxed),
    )
}

/// Times nothing; measures only allocation. `iters` is chosen per
/// workload so the cheap kernels still run enough to be stable and
/// the expensive ones do not dominate wall time.
fn profile(name: &str, iters: usize, mut op: impl FnMut()) {
    // Warm up: lazily-initialised constants (the AGM ln(2)/pi tables)
    // allocate on first use; exclude that from the per-op figure.
    for _ in 0..16 {
        op();
    }
    let (c0, b0) = snapshot();
    for _ in 0..iters {
        op();
    }
    let (c1, b1) = snapshot();
    let calls = c1 - c0;
    let bytes = b1 - b0;
    println!(
        "{name:<22} iters={iters:>5}  allocs/op={:>8.2}  bytes/op={:>10.1}",
        calls as f64 / iters as f64,
        bytes as f64 / iters as f64,
    );
}

fn bigfloat(s: &str, precision: u32) -> BigFloat {
    BigFloat::parse_str(s, precision, RoundingMode::NearestEven)
        .expect("literal parses")
        .0
}

fn main() {
    println!("pfloat slice 7f.0 allocation profile (representative kernels)\n");

    // Arithmetic mul at p=256: the multiplication path slice 7d also
    // calibrated; the foundation every kernel composes.
    let a = bigfloat("1.4142135623730950488016887242096980785696718753769", 256);
    let b = bigfloat("2.7182818284590452353602874713526624977572470937000", 256);
    profile("mul p=256", 5000, || {
        black_box(black_box(&a).mul(black_box(&b), RoundingMode::NearestEven));
    });

    // Transcendental exp at p=256: range reduction (needs the boosted
    // ln(2) constant) plus a Taylor compose, the elementary-function
    // allocation pattern.
    let x = bigfloat("1.2599210498948731647672106072782283505702514647015", 256);
    profile("exp p=256", 2000, || {
        black_box(black_box(&x).exp(RoundingMode::NearestEven));
    });

    // Special-function compose gamma at p=113: composes ln / ln(2*pi)
    // and several Stirling temporaries, the heaviest allocator class.
    let g = bigfloat("3.6", 113);
    profile("gamma p=113", 1000, || {
        black_box(black_box(&g).gamma(RoundingMode::NearestEven));
    });

    println!(
        "\nallocs/op is the ADR-0004 trigger metric: inline storage pays\n\
         off only if a hot path's per-op allocation count is high enough\n\
         that removing it materially moves the kernel."
    );
}
