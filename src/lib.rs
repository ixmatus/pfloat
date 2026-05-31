//! pfloat: pure Rust correctly-rounded arbitrary-precision floats.
//!
//! The public surface is stable under semver as of v1.0 (ADR-0054).
//! The arithmetic kernels, the elementary transcendental and special
//! function surface, and both precision profiles ([`BigFloat`] and
//! `FixedFloat<const PREC: u32>`) are implemented; the Phase 1
//! correctness sweep against the exhaustive binary32-input oracle
//! (every one of the `2^32` f32 values used as a test input; pfloat
//! computes at high working precision internally and the result
//! gets rounded to f32 for the bit-exact comparison) is complete
//! per ADR-0033. See `DESIGN.md` at the repository root for the
//! full design and `docs/decisions/` for the architecture decision
//! records that capture the load-bearing choices.
//!
//! # Quickstart
//!
//! ```
//! # #[cfg(feature = "big")] {
//! use pfloat::{BigFloat, RoundingMode};
//!
//! // BigFloat: runtime precision, heap-allocated mantissa.
//! // Square root of two at 200-bit precision, correctly rounded to
//! // nearest even. The call returns the result and a Status carrying
//! // any IEEE 754-2019 sticky exception flags raised by the operation.
//! let two = BigFloat::try_from_i64_exact(2, 200).unwrap();
//! let (_sqrt2, _status) = two.sqrt(RoundingMode::NearestEven);
//! # }
//! ```
//!
//! Five IEEE 754-2019 rounding modes are available
//! ([`RoundingMode::NearestEven`], [`NearestAway`](RoundingMode::NearestAway),
//! [`TowardZero`](RoundingMode::TowardZero),
//! [`TowardPositive`](RoundingMode::TowardPositive),
//! [`TowardNegative`](RoundingMode::TowardNegative)); every kernel
//! returns a `(value, Status)` pair so callers can inspect or
//! accumulate sticky flags without thread-local state. Under the
//! `std` feature, the same flags also accumulate into a thread-local
//! set accessible via [`flags`].
//!
//! # Scope target
//!
//! v1.0 ships an MPFR-equivalent surface: arithmetic with all five
//! IEEE 754-2019 rounding modes, sticky exception flags, correctly-
//! rounded transcendentals (exp, log, trig, hyperbolic, pow), and
//! special functions (gamma family, erf family, Bessel, zeta,
//! Ei/Si/Ci, Airy, AGM).
//!
//! Two precision profiles share the same operations:
//!
//! - [`BigFloat`]: runtime-determined precision. Requires the
//!   `alloc` feature.
//! - `FixedFloat<const PREC: u32>`: compile-time precision via a
//!   const generic. Stack-allocated, works without `alloc`. Lands
//!   in slice 1g.
//!
//! # Slices 1a–1c (currently shipped)
//!
//! 1a: [`BigFloat`] type, classification predicates, comparison
//! ([`partial_cmp`](BigFloat::partial_cmp),
//! [`total_cmp`](BigFloat::total_cmp), [`min`](BigFloat::min),
//! [`max`](BigFloat::max)), and exact integer construction
//! ([`try_from_i64_exact`](BigFloat::try_from_i64_exact)).
//!
//! 1b: full [`Status`] (all five IEEE flags), [`RoundingMode`], the
//! universal rounding pipeline, std-only thread-local flag accessors
//! ([`flags`](status::flags)), and rounding-required constructors
//! ([`try_from_i64_round`](BigFloat::try_from_i64_round) and
//! [`round_to_precision`](BigFloat::round_to_precision)).
//!
//! 1c: first arithmetic kernel.
//! [`add`](BigFloat::add) and [`sub`](BigFloat::sub) (plus
//! [`add_round`](BigFloat::add_round) / [`sub_round`](BigFloat::sub_round)
//! and `_with_flags` siblings) handle NaN propagation, signed
//! infinities, signed-zero arithmetic (including the IEEE sign rule
//! for `±0 ± ±0` under `TowardNegative`), and mantissa alignment by
//! `2^Δ` for arbitrary exponent gaps.
//!
//! 1d: [`mul`](BigFloat::mul) (plus
//! [`mul_round`](BigFloat::mul_round) and `_with_flags` siblings).
//! Schoolbook + Karatsuba multiplication via the shared
//! `ops::limbs` module; FFT (Schönhage-Strassen) deferred to 1.x
//! per ADR-0010. Handles `0 × ∞ → qNaN + INVALID` per IEEE 754
//! §7.2.
//!
//! 1e: [`div`](BigFloat::div) (plus
//! [`div_round`](BigFloat::div_round) and `_with_flags` siblings).
//! Bit-by-bit long division of the mantissas with sticky-bit
//! tracking from the remainder; routes through the rounding
//! pipeline. Raises `DIV_BY_ZERO` per IEEE 754 §7.3 for
//! `finite_nonzero / 0` and `INVALID` for `0 / 0` and `∞ / ∞`.
//!
//! 1f: [`sqrt`](BigFloat::sqrt) and
//! [`fma`](BigFloat::fma), the last two arithmetic primitives in
//! Phase 1 (plus `_round` and `_with_flags` siblings each). `sqrt`
//! uses bit-by-bit integer square root with parity-adjusted shift
//! so the result exponent splits cleanly; `fma` builds an exact
//! product `BigFloat` then re-rounds via `add_round` for the IEEE
//! 754 §9.4 single-rounding guarantee. The slice also fixes a
//! latent `addsub` bug exposed by FMA (the result-exponent formula
//! used `e_s - p_s + 1` instead of the genuine
//! `min(scale_l, scale_s)`; both happened to coincide whenever
//! `scale_s` was the minimum, which all of slice 1c's tests
//! exercised, but cross-precision FMA does not).
//!
//! 1g: [`FixedFloat<PREC>`](FixedFloat) (const-generic counterpart
//! to [`BigFloat`] with stack-allocated mantissa) plus `core::ops`
//! operator overloads (`Add`, `Sub`, `Mul`, `Div`, `Neg`, and the
//! `*Assign` siblings) for both types behind the `ops` feature.
//! Phase 1 arithmetic surface complete on both types.
//!
//! # Phase 2 (in progress)
//!
//! 2a: Decimal string parsing via
//! [`BigFloat::parse_str`](BigFloat::parse_str) and
//! [`FixedFloat::parse_str`](FixedFloat::parse_str). Accepts
//! signed decimal numbers with optional decimal point and `e/E`
//! exponent, plus case-insensitive `nan` / `inf` / `infinity`.
//! Routes through the rounding pipeline via a multi-precision
//! `m × 5^exp × 2^exp` decomposition.
//!
//! 2b: Decimal formatting via
//! [`BigFloat::to_decimal_string`](BigFloat::to_decimal_string)
//! and the standard [`Display`](core::fmt::Display) trait for both
//! types. `Display` uses
//! [`BigFloat::round_trip_digit_count`](BigFloat::round_trip_digit_count)
//! digits — enough that `parse_str` at the same precision recovers
//! the exact value. Output uses fixed-point notation for moderate
//! magnitudes and scientific notation otherwise. A finite value
//! whose magnitude exceeds the parse round-trip range saturates to
//! `inf` or `0` rather than rendering unbounded digits, and the
//! digit extraction is sub-quadratic (ADR-0051, ADR-0052). Shortest
//! round-trip output (Dragon4 / Steele-White) stays a deferred
//! follow-up under ADR-0029; the current output is correct and
//! round-trips, just not always minimal.
//!
//! # Phase 3 (in progress)
//!
//! 3a: First transcendental — [`BigFloat::exp`](BigFloat::exp) (and
//! `FixedFloat<PREC>::exp`) behind the `exp-log` cluster feature.
//! Algorithm: range-reduce by `ln(2)` (hardcoded 1024-bit constant)
//! so the residual `|r| ≤ ln(2)/2`, evaluate the Taylor series, and
//! compose with a free exponent shift for the `2^k` factor. Fixed
//! 64-bit guard above target precision; full Ziv-strategy retry
//! deferred. Lefèvre–Muller worst-case verification ships in
//! `tests/differential_lefevre_muller.rs`, asserting the kernel's
//! binary64 result matches an mpmath-derived oracle on a subset
//! of the CORE-MATH hard-to-round-case corpus.
//!
//! 3b: [`BigFloat::ln`](BigFloat::ln) (and `FixedFloat<PREC>::ln`).
//! Range-reduce by the binary exponent
//! (`x = m × 2^e` with `m ∈ [1, 2)`), then `ln(x) = ln(m) + e · ln(2)`.
//! The `ln(m)` part uses the atanh series
//! `ln(m) = 2·atanh((m-1)/(m+1))`, which converges roughly three
//! bits per term over `m ∈ [1, 2)`. `ln(2)` is shared with `exp`.
//! Same fixed 64-bit guard as 3a.
//!
//! 3c: [`BigFloat::pow`](BigFloat::pow) (and `FixedFloat<PREC>::pow`).
//! For positive finite `x` and finite `y`, evaluates
//! `exp(y · ln(x))` at working precision and rounds back to target.
//! Negative bases dispatch on integer parity of `y` (per IEEE 754-2019
//! §9.2.1): odd integer flips the sign, even integer leaves it
//! positive, non-integer raises `INVALID` and returns qNaN. Full
//! special-case table for `0`, `±∞`, `NaN`, and the `pow(±1, ±∞) = 1`
//! rule.
//!
//! 3d: thin wrappers around `exp` and `ln`.
//! [`expm1`](BigFloat::expm1) and [`log1p`](BigFloat::log1p) handle
//! the small-argument cancellation regime by boosting working
//! precision proportional to `-exponent(x)` before invoking the base
//! transcendental, then rounding back to target.
//! [`exp2`](BigFloat::exp2) / [`exp10`](BigFloat::exp10) compose as
//! `exp(x · ln(b))` for `b ∈ {2, 10}`; `ln(2)` reuses the hardcoded
//! 1024-bit constant, `ln(10)` is computed lazily per call.
//! [`log2`](BigFloat::log2) / [`log10`](BigFloat::log10) compose as
//! `ln(x) / ln(b)`. All six are mirrored on `FixedFloat`.
//!
//! 3f: forward trig family.
//! [`sin`](BigFloat::sin), [`cos`](BigFloat::cos), and
//! [`tan`](BigFloat::tan) on `BigFloat` and `FixedFloat`. Behind a
//! separate `trig` feature flag. Argument reduction multiplies `x`
//! by a hardcoded 4096-bit table of `2/π` (Payne-Hanek style),
//! yielding `q mod 4` (the quadrant index) and a reduced
//! `r ∈ [−π/4, π/4]`. The kernel then evaluates `sin(r)` or
//! `cos(r)` via Taylor and dispatches on the quadrant. The table
//! caps the supported input range at roughly `|x| < 2^3000`; past
//! that the kernel raises `INVALID` and returns qNaN. `tan` is
//! defined as `sin/cos` over the same reduction.
//!
//! 4a: erf family — opens Phase 4 (tier-1 special functions) with
//! [`erf`](BigFloat::erf) and [`erfc`](BigFloat::erfc) on `BigFloat`
//! and `FixedFloat`, behind a `specials` cluster feature. For
//! `|x|` inside a target-dependent threshold the kernel evaluates
//! the Maclaurin series of `erf` at working precision boosted by
//! roughly `x² · log₂ e` bits, then composes `1 − erf(x)` for
//! `erfc`. For larger `|x|` the kernel evaluates the divergent
//! asymptotic expansion of `erfc` to its smallest-term truncation,
//! then composes `1 − erfc(x)` for `erf`. A hardcoded 1024-bit
//! `2/sqrt(π)` constant supplies the leading coefficient.
//!
//! 4b: gamma family — adds [`gamma`](BigFloat::gamma) and
//! [`lgamma`](BigFloat::lgamma). `lgamma` uses Stirling's
//! asymptotic series at a target-dependent shift `z_min`, with
//! 17 hardcoded Bernoulli-derived coefficients (the largest pair
//! that fits in `i64`), and reflects negative `x` through
//! `ln Γ(x) = ln π − ln|sin πx| − lgamma(1 − x)`. `gamma`
//! composes `sign · exp(lgamma)`, with the sign for negative
//! non-integer `x` taken from `sin(πx)`. Pulls in `trig` because
//! the reflection routes through `sin`. Hardcoded `ln(2π)` at
//! 1024 bits supplies the Stirling constant.
//!
//! 4c: closes Phase 4 by adding
//! [`digamma`](BigFloat::digamma) and [`beta`](BigFloat::beta).
//! `digamma` reuses the lgamma Stirling shift, then evaluates the
//! derivative series `ψ(z) = ln(z) − 1/(2z) − Σ B_{2k}/(2k z^{2k})`
//! with coefficients derived from the same 17-pair table.
//! Negative non-integer `x` reflects through
//! `ψ(x) = ψ(1 − x) − π · cot(πx)`. `beta` evaluates
//! `exp(lgamma(a) + lgamma(b) − lgamma(a + b))` and currently
//! restricts to `a, b > 0` (negative/zero inputs return
//! `qNaN + INVALID`).
//!
//! 3g: inverse trig family.
//! [`atan`](BigFloat::atan) uses the identity
//! `atan(x) = π/2 − atan(1/x)` for `|x| > 1` to bring the argument
//! into `[0, 1]`, then applies the half-angle identity
//! `atan(y) = 2 · atan(y / (1 + sqrt(1 + y²)))` repeatedly until
//! `|y| < 1/16`, then sums the Taylor series. `asin(x)` and
//! `acos(x)` route through `atan` via cancellation-free
//! identities: `asin(x) = 2 · atan(x / (1 + sqrt(1 − x²)))`;
//! `acos(x) = 2 · atan(sqrt((1 − x)/(1 + x)))` for `x ≥ 0`, and
//! the reflected variant for `x < 0`. [`atan2`](BigFloat::atan2)
//! dispatches the full IEEE 754-2019 §9.2.1 special-case table for
//! `(y, x)` signs at zero, infinity, and the four quadrants, then
//! reduces to `atan(y/x)` plus a `±π` shift for the second and
//! third quadrants. All four are mirrored on `FixedFloat`.
//!
//! # Phase 5 — verification
//!
//! Phase 5 lands verification infrastructure on the Phase 1–4
//! surface. The work is internal: no public API changes. ADR-0012,
//! ADR-0013, and ADR-0014 record the architecture.
//!
//! 6a: scaffold slice. `src/verify/` ships Kani harnesses for
//! `add` (NaN propagation, ±∞ ± ±∞, ±0 ± ±0 sign rule, signaling-NaN
//! status). `fuzz/` ships a libfuzzer-sys subcrate with a `parse`
//! target (panic-freedom + `parse(fmt(x)) == x` round-trip).
//! `tests/differential/` ships the MPFR differential lane against
//! `rug`, with a canonical `differential_add` test that asserts
//! bit-for-bit agreement between `BigFloat::add` and `rug::Float`
//! addition for integer operands at precisions 53, 113, 256, 1024
//! and all five rounding modes. The Kani CI lane is advisory
//! (`continue-on-error: true`) per the `feedback_kani_ci_timeout_ok`
//! engineering memory; the fuzz and differential lanes are blocking.
//! Slices 6b–6g extend the pattern across the rest of the surface.
//!
//! 6h: differential lane local sweep + bit-exactness fixes. The
//! first local MPFR sweep surfaced three structural limitations
//! the original test infrastructure did not anticipate; all are
//! documented in ADR-0014's status update. At the time the
//! differential lane was made NearestEven-only (bit-exact conversion
//! via Display + `rug::Float::parse` is rounding-mode-aware and loses
//! up to 1 ULP under non-NE modes; full 5-mode coverage needs a
//! bit-exact converter), transcendental tests were precision-capped at
//! 256 bits via `TRANSCENDENTAL_PRECISIONS`, and the `pow` differential
//! used a 2 ULP tolerance. Later slices superseded the latter two:
//! slice 7b's AGM-based on-the-fly `π` and `ln(2)` lifted the 1024-bit
//! constant ceiling (so `TRANSCENDENTAL_PRECISIONS` now reaches 1024),
//! and slice 7c made `pow` correctly-rounded via the Ziv driver (the
//! differential now asserts bit-exact agreement, consistent with the
//! "Scope target" note above). 22 differential test files passed under
//! `--features=differential-mpfr` on macOS at the time.
//!
//! 6g: Phase 5 close. `fuzz/oss-fuzz/` ships the
//! upstream-submission scaffold (`Dockerfile`, `build.sh`,
//! `project.yaml`, `README.md`) for the OSS-Fuzz PR. ADR-0012,
//! ADR-0013, and ADR-0014 move to `accepted (Phase 5 complete)`
//! with status-update sections documenting what landed and what
//! was deferred. Phase 5 final counts: 196 Kani harnesses across
//! 38 op files, 7 fuzz targets, 22 differential test files.
//! 744 tests pass under
//! `--features=std,fmt,big,fixed,ops,exp-log,trig,specials`
//! unchanged across all of Phase 5.
//!
//! 6f: classification, comparison, parse, and fmt verification.
//! Four new `src/verify/<op>.rs` files cover the non-arithmetic
//! surface: classification totality on the canonical set (every
//! value is exactly one of NaN, infinite, zero, normal),
//! `signum` and `abs` sign-preservation rules, `total_cmp`
//! reflexivity, the `partial_cmp` NaN contract (returns `None`
//! for any NaN; quiet NaN raises no flag, signaling NaN raises
//! `INVALID`), `min` and `max` idempotence, parse of the
//! canonical literals (`nan`, `inf`, `-inf`, `0`), and Display
//! string contents for the canonical special values. One new
//! fuzz target (`fmt`) exercises the Display to parse round-trip
//! on i64 inputs. One new differential test
//! (`differential_parse`) compares pfloat's parser against rug's
//! on a hand-curated battery of decimal strings plus a
//! parse to Display to parse round-trip.
//!
//! 6e: tier-1 specials verification. Six new
//! `src/verify/<op>.rs` files (erf, erfc, gamma, lgamma, digamma,
//! beta) cover NaN propagation, sNaN INVALID, the gamma pole at
//! `±0` (returns `±∞ + DIV_BY_ZERO`), the gamma pole at every
//! non-positive integer (returns `qNaN + INVALID`), the lgamma
//! pole at `±0` (returns `+∞ + DIV_BY_ZERO`), the digamma pole
//! at `0` and at every non-positive integer (returns `−∞ +
//! DIV_BY_ZERO`), the erf saturation rules at `±∞` (returns
//! `±1`), the erfc saturation rules (`+∞ → +0` and `−∞ → +2`),
//! and beta's restricted domain (non-positive `a` or `b` returns
//! `qNaN + INVALID`; pfloat's current beta is positive-only per
//! slice 4c). One new fuzz target (specials) dispatches across
//! the six ops. Five new differential tests cover erf, gamma,
//! lgamma, digamma, and beta. The beta oracle is the lgamma
//! composition evaluated at higher working precision with 2 ULP
//! slack to absorb compounded rounding.
//!
//! 6d: trig + inverse + hyperbolic verification. Thirteen new
//! `src/verify/<op>.rs` files cover the circular and hyperbolic
//! transcendentals. NaN propagation, sNaN INVALID, domain checks
//! for asin/acos/atanh (qNaN + INVALID for `|x| > 1`), the
//! `atanh(±1) = ±∞ + DIV_BY_ZERO` endpoint, `acosh(x < 1) → qNaN +
//! INVALID`, and the `cos(±0) = +1` / `cosh(±0) = +1` /
//! `tanh(±∞) = ±1` identity rules. Two new fuzz targets (`trig`,
//! `hyperbolic`) dispatch across the seven trig ops and the six
//! hyperbolic ops respectively. Seven new differential tests
//! cover sin, cos, tan, atan2, sinh, cosh, asinh — one per
//! kernel whose MPFR counterpart admits straightforward bit-for-
//! bit agreement on integer inputs.
//!
//! 6c: exp/log family verification. Nine new `src/verify/<op>.rs`
//! files (`exp`, `expm1`, `exp2`, `exp10`, `ln`, `log1p`, `log2`,
//! `log10`, `pow`) cover NaN propagation, signaling-NaN INVALID,
//! the `ln(±0) = −∞ + DIV_BY_ZERO` form, the
//! `ln(negative_finite) → qNaN + INVALID` and `ln(−∞) → qNaN +
//! INVALID` forms, the `pow(NaN, ±0) = 1` and `pow(+1, NaN) = 1`
//! §9.2.1 trumping rules, `pow(±0, neg) = ±∞ + DIV_BY_ZERO`,
//! `pow(±0, pos) = ±0`, `pow(±∞, neg) = ±0`, `pow(±∞, pos) = ±∞`,
//! and the `pow(neg, non-integer) → qNaN + INVALID` rule. The
//! `exp_log_family` fuzz target dispatches across all nine ops.
//! Three new differential tests (`differential_{exp,ln,pow}.rs`)
//! cover one differential per mathematically distinct kernel
//! (the *2/*10/m1/1p variants share underlying kernels with `exp`
//! and `ln`).
//!
//! 6b: arithmetic-core verification. Kani harnesses for `sub`,
//! `mul`, `div`, `sqrt`, and `fma` cover NaN propagation,
//! signaling-NaN INVALID, the IEEE 754 invalid forms (`±0 × ±∞`,
//! `0/0`, `∞/∞`, `sqrt(negative_finite)`, `(0 × ∞) + finite`), the
//! `DIV_BY_ZERO` flag for `finite_nonzero / 0`, sign-of-product /
//! sign-of-quotient correctness, the `sqrt(−0) = −0` rule, and the
//! subtle §7.2 carve-out where `(0 × ∞) + NaN` propagates the NaN
//! without raising an extra INVALID. The fuzz arith target derives
//! an `(Op, a, b, c, mode)` tuple, exercises each of the six
//! arithmetic ops, and asserts the cheap identity invariants
//! `a + 0 ≡ a`, `a × 1 ≡ a`, and `a − a = 0`. Five new differential
//! tests (`differential_{sub,mul,div,sqrt,fma}.rs`) mirror the
//! integer-operand pattern from `differential_add.rs`. Slice 6b
//! closes DESIGN.md's verbatim Phase 5 properties for the
//! arithmetic core.
//!
//! 3e: hyperbolic family.
//! [`sinh`](BigFloat::sinh) uses
//! `(expm1(x) − expm1(−x))/2` to avoid the cancellation of
//! `exp(x) − exp(−x)` for small `|x|`.
//! [`cosh`](BigFloat::cosh) is the direct
//! `(exp(x) + exp(−x))/2` (both summands non-negative, no
//! cancellation).
//! [`tanh`](BigFloat::tanh) computes
//! `(1 − exp(−2|x|))/(1 + exp(−2|x|))` then flips the sign for
//! negative `x`; the form is well-behaved at `±∞` (returns ±1
//! directly).
//! [`asinh`](BigFloat::asinh) uses
//! `sign(x) · log1p(|x| + |x|²/(sqrt(|x|² + 1) + 1))`, avoiding the
//! catastrophic cancellation in `ln(x + sqrt(x² + 1))` for large
//! negative `x`.
//! [`acosh`](BigFloat::acosh) uses
//! `log1p((x − 1) + sqrt((x − 1)(x + 1)))` so the near-1 case
//! routes through `log1p`'s cancellation-aware path.
//! [`atanh`](BigFloat::atanh) uses
//! `(log1p(x) − log1p(−x))/2`. Domain check rejects `|x| > 1` with
//! `INVALID` and `|x| = 1` with `DIV_BY_ZERO`. All six are mirrored
//! on `FixedFloat`.

// pfloat depends on `feature(generic_const_exprs)` for the
// `FixedFloat<const PREC: u32>` storage spelling that lands in
// slice 1g. ADR-0011 records the trade-off (nightly toolchain
// required); the feature is `incomplete` upstream, so its lint is
// allowed at the crate root.
#![feature(generic_const_exprs)]
#![allow(incomplete_features)]
#![cfg_attr(not(feature = "std"), no_std)]
#![forbid(unsafe_code)]

#[cfg(feature = "alloc")]
extern crate alloc;

#[cfg(feature = "std")]
extern crate std;

#[cfg(feature = "big")]
mod big;
#[cfg(feature = "big")]
mod class;
#[cfg(feature = "big")]
mod classify;
#[cfg(feature = "big")]
mod cmp;
#[cfg(all(feature = "big", feature = "fixed"))]
mod fixed;
#[cfg(feature = "big")]
mod fmt;
mod mantissa;
#[cfg(any(feature = "exp-log", feature = "agm"))]
mod math;
#[cfg(feature = "big")]
mod ops;
#[cfg(all(feature = "big", feature = "ops"))]
mod ops_traits;
#[cfg(feature = "big")]
mod parse;
#[cfg(feature = "big")]
mod rounding;
mod sign;
mod status;
#[cfg(all(kani, feature = "big"))]
mod verify;

pub use sign::Sign;
pub use status::Status;

#[cfg(feature = "std")]
pub use status::flags;

#[cfg(feature = "big")]
pub use big::{BigFloat, BuildError, Parts};
#[cfg(feature = "big")]
pub use classify::IeeeClass;
#[cfg(all(feature = "big", feature = "fixed"))]
pub use fixed::{ClassFixed, FixedFloat};
#[cfg(feature = "big")]
pub use parse::ParseError;
#[cfg(feature = "big")]
pub use rounding::RoundingMode;

/// Per-kernel Ziv-driver calibration constants and the
/// thread-local trace consumer. Re-exported under the
/// `ziv-instrumented` feature so the pf-tqzz cross-check harness
/// (ADR-0039, slice p1g.3) can drain the most-recent Ziv trace and
/// look up the kernel's calibrated `error_guard`. Off in production
/// builds: the thread-local capture costs nothing without the
/// feature, and the constants stay `pub(crate)` to internal callers.
#[cfg(all(
    feature = "std",
    any(test, feature = "ziv-instrumented"),
    feature = "big",
    feature = "exp-log"
))]
pub mod ziv_instrumented {
    pub use crate::math::ziv::{take_last_trace, ZivTrace};
    pub use crate::math::ziv_calibration::*;
}
