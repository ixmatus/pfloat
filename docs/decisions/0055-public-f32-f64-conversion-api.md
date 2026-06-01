# ADR-0055: Public f32/f64 conversion API

Status: Accepted (2026-05-31)

## Context

pfloat shipped v1.0 with no public way to move between a `BigFloat`
and a hardware `f32` or `f64`. The conversions existed only inside the
Phase 1 oracle harness (`tests/oracle/convert.rs`), in two forms that
cannot ship in the library:

- the mode-aware `BigFloat → f32` path routes through `rug::Float`
  (MPFR), a C dependency the pure-Rust crate does not carry; and
- the pure-Rust `bf_to_f32_bits` rounds through a decimal `Display`
  plus Rust's `str::parse::<f32>`, which is `std` only (no `core`
  float parser) and `NearestEven` only.

The Phase 2 libm spinoff (`pfloat-libm`) is a thin shell that widens a
hardware float to a `BigFloat`, computes correctly rounded, and rounds
back to the hardware width on every call. It needs a public bridge that
is pure Rust, `no_std`, and mode aware. The `BigFloat → float` step is
a known double-rounding hazard: rounding to an intermediate width and
then to the format lands on the wrong neighbour in the subnormal range
(ferrodec flagged the analog as a real bug class).

## Decision

Add a public conversion API as methods on `BigFloat`, behind the `big`
feature, pure Rust and `no_std`:

- `from_f32(f32) -> BigFloat` and `from_f64(f64) -> BigFloat`: exact
  widening at precision 24 / 53. A hardware float is dyadic with at
  most that many significant bits, so the result is lossless.
- `to_f32_round(&self, RoundingMode) -> (f32, Status)` and
  `to_f64_round`: round to the format under the given mode, paired with
  the IEEE 754-2019 sticky flags raised (`INEXACT`, `OVERFLOW`,
  `UNDERFLOW`, and `INVALID` for a signaling-NaN operand).
- `to_f32(&self) -> f32` and `to_f64`: `NearestEven` convenience.

The narrowing rounds **once, straight to the format grid**, which is
what removes the double rounding. It reuses pfloat's verified
`round_to_precision` at the precision the format affords at the value's
magnitude: the full significand (24 / 53) for normals, and the reduced
significand the fixed exponent floor leaves for subnormals
(`exponent - sub_min_exp + 1`). The rounded value then sits exactly on
the grid, and its IEEE fields are read directly from the single,
left-aligned mantissa limb with a shift and a mask. No decimal round
trip, no `rug`, no second rounding. Overflow past the format maximum
yields `±∞` or `±MAX` per the mode's IEEE 754 §4.3 rule with
`OVERFLOW | INEXACT`; a magnitude below half the smallest subnormal
yields `0` or the smallest subnormal per mode with `UNDERFLOW |
INEXACT`.

The API is additive, so it lands as a semver-minor extension (1.1) of
the surface ADR-0054 froze; no existing signature changes.

## Consequences

- `pfloat-libm` widens and narrows through public, supported methods
  rather than test-only or `std`-only code, and its outer Ziv loop can
  build enclosures from `to_f32_round` directly.
- The conversion is itself oracle covered. A new differential lane
  (`tests/convert_oracle.rs`, `differential-mpfr`) cross-checks the
  pure-Rust result against MPFR over a broad off-grid sweep (p=53 from
  random f64s and p=128 ratios) across all five modes; the first run
  agreed bit-exact on every sampled value. Round-trip identity and the
  overflow, subnormal, and tie corners are unit tested in
  `src/convert.rs`.
- The harness keeps its `rug`-based `certified_round_bf_to_f32` as the
  independent oracle; the library code is the thing under test, not the
  oracle, so the two stay separate.
- `from_f32` returns a precision-24 value (exact). Callers that compute
  at a wider working precision widen with `round_to_precision`, which
  is exact when widening.
- `FixedFloat<PREC>` conversions are deferred; the shell uses
  `BigFloat`, and `FixedFloat` delegates to it.

## Related

- ADR-0054 (v1.0 public API freeze): this is the first additive 1.1
  extension of that surface.
- ADR-0035 (oracle protocol): the harness `convert.rs` whose shippable
  gap this fills; the rug path there becomes the cross-check oracle.
- The Phase 2 roadmap (libm spinoff) and the `pf-lm1a` slice.
