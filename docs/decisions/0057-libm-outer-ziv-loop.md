# ADR-0057: The pfloat-libm outer Ziv loop (the enclosure determines the float)

Status: Accepted (2026-06-01)

## Context

`pfloat-libm` is a thin shell over pfloat. A call widens its hardware
float to an exact `BigFloat` (ADR-0055 `from_f32` / `from_f64`),
evaluates one of pfloat's correctly rounded kernels, and rounds the
`BigFloat` result back to the hardware width. That second rounding is the
double rounding hazard: rounding an approximation to one width and then to
the format can land on the wrong neighbour, most visibly in the subnormal
range. ADR-0055 supplied the cure for the final step, `to_f32_round` /
`to_f64_round`, which round straight onto the format grid. This ADR
records how the outer loop uses that converter so no double rounding
survives, and settles the design points ADR-0056 left to the shell: the
enclosure construction, the precision schedule and fallback, and the
status composition.

The shell reimplements no IEEE special case. NaN, the infinities, signed
zero, and the domain errors (`ln` of a negative, `acosh` below one, the
reciprocal-trig poles) are already handled inside pfloat's kernels with
the correct value and `Status`; the shell rounds the kernel's result and
merges the flags.

## Decision

### 1. Enclose with the directed pair

At a working precision `w` ask the kernel for two values: `lo`, the
result rounded toward negative infinity to `w` bits (the largest `w`-bit
value not exceeding `f(x)`), and `hi`, the result rounded toward positive
infinity (the smallest `w`-bit value not below `f(x)`). The interval
`[lo, hi]` is a true enclosure of `f(x)`, at most one `ulp_w` wide. Round
both ends to the format under the requested mode. Every IEEE 754 rounding
mode is a monotone non-decreasing map, so when the two ends land on the
same hardware float the true value between them lands there too: commit
it. Otherwise grow `w`; the directed kernel calls tighten the enclosure
around `f(x)` from both sides.

### 2. Why the directed pair, not a single nearest-even call

The first design for this loop used a single nearest-even kernel call,
taking `r` correctly rounded to nearest at `w` bits and bracketing it with
the half-ULP error bound `[r - 2^(e - w), r + 2^(e - w)]`. That halves the
kernel evaluations, and a half-ULP bracket is a sound enclosure of the
true value, so the nearest modes round correctly. The directed modes do
not. A single nearest-even call discards the sign of the residual
`f(x) - r`, which is exactly what a directed mode needs.

The failure is concrete and not measure zero. Take `exp` of the smallest
`f64` subnormal `x = 2^-1074`, whose true value is `1 + 2^-1074`, a hair
above `1.0`. At every working precision below 1075 bits the kernel rounds
this to `1.0` and, because the residual `2^-1074` falls below its own
working precision, reports the result exact. A nearest-even shell then
sees an exact `1.0` and, under round toward positive, returns `1.0`, where
the correctly rounded answer is the float just above `1.0`. The symmetric
bracket around the on-grid `1.0` straddles the boundary and can never
decide which side `f(x)` is on, because the nearest-even call threw the
side away. The same shape appears for `sin(2^-149)` (just below the
smallest subnormal) and `cot(2^-126)` (just below `2^126`), each off by
one ULP under a directed mode.

The directed pair fixes it because each directed kernel call forces
pfloat's own Ziv loop to resolve the residual in that direction (growing
its internal precision up to `ZIV_GUARD_CAP`), so the enclosure carries
the side information directed rounding requires. It also dispenses with
the machinery the single-call path needed (a half-ULP construction, an
exactness short-circuit on the kernel's `INEXACT` flag, and a guard
against building `2^(e - w)` for a saturated exponent), and that exactness
short-circuit was itself the proximate bug above: the kernel's `INEXACT`
flag is not a reliable witness of exactness, since a residual below the
kernel's working precision reads as exact. Correctness decides this over
the doubled evaluation count; the exhaustive sweep that certifies the
surface is the cost that the doubling falls on, and it is shardable. A
single directed call suffices for a directed target mode (directed
double rounding in one direction is exact), so a future optimisation may
special-case the directed modes to one call; the uniform directed pair is
kept here for one correct code path.

### 3. Precision schedule and the hard-to-round fallback

The working precision is `prec + guard` for `guard` drawn in turn from
`{64, 128, 256, 512, 1024}`, five iterations that mirror pfloat's own Ziv
schedule (`ZIV_BASE_GUARD`, doubling, `ZIV_GUARD_CAP`, `ZIV_MAX_ITERS`;
those constants are private to pfloat, so the shell restates them). On the
measure zero hard to round input that exhausts the schedule the loop
returns the nearest-rounded value at the finest precision with `INEXACT`
set, the same caveat pfloat and MPFR document. No new `Status` flag is
added: the returned status stays the IEEE 754 sticky flags, not a
certification channel. The exhaustive `f32` sweep and `f64` differential
against an independent oracle (the following slice) are what certify the
grid.

### 4. Status composition

The returned status is the union of the two kernel statuses and the
statuses of converting the two ends. The kernel contributes `INVALID` and
`DIV_BY_ZERO` from its domain handling (those make a result non-finite, so
they are merged on the special-case path); the conversions contribute
`INEXACT`, `OVERFLOW`, and `UNDERFLOW`. The union reports `INEXACT`
exactly when `f(x)` differs from the committed float: it is clear only
when both ends equal the committed value exactly, which is to say `f(x)`
is exactly that format float.

One honest limitation remains in the flag. A result that is
mathematically exact but reached through a composed transcendental, such
as `log10(1000) = 3` (`ln(x) / ln(10)`) or `exp10(2) = 100`
(`exp(x * ln 10)`), may carry `INEXACT` because the kernel rounds
internally and cannot cheaply prove the result lands exactly on the grid.
Deciding it is the table maker's dilemma for the inexact flag, and neither
the kernel nor the shell settles it within bounded precision. The flag is
therefore conservative in the safe direction (it may over-report
inexactness, never under-report it) and the correctly rounded value is
unaffected. `log2` of a power of two clears `INEXACT` because it has an
exact exponent path; the composed kernels do not. This is recorded as a
pfloat-side refinement (filed discovered-from pf-lm2), not a shell defect.

## Consequences

- The shell is generic over the format through a small `Shell` trait
  carrying the precision, the exact widening, the grid-exact conversion,
  and bit identity. One driver serves both widths; a declarative macro
  emits the bare nearest-even entry and the mode-aware `_round` entry for
  each of the 25 unary functions, and `hypot` and `rootn` are written
  explicitly for their extra argument.
- The directed pair evaluates the kernel twice per outer iteration. For
  the overwhelming majority of inputs the enclosure agrees on the first
  iteration, so the cost is two kernel evaluations; the exhaustive sweep
  absorbs this and the directed-mode single-call optimisation above is
  available if it is ever needed.
- The returned per-call `Status` is authoritative. Under `std`, the
  conversions also raise into pfloat's thread-local flag set, so that
  thread-local is not meaningful across a shell call. This is documented
  rather than corrected: a save and restore would add machinery for a
  channel the per-call status already serves.
- The shell's certification is conditional on the kernel returning a
  correctly rounded directed result, which itself carries pfloat's measure
  zero `ZIV_MAX_ITERS` caveat. The shell fallback and the pfloat fallback
  are the same caveat, not a compounding one; the independent oracle sweep
  is where the `f32` grid is actually certified.

## References

- ADR-0055: the public `f32` / `f64` conversion API (the grid-exact
  converter the loop builds enclosures from).
- ADR-0056: the six direct libm kernels the shell wraps.
- ADR-0032: the direct-kernel policy that deferred the six to this phase.
- pfloat `DESIGN.md`, the Ziv strategy section; `src/math/ziv.rs` (the
  in-tree enclosure driver this loop mirrors); `src/convert.rs` (the
  grid-exact conversion).
- IEEE 754-2019 §4 (rounding), §7.6 (the inexact exception), §9.2
  (recommended functions).
