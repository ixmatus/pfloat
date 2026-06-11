# Algorithms in the pfloat family

This is a reading guide, not a derivation. Each entry names what an algorithm
does, the one idea that makes it work, the primary source it derives from, and
the architecture decision record that holds the full rationale. The ADRs under
`docs/decisions/` are the durable account; this page orients you into them and
points at the module that implements each piece. ADR numbers are stable; file
paths are a convenience and may drift. Every primary source named here has a
provenance record (verified citation, archived URL, license, rot risk) in the
reference registry at `docs/references/` (see its `INDEX.md`).

The three crates share one engine. `pfloat` is the scalar arithmetic and the
correctly rounded function kernels; `pfloat-ball` wraps a pfloat scalar in a
rigorous enclosure; `pfloat-libm` rounds a pfloat kernel result to a hardware
float. Read the `pfloat` section first, since the companions build on it.

## pfloat

### The Ziv correctly rounded driver

Every transcendental kernel runs the same loop: evaluate at a working precision
above the target, ask whether the rounded result is determined, and if not
widen the precision and retry. The soundness test is an interval test, not a
comparison of two adjacent guesses, because two adjacent guesses can agree on a
wrong digit. A cap of five iterations bounds the measure zero set of inputs no
finite precision resolves; near a zero of the function the driver boosts the
working precision to cover the cancellation. Source: Ziv (1991); the interval
test discipline. See `ADR-0022` and `ADR-0039`; `src/math/ziv.rs`.

### Multiplication: schoolbook, Karatsuba, Toom-3

The limb multiplier dispatches on operand size: schoolbook below 48 limbs,
Karatsuba above it, Toom-3 above 176 limbs (`KARATSUBA_THRESHOLD` and
`TOOM3_THRESHOLD`). The Toom-3 choice was made on an allocation argument rather
than an operation count: splitting by three allocates fewer intermediate limbs
than the recursion it replaces, which is what wins at large precision. The
helpers (the divide by three step, the Vandermonde interpolation) derive from
the published method, not from a transcription of GMP. Source: Brent and
Zimmermann, Modern Computer Arithmetic, section 1.3.3; Bodrato and Zanoni
(ISSAC 2007); Jebelean (1993) for exact division by three. See `ADR-0061` and
`ADR-0027`; `src/ops/limbs.rs`.

### The IEEE 754-2019 remainder kernel

`remainder(x, y) = x - n*y`, where `n` is the integer nearest to `x / y`, so the
result is exact and bounded by half of `y`. The quotient `n` can be enormous (an
exponent sized integer), so the kernel never forms it. It reduces by modular
exponentiation, computing the power of two that aligns the operands modulo the
divisor in a logarithmic number of multiplications, which keeps the operation
safe against an attacker supplied exponent. Source: IEEE 754-2019, section
5.3.1; verified differentially against MPFR. See `ADR-0069`;
`src/ops/remainder.rs`.

### Dragon4 shortest decimal formatting

The opt in decimal formatter emits the shortest digit string that round trips
back to the same value, built on the big integer routines the crate already
carries. `Display` keeps the round trip safe digit count. Source: Steele and
White (PLDI 1990); Burger and Dybvig (PLDI 1996). See `ADR-0071`, which
supersedes `ADR-0029` and builds on `ADR-0051`; `src/fmt.rs`.

### Elementary functions and the settled INEXACT flag

The elementary surface is argument reduction plus a series or an iteration,
every step wrapped in the Ziv driver: Payne and Hanek reduction for the large
argument trigonometric range, the arithmetic geometric mean for the logarithm
and the constants, Stirling and Spouge approximations for the gamma family. The
INEXACT flag is then settled by transcendence: a kernel pre dispatches the
finite set of exact inputs returning `Status::OK`, and forces INEXACT on the
transcendental fall through, because Lindemann Weierstrass and Gelfond Schneider
make the result irrational outside that exact set. Source: Payne and Hanek
(1983); Brent and Salamin for the AGM; the named transcendence theorems. See
`ADR-0015`, `ADR-0017`, `ADR-0060`, `ADR-0063`, `ADR-0064`;
`src/math/{exp,ln,sin,pow,gamma,agm}.rs`.

## pfloat-ball

The ball laws are stated as an in tree specification at `pfloat-ball/src/spec.rs`;
this section is the algorithmic view of those laws.

### Directed pair radius soundness

The midpoint of an arithmetic result is the nearest pfloat value, and the radius
bounds the kernel's own residual plus the propagated input error. The residual
bound comes free from a directed pair: compute the result toward minus infinity
and toward plus infinity, and half that spread (rounded up) bounds how far the
nearest result sits from the truth. This route never invokes the Ziv driver, so
it carries no separate Ziv residual term, which is the distinction the spec's
Law 2 draws. When the directed pair coincides the result is exact and the radius
is zero (Law 3). Source: the Fundamental Theorem of Interval Arithmetic. See
`ADR-0077`; `pfloat-ball/src/arith.rs`.

### Three enclosure shapes for elementary functions

A ball elementary function picks one of three shapes. A monotonic function is
enclosed by its endpoints: evaluate the kernel at the outward rounded interval
ends. A function with derivative bounded by one (`sin`, `cos`) is enclosed by
its midpoint value inflated by the input radius, since the value moves no faster
than the input. A composed function (`tan` as `sin / cos`) is built from the
other two by ball arithmetic, and a pole straddle returns the whole real line
with a divide by zero flag. Source: novel derivation from monotonicity and the
mean value bound. See `ADR-0082`; `pfloat-ball/src/elem.rs`.

### Mag: a radius that rounds up by construction

The radius type is a single limb binary float `m * 2^(e - 63)` with a `u64`
significand, no sign, and no NaN, so a negative or not a number radius cannot be
written. Every operation rounds the significand up, so a radius can only ever
over estimate. The 64 bit significand is wider than Arb's 30 bit magnitude type,
which buys tightness, and being a single limb with no heap it is a tractable
Kani target; the cost is a relative resolution floor near `2^-64`, which is
sound because the radius is only ever an upper bound. Source: novel derivation,
a deliberate divergence from Arb's `mag_t`. See `ADR-0074`;
`pfloat-ball/src/mag.rs`.

### The asymmetric conversion boundary

Ball to endpoints is exact: the lower and upper endpoints are the directed
differences `mid - rad` and `mid + rad`, the tightest representable bracket.
Endpoints to ball is sound but inflating: `from_interval` never assumes the
midpoint is centered, so it sets the radius to cover the farther endpoint, which
defeats a centered midpoint exclusion bug. The word lossless is reserved for the
exact direction only. Source: novel derivation. See `ADR-0076`;
`pfloat-ball/src/ball.rs`.

### The sealed scalar trait

The midpoint engine is a trait sealed to `BigFloat` and `FixedFloat`, so the
claim that a midpoint is always a correctly rounded pfloat scalar cannot be
broken from outside the crate. The seal moves an invariant from convention into
the type system. Source: the sealed trait pattern. See `ADR-0075`;
`pfloat-ball/src/scalar.rs`.

### The independent Arb range soundness backstop

A per release test only lane brackets each ball operation against Arb, reached
out of process through a python-flint subprocess so no interval library enters
the link graph. The point form checks that the ball admits Arb's enclosure of
the true value at sampled witnesses; the interval form (`BRACKETI`) brackets the
function over the whole input interval and asserts the ball is a superset of that
image, scoped to extremum straddles where the check is provably clean. Tightness
is measured per bucket, not asserted, because a correct ball can be tighter than
Arb's interval image away from the extrema. Source: the Fundamental Theorem of
Interval Arithmetic, range soundness. See `ADR-0082`, which extends `ADR-0078`;
`scripts/arb_oracle_worker.py` and `pfloat-ball/tests/differential_arb.rs`.

## pfloat-libm

### The directed pair outer Ziv loop

A correctly rounded hardware float needs more than widen, compute, and round:
the second rounding can double round. The shell brackets the true value with two
directed kernel calls, rounds both ends to the hardware format under the
requested mode, and commits only when the two agree, which proves the rounding.
A single nearest call would discard the residual sign and round a value such as
`exp(2^-1074)` the wrong way under a directed mode, which is the failure this
design closes. Source: IEEE 754-2019, sections 4 and 7.6. See `ADR-0057`;
`pfloat-libm/src/round.rs`.

### The guard bit schedule with a hard to round fallback

The outer loop widens through a fixed schedule of guard bits, `[64, 128, 256,
512, 1024]`, mirroring pfloat's own inner schedule. The measure zero set of
inputs that never resolve falls back to nearest even at the finest precision
with a forced INEXACT flag, the honest analogue of the scalar Ziv cap. Source:
the Ziv cap discipline. See `ADR-0057`; `pfloat-libm/src/round.rs`.

### Six direct kernels reused from pfloat

Six functions skip the high precision shell and call a pfloat kernel directly:
the reciprocal trigonometric functions (`cot`, `sec`, `csc`) at inflated
precision for a single rounding, the exact integer cube root, the integer root
`rootn` by Newton iteration with the standard special case table, and `hypot`
without Moler scaling because the saturating exponent makes the naive form safe.
Source: DLMF section 4.14 for the reciprocal definitions; IEEE 754-2019 section
9.2.1 for the `rootn` special cases. See `ADR-0056`; the kernels live in pfloat,
the shell wraps them at `pfloat-libm/src/f32.rs` and `f64.rs`.

### The saturation fast path

A transcendental whose argument is far past the format's overflow point cannot
change the answer by reducing the argument: it is already infinity or zero. A per
width threshold (1024 for `f32`, 2048 for `f64`) short circuits before the kernel
runs, returning the saturated value with the matching IEEE flags, bit for bit
identical to the kernel path. The bounded functions saturate to their limit (for
example `tanh` to plus or minus one) with INEXACT. A non finite input is never
fast pathed. Source: the pf-hzup sweep infeasibility finding that motivated it.
See `ADR-0063`; `pfloat-libm/src/saturate.rs`.

## See also

- `docs/references.md` for the standards, papers, and per function DLMF pointers.
- `DESIGN.md` for the cross cutting design narrative.
- `docs/rounding-status.md` for the per function rounding verification status.
- `docs/decisions/` for the full architecture decision records cited above.
