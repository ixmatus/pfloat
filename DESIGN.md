# pfloat: design

## Purpose

pfloat plants a flag for the canonical pure Rust answer to
arbitrary-precision floating-point arithmetic with correct rounding.
The point of reference is GNU MPFR. The point of differentiation
from existing pure Rust efforts is completeness against MPFR's
surface (especially special functions), shipped formal-verification
artifacts, and a dual API that serves both runtime-precision and
compile-time-precision users.

The horizon is decades. The crate ships under MIT or Apache 2.0,
operates in `no_std` environments, and follows the ferrodec
methodology for spec acquisition, conformance testing, formal
verification, and differential testing.

## Goals

1. IEEE 754-2019 binary arithmetic at any user-chosen precision,
   correctly rounded under all five rounding modes (RNE, RNA, RZ, RP,
   RM) with sticky exception flags (inexact, overflow, underflow,
   divide-by-zero, invalid).
2. Correctly-rounded elementary transcendentals: `exp` family, `log`
   family, trig and inverses, hyperbolic and inverses, `pow`.
3. Special functions matching MPFR's surface: `gamma`, `lgamma`,
   `digamma`, `beta`, `erf`, `erfc`, Bessel `J0/J1/Jn`,
   `Y0/Y1/Yn`, modified `I0/I1/In`, `K0/K1/Kn`, Riemann zeta on
   the real line, exponential integral `Ei`, sine and cosine
   integrals `Si` and `Ci`, logarithmic integral `li`, Airy `Ai`
   and `Bi`, arithmetic-geometric mean.
4. Two precision profiles in one crate. `BigFloat` carries a runtime
   precision and a heap-allocated mantissa (needs `alloc`).
   `FixedFloat<const PREC: u32>` carries a compile-time precision
   and a stack-allocated mantissa (works without `alloc`).
5. Verification artifacts shipped alongside the code: Kani-discharged
   safety properties, conformance test vectors from Lefèvre–Muller
   worst-case-rounding tables and IEEE-style binary-FP test corpora,
   differential testing against MPFR via `gmp-mpfr-sys` on a
   feature-gated CI lane.
6. `no_std`-first; `alloc`-only; `std` for ergonomic extras
   (thread-local sticky flags, `std::error::Error` impls).
   Cross-compiles to `thumbv6m-none-eabi` in CI for the embedded
   floor.

## Non-goals

- Complex arithmetic. Real-valued only; complex extensions, if any,
  ship in a sibling crate.
- Decimal floating-point. ferrodec covers IEEE 754 decimal.
- Interval arithmetic. `interval-1788` will sit on top of pfloat,
  not inside it.
- Arbitrary-precision integers as a first-class primary type.
  pfloat carries integers as a private intermediate; users wanting
  big integers pick a dedicated crate.
- Mirroring MPFR's C API. MPFR is the behavior oracle, not the API
  guide.

## Numeric representation

### Mantissa layout

Limbs are unsigned `u64` in little-endian order: limb 0 holds the
least significant 64 bits, limb `n - 1` holds the most significant.
The most significant bit of the most significant limb is `1` for
every normalized non-zero value (no implicit bit). Zero, infinity,
and NaN are not represented in the mantissa; they live in the
class tag.

Sign-magnitude. The sign rides as a separate field, never folded
into the mantissa.

`u64` is the standard limb across the literature (`Modern Computer
Arithmetic`, MPFR's source) and gives a clean width-doubling
identity for multiplication: `u64 × u64` fits in `u128`, and
`u128 / u64` is the standard division primitive. Targets where
`u128` is software-emulated (some embedded compilers) pay a small
constant; the alternative (`u32` limbs) doubles the limb count and
the per-operation overhead in the common case.

Decision recorded in ADR-0001.

### Precision granularity

Precision is tracked at the bit level. The minimum precision is one
bit, matching MPFR. The mantissa storage rounds up to whole limbs;
the high bits of the most-significant limb that exceed the working
precision must be zero in the canonical form, and the rounding
pipeline maintains that invariant.

Bit-level precision lets callers ask for IEEE 754 binary32 (24
bits), binary64 (53 bits), binary128 (113 bits), or any other
unusual width directly. 64-bit-aligned precision (astro-float's
choice) is faster but coarser and forces callers to pad; the cost
of bit-level granularity is a small mask on the boundary limb.

Decision recorded in ADR-0002.

### Special values

A tagged `Class` enum carries the value's kind explicitly:

```rust
enum Class {
    Zero { sign: Sign },
    Infinity { sign: Sign },
    Nan { quiet: bool, sign: Sign, payload: NanPayload },
    Normal { sign: Sign, exponent: i64, mantissa: Mantissa },
}
```

Reserved-exponent encoding (MPFR's `MPFR_EXP_NAN`, etc.) saves a
discriminant byte per value at the cost of fragility: every code
path has to remember which exponent values mean what. The tagged
enum costs one discriminant per value and rules out an entire
class of bugs at the type level.

Decision recorded in ADR-0005.

### Exponent

`i64`. Range vastly exceeds anything a real workload uses. Match
MPFR's choice for compatibility with the differential-testing
oracle. ADR-0006.

### Rounding modes and exception flags

Rounding mode is an enum passed at the call site:

```rust
pub enum RoundingMode {
    NearestEven,    // RNE (IEEE default)
    NearestAway,    // RNA
    TowardZero,     // RZ
    TowardPositive, // RP
    TowardNegative, // RM
}
```

A typestate wrapper (`Rounded<Mode, T>`) carries the rounding mode
in the type for callers who want compile-time discipline; the
runtime enum stays the primary interface for ergonomic reasons.

Exception flags are sticky and cumulative, matching IEEE 754. Two
storage strategies sit behind a feature flag:

- `std`: thread-local `Cell<Status>`, accessed via free functions
  (`pfloat::flags::test`, `pfloat::flags::clear`). The ergonomic
  default.
- `no_std`: `Status` argument passed explicitly through every
  arithmetic call. No global state. The honest default for embedded
  use.

ADR-0007.

## Type architecture

### `BigFloat`

Dynamic precision, heap-allocated mantissa.

```rust
pub struct BigFloat {
    class: Class,
    precision: u32, // bits
}

// where Mantissa = Vec<u64>, length = ceil(precision / 64)
```

Suitable for workloads that change precision at runtime, walk
through ranges of precisions for convergence-testing, or hold
values whose precision is determined by external data.

### `FixedFloat<const PREC: u32>`

Compile-time precision, stack-allocated mantissa.

```rust
pub struct FixedFloat<const PREC: u32> {
    class: ClassFixed<PREC>,
}

// where ClassFixed::Normal carries `[u64; ceil(PREC / 64)]`
```

Suitable for embedded use, hot loops at fixed precision, callers
who want the optimizer to see the precision as a constant, and
all `no_std` users without `alloc`.

`FixedFloat<53>` is binary64 with correct rounding under any
rounding mode. `FixedFloat<113>` is binary128. `FixedFloat<512>`
is the natural unit for some signal-processing workloads. The
const-generic surface is a strict superset of what fixed-precision
users would otherwise reach for.

### Conversions

`BigFloat::from(FixedFloat<PREC>)` is exact and infallible (the
fixed value's precision is exactly preserved).

`FixedFloat::<PREC>::try_from(BigFloat)` rounds under a chosen mode
and may set inexact / overflow / underflow flags. The rounding mode
threads through the conversion explicitly, matching the rest of the
API.

Both directions go through a shared `Mantissa` trait so the
arithmetic kernels avoid duplicating logic. The trait sits in a
private module; users see only the concrete types.

ADR-0003 records the dual API decision; ADR-0004 records the
storage choice.

## Arithmetic algorithms

References: Brent and Zimmermann, *Modern Computer Arithmetic* (MCA),
free PDF; Muller et al., *Handbook of Floating-Point Arithmetic*
(HFA, 2nd edition); Fousse et al. 2007 MPFR algorithms paper.

### Addition and subtraction (HFA §5, MCA §3.1.5)

Align the smaller exponent up to the larger; add or subtract limb
arrays; normalize; round. Three sub-cases by relative exponents
(close, far, swap-required) get separate code paths to keep the
common cases tight.

### Multiplication (MCA §1.3, §3.3)

Schoolbook for `n ≤ k1` limbs, Karatsuba between `k1` and `k2`,
Toom-Cook 3-way thereafter (deferred to 1.x; schoolbook plus
Karatsuba ships in 1.0).

`k1` and `k2` are tuned empirically against the bench harness. MPFR
defaults at roughly `k1 = 30` and `k2 = 100`; pfloat's Rust codegen
is likely to shift these. The thresholds are runtime constants
chosen at build time, not config.

Schönhage-Strassen FFT multiplication crosses over above ~10⁴ limbs.
Most users never reach the threshold. Defer to 1.x. ADR-0010.

### Division (MCA §1.4, §3.4)

Newton iteration on the reciprocal. Each step doubles the precision
of the approximation; the final residue corrects the rounding.
Schoolbook division below the Newton crossover.

### Square root (MCA §3.5)

Newton iteration on the reciprocal-square-root, then a final multiply
back. Avoids the integer square root path for high-precision
operands.

### FMA (HFA §5.5)

Compute `a × b` to double width without rounding, add `c` to the
double-width result, round once at the end. The naive implementation
is straightforward; correctness on the tie cases is the test target.

## String I/O

`parse_str(s, precision, mode)` parses a decimal string into a value
at the requested precision under the given rounding mode. The
algorithm is a generalization of Steele–White (Ryu does not extend
cleanly to arbitrary precision). `to_string(value, mode)` formats
round-trip-correct: `parse_str` at the same precision recovers the
exact value. It emits `round_trip_digit_count(p)` digits, enough to
guarantee round-trip, not the minimal number that round-trips;
shortest output (Dragon4 / Steele–White) is deferred to 1.x per
ADR-0029. Explicit-digit-count format is also exposed.

Both directions are differential-tested against MPFR's
`mpfr_set_str` and `mpfr_sprintf`.

## Transcendentals

References: HFA chapters 9–12, MCA §4.

### Range reduction

Trig functions: Cody-Waite for arguments where `|x| < 2^k` for some
small `k`; Payne-Hanek for arbitrary magnitude. Payne-Hanek requires
a high-precision representation of `2/π` covering enough bits to
absorb the worst-case argument; ferrodec ships such a table for
Decimal128 and pfloat will do the same for binary at user-selected
precision.

`exp` family: reduce `x = k · ln(2) + r`, evaluate the kernel on
`r` near zero, scale by `2^k`. Standard.

### Ziv's strategy

Compute the function at the target precision plus a guard. Detect
the rounding-boundary case (the result lies within an ULP of a
representable value at the target precision). If unsafe, double the
guard and retry. Cap the iteration count at a value documented in
the threat model (matches MPFR's posture).

For tabulated precisions where Lefèvre–Muller worst-case bounds are
known (binary64, binary128 and a few others), pre-cache the
guard-bit count so the first iteration succeeds. For arbitrary
precision, no such cache exists; the iteration cap is the honest
caveat.

### Coverage

Phase 3 (elementary): `exp`, `expm1`, `exp2`, `exp10`,
`ln`, `log1p`, `log2`, `log10`, `sin`, `cos`, `tan`, `cot`, `sec`,
`csc`, `asin`, `acos`, `atan`, `atan2`, `sinh`, `cosh`, `tanh`,
`asinh`, `acosh`, `atanh`, `pow`.

## Special functions

References: NIST DLMF, *Numerical Recipes* §6.

Each family dispatches by argument region.

- `gamma`, `lgamma`, `digamma`, `beta`. Stirling's series for large
  positive argument, Lanczos approximation in the transition region,
  reflection formula `Γ(x)Γ(1−x) = π/sin(πx)` for negative argument.
  `lgamma` carries its own dispatch to avoid overflow.
- `erf`, `erfc`. Taylor series near zero; asymptotic expansion for
  large argument; continued fraction in the middle. `erfc` has a
  separate path to avoid catastrophic cancellation in the tail.
- Bessel `J0/J1/Jn` (slice 6o, `bessel` feature; ADR-0023, shipped).
  Three regimes dispatch on the binary exponent of `|x|`: the
  convergent Maclaurin series for tiny `|x|` (DLMF 10.2.2), Miller
  backward recurrence normalized by the sum rule
  `1 = J0 + 2(J2+J4+⋯)` for moderate `|x|` (DLMF 10.6.1, 10.12.4),
  and the Hankel asymptotic summed to its smallest term for large
  `|x|` (DLMF 10.17.3). One descent yields every order, the
  recurrence and normalization machinery `Y`/`I`/`K` reuse. NaN
  propagates, `J0(±0)=1`, `Jn(±0)=0` for `n≠0`, `Jn(±∞)=+0` by the
  decaying-envelope convention. Bit-exact against MPFR `j0/j1/jn`
  under NearestEven.
- Bessel `Y0/Y1/Yn` (slice 6p, `bessel` feature; ADR-0024, shipped).
  `Y` is real only for `x>0`. The base pair `Y0`/`Y1` dispatches on
  the binary exponent of `x`, sharing the `J` asymptotic threshold:
  the DLMF 10.8.1 logarithmic series (digamma reduced to harmonic
  sums plus the in-tree `γ`, the `J_n` piece from slice 6o) below
  the cut, the DLMF 10.17.4 Hankel asymptotic above (the same
  `a_k(n)` coefficients as `J`, ADR-0023, with `Y`'s trig
  combination). `Yn` for `n≥2` climbs by the stable upward
  recurrence `Y_{k+1}=(2k/x)Y_k−Y_{k−1}` (DLMF 10.6.1), `Y` being
  the dominant solution (the opposite of `J`'s Miller descent; no
  recessive-normalization analog, so no cheap middle regime). Pole
  `Yn(+0)=−∞` raising `DIV_BY_ZERO`; `x<0`/`−0`/`−∞` give `NaN` +
  `INVALID` (complex, the `Ci`/`li` convention); `Yn(+∞)=+0`. The
  J/Y Wronskian `J_{n+1}Y_n−J_n Y_{n+1}=2/(πx)` (DLMF 10.5.2) is
  the re-enabled cross-tie. Bit-exact against MPFR `y0/y1/yn` under
  NearestEven.
- Modified Bessel `I0/I1/In`, `K0/K1/Kn` (slice 6q, `bessel`
  feature; ADR-0025, shipped). `I` is entire (the `J` domain shape);
  `K` is real only for `x>0` (the `Y`/`Ci`/`li` shape). `I` is the
  recessive solution in order, so it reuses the `J` template:
  three-regime dispatch on the binary exponent of `|x|` — the
  DLMF 10.25.2 all-positive Maclaurin (no cancellation, no
  `x·log₂e` boost) for tiny `|x|`, Miller backward recurrence
  normalized by the DLMF 10.35.5 sum rule `eˣ=I0+2·Σ_{k≥1}Iₖ`
  (every order, not `J`'s even-only rule) for moderate `|x|`, the
  DLMF 10.40.1 asymptotic above. `K` is the dominant solution, so it
  reuses the `Y` template: the DLMF 10.31.1 logarithmic series
  (digamma reduced to harmonic sums plus the in-tree `γ`, the `I_n`
  piece from `bessel_i`) below the cut, the DLMF 10.40.2 asymptotic
  above; `Kn` for `n≥2` climbs by upward recurrence
  `K_{k+1}=(2k/x)Kₖ+K_{k−1}` (DLMF 10.29.1 with the §10.25(ii)
  `e^{νπi}` convention flipping `K`'s sign relative to `I` — a
  derive-don't-recall catch, ADR-0025). The asymptotics reuse `J`'s
  `a_k(n)` (ADR-0023) with no trig (`I` alternating, `K`
  all-positive). Pole `Kn(+0)=+∞` (positive, opposite `Y`) raising
  `DIV_BY_ZERO`; `x<0`/`−0`/`−∞` give `NaN`+`INVALID` (complex);
  `Kn(+∞)=+0` is a genuine exponential-decay limit (not the
  decaying-envelope convention); `In(±∞)=±∞` a genuine infinite
  limit. Order parity is even with no sign (`I₋ₙ=Iₙ`, `K₋ₙ=Kₙ`).
  No MPFR `I`/`K` primitive exists, so the oracle is the Airy-style
  tiered one: an `mpmath` reference table (`p≤256`), the DLMF
  10.28.2 I/K cross-tie `I_{ν+1}K_ν+I_νK_{ν+1}=1/x`, and dyadic
  self-consistency, with `p=1024` pinned by in-module unit tests.
- Riemann zeta, real argument (slice 6r, `zeta` feature;
  ADR-0026, shipped). For `s>0` the Borwein / Cohen-Villegas-Zagier
  acceleration of the Dirichlet eta series (DLMF 25.2.3), whose
  weights carry a rigorous `2·(3+√8)^{−n}` error bound; for `s<0`
  the functional equation DLMF 25.4.2
  `ζ(s)=2·(2π)^{s−1}·sin(πs/2)·Γ(1−s)·ζ(1−s)` reflecting into the
  Borwein region. The roadmap originally specified Euler-Maclaurin
  reusing the `gamma_stirling` Bernoulli table; that was found to
  cap accuracy near 90 bits (`|B_{2k}/(2k)!|≈2/(2π)^{2k}`), short
  of the bit-exact `p=1024` lane, so the algorithm was changed and
  the rationale recorded in ADR-0026. `ζ(1)` is a pole
  (`+∞`+`DIV_BY_ZERO`); `ζ(0)=−1/2` and the trivial zeros
  `ζ(−2n)=0` are exact; `ζ(+∞)=1`; `ζ(−∞)` gives `NaN`+`INVALID`
  (an unbounded non-converging oscillation, explicitly not the
  decaying-envelope convention). Bit-exact against MPFR `zeta`
  under NearestEven across all precisions including `p=1024`.
  Complex arguments deferred (Riemann-Siegel out of scope for 1.0).
- Exponential integral `Ei`, sine and cosine integrals `Si`, `Ci`,
  logarithmic integral `li` (slice 6m, `integrals` feature). Each
  dispatches on the binary exponent of the argument like `erf`:
  convergent power series (DLMF 6.6.2/6.6.5/6.6.6) for small
  argument, divergent asymptotic summed to its smallest term
  (DLMF 6.12) for large. `Si`/`Ci` share the asymptotic auxiliaries
  `f`, `g`. `Ei`/`Ci`/`li` use the Euler-Mascheroni constant γ
  (slice 6m0, Brent-McMillan; ADR-0018). `li(x) = Ei(ln x)`.
  Real-only domain conventions: `Ei` for all `x ≠ 0` with
  `Ei(0) = −∞`; `Si` entire and odd; `Ci` and `li` require `x > 0`
  (`Ci(x<0)` and `li(x<0)` are complex, returning NaN + INVALID),
  with poles `Ci(0) = −∞` and `li(1) = −∞` and `li(0) = 0`. MPFR
  has only `eint`, so `Ei` differentials against MPFR directly,
  `li` via `eint(ln x)` for `x > 1`, and `Si`/`Ci` against a
  checked-in authoritative reference table plus self-consistency.
  ADR-0019 records the design.
- Airy `Ai`, `Bi`, `Ai′`, `Bi′` (slice 6n, shipped). Maclaurin
  series for small argument; sign-aware asymptotic for large (DLMF
  9.7, the optimally-truncated error scales as `e^{−2√ζ}`); the
  Wronskian `Ai·Bi′ − Ai′·Bi = 1/π` is the cross-check.
  ADR-0021 records the design.
- Arithmetic-geometric mean. Plain iteration; quadratic convergence.

## Verification

Four layers, each non-redundant.

### Conformance corpora

- Berkeley TestFloat-style vectors at multiple precisions for
  arithmetic.
- Lefèvre–Muller hardest-to-round arguments for transcendentals at
  binary64, binary128, and the precisions where their tables exist.
- DLMF reference values for special functions.

The corpus runs as integration tests. `cargo test --test conformance`
must pass at every commit.

### Differential testing (`gmp-mpfr-sys`)

A separate CI lane (`differential-mpfr` feature) runs random-input
differentials against MPFR for every primitive operation, every
transcendental, and every special function. Disagreements are bugs
in pfloat or in MPFR; either way they get surfaced.

The differential lane is Linux-only (CI cost: one extra runner). The
default lanes (Linux, macOS, embedded cross-compile) stay pure Rust.

ADR-0008.

### Kani harnesses

Phase 5 lands the harness layout copy-pasted from ferrodec, then
adapted. Initial properties:

- No panic on bounded-precision inputs for `+`, `−`, `×`, `÷`,
  `sqrt`, `fma`.
- Rounding direction is correct under each mode for fixed small
  precisions (Kani can tractably enumerate the relevant input
  domain at low PREC values via const-generic instantiation).
- Sign-of-zero correctness across all operations.
- NaN propagation matches IEEE 754-2019 §6.2.

ADR-0009.

### Fuzzing

`cargo-fuzz` harness for every parser entry. OSS-Fuzz integration
once the surface stabilizes.

## Feature gating and no_std

Build profiles that CI exercises:

| profile | flags | description |
|---|---|---|
| default | `--features=std,fmt,big` | desktop / server, dynamic precision |
| pure no_std | `--no-default-features --features=fixed` | embedded, fixed precision only |
| no_std + alloc | `--no-default-features --features=alloc,big` | embedded with allocator |
| full | `--all-features` | every flag, used for Kani and clippy |
| differential | `--features=...,differential-mpfr` | Linux differential lane |

Embedded targets (`thumbv6m-none-eabi`) build in CI under the pure
no_std and no_std-plus-alloc profiles.

## Phase plan

Realistic single-person calendar at the v1.0 = MPFR-equivalent bar:
six to nine months. Each phase ends with a tag that encodes its
exit criterion in CI.

| phase | scope | rough weeks |
|---|---|---|
| 0 | spec acquisition, ADRs, scaffolding | 1–2 |
| 1 | mantissa core, +, −, ×, ÷, √, fma, all rounding modes | 4–6 |
| 2 | string I/O, conformance harness | 2 |
| 3 | elementary transcendentals (exp, log, trig, hyperbolic, pow) | 6–8 |
| 4 | tier-1 specials (gamma family, erf family) | 4 |
| 5 | Kani harnesses, fuzz harnesses, OSS-Fuzz | 2 |
| 6 | tier-2 specials (Bessel, zeta, Ei/Si/Ci, Airy, AGM) | 6–8 |
| 7 | performance tuning, threshold calibration | 2–4 |
| 8 | docs, README conformance evidence, 1.0 tag | 1 |

Each phase merges to `main` via a signed commit; each ADR lands
in the same merge as the code that ratifies it.

Phases 5 and 6 swapped from the originally-tabulated order:
the verification surface (Kani, fuzz, OSS-Fuzz) landed before
the tier-2 specials. Git history reflects this directly. The
slices that built out verification are labelled `slice-6a`
through `slice-6k` (the leading `6` predates the renumbering
and is preserved verbatim as the git-historical record); the
ADRs they ratify (0012, 0013, 0014) carry status-update
sections that read against the renumbered phase plan.

## Caveats and open questions

- Ziv's strategy at unbounded precision is unproven to terminate
  in pathological cases. MPFR has the same caveat. Slice 7c
  (ADR-0022) implements the strategy for `pow` via the interval
  test, with the iteration cap fixed at 5 and stated in the
  `pow_round` doc comment; on the measure-zero exact-tie inputs that
  exhaust the cap the result may be 1 ULP off in directed modes.
  Kernels still on the fixed 64-bit guard (`exp`, `ln`, `sin`, …)
  carry the original caveat until a later slice extends the driver.
- `pow(x, y)` is correctly rounded under every IEEE rounding mode
  (subject to the Ziv cap above): an exact integer `y` takes a
  square-and-multiply fast path, every other case evaluates
  `exp(y · ln(x))` at working precision. It is the first
  transcendental off the NearestEven-only differential tier.
- Performance vs MPFR will not fully close in 1.0. MPFR carries
  decades of hand-tuned assembly via GMP. The target is "documented
  gap, never absurd"; principles forbid reaching for FFI to close
  the last percentages.
- `no_std`-without-`alloc` thread-safety of exception flags. The
  passed-context `Status` form is correct by construction but
  noisier in the API. Decided in ADR-0007; revisit if the noise
  hurts adoption.
- `property_jn`'s `self_consistent` check reconstructs its argument
  at two precisions and can spuriously fail when that argument is
  not a power of two: near a zero of `J_n` the amplification
  `|f'/f|` turns the small reconstruction mismatch between
  precision `p` and `p + 96` into a visible divergence.
  `property_yn` carried the identical shape and was given a dyadic
  argument at slice 6p.7. `property_jn` has not tripped in CI and
  ships as is for 1.0 (pf-ok9); the fix, when it lands, is the same
  one-line change to a power-of-two denominator.

## References

- IEEE 754-2019, *Standard for Floating-Point Arithmetic*.
- Brent, R. P., and Zimmermann, P. *Modern Computer Arithmetic*. Cambridge University Press. Free PDF at https://maths-people.anu.edu.au/~brent/pub/pub226.html.
- Muller, J.-M. et al. *Handbook of Floating-Point Arithmetic*, 2nd ed. Birkhäuser.
- Fousse, L. et al. "MPFR: A Multiple-Precision Binary Floating-Point Library With Correct Rounding." *ACM TOMS* 33:2 (2007).
- Lefèvre, V., and Muller, J.-M. "Worst Cases for Correct Rounding of the Elementary Functions in Double Precision." *ARITH-15* (2001) and follow-ups.
- NIST Digital Library of Mathematical Functions, https://dlmf.nist.gov.
