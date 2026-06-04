# ADR-0060: INEXACT flag fidelity for the transcendental exp/log and sin/cos kernels

Status: Accepted (2026-06-04)

## Context

IEEE 754-2019 §7.6 raises `INEXACT` exactly when the delivered result
differs from the exact result. pfloat's kernels carried two flaws against
that definition, both surfaced by the pfloat-libm shell (pf-lm2, pf-njs5).
The value was correctly rounded in every case; only the flag was wrong.

**Over-report.** `log10(1000) = 3`, `exp10(2) = 100`, `exp2(10) = 1024`
are exact, yet each kernel reached the answer through a composed
transcendental that rounds internally and set `INEXACT`. `log2` of a power
of two already cleared the flag through a pre-Ziv dispatch (ADR-0039); the
sibling kernels had none.

**Under-report.** `exp(2^-1074)` at a working precision below 1075 bits
returned `1.0` with `INEXACT` clear, because the true residual `2^-1074`
fell below the kernel's working precision so the rounding layer never
observed a rounding. The IEEE-correct flag is `INEXACT` set: `1 + 2^-1074`
is not `1`.

The bead framed the resolution as the table-maker's dilemma, needing the
kernel to grow precision until a residual becomes observable. It does not.
Transcendence settles the flag a priori.

## Decision

### 1. Transcendence makes INEXACT unconditional outside a decidable set

For each genuinely transcendental kernel the set of dyadic inputs whose
true result is exactly representable is decidable and bounded:

- `exp(x)` exact only at `x = 0`; `ln(x)` only at `x = 1`.
- `exp2(x)` exact iff `x` is an integer (then `2^x` is a single-bit power
  of two, exact at every precision).
- `exp10(x)` exact iff `x` is a non-negative integer and `10^x = 5^x·2^x`
  fits the target precision.
- `log2(x)` exact iff `x = 2^k`; `log10(x)` exact iff `x = 10^k`.
- `sin(x)`, `cos(x)` exact only at `x = 0`.

By Lindemann–Weierstrass (`e^α` transcendental for nonzero algebraic `α`;
`ln α` transcendental for algebraic `α ∉ {0, 1}`) and Gelfond–Schneider
(`2^α` transcendental for irrational algebraic `α`), every input outside
its kernel's exact set yields an irrational, hence non-representable,
result. `INEXACT` is therefore correct there without computing anything.
See Baker, *Transcendental Number Theory* (1975), Ch. 1-2.

### 2. Exact-input dispatch, then force on the fall-through

Each kernel gains a pre-Ziv dispatch that returns the exact value with
`Status::OK` (closing the over-report), and forces `INEXACT` on the
transcendental fall-through (closing the under-report). The detectors are
**sound**: a returned exact value is always exact, so the flag is never
wrongly cleared. Completeness is bounded by precision; a missed exact
input only leaves an over-report in place, never a wrong clear.

`exp2` reuses `integer_exponent` and constructs `2^k` directly. `exp10`
and `log10` share `ten_pow_if_fits`, which forms `10^k` by exact integer
exponentiation and accepts it only when its significant-bit count fits the
target. `log10`'s `power_of_ten_exponent` inverts the binary exponent to a
small candidate `k` and confirms by exact comparison. The kernels whose
non-finite cases flow through the composition (`log2`, `log10`, `sin`,
`cos`) force only on a finite normal result, so domain `qNaN + INVALID`
and pole `±∞ + DIV_BY_ZERO` results keep their status untouched.

### 3. The libm gate now checks INEXACT, against the oracle not the shell

`pfloat-libm`'s status gate hard-checks `INEXACT` for these eight
functions. The expectation is derived from the MPFR enclosure, not the
shell: the result is inexact unless the bracket collapses to a single
oracle-precision point that the committed hardware float reproduces
exactly. That rule is width-generic and also matches the conversion's
`OVERFLOW`/`UNDERFLOW` `INEXACT` at the format boundary, so the now-exact
pfloat values (`exp2(1024)` finite in pfloat's unbounded exponent range)
still gate correctly once rounded to `f32`/`f64`.

### 4. Scope

This ADR covers `exp`, `exp2`, `exp10`, `ln`, `log2`, `log10`, and the
`sin`/`cos` under-report the existing regression guard already asserted.
The same force-INEXACT rule extends to the rest of the trig, hyperbolic,
and special surface; each kernel needs its own exact-input table, deferred
to pf-uqd1.

## Consequences

- pfloat's `INEXACT` is reliable for the eight kernels, and the libm gate
  asserts it. The shell's composed-exact unit tests rise from value-only
  to flag-clear.
- No value changes: oracle `worst_ulp` is untouched; the fix is metadata.
- The cost is a cheap detector per call: one integer extraction, plus a
  bounded exponentiation that runs only on integer or power-of-ten
  candidate inputs.
- `exp2` of an integer now returns `2^k` exactly in pfloat's arbitrary
  exponent range rather than overflowing through `exp`; the hardware
  conversion still raises `OVERFLOW`/`INEXACT` when `2^k` exceeds the
  format.

## Related

- pf-njs5 (this work); pf-uqd1 (trig/hyperbolic/special follow-up).
- ADR-0057 records the directed-pair shell that surfaced the over-report;
  ADR-0039 the `log2` power-of-two dispatch this generalizes.
- Baker, *Transcendental Number Theory*, Cambridge University Press, 1975.
