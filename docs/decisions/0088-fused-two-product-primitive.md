# ADR-0088: the fused two-product primitive and its single-rounding proof

- **Status**: accepted
- **Date**: 2026-06-08

## Context

Slice C1 of the pfloat-complex node lands the one pfloat-core primitive the
complex kernels need that core did not already have. A complex product
`(a + bi)(c + di)` has real part `a·c − b·d` and imaginary part
`a·d + b·c`; complex division forms `a·c + b·d`, `b·c − a·d`, and the
denominator `c² + d²`. Every one of these is a fused two-product expression
`a·b ± c·d`, and a naive `mul` then `mul` then `add`/`sub` is not correctly
rounded (each `mul` rounds, so the difference inherits two roundings and, on
near-cancellation, can be wrong in the last place).

The Phase 4 scoping draft flagged two proof obligations as RECALLED-tier,
to be discharged before the complex node depends on them. This ADR records
the primitive and discharges both, and both resolve differently than the
draft framed them, in pfloat's favor. The discipline is the project's own:
derive from the primitives the implementation actually uses, do not recall.

## Decision

Add the fused two-product primitive to pfloat core as a pure additive
surface, and adopt the two proofs below as its correctness rationale.

### The primitive

`BigFloat::mul_add_mul(a, b, c, d, mode) = round(a·b + c·d)` and
`BigFloat::mul_sub_mul(a, b, c, d, mode) = round(a·b − c·d)`, each with the
`fma` four-variant shape (base, `_round` with explicit precision,
`_with_flags`, `_round_with_flags`) and a `FixedFloat` delegation. The
default result precision is `max(a, b, c, d)`. The surface is additive over
the frozen v1.0 API (ADR-0054), so it is a 1.x minor-version addition, not a
break.

### Proof 1: the primitive is correctly rounded, with one rounding and no Ziv loop

The kernel forms the exact product `c·d`, negates it for the difference
form, and hands it to `fma` as the addend:

```text
    cd_exact = c.mul_round(d, c.precision + d.precision)   // exact
    result   = a.fma_round(b, ±cd_exact, target, mode)     // one rounding
```

`cd_exact` is **exact**: the product of a `p`-bit significand and a `q`-bit
significand is an integer of at most `p + q` bits, so rounding it to
`p + q` bits is the identity (no INEXACT). `fma`'s contract (ADR-0054 method
conventions, IEEE 754-2019 §9.4) is a single rounding of the exact
`a·b + addend`: it forms the exact mantissa product `a·b` via the integer
`multiply_limbs` path and routes `(exact a·b) + addend` through `add_round`,
which rounds the exact exponent-aligned sum once. Composing the two, the
result is `round(a·b + c·d)` under `mode`, correctly rounded.

Catastrophic cancellation (`a·b ≈ ∓c·d`) does not threaten this, and this is
where the proof departs from the draft, which anticipated a "single Ziv
round" that "cancellation can demand more guard than 2p bits" to justify.
There is no Ziv round and no guard band. The exact `a·b + c·d` is a single
real value with a bounded representation: each product is finite precision,
so their exponent-aligned sum is exact in the arbitrary-precision
significand, and rounding an exact value once is correctly rounded by
definition. Cancellation reduces the number of significant bits in the exact
result; it introduces no error, because there is no rounding before the
final one for cancellation to corrupt. The draft's worry applies to a
*different* construction (products rounded to `2p` then subtracted at `2p`),
which pfloat does not use.

This makes the primitive's correctness a composition theorem over two
already-verified kernels (`mul`, `fma`), not a fresh numerical claim, and it
needs no Ziv driver, no error-guard calibration, and no differential lane to
establish (the unit tests, including exact and catastrophic cancellation and
a single-rounding-beats-round-then-subtract case, are confirmation, not the
argument).

### Proof 2: `hypot` is already correctly rounded by its Ziv kernel, not by a naive exact-square route

The complex magnitude `abs(a + bi) = hypot(a, b)` reuses pfloat's existing
`hypot` (ADR-0032, ADR-0056), so the complex node depends on `hypot` being
correctly rounded. The draft proposed proving this via "exact `a² + b²` then
one correctly-rounded `sqrt`." That construction is precisely the one
`hypot`'s implementation **rejects**: its doc comment records that squaring
at the target precision loses half the input precision, so
`sqrt(round(a²) + round(b²))` is not correctly rounded on hard-to-round
inputs. `hypot` instead evaluates `sqrt(a² + b²)` inside one `ziv_round`
closure at an inflated working precision and rounds once at the target.

The reason `hypot` needs the Ziv driver where `mul_add_mul` does not is
structural and worth stating: the sum of squares `a² + b²` is exact (it is a
two-product sum, exactly as in Proof 1), but its **square root is
irrational**, so it has no finite exact representation and the final
rounding cannot be read off an exact value. The Ziv interval test resolves
which representable value the irrational root rounds to. So `hypot` is
correctly rounded by its existing Ziv construction, certified by its oracle
status row (`CR` in `docs/rounding-status.md`). C1 *affirms* this and
depends on it; it does not re-derive `hypot` from the rejected route.

### Status fidelity

The returned `Status` ORs the `fma` status with the product status, so an
`INVALID` from a `0 × ∞` product or a NaN operand in `(c, d)` reaches the
result; the product is exact, so it never contributes a spurious `INEXACT`.
The complex-specific NaN and infinity conventions of C99 Annex G are layered
above this primitive in the complex kernels (slice C3/C4), not here.

## Consequences

- The complex `mul`/`div` kernels (C3) compose from `mul_add_mul` /
  `mul_sub_mul` and inherit correct rounding per component without a Ziv
  loop, error-guard calibration, or a new oracle lane: the correctness is a
  composition theorem over `mul`, `fma`, and (for `abs`) `hypot`.
- pfloat core gains a small additive surface; the v1.0 freeze (ADR-0054) is
  not touched (additive, 1.x).
- Two RECALLED-tier claims from the scoping draft are now grounded and
  corrected: the `ac − bd` "single Ziv round" framing is replaced by the
  stronger no-Ziv composition argument, and the `hypot` "exact squares then
  one sqrt" route is replaced by an affirmation of the existing Ziv kernel.

## Related

- ADR-0054: the v1.0 public API freeze (the `fma` / `mul` conventions this
  composes; additive surface stays inside the freeze).
- ADR-0056: `hypot` as a direct primary kernel (the kernel Proof 2 affirms).
- ADR-0089 (forthcoming): the `pfloat-complex` crate that consumes this.
- Plan: `~/.claude/plans/plan-tower-expansion-scope-goofy-raven.md` (slice
  C1) and `~/.claude/plans/pfloat-phase4-scoping.md` (the two RECALLED-tier
  obligations this discharges).
