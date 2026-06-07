# ADR-0077: ball arithmetic and the directed-pair radius-soundness law

- **Status**: accepted
- **Date**: 2026-06-07

## Context

This is the first slice that produces sound enclosures: `add`, `sub`,
`mul`, `div`, `sqrt` on `Ball<T>`. The midpoint is the existing
correctly-rounded pfloat kernel, unchanged; the work is the radius, and
the radius is where soundness is won or lost. An under-estimating radius
on any path turns every downstream enclosure into a falsehood, so the
radius-soundness law for this route must be written down and proven, not
left as a comment.

The scoping doc fixed the v1.0 route as the **directed-pair** route and
warned against conflating it with the (deferred) surfaced-half-width
route, which carries a different, stronger invariant.

## Decision

For each binary operation `f`, compute the midpoint as the
round-to-nearest pfloat result and bound the radius by

```text
    rad = rad_mid + propagated_input_error
```

- **`rad_mid`** is the midpoint's own rounding error, read off the
  directed pair: `lo = a.mid op↓ b.mid`, `hi = a.mid op↑ b.mid`, and
  `rad_mid = (hi − lo)/2` (an exact halving of the exact, small spread).
  Because nearest rounding lands on the closer of the two bracketing
  representables, `|mid − f(a.mid, b.mid)| ≤ (hi − lo)/2`.
- **`propagated_input_error`** bounds `|f(x, y) − f(a.mid, b.mid)|` over
  the input balls, per operation:
  - add / sub: `ra + rb`.
  - mul: `|a.mid|·rb + |b.mid|·ra + ra·rb`.
  - div: `(|b.mid|·ra + |a.mid|·rb) / (blo·|b.mid|)` where
    `blo = |b.mid| − rb` is a **lower** bound on `|y|`. If `blo ≤ 0` the
    divisor ball contains zero and the result is the entire real line
    with `DIV_BY_ZERO`.

Every radius scalar operation rounds **outward**: numerators and the
running radius toward `+∞`, the division denominator (`blo` and
`blo·|b.mid|`) toward `−∞`, so the quotient is an upper bound. The final
radius narrows up to `Mag` once at the end (the scoping doc's
narrow-upward-once strategy). `sqrt` is unary and monotonic, so it is
enclosed by evaluating the kernel at the outward-rounded input-interval
endpoints and building the ball from those; a ball dipping below zero is
sound over `[max(0, lo), hi]` with `INVALID`.

**The soundness law (route-specific):** with `mid` the nearest result
and the radius as above,

```text
    |f(x,y) − mid| ≤ |f(x,y) − f(a.mid,b.mid)| + |f(a.mid,b.mid) − mid|
                   ≤ propagated_input_error    + rad_mid    = rad
```

so `f(x, y) ∈ [mid − rad, mid + rad]` for every `x ∈ [a]`, `y ∈ [b]`
(FTIA). Soundness rests only on the directed kernels' correctness and the
`Mag` round-up; **there is no Ziv-residual term**, because this route
never calls the Ziv driver. This must not be conflated with the
surfaced-half-width route, whose invariant *does* add the Ziv residual.

Exactness in produces exactness out: when the directed pair coincides and
the inputs are exact, every term is zero, so `rad = 0`.

When an operation leaves the reals (division by a zero-containing ball,
sqrt of a wholly-negative ball) the result is the entire real line plus
the IEEE status flag, never a finite under-covering ball. `Status`
composes through the OR-monoid; on a ball `INEXACT` is the normal correct
outcome, so the radius is the primary accuracy channel.

## Consequences

- The enclosure law has a single written home and a one-line proof from
  directed-kernel correctness. The "no Ziv-residual term" caveat is
  recorded next to the route it applies to.
- The radius is computed in scalar arithmetic at midpoint precision and
  narrowed to `Mag` once, so the only `Mag`-resolution loss is a single
  final upward rounding.
- The law is property-tested by a randomized FTIA point-containment lane
  (sample points across the input balls, assert the high-precision true
  result is enclosed) and a corner-range lane (assert the ball contains
  the exact monotone-extreme range), plus boundary lenses
  (zero-straddling divisor, sqrt domain edges, entire propagation). The
  independent Arb containment backstop lands with the verification slice.
- `mul`/`div` are looser than the tightest possible enclosure (the
  propagation constants over-cover slightly), the standard ball tradeoff;
  tightness is logged, soundness is non-negotiable.

## Related

- Plan: `~/.claude/plans/melodic-knitting-ripple.md` (slice 6); design
  reference `~/.claude/plans/pfloat-phase4-scoping.md` (the directed-pair
  vs surfaced-half-width distinction, the highest-severity radius hazard).
- Beads: `pf-icgj.6` (under epic `pf-icgj`).
- Other ADRs: consumes `Mag` round-up (ADR-0074), `RealScalar` directed
  arithmetic (ADR-0075), and the `Ball` conversion boundary (ADR-0076);
  rests on the pfloat directed-kernel correctness
  (`feedback_directed_pair_for_libm`).
