# ADR-0063: INEXACT flag fidelity across the rest of the elementary transcendental surface

Status: Accepted (2026-06-04)

## Context

ADR-0060 corrected the IEEE 754-2019 §7.6 `INEXACT` flag for the exp/log
family and sin/cos: each transcendental kernel dispatches its decidable
exact-input set to `Status::OK` before the Ziv loop and forces `INEXACT`
on the transcendental fall-through, because by Lindemann-Weierstrass /
Gelfond-Schneider the result outside that set is irrational and so is
inexact even when it rounds onto a grid value. ADR-0060 scoped itself to
exp/log + sin/cos and filed the rest of the surface as pf-uqd1. This ADR
discharges pf-uqd1 for the elementary transcendentals.

## Decision

Apply the ADR-0060 pattern to `tan`, `cot`, `sec`, `csc`, `asin`, `acos`,
`atan`, `atan2`, `sinh`, `cosh`, `tanh`, `asinh`, `acosh`, `atanh`,
`expm1`, `log1p`, `erf`, and `erfc`: force `INEXACT` on the finite normal
fall-through (guarded on a `Class::Normal` result, so domain `qNaN` /
`INVALID`, the poles' `±∞` / `DIV_BY_ZERO`, and the exact non-finite
limits keep their status). The reciprocal trig kernels share
`reciprocal_via_ziv`, so one edit there covers `cot`/`sec`/`csc`.

### The audit was the load-bearing step

The force is sound only if every exact-input case is dispatched to `OK`
before the fall-through; forcing `INEXACT` on an undispatched exact result
would be a wrong over-report. A per-kernel audit confirmed every exact
case is already dispatched, and surfaced two facts that meant no
special-case arm needed changing:

- **The irrational-constant special cases already report INEXACT.**
  `asin(±1) = ±π/2`, `acos(0) = π/2`, `acos(−1) = π`, `atan(±∞) = ±π/2`,
  and the `atan2` axis cases return their constant through `pi_at_round` /
  `pi_over_2_at_round`, which compute the constant at `target + 128` and
  round to target. Rounding an irrational sets `INEXACT`, so those arms
  were already correct; they are not exact-input cases.

- **The tiny-argument fast paths already report INEXACT.** The ADR-0059
  small-argument paths (`sinh`/`tanh`/`asinh`/`atanh`/`expm1`/`log1p`)
  return through `round_with_infinitesimal`, whose final
  `round_to_precision` over the injected residue sets `INEXACT`.

So the genuinely exact inputs, all dispatched before the fall-through, are
the rational ones: `tan(0)=0`, `sec(0)=1`, `sinh(0)=0`, `cosh(0)=1`,
`tanh(0)=0`, `tanh(±∞)=±1`, `asinh(0)=0`, `acosh(1)=0`, `atanh(0)=0`,
`expm1(0)=0`, `expm1(−∞)=−1`, `log1p(0)=0`, `asin(0)=0`, `acos(1)=0`,
`atan(0)=0`, `atan2(+0, x>0)=0`, `erf(0)=0`, `erf(±∞)=±1`, `erfc(0)=1`,
`erfc(±∞)=0/2`. Everything else is transcendental.

### Scope

This covers the elementary transcendentals and `erf`/`erfc`. The heavier
special functions (`gamma`, `lgamma`, `digamma`, `beta`, `zeta`, the Airy
and Bessel families, the exponential/trigonometric integrals) are
deferred: their exact-input sets are subtler (factorial integer points and
poles for `gamma`, the trivial zeros and negative-integer rationals for
`zeta`, and so on), and a mistabulated exact set would clear `INEXACT`
wrongly. They are filed as a follow-up so the per-function exact-set
derivation gets its own careful pass.

### libm gate

`pfloat-libm`'s `inexact_is_gated` widens to the v0.1 libm surface members
of this set: `tan`, `cot`, `sec`, `csc`, `asin`, `acos`, `atan`, `sinh`,
`cosh`, `tanh`, `asinh`, `acosh`, `atanh`, `expm1`, `log1p`. (`atan2`,
`erf`, `erfc` are pfloat kernels not in the v0.1 libm surface.) The
expectation stays the enclosure-derived one from ADR-0060.

Gating these exposed and fixed a latent under-report in the libm
saturation fast path (`saturate.rs`): `tanh(±huge)` and `expm1(−huge)`
saturate to `±1`, which the fast path returned with `Status::OK`,
deliberately mirroring the kernel's pre-correction behavior (it cited
pf-njs5 as a deferred question). The corrected kernel now forces
`INEXACT` there — the true value is strictly inside the saturation limit —
so `pos_one`/`neg_one` return `INEXACT` to stay bit-for-bit identical to
the kernel, which the crate's `fast_path_matches_kernel` test reconfirms.

## Consequences

- `INEXACT` is reliable across the elementary transcendental surface, and
  the libm differential gates it for that surface. No value changes; the
  fix is metadata, and the oracle status is untouched.
- The reciprocal-trig edit in `reciprocal_via_ziv` is the only shared
  site; the rest are one local edit per kernel, all the same shape as
  sin/cos.
- The special-function exact-input tables remain the open follow-up.

## Related

- pf-uqd1 (this work); ADR-0060 (the pattern and the exp/log + sin/cos
  discharge); ADR-0059 (the tiny-argument fast paths that already carry
  INEXACT); the slice-p1.25 mode-aware irrational-constant returns.
- Baker, *Transcendental Number Theory*, Cambridge University Press, 1975.
