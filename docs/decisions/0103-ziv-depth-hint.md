# ADR-0103: the Ziv driver takes a lazy input-derived certification depth

- **Status**: accepted
- **Date**: 2026-06-12

## Context

The Ziv driver's guard schedule doubles from 64 to a fixed
`ZIV_GUARD_CAP = 1024` extra bits and then falls back with the last
attempt: a sound design for generic hard-to-round inputs (the
measure-zero caveat MPFR also documents), but wrong for inputs that
*encode* proximity deeper than the cap. pf-jl35's reproducers, both
re-verified bit-exactly by run before this slice:

- `ζ(1 − 2^-2000 @p2001) → 53` under TowardPositive returned
  `−2^2000` where the truth `−2^2000 + γ` demands
  `−(2^2000 − 2^1947)` — 1 ulp low, plain INEXACT, ~16 s burned in
  five futile exhaustion passes (release).
- `cos(RN2048(π)) → 53` under TowardPositive/TowardZero returned `−1`
  where the truth `−1 + 4.227e-1235` demands `nextUp(−1)`. The bead's
  claimed correct value (`−1 + 2^-53`) was itself misadjudicated by a
  binade — the representable just above −1 at p53 is `−(1 − 2^-54)` —
  the same correct-the-record pattern as R1's three reproducers.

The R1/R2.1 arcs fixed the *kernel-internal* side of this family:
ADR-0097/0098 grow kernel resolution to meet input-encoded depth, and
ADR-0102 dispatches around the driver entirely where an infinitesimal
rounding is provable. What remained was the *driver's certification*:
the charged half-width `|y|·2^-(working − error_guard)` is a
per-kernel claim about the eval closure, valid at **any** working
precision — the fixed cap was the only thing stopping a deeper,
legitimate certification. The disposition options the bead carried:
extend the driver with a depth hint, recalibrate the cap story and
document, or per-kernel cap scaling. The hint was chosen: it is the
one shape that also serves the residues R2.1 recorded (pf-fbjn, the
hypot/atan band holes), and per-kernel cap scaling is the same change
expressed with more duplication.

## Decision

**`ziv_round_with_depth(eval, target, mode, error_guard, depth_hint)`**
runs the unchanged legacy schedule; only on exhaustion does it call
`depth_hint()` (lazy — zero cost on every input the legacy schedule
settles), and if the hinted depth exceeds the legacy cap it runs one
further iteration at `target + depth + 64` under the identical
interval test. Certification there is legitimate for the same reason
it is at any other rung. A deep rung that still fails to certify
falls back as before — with the deep attempt as the better fallback,
unless it degenerated to NaN (a deeper working can exceed a kernel's
internal envelope), in which case the legacy attempt stands. The
hint must be input-proportional (the DoS-budget posture); both
adopters' hints are.

**zeta** derives the hint from the ADR-0098 conditioning probe:
`s − 1` Sterbenz-exact at the input's own precision, depth = its
negated exponent. Near the pole `ζ(1 ± ε) ≈ ±1/ε + γ` sits `|e(ε)|`
relative bits from the on-grid `±2^k`.

**The trig reduction family** (sin, cos, tan, cot, sec, csc) derives
the hint from the reduction residual: `reduction_depth_hint` runs
`reduce` once at the pre-check width — the ADR-0098 growth inside
`reduce` resolves the residual's position at any working width — and
returns `2|e_r|`. The cos-shape arms evaluate to `±(1 − r²/2 + …)`,
the sin/tan-shape arms carry their `r³` corrections, and the
reciprocal pole arms their `1/r` series, all parked at that depth.

**The reduction range cap is decoupled from the working precision.**
`reduce`'s up-front check (`e_x + working + 64 < 4096`) predated
ADR-0098's live `two_over_pi_at` and conflated "x too large" with
"working too deep": at deep working precisions it refused even
`e_x = 1`, which surfaced as cos turning NaN + INVALID at the deep
rung (and, independently, had been refusing mid-range exponents at
high *targets*: `sin(2^3000)` at target 1000 was NaN + INVALID). The
cap is now on the input alone (`e_x + 64 < 4096`, i.e. `|x| < 2^4032`
— the documented range policy and the live-computation cost bound);
the lane pins both the deep rung and the formerly-refused
high-target case with a refinement-consistency oracle.

## Consequences

- The lane pins: zeta's named reproducer (TP corrected, NE control)
  in release builds only — `#[cfg(not(debug_assertions))]`, running
  in the MPFR full-union release job at ~1 minute where the debug
  matrix would pay ~15. There is no cheap debug instance of the
  zeta wiring: a shallower depth survives by *uncertified fallback*
  (at depth ≤ ~1080 the capped eval still carries γ, so the
  fall-through rounds correctly — verified at e8b1284 with depth
  1060, which is therefore a control, not a guard); red requires
  the full collapse, and the collapse requires the cost. The trig
  rows guard the shared driver mechanism in the debug matrix at
  millisecond cost,
  cos's named reproducer (TP/TZ = `nextUp(−1)`, NE control), sin's
  exposed cos-shape arm at `RN2048(π)/2` (inward modes = `pred(1)`),
  sec through the shared reciprocal helper (TN = `−(1 + 2^-52)`),
  and `sin(2^3000)` at target 1000 computing instead of refusing.
  All red at the pre-fix commit by run (sin/sec/range verified in a
  worktree at e8b1284), all directions mpmath-pinned first.
- Cost: the legacy burn before the deep rung remains by design (the
  early rungs settle every shallow input; replacing them with a jump
  to the hint would regress inputs whose realized boundary distance
  is far shallower than the encoded depth). The zeta named case pays
  ~61 s release total (the ~16 s legacy burn plus a ~45 s certified
  deep evaluation at ~2100 working bits) — the honest
  input-proportional price, paid only by inputs that encode the
  depth. The trig deep cases are milliseconds (the reduction was
  already input-proportional).
- `ZIV_MAX_ITERS` keeps its meaning (the legacy schedule length);
  the loop now counts it explicitly and the trace/capture shapes are
  unchanged, so the pf-tqzz cross-check contract is untouched.
- Failure modes considered (inverted):
  1. **A wrong hint cannot create a wrong answer** — too small and
     the deep rung is skipped or fails to certify (the legacy
     fallback path, today's behavior); too large and the only cost
     is an over-provisioned evaluation. Certification itself never
     rests on the hint's value, only on the interval test.
  2. **The deep attempt as fallback displaced a sound legacy
     attempt** — found while wiring cos: the deep eval came back
     NaN through the (then working-coupled) reduction cap and the
     first draft of the loop returned it, turning a 1-ulp INEXACT
     into NaN + INVALID. Fixed twice over: the NaN-guarded fallback
     preference, and the range-cap decoupling that removes the NaN
     at its source.
  3. **Diophantine residue**: an input whose truth sits within the
     deep rung's half-width of a boundary still exhausts — 1 ulp,
     INEXACT, the measure-zero caveat one layer deeper. This is the
     honest floor; no effective bound exists without per-value
     irrationality measures.
  3b. **The hint cannot repair a violated error model.** The
     half-width's validity at any working precision is the KERNEL's
     obligation; the slice's adversarial verification found
     `zeta_fe` violating it near the trivial zeros (no conditioning
     on s's proximity to the negative integers: `ζ(−2 − 2^-1200)`
     certified ~296 orders of magnitude wrong at the FIRST rung,
     pre-existing, byte-identical at baseline) — the driver is not
     even exhausted there, so no hint fires. Filed as pf-hkoj with
     the ADR-0098 `pole_proximity_depth` fix shape; recorded here
     because it bounds what this mechanism can ever fix.
  4. **The widened trig range** (mid-range exponents at high
     targets now computed instead of refused) strictly converts
     NaN + INVALID into correct values; no caller could have relied
     on the refusal as a *value*, and the range policy on `|x|`
     itself is unchanged.
- Adoption guidance recorded for the rest of the family: pf-fbjn
  (the ADR-0059 tiny-x fast paths) and the hypot/atan band residues
  from ADR-0102 can now pass `max(2|e|, p) + 64`-shaped lazy hints
  and re-arm their triggers; that lands with pf-fbjn's own lane
  rows, not silently here.

## Related

- Issues: pf-jl35 (closed by this ADR), epic pf-8iji; pf-fbjn (the
  family adoption, next), pf-7nnw (deep-tiny directed modes, Opus
  arc); pf-t6ht/pf-2thy (guard-model probes, adjacent).
- Other ADRs: ADR-0038 (the shared driver), ADR-0097/0098 (the
  kernel-side depth pattern this completes), ADR-0102 (the dispatch
  alternative when an infinitesimal rounding is provable), ADR-0096
  (the certified-floor analogue), ADR-0080 (directed-mode posture).
- References: `docs/references/zeilberger-zudilin-pi-2020.md` (the
  reduction-growth termination bound this slice's trig hints lean
  on, unchanged).
