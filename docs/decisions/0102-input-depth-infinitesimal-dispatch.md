# ADR-0102: input-encoded depth dispatches through the infinitesimal rounding — hypot, atan, atan2

- **Status**: accepted
- **Date**: 2026-06-12

## Context

The R1 verifiers filed two members of the certified-wrong-answer family
(epic pf-8iji) that R2 opens with:

1. **hypot falsely exact** (pf-71u2): `hypot_round(1, 2^-2000, 128,
   TowardPositive)` returned `(exactly 1, Status::OK)` — wrong under
   TowardPositive (the truth `sqrt(1 + 2^-4000) = 1 + 2^-4001 − …` is
   strictly above 1) and falsely exact in every mode. With exponent gap
   `g = e_big − e_small`, the Ziv eval's sum absorbs `small²` whenever
   `2g ≥ working`; the value collapses onto `|big|`, exactly on-grid,
   the interval test never converges, and the `ZIV_MAX_ITERS`
   fall-through returns the collapsed candidate, whose rounding is
   exact — Status OK on a wrong value, observable by every hypot
   caller (pfloat-complex's clog defends itself since ADR-0100; nothing
   else does).
2. **atan2 returns the argument** (pf-e2ow): `atan2_round(2^-d, 1, 64,
   TowardZero/TowardNegative)` for `d ≥ ~545` returned `2^-d` where the
   truth `atan(ε) = ε − ε³/3 + … < ε` demands `pred(2^-d)` under the
   inward modes (INEXACT — flag honest, value wrong). The boundary at
   `d ≈ 545` is exactly where `2d` outruns `target + ZIV_GUARD_CAP`.

Re-deriving the second against the grid located the root one layer
down: **atan itself** has no tiny-x dispatch at all, so
`atan_round(2^-545, 64, TowardZero)` fails identically, at any depth
up to and past every feasible working precision (`atan(2^-2^40)`
returns the argument with INEXACT; this ADR's first draft claimed
"certified zero" there from a Display-truncated read of the failing
test's output — the same misadjudication family as the review's asin
oracle, caught by the slice's verifier). atan2 composes atan, so the
root fix lands in atan and atan2 forwards to it.

The mechanism is ADR-0097/0098/0100's one more time: the input encodes
a proximity (here, of the result to an on-grid value) deeper than any
resolution the driver visits, and `half_width`-based certification
cannot speak about a value that collapses exactly onto the grid. A
precision boost inside the eval closure cannot repair certification
either: the driver charges its half-width against the *driver's*
working precision, so a boosted closure resolves the value but the
charged interval still straddles — the answer then arrives only as an
uncertified cap-exhaustion fallback. The honest dispatch is *before*
the driver.

## Decision

When the input-encoded depth puts the truth strictly inside the
boundary-free zone next to an exactly known base, round the base with
a directional infinitesimal (`round_with_infinitesimal`, the ADR-0059
machinery) instead of entering the driver.

**The boundary-gap argument** (shared by all three kernels). Let the
base `b` be a `p`-bit value with exponent `e`, and let the truth be
`b ± c` with the correction window `c ∈ (2^(e−D−k₁), 2^(e−D+k₂))` for
input-computable depth `D` and small constants. The target's
rounding-change points (grid points for the directed modes, midpoints
for the nearest modes) near `b` lie on a grid that, intersected with
`b`'s own representation, leaves a boundary-free open zone on each
side of `b` of width at least `2^(e − max(p, target) − 1)` (one input
ulp when `b` sits off the change grid; half a change step when on it;
the worst case is the shrink direction from a power-of-two base under
the nearest modes, where the lower binade's grid is finer). The
infinitesimal residue sits at `e − max(p, target) − 2`
(`round_with_infinitesimal` uses `wide = max(p, target) + 3`) — half
the worst-case zone, a one-bit margin, with the verifier's probes
covering the tightest geometries (power-of-two bases both directions,
all-ones bases at binade tops, bases exactly on midpoints). If
`D ≥ max(p, target) + 6`, the correction window tops out at
`2^(e − max − 4)`, a further two bits inside the residue's position,
so correction and residue sit strictly inside the zone on the same
side and `round(b ± ε) = round(truth)` under every mode; INEXACT is
honest (the truth is irrational, or in hypot's case its square's
bit-span exceeds the target).

**The rim guard** (added after the slice's adversarial verification —
see the inversion below). All three triggers additionally require the
base exponent to clear `i64::MIN + max(p, target) + 5` (for atan2,
`e_t_hi ≥ i64::MIN + w2 + 6` with the quotient's `wide = w2 + 3`):
within that reach of the rim, `round_with_infinitesimal`'s residue
placement `e − wide + 1` saturates (pf-a77o's named site,
`src/rounding.rs:447`), `base ± ε` becomes exactly representable, and
the dispatch would certify a wrong value with Status OK. Refused
inputs fall through to the driver — the pre-existing rim behavior —
until pf-a77o fixes the placement at the root, at which point the
guards can be relaxed with rim lane tests of their own.

**hypot** (`src/ops/hypot.rs`): for two finite nonzero operands the
depth is `2g`, with `δ = small²/(hypot + |big|) ∈
(2^(e_b − 2g − 3), 2^(e_b − 2g + 1)]`, always positive. Trigger
`2g ≥ max(p_big, target) + 6`, dispatch
`round_with_infinitesimal(|big|, +, grows)`. Inside the band the
driver is untouched; it certifies wherever its cap clears the depth
(input precision up to roughly `target + 1030`), and the
higher-precision band residue is pf-fbjn/pf-jl35 territory (see the
inversion below) — exact Pythagorean triples keep Status OK
bit-for-bit.

**atan** (`src/math/atan.rs`): the correction `c = |x| − |atan x| ∈
(2^(3e−2), 2^(3e+2))` gives depth `2|e|`. Trigger `e < 0` and
`2|e| ≥ max(p, target) + 6`, dispatch shrink-direction with `x`'s
sign. Unlike the ADR-0059 siblings the trigger carries the
**input-precision arm** from birth (see the inversion below).

**atan2** (`src/math/atan2.rs`): in the right half-plane with
`e_y − e_x ≤ −(target + 5)`, compute the quotient once at
`2·target + 2` bits; if it is **exact**, forward
`atan_round(q, target, mode)` — atan's dispatch is then guaranteed to
fire (`2|e_q| ≥ 2·target + 10 > max(2·target + 2, target) + 6`), and
values, flags, and INEXACT arrive whole. An inexact quotient carries
the truth's grid position in its own expansion (the driver's
fall-through rounds it correctly outside a measure-zero proximity
class) and stays with the driver; a rim-saturated quotient flags
OVERFLOW/UNDERFLOW, fails `is_ok()`, and is never forwarded. `x < 0`
keeps the quadrant shift: the result there is `≈ ±π` with no
tiny-result collapse (its only failure shape needs π itself parked at
a target boundary — the Diophantine-thin Ziv-cap class).

## Consequences

- The regression lane pins all named reproducers plus controls:
  hypot's deep gap mode-by-mode (TowardPositive gets `nextUp(1)`,
  everything else exactly 1, INEXACT always — OK was the defect), the
  2000-bit-target variant, an exact Pythagorean band control
  (`(2^62−1, 2^32) → 2^62+1`, Status OK), atan at depth 545 and at
  depth `2^40` (where the baseline certified **zero**), atan2's named
  reproducer both signs, an exact ratio past every feasible precision,
  and an inexact-ratio control (`atan2(2^-600, 3)`) pinned bit-exact
  against mpmath 1.4.1. Directions pinned with mpmath at 4500–8000
  bits before any code was written.
- Deep exact-ratio atan2 and deep hypot now cost O(target) instead of
  five futile cap-iterations; the certifiable bands are bit-for-bit
  untouched, still burning their retries (a cost note, not a
  correctness one — the depth-hint question below owns it).
- Failure modes considered (inverted):
  1. **A closure-level precision boost was the first design and is
     wrong twice over**: it cannot certify (the charged half-width is
     driver-w-based, so the straddle persists and the value arrives
     only as uncertified fallback), and for exponent-encoded depths
     (`2^40` and beyond) the boost saturates u32 — infeasible by
     memory alone. The infinitesimal dispatch costs O(target) at any
     depth. This inversion is why the dispatch sits before the driver.
  2. **The trigger could fire on a zero/special operand** — the
     dispatches match `Class::Normal` on both operands (hypot) or the
     argument (atan) first; `hypot(x, 0) = |x|` keeps its exact OK.
  3. **The same derivation applied to the existing ADR-0059 family
     found a real shipped defect** (filed as pf-fbjn, P1,
     run-verified on atanh and asinh): the six tiny-x triggers test
     only `e ≤ −(target+2)`, so a high-precision input parked next to
     a change point (`gap < correction`) rounds the wrong side. A
     trigger-arm fix was drafted and deliberately **reverted**: the
     Ziv fall-through collapses the input onto the same midpoint
     (working caps below the input precision) and returns the
     identical wrong side — re-routing, not repair. The genuine fix
     (exact tail-vs-correction comparison, or a driver depth-hint
     that makes a boosted closure certifiable) belongs to the
     pf-jl35 disposition slice, which also owns hypot's and atan's
     analogous band/hole residues (input precision > target + ~975
     with mid-band depth).
  4. **An exact target-precision hypot in the deep zone would make
     INEXACT a lie** — impossible: `D² = big² + small²` forces `D`'s
     bit-span past `g − p`-ish, beyond the target whenever the
     trigger holds.
  5. **The atan2 forward could double-raise or lose flags** — the
     forward returns atan's `(value, status)` whole, and the
     quotient's own rounding never reaches the caller (exactness is
     the gate).
  6. **The slice's adversarial verification refuted the first draft
     at the exponent rim**: without the rim guard, all three
     dispatches composed through `round_with_infinitesimal`'s
     saturating residue placement and returned certified Status-OK
     wrong values within `max(p, target) + 2` of `i64::MIN` —
     `atan(2^(i64::MIN+1))` came back **exactly 0 with OK** where
     the baseline failed loudly (debug panic in the driver's
     half-width, garbage with INEXACT in release). A new silent
     surface of the very family under remediation, converted back
     to the loud pre-existing behavior by the rim guard
     (run-verified at the measured boundary: wrong at
     `i64::MIN + max + 1`, correct from `i64::MIN + max + 2`; the
     guard sits at `+5`). The root and its mode-aware repair stay
     with pf-a77o; the verifier's rim reproducers are recorded
     there.
  7. **Same verifier, the record corrected twice**: the first
     draft's claim that the baseline certified zero at depth `2^40`
     was a Display-truncation misread (it returns the argument,
     INEXACT), and the first draft's prose overstated the
     residue-to-zone margin (three bits claimed, one bit real —
     the code was sound, the argument as written was not).
- atan2's inexact-ratio deep zone keeps today's behavior knowingly:
  value correct outside the thin proximity class, INEXACT honest,
  certification absent. Recorded for pf-jl35 rather than silently
  declared fixed.

## Related

- Issues: pf-71u2, pf-e2ow (closed by this ADR), epic pf-8iji;
  pf-fbjn (the family hole, filed from this slice); pf-jl35 (the
  driver depth-hint disposition, R2.2); pf-7nnw (deep-tiny directed
  modes, Opus arc).
- Other ADRs: ADR-0059 (the infinitesimal machinery and its original
  fast paths), ADR-0097/0098 (input-encoded depth, scalar side),
  ADR-0100 (clog's consumer-side defense, now redundant for this
  defect but kept as depth), ADR-0080 (directed-mode saturation
  posture).
