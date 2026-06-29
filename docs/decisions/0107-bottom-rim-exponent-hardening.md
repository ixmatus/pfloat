# ADR-0107: exponent arithmetic at the bottom rim — the lift, the Taylor floor, and the deep-rung ceiling

- **Status**: accepted
- **Date**: 2026-06-12

## Context

pf-a77o collected the i64::MIN-adjacent exponent-arithmetic sites the
2026-06 review found, and the R2 arc's own verifiers kept adding
reproducers to it:

- **(a) `round_with_infinitesimal`'s residue placement**
  (`src/rounding.rs`): the residue's position `e − wide + 1`
  saturated near the rim, parking the residue at `i64::MIN` instead
  of below the base's tail; `base ± ε` then became exactly
  representable and the rounding certified wrong values **with
  Status OK** through every tiny-x dispatch
  (`sinh(2^(i64::MIN+10))`, `atanh` likewise; the R2.1 verifier found
  the same composition turning the new ADR-0102 dispatches into
  fresh OK-lies — `atan(2^(i64::MIN+1))` returned exact 0 with OK —
  which is why those three dispatches carried interim rim guards).
- **(b) the Taylor loops' saturated terms** (`sin_taylor`,
  `cos_taylor`, `atan_taylor`): for deep-tiny `r` the `r²`/`r³`
  terms' TRUE positions lie below `i64::MIN`, but the mul saturates
  the exponent AT the rim — the term then sits a few bits below `r`
  instead of `2|e_r|` below, and adding it corrupted the sum
  (`sin(2^(i64::MIN+10))` returned `0.9987·x`, ~2^43 ulps, certified;
  the rim-saturated atan2 quotient path returned `−0` for positive
  arguments through the same arithmetic).
- **(c) the driver's half-width** (`src/math/ziv.rs`): the raw
  `exponent − shift` overflowed for values near the rim — a debug
  panic, wrapped garbage in release.
- **(d)** found while deriving (b)'s fix: once the Taylor floor makes
  deep-tiny evaluations exact-on-grid, directed modes exhaust the
  driver and the ADR-0103 hints — saturated to `u32::MAX` by
  exponent-encoded depths — would request a four-billion-bit
  evaluation: a de facto hang from a 16-byte input.

## Decision

**The lift** (site a). When the base exponent sits within
`wide + 2` of `i64::MIN`, `round_with_infinitesimal` scales the
computation up by an exact power of two (`wide + 4` bits), rounds in
the lifted frame, and scales back — rounding commutes with exact
scaling, and the lifted exponent stays far from the top rim. The
un-lift can saturate in exactly one shape: the base was `±2^(i64::MIN)`
and a magnitude-shrinking round borrowed below it, putting the truth
strictly inside `(0, ±MinPos)` where the no-subnormal grid has
nothing — the inward modes give `±0` with `UNDERFLOW|INEXACT` (the
ADR-0096 sliver shape); the nearest and outward modes never borrow
(an infinitesimal cannot cross the to-nearest midpoint or pull an
outward round below the base), so their arms are defensive. With the
root fixed, the three ADR-0102 interim rim guards are removed; the
dispatches are sound at every Normal exponent.

**The Taylor floor** (site b). Each series loop returns its exact
leading term before forming `r²` when the next term cannot reach the
working window (`2·e_r < −(working + 4)`, saturating): `r` (sin,
atan), `1` (cos) IS the correct w-bit evaluation there, and the
saturated garbage term is never formed.

**The half-width** (site c) uses `saturating_sub`: the clamp at
`i64::MIN` OVERSTATES the half-width, which only refuses
certification — the sound direction.

**The deep-rung ceiling** (site d). `ziv_round_with_depth` clamps
the deep rung's guard to `2^24` extra bits (one clamped rung still
runs — cheap for the tiny-x evals, whose Taylor floors return
immediately). Exponent-encoded depths (an i64 an adversary writes
for free) are bounded by the ceiling; a depth genuinely past it
falls through with the documented 1-ulp INEXACT caveat (the slice
verifier pinned the caveat side: `sin(2^−(2^30))` under TowardZero
returns the argument, not pred, in milliseconds).
Input-precision-proportional hints below the ceiling (a 2 MB operand)
still certify. The deep-tiny evaluations themselves stay cheap at any
working precision (the Taylor floor returns immediately; near-tie
moderate depths converge in a couple of terms).

## Consequences

- The lane pins: atan at `i64::MIN + {1, 10, 54}` (inward = pred,
  nearest = the argument, INEXACT — exact-0-with-OK and
  value-with-OK were the recorded defects), the grow-direction family
  (sinh/atanh) at the rim, hypot's dispatch at the rim (falsely-exact
  was the defect), the bottom-most borrow (`atan(2^MIN)`: inward
  `+0` with `UNDERFLOW|INEXACT`, nearest `MinPos` with no underflow —
  after-rounding tininess, the exp-window convention), the
  `sin(2^(MIN+10))` corruption row, and the rim-saturated atan2
  quotient's sign. All red at their recorded baselines by run
  (debug panics and OK-lies), all green in both profiles now.
- The pre-existing debug panic at `ziv.rs` half-width and the
  `−0`-for-positive atan2 rim case (both recorded on pf-a77o by the
  R2.1/R2.2 verifiers) are gone.
- Failure modes considered (inverted):
  1. **The lift could itself saturate at the top** — it engages only
     for `e ≤ i64::MIN + wide + 2` and adds `wide + 4 ≪ 2^62`.
  2. **The borrow dispatch could mis-flag** — the nearest-mode
     MinPos result carries no UNDERFLOW (after-rounding tininess:
     the unbounded-exponent rounding of `MinPos − ε` is MinPos,
     not tiny), matching exp's window convention; the inward zeros
     flag `UNDERFLOW|INEXACT` (genuinely below every representable).
  3. **The Taylor floor could fire where the term matters** — the
     trigger demands the FIRST correction term below the window's
     resolution plus slack; at the trigger boundary the ordinary
     loop still runs (probed by the slice verifier across the
     boundary).
  4. **The ceiling could refuse a certifiable depth** — only above
     16 Mbit of needed working precision, where the evaluation cost
     is the binding constraint anyway; refusing certification is
     never unsound (the fall-through is the documented caveat).

## Related

- Issues: pf-a77o (closed by this ADR), epic pf-8iji; pf-kh3z (the
  add/sub silent-saturation flags, separate arc); pf-zmk3 (the
  rim-saturated-zero underflow flags, filed by R2.3's verifier).
- Other ADRs: ADR-0102 (the interim rim guards this removes),
  ADR-0103 (the deep rung this bounds), ADR-0096 (the sliver shape
  the borrow dispatch reuses), ADR-0059 (the infinitesimal
  machinery).
