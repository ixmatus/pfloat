# ADR-0096: exp at the exponent rim — certified dispatch and pinned reduction

- **Status**: accepted
- **Date**: 2026-06-11

## Context

The 2026-06-10 workspace review (epic pf-8iji, finding pf-7z66) confirmed
three distinct exp failures at the i64 exponent rim, all reproducer-verified:

1. `exp(-1e300)` returned a garbage Normal near exponent `i64::MIN` with
   INEXACT only: `k = round(x/ln2)` saturated to `i64::MIN`, leaving
   `r = x − k·ln2 ≈ -1e300` — far outside the Taylor series' domain — and
   the compose shift saturated silently.
2. In the wrap window (`x/ln2 ∈ [2^63 − 0.5, 2^63)`), the rounded magnitude
   was exactly `2^63` and `magnitude as i64` wrapped to `i64::MIN`,
   returning a tiny garbage Normal where the truth is `≈ 2^(i64::MAX)`.
3. Just below the window, `exp(6393154322601327829.5)` returned `+∞` where
   the truth is the representable finite `1.34828…·2^(i64::MAX)` (mpmath
   1.4.1 @400 bits).

The module doc promised overflow to `+∞` and underflow to `+0`; none of the
three paths honored it. The Ziv driver cannot catch any of these: the wrapped
or saturated `k` is independent of the working precision, so every iteration
returns the same wrong value and the interval test certifies it.

## Decision

Triage on `e_x = exponent(x)` before the generic reduction, keeping the
established path bit-for-bit untouched for `e_x ≤ 61` (there
`|x|/ln2 < 2^62.6`, clear of both rims).

**Certain regions (`e_x ≥ 63`, or a certified floor past the rims).**
`n = floor(x/ln2)` is the result's binary exponent. `n ≥ 2^63` means the
truth is at least `2^(i64::MAX+1)`, more than an ulp past MaxFinite, so the
IEEE §7.4 overflow shape is forced: `+∞` under the nearest and upward modes,
MaxFinite under TowardZero/TowardNegative, `OVERFLOW|INEXACT` always. This
deliberately diverges from the ops' saturate-to-finite contract (mul/div/fma
keep a meaningful mantissa under a clamped exponent; a deep exp overflow has
no approximation within bounded relative error). `n ≤ i64::MIN − 2` puts the
truth below `MinPos/2`: `+0` everywhere except TowardPositive's MinPos,
`UNDERFLOW|INEXACT`. `n = i64::MIN − 1` is the sliver `[MinPos/2, MinPos)`:
strictly above the to-nearest midpoint (the mantissa `2^{x/ln2 − n}` is
strictly inside `(1, 2)` because `x/ln2` is irrational), so the nearest
modes and TowardPositive give MinPos, the inward modes `+0`,
`UNDERFLOW|INEXACT`.

**Certified floor (`e_x = 62` only).** The band where `n` can fall on either
side of either rim. `ln 2` is bracketed by the neighbours of its q-bit
correct rounding, `x` by directed rounding, and the eight directed quotients
bound `t = x/ln2`; both ends flooring to the same integer certifies `n`,
else `q` doubles up to a cap that scales with `x`'s precision:
`4·(precision(x) + 64) + 1024`. By `μ(ln 2) ≤ 3.57455` (Marcovecchio 2009,
`docs/references/marcovecchio-log2-2009.md`), a `px`-bit dyadic `x` cannot
place `t` closer to an integer than `~2^(64 − μ·(px + 64))`, so the scaled
cap always certifies (`4 > μ` with slack — the same derivation as the
pinned reduction's retry cap). This ADR's first draft capped `q` at a fixed
1024 under a "≤64-bit-span dyadic" premise; the slice's adversarial
verification refuted it with an 1100-bit `x` one part in `2^1037` from
`i64::MIN·ln2`, where the fall-through crossed a DISPATCH rim (sliver
instead of window: wrong TowardZero value, spurious UNDERFLOW) — unlike a
one-off `k` inside the window, which is self-correcting because any `k`
with bounded `r` reproduces `e^x` exactly through the compose. The
fall-through return of the lower floor remains as the defensive total
fallback, now genuinely unreachable by the measure bound.

**Representable window.** `k = clamp(n, i64::MIN+1, i64::MAX−1)` pins the
reduction; the Ziv driver certifies the *unscaled* `s = exp(x − k·ln2)` at
the target precision and mode, and the exact `scale_by_pow2(k)` composes
afterward. Exact power-of-two scaling commutes with rounding, so rounding
`s` rounds the result; the clamp keeps `exponent(s) ∈ {−1, 0, 1}` and the
composed exponent inside i64 for every in-range truth. A target-rounding
carry (certified `s → 4.0` at `k = i64::MAX − 1`) means the
unbounded-exponent rounding of `exp(x)` under the mode lands at
`2^(i64::MAX+1)` — a genuine §7.4 overflow — and dispatches to the
mode-aware overflow result. Letting it reach `scale_by_pow2`'s clamp
instead (this ADR's first draft) returned a non-monotone `1.0·2^(i64::MAX)`
after the carry had replaced the mantissa: below the same input's
TowardZero answer and about half the truth; refuted by the adversarial
verification with `x = 2^63·RD_130(ln2)`. Only upward-rounding modes can
carry (an inward certified rounding never exceeds the true `s < 4`), and
TowardZero correctly stays at MaxFinite WITHOUT the OVERFLOW flag there:
its unbounded-exponent rounding is exactly MaxFinite, which §7.4 does not
class as overflow. The bottom compose cannot saturate at all
(`k ≥ i64::MIN + 1`, `exponent(s) ≥ −1`). The reduction runs at
`wr = w + 256` (covers any `x` up to ~70
bits of precision outright) and grows on realized collapse — `r` must clear
the reduction noise floor `2^(e(k·ln2)+1−wr)` by `w + 8` bits — up to the
measure-derived cap `w + 4·(precision(x) + 64) + 1024` (`4 > μ`). Past the
cap, unreachable by the bound, a collapsed `r` reproduces the Ziv-cap
measure-zero caveat (possible final-ulp directed-mode error), not certified
garbage.

`round_bigfloat_to_i64` also saturates instead of `as i64`-wrapping; the
triage makes the wrap unreachable, but the helper's contract should not
depend on its callers.

## Consequences

- The three reproducers and the full window/sliver/deep classification now
  return mpmath-pinned values with honest flags
  (`tests/regression_review_2026_06_10.rs`, six tests across all five
  modes). Downstream kernels composing through `exp_round` (cosh, sinh,
  tanh, pow chains) inherit the rim behavior.
- The legacy path is untouched for `e_x ≤ 61`: no cost and no behavior
  change for every non-rim input. The Taylor block is extracted to
  `exp_taylor` verbatim for reuse.
- The pf-t6ht probe (flat 24-bit error guard at `e_x ≈ 62` claimed against
  the *legacy* reduction's cancellation) is adjacent but separate: the
  pinned path's deep reduction makes the new band immune, and the legacy
  band below `2^62` stays as probed.
- Failure modes considered (inverted): (1) a mis-certified floor near a rim
  boundary dispatches the wrong region — this ADR's first draft HAD this
  failure (the fixed q cap), found by the adversarial verifier and closed
  by scaling the cap with input precision; the regression lane now pins the
  1100-bit reproducer. (2) The sliver's NE direction depends on the
  mantissa being strictly inside `(1, 2)`; equality would need
  `x/ln2 ∈ ℤ`, impossible for dyadic `x ≠ 0`. (3) An adversarial
  high-precision `x` agreeing with `k·ln2` beyond the retry cap would
  collapse `r`; the cap is derived unreachable from the irrationality
  measure, with factor-4 headroom over the published bound (run-verified
  by the verifier at a 300-bit construction needing `wr ≈ 1240`). (4) A
  compose carry mishandled as saturation — also a first-draft failure,
  closed by the §7.4 overflow dispatch and pinned in the lane with the
  `RD_130(ln2)·2^63` reproducer plus a monotonicity guard.
- Adjacent finding (not this slice): cosh/sinh discard the inner exp
  statuses inside their Ziv closures (`cosh(1e300) → (+∞, Status::OK)`
  observed while the thread-local carries OVERFLOW|UNDERFLOW|INEXACT);
  filed as a discovered-from bead into the flag-fidelity arc.

## Related

- Issues: pf-7z66 (closed by this ADR), epic pf-8iji; pf-t6ht (adjacent
  probe, separate); pf-lkno (the libm-shell saturation analogue, separate
  arc).
- Review: `~/.claude/plans/pfloat-workspace-review-2026-06-10.md` Theme 1
  item 9; reproducer checks E1/E2 in `~/.claude/plans/pfverify-harness/`.
- References: `docs/references/marcovecchio-log2-2009.md`.
- Other ADRs: ADR-0022 (exp under Ziv), ADR-0060 (unconditional INEXACT),
  ADR-0080 (directed-mode saturation posture), ADR-0095 (this arc's agm
  slice).
