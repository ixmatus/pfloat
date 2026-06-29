# ADR-0106: agm normalizes its operands toward exponent 0

- **Status**: accepted
- **Date**: 2026-06-12

## Context

pf-06lk (filed by the R1 agm slice's verifier, re-verified by run at
the R2.4 baseline): operand exponents within ~2^62 of the i64 rim made
`agm` certify corrupted iterates. The Gauss iteration's `mul`/`sqrt`
saturate their result exponents per the no-emax contract and flag it,
but the eval closure discarded every per-op Status
(`let (prod, _) = a_n.mul(...)`), and the clamp is independent of the
working precision, so the Ziv interval test certified the corrupted
value. Run-verified: `agm(2^(2^62), 3·2^(2^62))` at 53 NE was 47% low
with INEXACT only; `agm(2^(2^62+1000), 3·2^(2^62+1000))` was ~10^31
wrong **with Status OK**; the negative-exponent mirror ~10^301 wrong.
The boundary is sharp at `k = 2^62` (`k = 2^62 − 30` was correct).
ADR-0095 recorded the exclusion when its verifier found the class.

## Decision

Two regimes, split at exponent spread 2^33.

**Huge spreads take a closed form.** The Gauss iterates *converge*
toward the AGM's exponent ≈ `s − log₂(s)` (`s` = half the spread), so
near-convergence products sit near `2s` — no static normalization
keeps a spread ≳ 2^63 loop off the rim (the slice verifier's
refutation 3: the opposite-rim pair certified a ~2^42-wrong value
with Status OK, pre-existing behavior the first draft's ADR claimed
discharged). For `a ≥ b`, `AGM(a, b) = π·a / (2·ln(4a/b))` with
relative error `O((b/a)²)` (the `K(k) → ln(4/k′)` limit; pinned
numerically against mpmath, error shrinking like `(b/a)²`), so once
`2·spread` clears every expressible target-plus-guard the closed form
is the correctly roundable value: `spread ≥ 2^33` gives margin over
`u32::MAX + 1024`. The Ziv driver certifies the composition as usual;
`π/(2L)` multiplies before `·a` so no intermediate can saturate.

**Everything else normalizes toward exponent 0.** Degree-1
homogeneity — `agm(2^m·a', 2^m·b') = 2^m·agm(a', b')` — with
`m = max(⌊(e_a + e_b)/2⌋, i64::MIN + 1)`: the floor midpoint is pure
shift-and-mask arithmetic (`(e_a >> 1) + (e_b >> 1) + (e_a & e_b &
1)`), and the clamp keeps `−m` representable (the verifier's
refutation 1: both operands at the bottom rim made `m = i64::MIN`,
whose wrapped negation scaled both operands onto the same saturated
value and certified an equal-pair Status-OK lie). With the asymptotic
branch owning every spread ≥ 2^33, the normalized exponents stay
within ±2^32 of 0 and nothing in the loop can reach the rim.

**The scale-back surfaces its saturation.** Rounding the normalized
AGM at a target coarser than the operands can carry past `max(a, b)`
into the next binade (the verifier's refutation 2 constructed it
deterministically: operands just below `2^(i64::MAX + 1)` whose AGM
rounds up at 53 bits) — at the top rim that binade is
unrepresentable, `scale_by_pow2` applies the documented saturation
contract, and its flag is merged into the returned Status. The first
draft debug-asserted the scale-back exact: a new panic in debug and a
silently dropped OVERFLOW in release.

## Consequences

- The lane pins, all bit-exact against integer-rounded mpmath
  oracles (structural mantissa/exponent assertions — Display
  saturates at astronomical exponents): the three filed reproducers
  (the 47%-low boundary, the Status-OK worst case, the negative
  mirror) plus the in-range control; the bottom-rim corner
  `agm(2^MIN, 3·2^MIN)`; the top-rim carry (NE takes the saturation
  value WITH the overflow flag, TowardZero stays exact-scaled below
  the binade); the opposite-rim pair and the `s = 2^62 + 100` band
  through the asymptotic branch; and the loop/asymptotic seam at
  spread `2^33 ∓ 2`. ADR-0095's recorded exclusion is discharged —
  by two mechanisms, not one.
- The convergence criterion and iteration budget (ADR-0095) are
  untouched; mid-range operands run bit-identically (the
  normalization scaling is exact and `m = 0` in the common case is
  not special-cased — the scaling is a no-op shift).
- Failure modes considered (inverted):
  1. **All three of this ADR's first-draft claims were refuted by
     the slice's adversarial verification, before commit**: the
     midpoint's negation overflowed at the bottom rim
     (wrong-with-OK, worse than the baseline corruption it
     replaced); the "scale-back cannot saturate" impossibility
     argument missed that the rounding TARGET can be coarser than
     the operands; and "nothing can reach the rim" ignored that the
     iterates converge to the AGM's exponent, not the midpoint.
     The repairs: the `i64::MIN + 1` clamp, the surfaced
     saturation status, and the asymptotic branch.
  2. **The asymptotic threshold could under-cover** — validity
     needs `2·spread > target + guard + margin`; `spread ≥ 2^33`
     dominates `u32::MAX + 1024` for every expressible target, and
     the seam rows pin both sides against the same oracle.
  3. **The closed form could be mis-derived** — its error model is
     pinned numerically (mpmath agreement improving like `(b/a)²`
     at computable spreads).
  4. **The asymptotic branch's first draft was refuted twice by a
     second verification round**: (a) it called the trig-gated
     `pi_at`, breaking the `big,agm` CI combo build — π now comes
     from the agm feature's own Brent–Salamin iteration; (b) it
     computed `ln(big)` and `ln(small)` whole, and for SAME-SIGN
     rim exponents their near-cancellation amplified the logs'
     absolute error by up to `(|e_a| + |e_b|)/spread ≈ 2^31` past
     the charged half-width — 40 constructed near-tie reproducers
     certified the wrong side (1 ulp, INEXACT). The first
     derivation covered only the symmetric and opposite-sign cases:
     incomplete case analysis on the exponent signs. The exact
     exponent split (`L = (e_big − e_small + 2)·ln2 + ln(m_big) −
     ln(m_small)`, integer part exact in two i64 halves, mantissa
     logs O(1)) removes the amplification uniformly; the lane pins
     the verifier's same-sign near-tie (mpmath @2600 bits, NE
     resolving the 0.5000000000009 tie fraction upward).

## Related

- Issues: pf-06lk (closed by this ADR), epic pf-8iji; pf-kh3z (the
  scalar add/sub silent-saturation arc, now the only consumer-side
  reason left to plumb statuses here); pf-a77o (the rim family,
  R2.5).
- Other ADRs: ADR-0095 (the convergence criterion and the recorded
  exclusion this discharges), ADR-0015/0038/0039 (the agm lineage),
  ADR-0099 (the Ball-side saturation soundness analogue).
