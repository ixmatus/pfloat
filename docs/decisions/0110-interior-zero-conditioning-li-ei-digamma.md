# ADR-0110: interior-zero conditioning for li, Ei, and digamma — geometric cancellation growth, an Ei zero-window, and a Spouge digamma

- **Status**: accepted
- **Date**: 2026-06-29

## Context

The R2.3 (ADR-0105) adversarial verifier surfaced the value-side
siblings of the certified-wrong-answer family: three kernels returned
hundreds of orders of magnitude wrong — two with the **wrong sign** —
at inputs encoding deep proximity to an *interior zero* (pf-0r1l, all
re-verified bit-exact vs mpmath 1.4.1 at the bit-identical p-bit
roundings of each zero). All three certified the corrupted value.
Three distinct mechanisms, one family (the input-encoded proximity to
an irrational zero is rounded away before the kernel resolves it):

1. **Ei at its zero** `x₀ ≈ 0.3725`. `ei_series` rounds `x` to a flat
   `target + 64 + extra` working width and had **no near-zero
   cancellation handling**, so every input deeper than ~245 bits
   SATURATED to the same `~-1.76e-74`, and an input landing the far
   side of the zero came back the **wrong sign**
   (`Ei(RN1500(x₀)) → 53` returned `-1.76e-74` where the truth is
   `+1.745e-452`).
2. **li at the Ramanujan–Soldner constant** `μ ≈ 1.4514`. `li(x) =
   Ei(ln x)` already wraps its composition in `cancellation_boosted`,
   but that probe's **12-iteration linear crawl** (`needed = working +
   cancel + 8`, ~`working` bits per step) saturated at `~1492` bits:
   for an input deeper than that the probe never reached `w ≥
   precision(x)`, where `round(x, w) = x` finally exposes the true
   input and cancellation (`li(RN2000(μ)) → 53` was ~139 orders wrong,
   NE and TowardZero identical — the collapse tell).
3. **digamma at its positive root** `r ≈ 1.4616`. The same crawl
   saturation (wrong sign at `RN1500(r)`), **plus** a cost defect: once
   the Ziv working precision climbed past the 17-pair Stirling table's
   reach, `digamma_positive_at_w` shifted `x` up to `z_min`, capped at
   `2^28` — a ~268M-term recurrence sum costing **~28 minutes** per
   evaluation.

## Decision

**`cancellation_boosted` grows geometrically.** The probe rounds the
input inside the eval closure, so below the true depth `D` it sees a
`w`-bit rounding artifact and reports `cancel ≈ w` — a crawl of
`working` bits per iteration that the legacy 12-iteration cap truncated
at `~1492`. Doubling `w` each iteration (`w = needed.max(2·w)`) reaches
any `D` in `~log₂(D/working)` steps regardless of the Ziv rung, then
converges within one more step once `w ≥ D` fixes `needed`. The
overshoot is at most `2·D` bits — input-proportional (the DoS-budget
posture) — and only ever *adds* accuracy, which the outer Ziv driver
rounds away. The cap is raised to 32 as a pure backstop (real
termination is convergence at `w ≥ precision(x)`). This fixes li and
digamma's conditioning at once; **li needs no other change** (its
existing wrapper inherits the growth).

**Ei wraps its near-zero series in `cancellation_boosted`.** A fixed
dyadic window `[11/32, 13/32]` brackets `x₀` (the lgamma/digamma
root-window precedent, ADR-0097); inside it, `ei_series` runs through
the realised-cancellation boost (returning the O(1) leading magnitude
`γ + ln|x|` as the operand scale), exactly as li already wraps its Ei
composition. Outside the window Ei pays nothing. The boost preserves
the input-encoded proximity (the geometric climb reaches
`precision(x)`), so the sign and magnitude are both right.

**digamma dispatches to a Spouge-derivative past the Stirling reach.**
`spouge_digamma_scaled` is the analytic derivative of Spouge's `ln Γ`
(the [`spouge_lgamma_scaled`] companion):

```text
ψ(z) = ln(z+a) + (z+1/2)/(z+a) − 1 + S'(z,a)/S(z,a) − 1/z
  S(z,a)  = √(2π) + Σ_{k=1}^{a−1} c_k/(z+k)
  S'(z,a) = −Σ_{k=1}^{a−1} c_k/(z+k)²
```

reusing the **same** memoized Spouge coefficients `c_k`. It sums `a − 1`
terms (`a ∝ working_prec`) rather than shifting `z_min = 2^28` times, so
the ~28-minute shift bomb is gone; `digamma_positive_at_w` dispatches to
it past `STIRLING_REACH_THRESHOLD = 600` (lgamma's value; the shared
17-Bernoulli table caps both at ~895 bits), keeping the
`differential_digamma` lane (p ≤ 256) and the shallow Ziv rungs on
Stirling. The formula was **derived from Spouge**, not recalled, and
cross-checked against mpmath. (Cost note: the per-coefficient `exp/ln`
make a single high-precision evaluation roughly `O(w² log w)`, paid
once per working precision via the coefficient memoization and bounded
by the caller-supplied precision; a high-precision digamma is therefore
seconds, not the prior ~28 minutes.)

**The whole Spouge regime re-drives through `cancellation_boosted`.**
`S` and `S'` alternate, so their largest term exceeds the sum by an
internal cancellation hidden behind the `S'/S` division (the pf-wmv7
Spouge-sum lesson, ADR-0097). The first draft tried to self-absorb it
with a fixed `working/8 + 96` margin and a round-back — but this
slice's adversarial verifier refuted that: **the cancellation grows
with the argument**, ~0.1·w at `z≈2.5` (where it passed) but ~0.4·w at
`z≈1e6`, so a fixed margin under-covers large `z` and certified
`digamma(1000000)@1024` ~190 bits short (a silent wrong value — the
"lying inner kernel" shape). The sound recovery is the one
`spouge_lgamma_scaled` already uses: charge the **measured** depth
(`max_term − sum_exponent` for both sums) into the returned scale, and
re-drive through `cancellation_boosted`, which iterates to convergence
— geometric growth reaches `working + depth` for *any* rate. So
`digamma_positive_at_w`'s positive branch routes through
`cancellation_boosted` for the **whole** Spouge regime
(`working > 600`), not only the root window; below the threshold the
shift-Stirling path has no sum cancellation and runs directly.

## Consequences

- Lane rows (mpmath 1.4.1, two-precision cross-checked, bit-identical
  p-bit-rounded-zero inputs): `Ei(RN256(x₀))` resolves the saturation
  (debug); `Ei(RN1500(x₀))` keeps the sign (debug, ~1.5 s);
  `li(RN256(μ))` shallow control stays correct (debug, guards the
  geometric change); `li(RN2000(μ))` certifies (release-gated, ~2.8 s);
  `digamma(2.5)@1024` exercises the Spouge path's correctness and cost
  in the **debug** matrix (a 330-digit reference pins it to 1024 bits,
  catching the sum-cancellation loss; ~7 s — NOT release-gated, so a
  Spouge regression is caught fast); `digamma(RN1500(r))` keeps the
  sign (release-gated, ~2 s — was wrong-sign + ~28 min).
- The cost collapse: `digamma(2.5)@1024` and the deep root went from
  ~28 minutes to sub-second/seconds; the deep li/Ei/digamma roots are
  now release-gated at seconds rather than infeasible.
- `differential_digamma` (p ≤ 256) and the existing `digamma`/`lgamma`
  root-window lane rows stay on their prior paths (Stirling, below the
  threshold) and are unchanged; the geometric growth only converges the
  existing in-window calls faster (more internal precision, same
  correctly-rounded result).
- Failure modes considered (inverted):
  1. **The geometric overshoot could be a DoS** — it is bounded by
     `2·precision(x)`, the caller-supplied operand size; convergence
     (not the cap) terminates it once `w ≥ precision(x)` stabilises the
     probe.
  2. **The Spouge formula could be mis-derived** — cross-checked
     against mpmath at several `z` and the production `spouge_a_for`,
     and pinned by the `digamma(2.5)@1024` bit-exact row and the deep
     root row.
  3. **The Spouge sum cancellation could leave the value short** — the
     first draft did (`digamma(2.5)@1024` ~979/1024 bits, found by
     this slice's own pre-commit verification); fixed by internal
     absorption, mpmath-verified to `≥ working` bits, and a too-short
     50-digit lane reference that masked it was corrected to 330
     digits.
  4. **The Ei window could be too narrow** — any dyadic within `2^-k`
     of `x₀` for large `k` sits well inside `[11/32, 13/32]`; a dyadic
     outside the window has `Ei` bounded away from 0, so no deep
     cancellation regardless of input precision.
  5. **The same Spouge-sum cancellation bug exists in
     `spouge_lgamma_scaled`'s non-window high-precision path**
     (run-verified `lgamma(2.5)@1024` ~73 bits) — pre-existing, untested
     (the differential lane caps at p256, below the Spouge threshold),
     a certified-wrong VALUE. Its blast radius reaches gamma/beta/zeta,
     so it is **filed as pf-rlrb** (the same internal-absorption fix,
     its own verified slice) rather than folded here.

## Related

- Issues: pf-0r1l (closed by this ADR), epic pf-8iji; pf-rlrb (the
  lgamma sibling, filed); recorded as the value-side siblings in
  ADR-0105 §5.
- Other ADRs: ADR-0097 (the `cancellation_boosted` probe and the
  Spouge-sum-cancellation lesson), ADR-0098 (the input-structure
  conditioning family), ADR-0103 (the depth-hint posture), ADR-0105
  (the conditional INEXACT posture — this fix keeps the value right,
  not just the flag), ADR-0026/0030 (the zeta/beta Spouge consumers).
- References: `docs/references/spouge-1994.md` (the Γ/ψ/trigamma
  approximation this derivative is taken from), `docs/references/dlmf.md`
  (6.6.2 Ei series, 6.2.8 li = Ei∘ln).
