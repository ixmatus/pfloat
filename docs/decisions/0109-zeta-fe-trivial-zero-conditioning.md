# ADR-0109: zeta_fe pre-boosts its working precision by the input's proximity to the negative integers

- **Status**: accepted
- **Date**: 2026-06-29

## Context

The 2026-06 review-remediation arc R2 (ADR-0103) recorded, as the
bound on what its driver-side depth hint could ever fix, a value-side
defect its adversarial verifier surfaced in the zeta functional
equation: `ζ(s)` for `s` just past a *trivial zero* certified a value
hundreds of orders of magnitude too small (pf-hkoj, run-verified
bit-exact vs mpmath 1.4.1, byte-identical at the R2 baseline e8b1284).

`ζ(s) = 2·(2π)^{s−1}·sin(πs/2)·Γ(1−s)·ζ(1−s)` (DLMF 25.4.2). The
trivial zeros `ζ(−2n) = +0` are dispatched exactly upstream, but the
*neighbourhood* `s = −2n − ε` feeds `zeta_fe`, and there the entire
result rides on `sin(πs/2)`: near an even integer `sin(∓nπ) = 0`, so
`ζ(−2n − ε) ≈ −ε·ζ'(−2n)` has magnitude `|ε|`. `zeta_fe` computed at a
flat `working = target + 96` and rounded `s` to that width *before*
forming `sin`, so for `ε` deeper than the working precision `s`
collapsed to exactly `−2n`, `sin(πs/2)` collapsed to (working-precision
noise around) zero, and the product came out `|e(ε)|` orders too small.
The Ziv interval test then certified it at the **first rung**: the
half-width model `|y|·2^-(working − guard)` claims the eval is accurate
to `working − guard` bits relative to `|y|`, but the flat boost only
delivers `(working + 96) − depth` bits, which is negative once the
proximity depth exceeds `96 + guard` (ADR-0103 §3b — the kernel's
obligation the driver hint cannot repair).

Reproducers, re-derived from mpmath at the bit-identical exact inputs:

- `ζ(−2 − 2^-1200) → 53` returned `9.944e-67` where the truth is
  `+1.76836e-363` — ~296 orders wrong, NE and TowardZero identical.
- the shallow `ζ(−2 − 2^-200) → 53` returned a value ~5.5e-5 relative
  (~39 wrong bits), all modes, certified.

This is the pf-gg96 mechanism (ADR-0098) one branch over: there the
proximity was to the *pole* at 1 on the positive axis; here it is to
the *zeros* on the negative axis, and the cancelling factor is `sin`
rather than `1 − 2^{1−s}`.

## Decision

`zeta_fe` pre-boosts its working precision by the exact proximity of
`s` to the nearest integer, the **same `pole_proximity_depth` pattern
the lgamma and digamma reflections already use** (ADR-0098): the gap to
the nearest integer is Sterbenz-exact on `s`'s own grid, so its negated
exponent is an exact, input-encoded depth. With
`working = target + 96 + pole_proximity_depth(s)`, the round of `s`
preserves the `ε` offset, `sin(πs/2)` carries it, and the eval is
accurate to `target + 96` bits relative to `|y|` at every Ziv working
precision — restoring the half-width premise.

The boost is bounded by the input precision (the gap exponent cannot
exceed `precision(s) − e_s`), so the cost is proportional to what the
caller supplied (the bignum DoS-budget posture); the working-precision
cap gains the matching `+ precision(s)` term, exactly as `zeta_borwein`
already carries for the pole side. There is **no new kernel and no new
cancellation probe**: the factors of the functional equation still
multiply without cancellation (`Γ(1−s)` grows, `(2π)^{s−1}` decays,
`sin` carries the small factor), so a flat pre-boost — not
`cancellation_boosted` — is the right shape.

The trivial-zero dispatch (`ζ(−2n) = +0` exact) is unchanged; this
fix serves the near-but-not-at neighbourhood the dispatch feeds.

## Consequences

- The lane pins three rows (mpmath 1.4.1, two-precision cross-checked,
  exact inputs): the shallow `ζ(−2 − 2^-200)` positive and bit-exact,
  NE and TowardZero (debug, ~4 s); the deep `ζ(−2 − 2^-1200)` positive
  and bit-exact, NE and TowardZero (release-gated, ~9 s release / ~120 s
  debug, the pf-jl35 precedent); and an **odd-integer control**
  `ζ(−3 − 2^-200)` where the boost fires (depth 200) but `sin` *peaks*
  rather than vanishing, so the value must stay the smooth `RN53(1/120)`
  — guarding against the boost perturbing the non-cancelling case.
- Generic `s` (not within `2^-1` of an integer) gets
  `pole_proximity_depth = 0`, so `working` is the unchanged `target + 96`
  and there is no cost regression off the integer neighbourhoods.
- Near the *odd* integers the boost over-provisions (the result is a
  smooth nonzero rational, no cancellation), but only by the
  input-bounded depth, and only inside `2^-1` of an integer — a bounded,
  rare cost, not a correctness issue (the control row proves the value
  is right there).
- Failure modes considered (inverted):
  1. **The boost could mask a genuine collapse** — it cannot: the
     trivial zeros are dispatched exactly upstream (`is_negative_even_integer`),
     so a zero `sin` factor is only reachable through unresolved
     proximity, which the boost resolves rather than hides.
  2. **The depth could be computed on the wrong grid** —
     `pole_proximity_depth` takes the gap on `s`'s own precision
     (Sterbenz-exact), the same helper the lgamma/digamma reflections
     are already verified against; the lane's deep row exercises the
     1200-bit depth end-to-end.
  3. **Incomplete case analysis on parity** — both even (cancelling,
     boost load-bearing) and odd (peaking, boost over-provisioning)
     integer neighbourhoods are covered by lane rows; the pole at `s = 1`
     remains on the positive-axis Borwein path with its own
     conditioning (ADR-0098) and depth hint (ADR-0103), untouched here.
  4. **The `s → 0⁻` neighbourhood is a residual, out of scope, improved
     not regressed** (found by this slice's adversarial verification).
     `ζ(0) = −1/2` is dispatched exactly for `±0`, but `s = −2^-k`
     (tiny negative) reaches `zeta_fe`, and there the functional
     equation is a `0 × ∞` form: `sin(πs/2) → 0` *and* `ζ(1 − s)`
     approaches the *pole at 1*. The `sin`-side boost this ADR adds
     and `zeta_borwein`'s pole-side conditioning (ADR-0098) then
     compound, and for `k` past the working-precision cap `1 − s`
     rounds onto the pole, where `zeta_borwein`'s defensive belt
     returns an **honest NaN** the driver cannot certify — never a
     certified-wrong finite (the class this arc targets). Verified:
     `ζ(−2^-500 @p564) = −0.5` correct (baseline `35c4415` NaN'd past
     depth ~1172, so the fix strictly *improves* this region);
     `ζ(−2^-8000 @p264)` is honest NaN; deep cases under directed modes
     are slow-but-correct (the depth hint targets the pole at 1, not
     the output boundary). Filed as pf-qt7v (the `0 × ∞` near-zero
     conditioning) rather than expanded here — it needs the
     pole-and-zero double conditioning, a distinct mechanism from the
     trivial-zero proximity pf-hkoj names.

## Related

- Issues: pf-hkoj (closed by this ADR), epic pf-8iji; pf-qt7v (the
  `s → 0⁻` `0 × ∞` residual this slice's verifier surfaced, filed not
  fixed); recorded as the bound on ADR-0103's mechanism (§3b) and the
  value-side sibling in ADR-0105 §5.
- Other ADRs: ADR-0098 (the `pole_proximity_depth` pattern and the
  pole-at-1 conditioning), ADR-0103 (the driver depth hint and the
  half-width-is-the-kernel's-obligation framing), ADR-0026 (the zeta
  algorithm and functional-equation branch), ADR-0105 (the conditional
  INEXACT posture; this fix keeps the value right, not just the flag).
- References: `docs/references/dlmf.md` (25.4.2, 25.6.4).
