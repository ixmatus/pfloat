# ADR-0105: acos's gap boost, and the INEXACT force covers the zero-result gap

- **Status**: accepted
- **Date**: 2026-06-12

## Context

pf-9761 (filed by the R1 cancellation slice's verifier, run-verified
pre-existing): `acos(1 − RN150(π)·2^-202 @p400) → 53` returned **+0
with Status OK** where the truth is ≈ 9.887e-31. acos does not compose
asin — it evaluates `2·atan(sqrt((1∓x)/(1±x)))` with `1 − x_w` formed
at the Ziv working precision and no input-derived boost — so
input-encoded proximity to 1 collapses `x_w` to exactly 1, the
evaluation returns exact 0, `half_width(0) = 0` certifies on the first
rung, and the defensive INEXACT force fired only on `Class::Normal`
results: the wrong zero also claimed exactness. Two defects in one
reproducer: the collapse (acos-local) and the posture hole
(class-wide — ADR-0101 had already hit it again in expm1/sinh/cosh).

## Decision

**acos: the ADR-0097 gap boost, both branches.** `1 − |x|` is
Sterbenz-exact at the input's precision; when its exponent is ≤ −2,
every Ziv evaluation runs boosted past the proximity depth
(`−exponent(gap) − 1 + 8`, the asin constant). Both branches share the
gap: the negative branch's `1 + x_w` cancels at the same depth for
`x` near −1. acos near ±1 then evaluates `≈ sqrt(2δ)` — an ordinary
normal value with generic boundary distance — so the driver certifies
normally and no depth hint is needed (unlike the on-grid-adjacent
ADR-0102/0104 family).

**The posture: the INEXACT force covers Zero results, class-wide.**
The 33 per-kernel `matches!(result.class, Class::Normal { .. })`
INEXACT-force sites (ADR-0063's family) now route through one shared
helper, `force_transcendental_inexact`, which forces INEXACT on
`Normal` **and `Zero`** results. Soundness rests on a per-kernel
audit: a Zero reaching the post-driver force can only be collapse,
because every kernel in the family either dispatches its exact-zero
inputs before the driver or has no exact zeros at dyadic inputs:

- Dispatched exact zeros: sin/asin/atan/atanh/asinh/sinh/tanh/expm1/
  log1p/erf/Si at x = 0; acos at x = 1; acosh at x = 1; lgamma at
  x ∈ {1, 2}; zeta's trivial zeros ζ(−2n); beta's Γ-pole zeros via
  the ADR-0098 exact-sum classification; atan2's axis zeros.
- No dyadic zeros (Lindemann–Weierstrass or kernel range): cos, tan,
  cot at irrational multiples of π/2; sec/csc/cosh/gamma/erfc never
  zero; digamma's root irrational; real ζ has only the trivial
  zeros; Bessel J/Y/I/K and Airy zeros transcendental (Siegel).
- Conditional (open problems, the ADR-0066 posture): li's zero at
  the Ramanujan–Soldner constant μ, and Ei's and Ci's real zeros —
  not proven irrational; the force is honest conditional on those
  constants not being dyadic. li's in-code comment claimed μ
  transcendental, overstating the literature; corrected to the
  conditional phrasing in this slice (the map-versus-territory
  rule).

Infinity results are NOT folded in: a post-driver Infinity either
came from an honest rim dispatch carrying its own flags (ADR-0096/
0101) or from discarded internal statuses — the latter is pf-kh3z's
flag-fidelity arc, not a posture this helper can repair blindly
(forcing INEXACT on a legitimately-flagged overflow would be noise,
and the zeta certified-Inf case was fixed structurally in ADR-0098).

## Consequences

- The lane pins the named reproducer bit-exactly (NE and TowardZero,
  plus the negative branch through π − acos and a shallow control).
  The R1 asin dense-delta construction is reused, so the two kernels'
  guards stay aligned.
- The 33-site rewrite is behavior-preserving everywhere a result is
  Normal (the overwhelming case) and changes only the collapsed-zero
  statuses; the full lib suite (872) and lane (49) pass unchanged,
  which is itself the audit's run-verification: no kernel's
  legitimate exact zero reaches the force.
- Failure modes considered (inverted):
  1. **A kernel with a genuine post-driver exact zero would now lie
     INEXACT** — the audit above is the guard, and it is recorded
     here precisely so the next kernel added to the family checks
     itself against it (the helper's doc points back).
  2. **The boost could under-cover acos's composition** — the asin
     derivation transfers: the gap-relative error of `1 ∓ x_w` is
     `2^-(w+boost−d)`, halves through sqrt, passes atan with
     sensitivity ≤ 1, and the `π −` shift on the negative branch is
     cancellation-free (result ≈ π there); the verifier probes both
     branches deeper than the lane.
  3. **The conditional zeros** (li/Ei/Ci) could in principle make
     the force a lie if μ etc. turned out dyadic — the same
     conditional soundness ADR-0066 already ships for γ and ζ(5),
     now stated rather than implied.
  4. **The slice's adversarial verification refuted the first draft
     three ways**, all fixed pre-commit: (a) the new helper used a
     `Status` import cfg-gated behind `trig`, breaking the
     `--no-default-features --features=big,agm` CI combo outright
     (now fully qualified, no gated import); (b) the helper had been
     spliced into the middle of `ln_2_at`'s rustdoc, leaving it a
     false opening sentence and `ln_2_at` an orphaned doc fragment;
     (c) the comment correction stopped at li.rs while ei.rs and
     ci.rs carried the identical transcendence overclaim — and the
     verifier further showed the three kernels' BLANKET
     transcendence sentences ("takes transcendental values at
     algebraic arguments") are themselves not theorems: li/Ei/Ci
     values at algebraic points are γ-entangled, and only the
     E-function component carries a Shidlovskii-class result. All
     three comments now state the ADR-0066 conditional posture,
     which ADR-0064's scope split had already recorded as the
     honest basis.
  5. **The verifier's deep-zero probes surfaced the value-side
     siblings** (pre-existing, bit-identical at baseline, all
     INEXACT-flagged — the posture fix removes only the exactness
     lie): li at the p2000-parsed Soldner constant is ~139 orders
     wrong; Ei at its p1500-parsed zero has the WRONG SIGN (378
     orders); digamma at its p1500-parsed root has the wrong sign
     (and costs ~28 minutes for the single evaluation); ζ near its
     trivial zeros returns sign-independent wrong values (pf-hkoj's
     scope, confirmed). Ci at p1500 and lgamma at depth 800 are
     bit-exact, so the cancellation machinery exists — li, Ei, and
     digamma's deep-root paths lack it. Filed as a discovered-from
     bead rather than silently absorbed here.

## Related

- Issues: pf-9761 (closed by this ADR), epic pf-8iji; pf-kh3z (the
  Infinity-side flag fidelity, separate arc); pf-hkoj (zeta_fe near
  trivial zeros — the value-side sibling, this arc).
- Other ADRs: ADR-0097 (the gap-boost pattern), ADR-0063/0064/0066
  (the INEXACT posture lineage), ADR-0101 (the posture hole's
  previous sighting), ADR-0098 (beta's exact-sum zero dispatch).
