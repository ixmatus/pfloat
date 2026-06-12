# ADR-0101: the exponent-rim dispatch composed through the shell and the exp-family siblings

- **Status**: accepted
- **Date**: 2026-06-11

## Context

The R1 merge (`dfcfd48`, ADRs 0095–0100) turned main red: pfloat-libm's
consistency lane caught the shell's `drive()` returning `+0` where binary32
demands the smallest subnormal under TowardPositive. The root pattern, one
more instance of ADR-0099/0100's "lower layers compose wrongly with the rim":

1. **The shell's drive.** `drive()` converted only `lo` whenever either
   directed bracket end was non-normal, on the assumption "both directed
   roundings agree" — true while pfloat's exp returned garbage Normals at
   the rim, false once ADR-0096 made it honest: the deep-underflow pair is
   now legitimately MIXED (`[+0, MinPos]`), and converting `lo` alone is
   mode-blind (pf-lkno's named defect, observable at last).
2. **The exp-family siblings.** exp2, exp10, expm1, cosh, and sinh compose
   `exp` inside their Ziv closures and DISCARD its Status; exp's mode-aware
   rim dispatch arrived as a bare `+inf`/`+0` that `half_width(non-Normal)
   = 0` certified. exp2/exp10 emitted INEXACT-only; expm1/sinh/cosh emitted
   **Status::OK** on transcendental results (their defensive INEXACT force
   fires only on `Class::Normal` — the pf-9761 posture hole, found again).
   Pre-R1 these paths returned garbage Normals instead; the R1 fix made the
   values honest and the flag losses visible.
3. **The shell's own test oracle.** `exp_matches_single_rounding`'s gold
   standard was `kernel(x, 2048, NearestEven)` then a directed conversion —
   mode-blind at the rim by construction (converting NE's `+inf` under
   TowardZero cannot recover MaxFinite). The test was pinning the defect.

## Decision

**drive(): nudge the special end of a mixed bracket to its adjacent
representable.** NaN pairs keep the early return. Otherwise a zero end of a
mixed bracket nudges via `next_up`/`next_down` (to ±MinPos) and an infinite
end of a mixed bracket nudges to ±MaxFinite, after which the single generic
both-ends-convert-and-compare path decides. Soundness: a directed pair
`[+0, positive]` proves the truth strictly positive (a zero truth collapses
both ends), a mixed `[finite, +inf]` proves it finite (an infinite truth
collapses both ends), and pfloat's MinPos/MaxFinite sit so far outside
every hardware format that the nudged end converts identically to the open
interval end it stands for, under every mode. Agreeing special pairs
(`[+0,+0]`, `[+inf,+inf]`) convert bits-equal and settle unchanged.

**The siblings forward exp's rim instead of discarding it.**
- **exp2** gets the full ADR-0096 triage, simpler than exp's: the result
  exponent IS `floor(x)`, exact integer arithmetic on the input's grid —
  no certified division, no irrationality measure. Certain regions and the
  sliver reuse exp's mode-aware dispatch constructors; the representable
  window Ziv-certifies the unscaled `2^frac` (`frac = x − k` exact) and
  composes by exact scaling with the carry routed to the overflow dispatch.
- **exp10** forwards: the product `x·ln10` is computed once at
  `max(precision(x), target) + 128` bits and handed to `exp_round(…, mode)`
  whole — values, flags, and the certified band classification inherited.
  The forwarding perturbs the argument by `< 2^-(target+126)` relative,
  which the result inherits as a `~2^-(target+65)` relative error (the
  argument's absolute error `|t|·2^-(target+128)` IS the result's relative
  error; the verifier corrected the first draft's constant), so the result
  differs from the true `10^x` only on inputs within that distance of a
  rounding boundary: the documented measure-zero caveat class. The same
  correction applies to the cosh/sinh `+128` rationale below.
- **expm1** forwards verbatim for positive `x` at the rim: `e^x − 1` and
  `e^x` round identically at every expressible target there (the `−1` sits
  `≥ 2^(2^62 − 2^32)` below the ulp) — no caveat at all. The negative side
  already had the ADR-0080 mode-aware `−1 + ε` dispatch.
- **cosh/sinh** forward `exp(|x| − ln2)` (the `e^{−|x|}` term is below
  `2^-(2^62)` relative), with the argument at `+128` bits as for exp10;
  sinh's negative side computes under the mirrored mode and negates
  (`mirror_mode_for_negation`, the pf-l38k class avoided). This closes
  pf-6nn5: the returned Status now matches the forwarded flags instead of
  reporting OK against a raised thread-local.

**The oracle takes the mode.** `kr` in the single-rounding tests is now
`fn(&BigFloat, RoundingMode) -> BigFloat`; exp's closures pass the mode
through (its kernel is rim-mode-aware), the other five single-rounding
tests (ln, sqrt, sin, cbrt, cot) keep an explicit NearestEven — their
original semantics, none of them rim-mode-aware yet.

## Consequences

- The CI lane (`cargo test -p pfloat-libm --release
  --features=differential-mpfr`) is green again, and that lane joins the
  local gate battery permanently — its absence from the R1 batteries is
  how this shipped (the local gate must match CI EXACTLY; the lesson's
  third strike).
- The pfloat lane pins the kernel-level behavior (exp2 rim both signs +
  window + mode-aware values; the four composers' flags incl. sinh's
  mirrored negation), independent of the differential-mpfr feature.
- pow still composes exp internally and was not surfaced by the
  consistency lane (no `sat_pow` table entry); its rim behavior goes with
  pf-vzim/pf-l38k in the Opus arc — recorded here so it is not
  rediscovered.
- Second-round verifier refutations, fixed before commit: (1) the
  cosh/sinh forward landed `|x| ∈ [2^62, 2^62+ln2)` on exp's LEGACY path,
  whose reduction cancels ~`e_x` bits against the flat 24-bit guard — a
  1-ulp NE value regression at percent-level density (`cosh(2^62+0.5)`
  certified the wrong midpoint side). Closed at the root: the legacy
  reduction now carries `e_x + 8` extra bits, repairing the pre-existing
  pf-t6ht band for EVERY exp caller (exp10's forward inherited the same
  hazard). (2) exp2's exact integers past `integer_exponent`'s i64
  magnitude cap: `x = −2^63` (a representable exact power; the window's
  on-grid value would defeat the Ziv interval test into a spurious
  INEXACT) and `x = −2^63 − 1` (the truth EXACTLY at the MinPos/2 tie,
  where the sliver's strict-interior justification fails) get explicit
  rows; the NE tie resolves to +0 (the zero-significand-even convention,
  recorded as a convention since pfloat's no-subnormal grid leaves both
  candidates degenerate).
- Failure modes considered (inverted): (1) the nudge could fire on an
  agreeing-zero pair and manufacture a phantom infinitesimal — the
  mixedness conditions require the other end nonzero/finite-opposed;
  (2) exp10/cosh/sinh forwarding could double-count INEXACT or lose
  INVALID — the forwarded status is returned whole, and the specials were
  dispatched before the forward; (3) expm1's "identical rounding" claim
  fails if a target ≥ 2^62 were expressible — targets are u32, bounded
  three decades below the gap.

## Related

- Issues: pf-lkno, pf-qm0h, pf-6nn5 (closed by this ADR), epic pf-8iji;
  pf-9761 (the same INEXACT-force posture hole in acos, R2 arc); pf-l38k,
  pf-vzim (pow's siblings, Opus arc).
- Review: the R1-merge CI failure (run 27379891399); ADR-0057 (the shell's
  directed-pair architecture), ADR-0080 (the saturation class), ADR-0096
  (the rim dispatch being composed).
- Other ADRs: ADR-0095..0100 (the R1 arc this completes the composition
  story for).
