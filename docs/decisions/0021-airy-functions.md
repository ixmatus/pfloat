# ADR-0021: Airy functions Ai, Bi, Ai′, Bi′

- **Status**: accepted
- **Date**: 2026-05-16

## Context

Roadmap slice 6n adds the Airy functions to the tier-2 specials. It
is the first kernel to carry a Bessel-flavoured asymptotic expansion
without the recurrence machinery the later Bessel slices need, so it
sets the asymptotic-plus-Taylor dispatch shape those slices reuse.
All four functions are entire on the real line, which keeps the
special-value handling simpler than the slice-6m integral cluster
while the regime dispatch is harder than `erf`.

Scope is all four as first-class public methods (`ai`, `bi`,
`ai_prime`, `bi_prime`, each with the `_round` variant and a
`FixedFloat` delegation): DLMF Chapter 9 treats Ai/Ai′/Bi/Bi′ as one
coherent set, the derivative series and asymptotic fall out of the
same machinery, and the Wronskian test then exercises real public
API.

## Decision

**Boundary constants composed at runtime, no new table.** DLMF
9.2.3 to 9.2.6 give `Ai(0)=1/(3^{2/3}Γ(2/3))`,
`Ai′(0)=−1/(3^{1/3}Γ(1/3))`, `Bi(0)=1/(3^{1/6}Γ(2/3))`,
`Bi′(0)=3^{1/6}/Γ(1/3)`. These are built at working precision from
the in-crate gamma kernel and `exp`/`ln`. `3^{k/6}` is formed as
`exp((k/6)·ln 3)` rather than via `pow`: `pow` composes the same
exp/ln and slice 7c (unshipped) is what tightens it to 1 ULP, so the
direct form keeps the error budget explicit. Memoization of the four
boundary constants (the slice-7b1 cache infrastructure) is
deliberately deferred: no perf machinery without a measurement.

**Three-regime sign-aware dispatch.** Small `|x|` uses the
convergent Maclaurin series (DLMF 9.4.1 to 9.4.6) in the two entire
solutions `f`, `g` and their term-by-term derivatives, combined with
the boundary constants and `√3`. Large positive `x` uses the
exponential asymptotic (DLMF 9.7.5 to 9.7.8) in `ζ=(2/3)x^{3/2}`.
Large negative `x` uses the oscillatory asymptotic (DLMF 9.7.9 to
9.7.12) with phase `φ=ζ−π/4`. All asymptotic sums truncate at the
smallest term.

**Provenance correction worth recording.** The asymptotic
coefficient recurrence was first drafted from memory as
`u_k=(6k−5)(6k−3)(6k−1)/(216k)·u_{k−1}`, which is wrong: it drops
the `(2k−1)` divisor. The error was invisible at `k=1` (where
`2k−1=1`) and inflated every later term, costing all precision past
a few digits. Fetching the DLMF primary source pinned the correct
recurrence
`u_k=(6k−5)(6k−3)(6k−1)/((2k−1)·216·k)·u_{k−1}`, with closed form
`u_1=5/72`, `u_2=3465/93312` and `v_k=(6k+1)/(1−6k)·u_k`. This is
exactly the recalled-implementation failure mode the provenance
discipline guards against; the lesson is that a coefficient
recurrence must be derived from the spec, not recalled, and pinned
against an independent reproduction.

**The `e^{−2√ζ}` accuracy law.** The Airy coefficient ratio grows
like `k²`, so the asymptotic series minimises near `k≈√ζ` and its
optimally-truncated error is `≈ e^{−2√ζ}`, not `e^{−2ζ}`. The
regime threshold therefore solves `|x|³ ≥ (9/4)((p+32)/(2 log₂e))⁴`
in integer arithmetic; the fourth-power growth means the asymptotic
only takes over for genuinely large `|x|`.

**Uncapped Maclaurin guard.** The Maclaurin path is the correctness
backstop: it is valid at any precision and any `|x|`, with a working
precision boosted by `≈ (2/3)|x|^{3/2}·log₂e` for the peak term and
the `c1·f − c2·g` cancellation. Unlike `erf`/`si`, this guard is not
capped, matching the library's no-caps ethos; the threshold hands
large `|x|` to the asymptotic, so the series is only ever used where
the guard is bounded.

**Limit conventions.** `Ai(+∞)=+0`, `Ai′(+∞)=−0`,
`Bi(+∞)=Bi′(+∞)=+∞` with `Status::OK` (the exact limits at an
infinite argument, the `exp(+∞)`/`gamma(+∞)` convention). At `−∞`
all four return `+0` with `Status::OK` by the decaying-envelope
convention: the true behaviour is a bounded oscillation with no
limit, so the conservative total result keeps the functions total.

**Tiered oracle.** `rug` 1.27 exposes MPFR's `mpfr_ai` (Airy Ai,
all real arguments) but no Bi/Ai′/Bi′. So Ai gets a true MPFR
differential lane; Bi/Ai′/Bi′ are pinned to a checked-in
authoritative reference table (mpmath, the slice-6m recipe) plus the
Wronskian cross-tie `Ai·Bi′ − Ai′·Bi = 1/π` (DLMF 9.2.7) and a
self-consistency sweep. The differential lanes cap at p=256 because
each evaluation drives the gamma kernel twice and MPFR's own
`mpfr_ai` is costly for large `|x|`; the p=1024 path is pinned
independently by the boundary-constant and Wronskian unit tests.

## Consequences

- All four Airy functions are correct at any precision and any real
  argument (Maclaurin backstop), with the asymptotic as the fast
  path where it is accurate. The Wronskian holds to p−8 over the
  full point table at p up to 256.
- The asymptotic accuracy law is documented, so the regime boundary
  and the test point selection are principled rather than tuned.
- Cost: large `|x|` at high precision routes through the uncapped
  Maclaurin guard, which is slow (a deliberate correctness-over-speed
  trade); the asymptotic optimises the common large-argument case.
  A later slice can add a faster guaranteed-precision large-argument
  method if a workload needs it.

## Related

- Plan: `plans/resume-pfloat-work-slices-sprightly-goose.md`
- Other ADRs: related to ADR-0019 (integral specials, the same
  dispatch and tiered-oracle shape).
- Primary source: DLMF Chapter 9 (§9.2, §9.4, §9.7).
