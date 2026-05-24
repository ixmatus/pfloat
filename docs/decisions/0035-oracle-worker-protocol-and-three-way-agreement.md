# ADR-0035: Oracle worker reports certified `f32` directly; three-way agreement architecture

- **Status**: proposed
- **Date**: 2026-05-24

## Context

ADR-0034 set up the Oracle layer with an `Enclosure { lo: rug::Float,
hi: rug::Float }` bracket type and an `OracleBackend` trait that returns a
proven enclosure of the true function value. The MPFR backend (in-process
via `rug`) and the Arb backend (out-of-process via `python-flint`) both
implement that contract. The Rust-side verifier then runs the Ziv-at-oracle
loop: at increasing working precision, the oracle returns a bracket, the
verifier rounds both endpoints to `f32` under the requested mode, and
certifies a unique `f32` when both endpoints round to the same value.

Slice p1.7's diagnostic pass on pf-6a4e (BesselI1's 14030 reported
mismatches on f32 subnormal inputs) uncovered three layered defects in
the Arb backend specifically and surfaced an architectural concern about
the bracket protocol generally:

1. **f32 input precision loss in the Arb worker.** The worker decoded
   the f32 bit pattern to a Python `float`, then invoked
   `arb(repr(python_float))`. Python's `repr` is the shortest decimal
   round-tripping through f64, which for f32 subnormals is at most 17
   significant digits even though the exact decimal of a subnormal
   carries ~105 digits. Arb parses the `repr` decimal **literally** at
   `ctx.prec` bits, so the input is perturbed by up to ~7e-62 at the
   f32-subnormal scale. For `BesselI1(x) ≈ x/2` on f32-subnormal
   inputs, the output sits exactly on the f32-subnormal-grid midpoint
   between two neighbors; any sub-LSB perturbation of the input flips
   the NE rounding direction of the output, so the oracle systematically
   certified the WRONG f32 neighbor for 14030 inputs. **pfloat's
   `BesselI1` kernel was correct; the oracle was reporting wrong
   answers.** The probe ran two independent ball-arithmetic computations
   (`arb(repr(v)).bessel_i(arb(1))` vs `arb(exact_bits).bessel_i(arb(1))`)
   and showed the decimals differing at the 17th significant figure,
   exactly where Python `repr`'s precision floor sits.

2. **Decimal bracket-collapse in the verifier's parse step.** Even with
   the input lifted exactly, the Arb worker emitted a decimal bracket
   `[lo_decimal, hi_decimal]` derived from `arb.mid_rad_10exp(n)`. The
   Rust verifier parsed each endpoint via `rug::Float::parse` at the
   Ziv-at-oracle `working_prec` (starting at 64 bits and doubling). At
   low working_prec, the bracket's decimal width was many orders of
   magnitude **below the binary ULP of the parse target precision**.
   The parser's nearest-rounding then collapsed both endpoints to the
   same binary value (the f32-grid midpoint), `to_f32_round(Nearest)`
   tied to even, and `certified_round_f32` reported `Some(wrong_f32)`
   when in fact the bracket spanned a rounding boundary. The Ziv loop,
   seeing a certified answer, exited at low precision without ever
   probing the higher precisions where the same bracket would have
   correctly resolved. The defect manifested as silent false
   certification.

3. **Oracle independence depth.** The Phase 1 plan layered two
   oracles (MPFR and Arb) for redundancy on functions MPFR could cover.
   For the 12 functions that only Arb covers (Si, Ci, li, Bi, Ai_prime,
   Bi_prime, BesselI{0,1,n}, BesselK{0,1,n}), the oracle layer reduces
   to a single backend. The slice p1.5/p1.6 sweeps reported 9 of 10
   non-parametric Arb-primary rows correctly-rounded. Defects 1 and 2
   above prove that the Arb backend alone cannot substantiate the
   correctly-rounded claim with confidence: any silent bug in the
   worker's pipeline (input encoding, ball arithmetic, decimal bracket,
   verifier parse) becomes a silent bug in the v1.0 status table. The
   correctness lift ADR-0034 designed for the layer was real and the
   trait shape was right, but the chain from "Arb's ball arithmetic on
   the true input" to "Rust-side certified f32 bit pattern" passes
   through too many independently-buggable steps to support v1.0's
   correctness claim without independent corroboration.

The probe artifacts (`tests/oracle_i1_probe.rs`,
`tests/oracle_rug_f32_probe.rs`) capture the diagnostic chain that
surfaced both defects. The bucket-by-exponent corpus analysis showed
the 14030 mismatches distributed across every f32 subnormal exponent,
not concentrated at a single magnitude, ruling out a localized kernel
threshold bug; the multi-precision kernel probe showed pfloat's I1
producing the correct decimal at every precision; the multi-precision
Arb probe showed the worker's bracket parse-collapsing at low Ziv
precision. The slice p1.6 worker `rad=0` fix for the inconclusive cases
turned out to be **correct for the cases it covered** (li/si/i1 at +0,
where the rad WAS genuinely zero) but the same code path was hiding the
parse-collapse bug for other inputs whose rad was nonzero-but-small.
Slice p1.6's li closure was incidentally safe (li's output is not at an
f32 midpoint at the tested inputs), but the safety was coincidence
rather than design.

## Decision

This ADR refines ADR-0034's worker protocol. ADR-0034's framing of the
Oracle layer (enclosure posture, two-backend layering, Ziv-at-oracle
correctness argument) stays in force; this ADR replaces the
**bracket-as-decimal-pair** wire format with **certified-f32-directly**
and introduces the architecture that defends the correctness claim end
to end.

### 1. Worker reports the certified f32 directly

The worker's wire-format response changes from a decimal bracket
`[lo_decimal, hi_decimal]` to a `f32` bit pattern (or `INC` if the
worker's ball arithmetic cannot certify a unique `f32` even at its
maximum precision):

```text
request:   <fn_id> <order_or_dash> <input_bits_hex> <mode>
response:  OK <f32_bits_hex>     | INC | ERR <message>
```

The worker's ball arithmetic runs in-process (Arb / mpmath / Maxima all
have native arbitrary-precision interval or ball types), so the worker
can decide certification locally without the decimal bridge. The
worker's loop:

```text
prec = 64
while prec <= 8192:
    ctx.prec = prec
    result_ball = dispatch(fn_id, order, x_arb)
    lower, upper = ball_to_exact_rational_bounds(result_ball)
    cert = certified_round_f32(lower, upper, mode)
    if cert is not None:
        return cert
    prec *= 2
return INC
```

The Ziv-at-oracle loop migrates from the Rust verifier into the worker.
The worker owns it because the worker has access to the ball arithmetic
that tells it whether further precision will tighten the ball; the
verifier lost that information at the decimal step. The
`MAX_PREC = 1024` cap from `tests/oracle/verify.rs` lifts to 8192 inside
the worker (the worker pays the ball-arithmetic cost, not the decimal
bridge cost, so higher precision is much cheaper).

### 2. Shared certified-rounding routine, library-agnostic

The `certified_round_f32(lower, upper, mode)` routine takes exact
rationals `lower <= upper` and returns the f32 bit pattern every value
in `[lower, upper]` rounds to, or `None`. The routine has no dependency
on Arb, mpmath, or Maxima; it operates purely on `fractions.Fraction`
and computes f32 rounding from first principles on the rational
representation. The routine lives at
`scripts/oracle_workers/certified_rounding.py` and is imported by every
worker.

The contract:

- Some(`f`) iff for all `x` in `[lower, upper]`: `round_f32(x, mode) == f`
- None iff there exist `x, y` in `[lower, upper]`: `round_f32(x, mode) != round_f32(y, mode)`

The routine's correctness can be **property-tested exhaustively** in
isolation from any function computation, since it operates purely on
rationals and IEEE 754 grid logic. This separation is the load-bearing
verification handle: a bug in Arb's `bessel_i` or mpmath's `besseli`
will be caught by the three-way cross-check (#3 below); a bug in the
rounding routine itself is caught at the routine's property tests
before any function computation ever runs.

### 3. Three independent oracles with three-way agreement

For every Arb-primary FnId in scope, the verification harness runs
**three independent oracle workers**:

- **Arb** via `python-flint` (already in place; per-call cost lowest)
- **mpmath** via the `mpmath` package (Python; pip-installable; uses
  its own arbitrary-precision arithmetic implemented in pure Python,
  totally independent of FLINT/Arb)
- **Maxima** via `wolframscript`-equivalent CLI (the `maxima --batch`
  command); Maxima's bigfloat arithmetic lineage traces to MIT Macsyma
  in 1968 and shares no code with FLINT/Arb or mpmath; covers every
  Arb-primary FnId at arbitrary precision

The three workers MUST agree on the certified f32 bit pattern for every
input where they are queried. Disagreement halts the slice and triggers
investigation; this is the architectural escalation that catches silent
oracle bugs the slice p1.5/p1.6 episode demonstrated were possible with
a single oracle.

Use-mode by tier:

- **Arb (primary)**: full f32 sweep (65536 inputs per FnId per release).
  Lowest per-call cost; carries the bulk of the verification work.
- **mpmath (cross-check)**: full f32 sweep per release. Higher per-call
  cost than Arb but still manageable. Runs alongside Arb in the
  release-cadence sweep, and its certified f32 must match Arb's for
  every input; any divergence halts.
- **Maxima (sampling)**: hand-derived tricky-case corpus + tie-breaker
  for any Arb/mpmath disagreement + random N-sample per FnId per
  release. Highest per-call cost; sampling layer rather than full
  sweep, because three-way agreement at the high-value cases (corpus +
  disagreements) is the load-bearing signal.

The per-push CI gate stays exactly as before slice p1.7: MPFR-only,
Python-free, fast. The three-way oracle work runs at slice-close
cadence as part of the release sweep, not on every commit.

### 4. Pinned worker-output corpus in tree

The output of the certified Arb worker for a defined input corpus
(every Arb-primary FnId at the L-M-style hard-to-round inputs we have
derived plus the slice p1.7 hand-derived corpus) is **pinned in
`tests/oracle/pinned/`** as a TOML table: input bit pattern, rounding
mode, certified f32 bit pattern, the date the pin was generated, the
oracle (Arb/mpmath/Maxima) that produced the answer, and a short
provenance note for hand-derived cases.

The per-push CI gate diffs the live worker's output against the pinned
corpus and fails on any mismatch. Generated per-slice; updated only by
explicit pin-regeneration commits that explain why the pin moved (a
fix? a discovered defect? a Maxima/mpmath disagreement resolution?).

The pinned corpus is the durable contract. Any future worker change
that intends to alter behavior on the pinned inputs must update the pin
with an explanation, and the explanation goes in the commit message and
the corpus's provenance column. This converts the "I trust the oracle"
question into "look at the diff and the commit log."

### 5. Re-audit slice p1.5 and p1.6 closures under the new architecture

The slice p1.5 and p1.6 sweeps reported 33 + 9 of 10 = 42 Arb-primary +
MPFR-primary rows as correctly-rounded. The pf-6a4e diagnostic showed
that at least one of those rows (BesselI1, now reclassified as oracle
defect not pfloat defect) had a silent disagreement between the oracle
and pfloat. Other rows could have analogous silent issues that the
single-oracle protocol could not detect.

Slice p1.8 onwards re-runs the full Arb-primary sweep under the new
three-way-agreement protocol. Each row's `oracle` field gains a
`oracles_concur` companion field recording which oracles (Arb, mpmath,
Maxima) agreed on the certified f32 for the sampled inputs. Any row
that flips status under the new protocol gets a dedicated investigation
commit; the slice p1.5/p1.6 status tables are not rewritten in place
but supplemented with new rows reflecting the re-audited findings.

## Consequences

### Wins

- **Decimal-bridge bug class eliminated.** The class that includes
  ADR-0034's slice p1.5/p1.6 worker `rad=0` interaction and the slice
  p1.7 bracket-parse-collapse goes away: there is no decimal
  representation crossing the subprocess boundary in either direction
  (only f32 bit patterns and exact-bits input encodings).
- **Verifier surface shrinks.** The Rust verifier's
  `certified_round_f32` and Ziv-at-oracle loop become unused for
  Arb-primary FnIds; only the MPFR-primary path retains them (and MPFR
  is in-process so the decimal bridge never appears there either). The
  Rust code path for the Arb backend reduces to "send hex, parse 4
  hex chars, compare to pfloat output."
- **Three-way agreement defends the correctness claim.** Three
  independent libraries agreeing on every sampled f32 is much stronger
  evidence for v1.0's "correctly-rounded under the status table"
  claim than one library's answer. Silent bugs in any one oracle now
  surface as visible disagreement.
- **Pinned corpus converts trust into auditable diffs.** The pin file
  is the contract; review is reading a diff plus a commit message
  rather than reading an oracle.
- **Property-tested correctness of the rounding routine.** The shared
  routine is small (~150 lines) and operates on exact rationals; it
  can be exhaustively property-tested with `hypothesis` over every f32
  boundary class (normals, subnormals, midpoints, ties, inf/NaN). The
  proof-by-tests is far stronger than what was practical with the
  decimal-bridge protocol.

### Costs

- **Two more dev-only dependencies.** `mpmath` (Python, pip; GPL-free
  BSD) joins `python-flint` in the venv. Maxima joins as a system
  dependency (nixpkgs `maxima`; GPL, subprocess-only so does not
  affect pfloat's MIT/Apache licensing per the "mere aggregation"
  clause and the precedent ADR-0034 set for Arb's LGPL).
  `scripts/setup_oracle_workers.sh` supersedes
  `scripts/setup_arb_oracle.sh` and handles all three.
- **Worker complexity grows.** Each worker now carries the Ziv loop
  and the certified-rounding call. The certified-rounding routine
  itself is non-trivial: rational-based f32 rounding across all five
  IEEE modes, with NE basin logic at midpoint ties.
- **Slice p1.5/p1.6 status table re-audit.** The status rows in
  `tests/oracle/status/` that recorded "correctly-rounded" under the
  old single-oracle protocol must be re-verified under the new
  three-way protocol. Any row whose status changes gets a dedicated
  investigation commit. This adds re-verification cost to slice p1.8
  but does not invalidate the slice p1.5/p1.6 sweeps as a body of
  evidence; it strengthens them.
- **Per-release sweep wall time grows.** mpmath at full sweep is
  slower than Arb at full sweep (rough estimate 3-5× slower per FnId).
  Maxima at sampling is acceptable. Total per-release sweep budget
  doubles relative to slice p1.5/p1.6, but stays bounded (single-digit
  hours, not days). Per-push gate cost unchanged.
- **The certified-rounding routine carries the rigor.** A bug in the
  routine becomes a bug in every oracle's output; the property tests
  and three-way agreement together carry the defense. A bug present
  identically in Arb and mpmath and Maxima would slip through; this is
  the irreducible trust we extend, and it is much smaller than the
  trust extended to a single oracle.

### Out of scope at this ADR

- **Formal proof of the certified-rounding routine.** Kani / Creusot
  do not reach Python; the property tests carry the rigor instead.
  A future v2 might port the routine to Rust and prove it under Kani;
  this ADR does not commit to that.
- **Worker output performance optimization.** The worker's Ziv loop
  may issue many ball-arithmetic calls per f32 input at hard-to-round
  cases. Optimization (early-out heuristics, precision-doubling
  cadence tuning) is deferred until the protocol lands and a measured
  baseline exists.
- **Sharding across cores.** ADR-0034 deferred per-FnId parallelism;
  this ADR carries that deferral forward. Each oracle worker is a
  single Python subprocess at a time.

## Related

- **Refines ADR-0034.** ADR-0034's enclosure posture, two-backend
  layering, and Ziv-at-oracle correctness argument stay in force; this
  ADR replaces the decimal bracket wire format with certified f32
  bits and adds three-way agreement.
- **Reuses ADR-0022.** The Ziv interval test moves from the Rust
  verifier into the Python worker; the underlying argument (precision
  doubling until both endpoints round to the same target) is
  unchanged.
- **Continues ADR-0033.** Phase 1 still runs to completion before v1.0;
  this ADR strengthens the correctness claim by addressing a silent
  failure mode that would have shipped under the slice p1.5/p1.6
  protocol.
- **Plan**: `plans/phase-1-correctness-sweep.md` ("Slice p1.7 +
  follow-on slices" section added at slice p1.7 close).
- **Issues**: pf-tejz (worker exact-bits input, now subsumed by slice
  p1.8 worker rewrite); pf-6a4e (closed at slice p1.7 with the
  reclassification: the 14030 reported mismatches were oracle defects,
  not pfloat kernel defects); pf-1x1b (re-sweep under the new
  protocol); pf-ibu4 (re-verify slice p1.6 li closure under the new
  protocol); pf-izn0 (prose corrections).
