# ADR-0092: the pfloat-complex verification posture

- **Status**: accepted
- **Date**: 2026-06-08

## Context

Slice C4 (ADR-0091) landed the complex magnitude, phase, and elementary core
(`csqrt`, `cexp`, `clog`) with C99/C11 Annex G branch cuts, plus the §G.5.1
complex-infinity recovery for `div` and `mul`. Each kernel was adversarially
verified at landing time. C5 is the standing verification pass that pins that
surface for the 1.0 cut and keeps it pinned.

Two hazards shape the design. First, the **branch-cut claim is semantic, not
numeric**: which branch, which signed zero. A numeric oracle cannot see it,
because the rigorous-arithmetic engines (Arb) have no signed zero. Second, the
enumerated tables and the algebraic identities are written by a human reading
the standard and the kernels, so a shared transcription error could let a wrong
test agree with a wrong kernel (the self-oracle circularity). The verification
must therefore combine an INDEPENDENT engine for the numeric claim with
human-readable enumeration for the semantic claim, and the enumeration must be
cross-checked against a primary-source re-derivation that never saw the kernels.

## Decision

The verification is five lanes, each covering what the others structurally
cannot. The split is the point: no single lane is trusted alone.

### 1. Enumerated Annex G special-value tables (the semantic claim)

`tests/annex_g_special_values.rs` enumerates every Annex G special-value row for
`csqrt` (§G.6.4.2), `cexp` (§G.6.3.1), `clog` (§G.6.3.2), and the §G.5.1
`div`/`mul` infinity recovery, through the PUBLIC `Complex` API, across
precisions `{53, 113}`. Every signed-zero row is asserted under all five
rounding modes: that is the regression guard for the C4 `resolve`-mode-sign
defect class (a result zero whose sign is input-fixed was wrong in four of five
modes before the `copysign(0, input)` fix; ADR-0091). Rows the standard leaves a
sign genuinely unspecified are pinned to the crate's chosen representative with
an explicit `[REP]` note. This lane carries the branch-cut and signed-zero
convention that no numeric oracle can.

### 2. Dispatch totality by exhaustive finite-grid enumeration

`tests/dispatch_totality.rs` enumerates the full IEEE class grid (8×8 unary
inputs, 64×64 `mul`/`div` operand pairs). Because the special-value dispatch
branches on a finite grid of component classes, an exhaustive enumeration is a
COMPLETE totality proof: every class combination returns without panicking, in
every mode. This is stronger than a sampled model-checker run would be, and it
runs on every push. It also asserts three invariants that do not re-encode the
row tables (so a shared table error cannot mask them): signaling-NaN ⇒ INVALID,
quiet-NaN-without-infinity ⇒ `(NaN, NaN)`, and the Annex G conjugation symmetry
`f(conj z) = conj(f z)`. The symmetry is asserted bit-exactly only under the
sign-symmetric modes `{NE, NA, TZ}`: under directed rounding the even real part
wants the same mode while the odd imaginary part wants the mirror mode, so the
single-mode identity holds only to ~1 ULP off the symmetric modes (correct
behavior, not a defect).

### 3. Algebraic identities (oracle-free consistency)

`tests/identities.rs` cross-ties the kernels to one another to a few ULPs:
`csqrt(z)² = z`, `cexp(clog z) = z`, `clog(cexp z) = z` on the principal strip,
`cexp(z+w) = cexp(z)cexp(w)`, and `|cexp(x+iy)| = e^x`, over a grid that includes
the named cancellation regimes (`clog` near `|z|=1`, `csqrt` near the cut, `cexp`
near `y = kπ/2`) and Lefèvre-Muller hard-to-round seeds. A defect in any one
kernel breaks an identity the others witness, without an external reference.

### 4. The acb componentwise certified-rounding differential (the numeric claim)

`tests/differential_acb.rs` breaks the self-oracle circularity. It computes each
operation's true value in an INDEPENDENT engine, python-flint's rigorous `acb`
(complex Arb) ball arithmetic, and checks pfloat-complex's output is the correct
rounding bit-for-bit, per component. The check is the certified-rounding test:
Arb encloses the true value `v ∈ [lo, hi]`; when `round(lo, p, mode) ==
round(hi, p, mode)` (value AND sign), that is the unique correct rounding `cr`
and the component must equal it; the oracle precision tightens until it
certifies, and an uncertifiable case is a counted Ziv-cap residual, not a
failure. Finite nonzero components only: Arb has no signed zero, so the
signed-zero rows and the inf/NaN specials are lane 1's job. A negative control
proves teeth (a one-ULP-off value is rejected), and a probe stresses the named
hard regimes. Measured at landing: 2940 componentwise certified-rounding checks
(`csqrt`/`cexp`/`clog` 480 each, `cmul`/`cdiv` 840, Lefèvre-Muller 180, residual
probe 480), zero mismatches, zero Ziv-cap residuals.

The oracle is `acb`, not `rug`/MPC. `rug` is not built with the MPC "complex"
feature in this workspace, and the out-of-process Python subprocess keeps
FLINT/Arb (LGPL) out of the link graph at rest and under test (the ADR-0034 /
ADR-0035 posture, reused). The lane is per-release and venv-gated
(`scripts/setup_arb_oracle.sh`); `PFLOAT_ARB_REQUIRED=1` turns a missing venv
into a hard failure so the backstop cannot silently no-op.

### 5. Kani: advisory, the Status merge only

`src/kani_harness.rs` proves the componentwise `Status` merge contract -- every
operation returns `s_re | s_im`, which must be exactly the union of the two
component statuses -- as a union monoid (exact union, OK identity, commutativity,
associativity, idempotence) for all flag combinations, lifting it from pfloat's
example tests to proof tier. That is the whole Kani surface: every `Complex`
operation runs in `BigFloat`, which is `Vec`-backed and CBMC-hostile (ADR-0062),
so the kernels themselves are unverifiable by model checking and rest on lanes
1-4. The harness discharges with `cargo kani --no-default-features` (5/5). The
`--no-default-features` invocation is load-bearing: the default profile pulls
`pfloat/big`, which activates pfloat's own `cfg(all(kani, feature = "big"))`
verify suite needing math features this crate does not forward; the merge
harness needs only `pfloat::Status`, available in the bare build.

### The independent primary-source re-derivation (verify-the-verdict)

Before the enumerated tables were trusted, an out-of-band multi-agent workflow
re-derived all four tables (`csqrt`, `cexp`, `clog`, §G.5.1) from PRIMARY sources
(N1570 Annex G, cppreference, POSIX, the compiler-rt `__muldc3`/`__divdc3`
reference) WITHOUT reading the pfloat kernels, then adversarially refuted each
row. The verdict was SOUND with zero value or sign discrepancies against
ADR-0091 and the kernels (csqrt 25, cexp 19, clog 24, div/mul 17 rows confirmed).
Two flag-only notes were raised, both standard-permitted and already documented
as best-effort (§G.5.1p5): the `z/0` direction path raises INVALID via the
internal `∞·0` rather than a mandated `DIV_BY_ZERO`, and `cexp`'s mandatory-vs-
optional INVALID split (`x + ∞i` shall, `x + NaN·i` / `NaN + iy` may) is
implemented as the shall/may distinction. Its six suggested coverage ADDs were
each re-derived before adding (the verify-the-verdict discipline,
`feedback_adversarial_review_verify_the_verdict`); five were genuine and added,
and one was itself wrong -- a `(2+3i)·(∞+∞i)` "boxing asymmetry" test that never
reaches the recovery path (naive `im = +∞`, not NaN), so the table pins the
recovery-firing CONDITION and the partial complex-infinity `(NaN, +∞)` instead.

## Consequences

- The numeric correct-rounding claim and the semantic branch-cut claim are
  pinned by DIFFERENT oracles (acb for the value, primary-source enumeration for
  the branch), so neither a kernel bug nor a transcription bug can pass both.
- The dispatch totality is a complete proof over the finite class grid, not a
  sample; the Annex G conjugation symmetry is a standard-mandated invariant
  independent of the row values, so it cross-checks the signed-zero rows.
- The Ziv-cap residual is a measured quantity (zero at landing, including the
  near-tangent and near-unit-circle probe), not an assumed one; a future
  enclosure change that introduced a residual would surface it.
- The always-on lanes (enumerated, totality, identities) ride the per-push
  `--features=fixed,trig` CI job. The acb differential and the Kani harness are
  per-release / manual, dependency-light by construction (no MPFR, no Arb in the
  link graph).
- Kani's reach here is genuinely narrow and the ADR says so rather than
  overclaiming: the soundness of the elementary kernels is an enumeration +
  differential property, not a discharged model-checking proof.

## Related

- Plan: `plans/magical-skipping-lagoon.md` (C5)
- Builds on: ADR-0091 (the C4 surface and its tables), ADR-0090 (mul/div and the
  directed-pair enclosure), ADR-0078 (the ball's independent Arb backstop, whose
  worker and codec this lane mirrors), ADR-0062 (BigFloat is CBMC-hostile),
  ADR-0034 / ADR-0035 (LGPL-out-of-link-graph subprocess oracle posture)
- Tests: `pfloat-complex/tests/{annex_g_special_values,dispatch_totality,
  identities,differential_acb}.rs`, `pfloat-complex/src/kani_harness.rs`,
  `scripts/acb_complex_worker.py`
