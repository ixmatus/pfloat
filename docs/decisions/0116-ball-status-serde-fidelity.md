# ADR-0116: Ball kernel-status fidelity and validating serde deserialize

- **Status**: accepted
- **Date**: 2026-07-03

## Context

The 2026-06-10 workspace review (epic pf-8iji) filed two pfloat-ball
defects that this slice (review remediation R4) closes. Both are
behavior-only fixes under the frozen 1.0 API (ADR-0086), landing as 1.0.1.

1. **pf-t9qq — the elementary functions and `sqrt`/`cbrt` discarded the
   kernel `Status`.** The scalar kernels the ball wrappers call return a
   `(value, Status)` pair, but the wrappers kept only the value and
   hard-coded `Status::OK`. The arithmetic ops `add`/`sub`/`mul`/`div`
   already thread the status (ADR-0099); the elementary surface did not.
   Two observable consequences, both reproduced red before the fix:
   - `Ball::sqrt(point(2))` returned a correct nonzero-radius ball but
     `Status::OK`, when the operation was inexact. Law 5 names `INEXACT`
     the *normal* correct outcome for a ball op, so dropping it left the
     secondary IEEE-flag channel lying.
   - `Ball::exp(point(2^63))` overflowed pfloat's `i64` exponent ceiling
     and returned the *entire* line (a sound but degenerate enclosure)
     with `Status::OK` and no `OVERFLOW`. A caller reading the flags alone
     could not tell an unbounded enclosure from an ordinary one.

   A latent soundness edge sharpens this beyond cosmetics. A directed
   endpoint kernel raises `OVERFLOW` only at the exponent rim (pfloat's
   no-emax contract), where it returns a *finite* saturated value whose
   true image is larger. A monotonic wrapper that builds a bounded ball
   from such an endpoint would produce an interval that EXCLUDES the truth
   — exactly the Law-1 family ADR-0099 closed for `mul`/`div`, now
   possible on the elementary/`sqrt`/`cbrt` surface too.

2. **pf-ufb8 — `Mag` and `Ball` derived `Deserialize`, bypassing their
   construction invariants (code-read, not run).** `Mag::Finite` documents
   a canonical-mantissa precondition (top bit set, `mantissa >= 2^63`) that
   makes the derived `Ord` equal to value order (mag.rs). A hand-crafted or
   tampered wire form such as `Finite{exponent: 5, mantissa: 1}` denotes
   the value `2^-58` yet sorts by `(5, 1)`, which the derived `Ord` ranks
   ABOVE a genuine `2^0` stored as `(0, 2^63)`. That inversion poisons the
   total order the radius pipeline relies on for `Mag::max` and radius
   comparisons; a poisoned `max` can pick the smaller magnitude and
   under-size a radius, a Law-1 break. Likewise a derived `Ball`
   deserialize would admit a NaN or `±∞` midpoint, violating the
   finite-midpoint invariant every constructor enforces.

## Decision

**pf-t9qq — thread the kernel status, mirroring the arithmetic ops, and
flag a degenerate enclosure.** The monotonic enclosure helpers
(`enclose_increasing`/`enclose_decreasing`) now capture both directed
endpoint statuses and route through a shared `finish_enclosure`, which
returns `(Ball, Status)`. The 1-Lipschitz route (`sin`/`cos`) threads the
round-to-nearest midpoint status, exactly as `add`/`sub`/`mul`/`div` thread
theirs. `sqrt`/`cbrt` `OR` their two directed kernel statuses into the
returned flag. The domain functions (`ln`, `log1p`, `cosh`, `acosh`,
`atanh`, `asin`, `acos`) `OR` the kernel status onto their existing
domain-driven `INVALID`/`OK`.

**The degenerate-enclosure flag: raise `OVERFLOW` when a finite input's
enclosure degenerates to unbounded.** `finish_enclosure` raises
`Status::OVERFLOW` in exactly two cases: (a) a directed endpoint kernel
reports `OVERFLOW` (its finite rim saturation would under-cover the true
image, so the sound and only representable response is to widen to the
entire line); or (b) `from_interval` rejects a non-finite (`±∞`) endpoint,
which reaches the same unbounded conclusion. This is the same widening
ADR-0099 chose for `mul`/`div`, applied to the elementary/`sqrt`/`cbrt`
surface. Deliberate scoping:
   - Domain-driven `entire` results (out-of-domain inputs: `ln` of a
     zero-crossing ball, `asin` past `±1`) keep their `INVALID` flag and
     are NOT reclassified as `OVERFLOW`. They return before reaching
     `finish_enclosure`.
   - `INEXACT` and `UNDERFLOW` pass through the OR-monoid as the secondary
     channel; they are never used to widen. `UNDERFLOW` at the rim clamps
     the exponent UP (over-stating magnitude), which over-covers and stays
     sound; only `OVERFLOW`'s clamp-down under-covers, so only it widens.
   - `OVERFLOW` is a signal the scalar kernels raise only at the rim, so
     ordinary results are untouched and tightness is unaffected (the
     truth at the rim is genuinely not representable in a bounded ball,
     so `entire` is the honest answer, not a tightness loss).

**pf-ufb8 — replace the derived `Deserialize` with a validating impl, for
both types, gated behind `serde` exactly as the derive was.** Each impl
deserializes into a raw shadow (a mirror enum for `Mag`, a mirror struct
for `Ball`) whose wire form matches the still-derived `Serialize`
byte-for-byte, then revalidates:
   - `Mag::Finite` with the top bit clear (`mantissa < 2^63`, which also
     catches `mantissa == 0`) is REJECTED with a serde error, not silently
     coerced. A serialized `Mag` is always canonical by construction, so a
     non-canonical wire form is corrupt or adversarial; rejecting surfaces
     it rather than repairing it silently, honoring the honesty-of-the-
     record precedence.
   - `Ball` routes its deserialized parts through `Ball::new`, which
     rejects a non-finite midpoint; its `rad` field is a `Mag` whose own
     validating deserialize guards the radius, so the composition validates
     both fields.

   This mirrors pfloat's `BigFloat`/`FixedFloat` deserialize discipline
   (ADR-0068): deserialize is a trust boundary that revalidates the
   canonical form rather than trusting the bytes. `BallError` keeps its
   derived serde impls (it carries no invariants).

## Consequences

- `pfloat-ball/tests/regression_pf8iji.rs` pins both classes: sqrt/cbrt/
  exp/ln/sin/atan surface `INEXACT` where they were `OK`; `exp`/`cosh` at
  `2^63` surface `OVERFLOW` and stay unbounded (sound); exact results
  (`sqrt(9)`, `cbrt(27)`) stay `OK` and exact (Law 3 preserved); canonical
  `Mag`/`Ball` round-trip unchanged while a non-canonical mantissa, a NaN
  or `±∞` midpoint, and a non-canonical radius are each rejected. Eleven of
  the file's tests were confirmed red against the pre-fix source.
- Two elementary-module unit tests in `src/elem.rs` asserted the buggy
  always-`OK` status (`ln(10)` and `atan(1)`, both irrational) and were
  corrected to assert `INEXACT`; they encoded the defect.
- `serde_json` is added as a dev-dependency (dev-only, not in the shipped
  graph), matching pfloat's serde test setup.
- The soundness self-consistency lane (`property_ftia`) and the ADR-0099
  rim-regression lane both stay green: threading a status does not move any
  interval, and the `differential-arb` Arb containment lane was NOT run in
  this environment (Arb is not installed here); it is a per-release lane
  and this change does not alter any radius, so it is not a gate here.
- The public API is unchanged: every touched op already returned
  `(Ball, Status)`; only the status *value* and the serde *impl kind*
  changed. This is a 1.0.1 behavior fix under the ADR-0086 freeze.

### Inversion (the failure paragraph)

Three ways this fix could be wrong, harmful downstream, or stupid in
hindsight, each checked:

1. **Widening on `OVERFLOW` could destroy tightness if the scalar kernels
   raised it spuriously.** Checked: `OVERFLOW` is raised only by the
   documented rim-saturation paths (ADR-0099 established this for the same
   scalar engine, probed directly there); the regression lane confirms
   ordinary results (`sqrt(2)`, `exp(1)`, `sin(1)`) keep a bounded ball and
   carry only `INEXACT`. The direction analysis is the guard: only the
   clamp-down (`OVERFLOW`) under-covers, so only it widens; `UNDERFLOW`
   clamps up and is left to pass through.
2. **Rejecting a non-canonical `Mag` rather than canonicalizing it could
   break a legitimate round-trip.** Checked: a serialized `Mag` is
   canonical by construction (every constructor produces a top-bit-set
   mantissa), so the only inputs rejection turns away are corrupt or
   adversarial. Silently repairing radius data at a trust boundary would
   mask tampering; rejecting is the honest choice and matches ADR-0068. The
   round-trip test confirms canonical values are untouched.
3. **The shadow enum/struct could drift from the wire form and silently
   fail every deserialize.** Checked: the round-trip tests serialize a real
   value with the still-derived `Serialize` and deserialize it back through
   the new impl; a field or variant-name mismatch would fail them. They
   pass, so the shadow and the derive share one wire form.

## Related

- Issues: pf-t9qq, pf-ufb8 (closed by this ADR), epic pf-8iji.
- Files: `pfloat-ball/src/elem.rs`, `pfloat-ball/src/arith.rs`
  (`sqrt`/`cbrt`), `pfloat-ball/src/mag.rs`, `pfloat-ball/src/ball.rs`,
  `pfloat-ball/Cargo.toml`; lane `pfloat-ball/tests/regression_pf8iji.rs`.
- Other ADRs: ADR-0099 (the `mul`/`div` rim-saturation widening this
  extends to the elementary/`sqrt`/`cbrt` surface), ADR-0068 (pfloat's
  validating serde deserialize, the precedent mirrored here), ADR-0074
  (`Mag` and its canonical-mantissa invariant), ADR-0076/0077 (the ball
  spec and arithmetic soundness laws, including Law 5), ADR-0086/0087 (the
  1.0 freeze and accuracy posture: this is an internal-behavior fix).
```
