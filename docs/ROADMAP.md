# pfloat ecosystem roadmap

This document records the direction pfloat is built toward, both
inside the crate and as part of the wider pure Rust numerics
ecosystem the project plants flags in. State that moves (which
slice is in flight, which beads are claimed) lives in the project's
ADRs, the per-slice plans under `docs/decisions/plans/`, and the
out-of-tree work tracker. This file is direction, with rejected
directions noted so a declined idea is not relitigated from zero.

Written 2026-05-22, from a conversation with the project's lead
consumer about where pfloat goes after the v1.0 tag.

## The umbrella

Pure Rust, `no_std` real, permissively licensed (MIT or
Apache 2.0), correctness first numerics, built from arbitrary
precision float upward. Each crate in the line is justified on the
day it ships; the philosophy is the brand. The horizon is decades.

## Phase 0: finish what is on the bench

ferrodec (decimal, shipped): close the 1.13 medium tier debt.
Optionally add a Python decimal differential lane, metamorphic
transcendental identities, cross precision consistency, and a
quarterly mutation pass. Hardening, not new construction.

This phase belongs to ferrodec, not to pfloat. It is recorded here
only so pfloat's own roadmap stays inside the wider sequence.

## Phase 1: pfloat correctness sweep

Build an exhaustive `f32` verification harness pointed at pfloat
directly. One artifact, three jobs at once: the correctness audit,
the per function correct vs faithful status table, and the
verifier the eventual libm reuses verbatim.

The oracle layer is built behind an abstraction. MPFR (via `rug`)
is the primary oracle for functions MPFR has a primitive for. Arb
(via `python-flint`, run as a long lived subprocess) is the
primary oracle for functions with no MPFR primitive (modified
Bessel I and K, Si, Ci, li, Airy, future Lambert W and incomplete
gamma or beta) and the spot check for the MPFR primary subset.
Both backends return a proven bracket of the true value, not a
rounded scalar; "the bracket determined the rounded `f32`" and
"the bracket did not determine the rounded `f32`" are first class
outcomes of the verifier. Lefèvre and Muller worst case rounding
vectors layer onto the standard elementary functions, so those
functions earn a stronger (not merely statistical) correctly
rounded claim.

Scope is the unary surface only. Multi argument functions (`pow`,
`atan2`, `beta`, `fma`, `agm`) cannot be exhausted over `f32`
(2⁶⁴ plus input space) and stay on differential plus worst case
vectors; their rounding rigor is a later track.

The full work breakdown lives at
`docs/decisions/plans/phase-1-correctness-sweep.md`. The one
discrete decision discharged from Phase 1 so far is ADR-0032: cot,
sec, csc, cbrt, hypot, rootn ship as direct primary kernels in
Phase 2, not as aliases over the existing tan, cos, sin, sqrt,
pow surface; the original Phase 1 plan proposed adding them as
"trivial aliases" pre-freeze, and that step was removed because a
correctly rounded `tan` followed by a correctly rounded reciprocal
is not a correctly rounded `cot`, so the alias would have enabled
a status table overclaim.

## Phase 2: libm spinoff

A thin shell over verified pfloat: widen, compute, round, return,
with the Ziv loop. The Phase 1 harness is re-pointed at the shell
to catch the `BigFloat → f32` double rounding step (ferrodec
flagged the analog as a real bug class).

The claim the spinoff earns: "correctly rounded unary elementary
and special functions, exhaustively verified over `f32`," sharper
than anything the existing pure Rust contenders can say and the
wedge against the mature incumbents.

cot, sec, csc, cbrt, hypot, rootn land in this phase as direct
primary kernels with their own derivations from cited primary
sources (DLMF §4.14 for the trig reciprocals' range reduction,
IEEE 754-2019 §9.2.1 for hypot, §9.2 for rootn; CRlibm and Sun
fdlibm as the state of the art reference implementations). The
libm phase's first commit is a per function kernel list document
that cites ADR-0032 against each of the six as `direct kernel
required, not aliased`.

## Phase 3: pfloat adoption polish

The polish that earns 1.x adoption. Sequenced after the in repo
Phase 8 v1.0 tag and after the Phase 1 status table.

- The per function rounding status table (Phase 1's output)
  published in README and docs.
- `serde` impls behind a feature.
- `num-traits` impls behind a feature.
- Shortest round trip formatter (Dragon4 or a Ryū generalization).
  ADR-0029 deferred this from v1.0; this phase picks it up.
- A constants at precision API (so users can request π or ln 2 at
  any precision without learning the AGM module directly).
- README catches up to the code. An honest "pfloat vs `rug`,
  choose this when" section stating the no C toolchain,
  permissive license, real `no_std` wins, without overclaiming on
  speed or maturity.

## Phase 4 and beyond: the tower

One verified star at a time. Each waits for its dependency to
ship 1.0.

- Complex arithmetic over pfloat. The MPC analog; a genuine gap.
- Ball or interval arithmetic over pfloat. The Arb analog; an
  `interval-1788` flavored crate. Also the rigorous oracle for
  Phase 1, so it closes a loop: pfloat verifies pfloat.
- Verified numerics on the ball layer. Taylor models, verified
  root finding, verified quadrature, rigorous ODE. Smallest
  audience, deepest gap, most distinctively in the project's
  taste.
- Surface pfloat's existing special functions cleanly throughout
  the tower (gamma, zeta, Bessel, Airy, the integrals). Already
  implemented; a head start the dependent crates inherit.
- Bonus star if it calls: exact or computable reals.

Shared traits and types are extracted from concrete crates after
they exist, never pre designed. The `ferrodec-ieee` extraction
pattern, generalized. Rust's abstract algebra trait graveyard is
the warning against the opposite approach.

## Governing rules

- Own the gaps; depend on the incumbents. Build float and up;
  lean on `malachite` or `dashu` for integers and rationals. Do
  not reimplement bignum.
- Stay out, permanently: numerical (inexact) linear algebra
  (`faer`, `nalgebra` own this), random number generation
  (`rand`), any computer algebra system.
- The DAG order is the build order. Strictly. Never start a node
  whose dependency is not at 1.0. Ball arithmetic waits for
  pfloat 1.0; verified ODE waits for ball arithmetic.
- Ship before next. Each node stands alone, worth having if the
  project stops at it. No multi crate vaporware.
- One active build at a time. The scarce resource is attention,
  not code. The queue is not the program; the queue is the
  shaping force that keeps the program one node deep.

## Currently in flight

pfloat in repo Phase 8 (v1.0 tag preparation). Per slice detail
in the project's per slice plans. Phase 1 and everything past it
is queued, not in flight. ADR-0032 is the only discrete decision
that has already been discharged from a future phase, captured
early because the discussion that produced it happened mid Phase
8 and the reminder needed to be durable before the libm phase
opens.

## Related

- `docs/decisions/plans/phase-1-correctness-sweep.md` — the
  detailed work breakdown for Phase 1.
- `docs/decisions/0032-libm-reciprocal-and-root-kernels-direct.md`
  — the one Phase 1 / Phase 2 decision already discharged.
- `DESIGN.md` — current pfloat architecture (the thing Phase 1
  audits and Phase 2 builds on).
