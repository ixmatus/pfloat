# Kickoff prompt — Phase 2 FFT/Schönhage–Strassen multiplication

Paste this into a fresh session to start `/plan`-ing the FFT
multiplication work for pfloat. The bead `pf-rh4c` is the
deliverable; the bead body has the acceptance criteria.

---

I'm starting the Phase 2 perf workstream that lands FFT-based
multiplication ahead of v1.0 ship. Phase 1g (verification
architecture closure) merged into main at `7ed19a4` (2026-05-26);
the four Phase 1g beads (`pf-kk16`, `pf-yupm`, `pf-tqzz`, `pf-hdh8`)
all closed. The user chose to sequence perf work BEFORE the
6–12-hour pf-tqzz full per-release sweep — the sequencing reasoning
lives in MEMORY entry `project_perf_before_full_sweep`.

**This reverses ADR-0010** ("Schönhage–Strassen FFT multiplication
deferred to 1.x", accepted 2026-05-10). The new sequencing folds
FFT into the pre-v1.0 scope. The plan needs to handle the ADR
supersession alongside the implementation work.

Please `/plan` the FFT phase. The bead is **`pf-rh4c`**:

- `bd show pf-rh4c` for the full description, acceptance criteria,
  and design notes.
- `pf-rh4c` blocks `pf-g8h` (v1.0 version bump). The sibling
  workstream `pf-6fvx` (kernel-specific algorithmic improvements
  for Spouge / asymptotic-series cutoffs / Bessel Miller depth)
  also blocks `pf-g8h` and can run in parallel or sequenced after
  FFT — that's a plan-phase decision.

**Durable references** (anchors for the plan agent):

- `docs/decisions/0010-fft-deferred.md` — the prior decision the new
  work reverses. The plan should produce an ADR amendment or a
  successor ADR (ADR-0040 candidate).
- `docs/decisions/0027-karatsuba-threshold-calibration.md` — slice
  7d's empirical-threshold-calibration pattern. The FFT slice should
  follow the same bench-first discipline.
- `benches/mul_thresholds.rs` — the existing Karatsuba threshold
  bench harness. The FFT crossover measurement extends it.
- `docs/decisions/0037-mantissa-payload-inline-storage-rejected.md`
  (pf-cvs) — the precedent for "measure before shipping a perf
  change; reject if the measurement doesn't justify the win".
- `docs/decisions/plans/project_perf_before_full_sweep.md` is NOT a
  doc; the sequencing memory is in the MEMORY system at
  `~/.claude/projects/-Users-parnell-Development-pfloat/memory/project_perf_before_full_sweep.md`.
- ADR-0001 (mantissa storage layout), ADR-0002 (top-bit-set rule),
  ADR-0006 (i64 exponent), ADR-0011 (nightly toolchain pin) are the
  arithmetic-core constraints the FFT implementation must respect.
- The `src/ops/` directory holds the existing arithmetic ops (add,
  sub, mul, div, sqrt, fma); the FFT module lands here.

**Plan-time decisions the plan agent must surface for explicit
user direction:**

1. **NTT vs floating-point FFT.** ADR-0010 §Future-work names the
   negacyclic-transform-in-Z/(2^N + 1) path (Schönhage–Strassen
   original). Modern alternatives exist (Harvey–van der Hoeven
   "fastest known", FP-FFT variants with rigorous round-off
   bounds). Each has trade-offs around correctness story complexity
   vs constant factor. The plan should present the options with
   pros/cons before picking.

2. **Bench-first scope.** Before any FFT implementation lands,
   `benches/mul_thresholds.rs` measures where pfloat's current
   tuned Karatsuba actually loses to a candidate FFT implementation
   on real workloads (not the GMP-literature numbers which may not
   apply). The plan should include the measurement slice as the
   decision gate: if pfloat's Karatsuba covers the precision tail
   end-users hit, the FFT work ships as documentation-tier "we did
   the measurement, FFT is not the win". This is the pf-cvs
   precedent — measure before committing.

3. **Slice cadence.** Long-arc Phase 2 branch (the Phase 1f / Phase
   1g precedent) vs per-slice signed-merge cadence. The work has
   multiple natural sub-slices (bench measurement → variant
   selection → implementation → threshold tuning → ADR landing) and
   each could ship independently. The user's prior preference
   (Phase 1g closure memory) is the long-arc for tightly coupled
   work.

4. **Correctness preservation gate.** Per the sequencing memory
   caveat, perf changes that affect rounding ordering (FMA usage,
   accumulator reordering) can shift "correctly-rounded" results at
   exact-tie inputs. The 47 status TOMLs in `tests/oracle/status/`
   are the cheap-to-rerun sanity gate after each FFT slice; the
   pf-tqzz cross-check smoke is the secondary backstop. The plan
   should make these gates explicit per slice.

5. **Sibling workstream coordination.** `pf-6fvx` (kernel-specific
   algorithmic improvements) is the parallel Phase 2 workstream.
   Does the FFT plan name a hand-off point with `pf-6fvx`, or stay
   purely independent? The kernel-specific work could potentially
   benefit from FFT being landed first (Spouge precision-pegging
   measurements would benefit from FFT-accelerated mul); or it
   could ship first to deliver per-kernel wins faster while FFT
   churns. Plan-time decision.

**Starting moves:**

- `bd show pf-rh4c` to read the bead description in full
- `bd update pf-rh4c --claim` to claim the bead
- `bd show pf-6fvx` for the sibling workstream context
- `cat ~/.claude/projects/-Users-parnell-Development-pfloat/memory/project_perf_before_full_sweep.md`
  for the sequencing context
- Then `/plan` — the plan agent should surface (1)-(5) above as
  AskUserQuestion checkpoints before writing the final plan

The branch should be off `7ed19a4` (current main, post-Phase-1g
merge). The phase name "Phase 2 perf" is reasonable; the FFT slice
within it is "Phase 2a" or analogous (the kernel-specific work
becomes Phase 2b). Plan-phase naming decision.

The session-end signed merge boundary is the YubiKey-gated workflow
per `feedback_phase1f_long_arc_workflow` + `feedback_slice_landing_workflow`;
prompt before the merge command.
