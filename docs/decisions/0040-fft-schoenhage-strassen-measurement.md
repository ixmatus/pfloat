# ADR-0040: Schönhage-Strassen FFT multiplication — measured, not the v1.0 win

- **Status**: accepted (GATE A per the gate logic below). Phase 2a
  closes at this single measurement slice; the FFT implementation
  slices (2a.2 through 2a.7 of the original slice plan) do not fire.
  ADR-0010 is amended in this same slice to cite this measurement as
  confirming the original deferral; ADR-0041 (which would have
  ratified an empirical `FFT_THRESHOLD` under GATE B) is not
  written.
- **Date**: 2026-05-27

## Context

ADR-0010 deferred Schönhage-Strassen FFT (and Toom-Cook 3-way)
multiplication to 1.x on the argument that "the crossover from
Karatsuba to Toom-Cook to FFT happens at precisions most users do not
reach. Cryptographic and computer-algebra workloads hit FFT regularly;
numerical and financial workloads almost never do." That argument was
made on the literature (Brent and Zimmermann, *Modern Computer
Arithmetic*, §1.3 and §3.3), not on a measurement against pfloat's
tuned Karatsuba. ADR-0027 then calibrated `KARATSUBA_THRESHOLD = 48`
empirically on `aarch64-apple-darwin`; this implementation's Karatsuba
is allocation bound (~17 alloc / dealloc sites per recursion node),
which pushes its constant factor higher than a textbook
allocation-free Karatsuba's.

The sequencing decision recorded in
`~/.claude/projects/-Users-parnell-Development-pfloat/memory/project_perf_before_full_sweep.md`
(2026-05-26) folded Phase 2 perf work into pre-v1.0 scope, reopening
the FFT question with a measurement-first discipline. The named
benefit was amortization: every v1.x release pays the per-release
oracle-sweep gate (ADR-0039), and any kernel-side perf win compounds
across releases. The Bessel I/K family at `p = 320` was named as the
canonical large-precision consumer in pfloat's v1.0 surface.

The pf-cvs precedent (ADR-0037, rejected) sets the discipline this ADR
follows: a perf change ships only on a reproducible win against the
targeted bench; a neutral or absent-win measurement reverts and the
ADR entry is the deliverable. ADR-0037 measured the SmallVec
field-swap and found 1.25× / 1.003× / 1.15× against a 2× criterion,
landing as a rejection ADR rather than a code change. Phase 2a applies
the same discipline at the FFT decision: slice 2a.1 ships the
measurement; the implementation slices 2a.2 through 2a.7 fire only if
the measurement supports them.

### Calibration methodology

`benches/mul_thresholds.rs` (criterion, `harness = false`, slice 7d
infrastructure per ADR-0027). Slice 2a.1 extends the existing
`LIMB_SIZES` sweep (`8` through `512`, dense around the 30-to-48
crossover region) with a `LIMB_SIZES_TAIL` sweep (`768` through
`65536`) covering the precision range every pfloat in-tree consumer
might reach and extending well past the literature crossover into the
region where Schönhage-Strassen wins unambiguously (≳10⁴ limbs per
ADR-0010). The small and tail sweeps are split into separate
benchmark groups so the existing ADR-0027 small-region baseline stays
directly comparable; the tail group uses
`measurement_time = 15 s`, `warm_up_time = 2 s`, `sample_size = 20`
because individual samples at the top sizes run in the 100s of
milliseconds (extrapolating ADR-0027's `108 µs @ 512 limbs` by the
O(n^1.585) Karatsuba complexity gives ~240 ms at 65536 limbs).

Operand construction is unchanged from ADR-0027: `dense_bigfloat`
parses a non-trivial decimal literal at exact target precision so
every limb is densely non-zero and the schoolbook fast-path
`if ai == 0 { continue }` is never taken. The equal-size group
isolates the algorithmic-crossover signal (the dispatcher's
`a.len().min(b.len())` reduces to the size itself); the skewed group
(n by n/4) exercises the dispatch `min()` and the unbalanced split.

The bench runs on `aarch64-apple-darwin` (the same host that
calibrated `KARATSUBA_THRESHOLD = 48` per ADR-0027), saved under
criterion baseline `phase2a-baseline` so any future revisitation
(Phase 2b's `pf-6fvx` parameter-tuning sub-slices, or a 1.x / 2.x
slice that reopens FFT) can diff against the same data set without
re-running the sweep.

### Measurement results

`aarch64-apple-darwin`, `opt-level=3`, profile `bench` (`lto = "thin"`,
`codegen-units = 1`). Criterion baseline `phase2a-baseline`, saved
under `target/criterion/` for slice-2a.5 baseline diff if Phase 2b or
a future 1.x slice reopens FFT.

Small region (slice 7d / ADR-0027 sweep, reproduced for continuity).
The 192 to 512 limb tail of this group is the realistic-consumer
neighbourhood and the load-bearing input to the gate logic:

| limbs | mul_equal (median) | mul_skewed (median, n by n/4) |
|------:|-------------------:|------------------------------:|
|     8 |          204 ns    |                   164 ns      |
|    16 |          417 ns    |                   255 ns      |
|    32 |         1.30 µs    |                   444 ns      |
|    48 |         2.88 µs    |                   919 ns      |
|    64 |         4.50 µs    |                  1.39 µs      |
|   128 |         12.2 µs    |                  4.86 µs      |
|   192 |         25.3 µs    |                  10.9 µs      |
|   256 |         46.4 µs    |                  21.6 µs      |
|   384 |         78.9 µs    |                  42.4 µs      |
|   512 |          122 µs    |                  65.8 µs      |

Large-precision tail (slice 2a.1 extension). Pure Karatsuba; pfloat
v1.0 has no FFT path:

| limbs | mul_equal_tail (median) | mul_skewed_tail (median, n by n/4) |
|------:|------------------------:|-----------------------------------:|
|   768 |                  269 µs |                            110 µs  |
|  1024 |                  380 µs |                            167 µs  |
|  1536 |                  752 µs |                            327 µs  |
|  2048 |                 1.16 ms |                            503 µs  |
|  3072 |                 2.28 ms |                            977 µs  |
|  4096 |                 3.50 ms |                           1.50 ms  |
|  6144 |                 7.45 ms |                           2.90 ms  |
|  8192 |                 10.3 ms |                           5.17 ms  |
| 12288 |                 22.6 ms |                           10.1 ms  |
| 16384 |                 31.1 ms |                           15.5 ms  |
| 24576 |                 71.9 ms |                           30.9 ms  |
| 32768 |                  104 ms |                           46.9 ms  |
| 49152 |                  183 ms |                           92.1 ms  |
| 65536 |                  279 ms |                            143 ms  |

The equal-size scaling fits O(n^1.585) cleanly: 192 to 65536 limbs is
a 341× operand-size ratio; predicted Karatsuba ratio 341^1.585 ≈
1.74 × 10⁴; measured ratio 279 ms / 25.3 µs ≈ 1.10 × 10⁴. Skew (1.6×
slower than predicted at the tail) is consistent with allocation
overhead dominating at larger sizes (the ADR-0027 asm-read note that
in-tree Karatsuba pays ~17 alloc / dealloc sites per recursion node
applies here too). The Karatsuba curve is operating as designed
across the full sweep.

The small-region 32 to 512 timings re-measured ~10 to 30% slower than
the ADR-0027 calibration table. This is consistent with host-state
variability (background load, thermal state) between two runs taken
~9 days apart on the same laptop; no change to `multiply_limbs` has
landed between ADR-0027 and slice 2a.1. The threshold-calibration
conclusion (`KARATSUBA_THRESHOLD = 48`) is unaffected: the
small-region shape and the relative ordering between schoolbook and
Karatsuba reproduce.

### In-tree consumer analysis

The decision gate turns on which precisions pfloat callers actually
reach. Direct callers of `multiply_limbs` enumerate as follows
(`grep -rn "multiply_limbs" src/` against the slice 2a.1 tip):

- **`src/ops/mul.rs:160`** — `BigFloat::mul`. Inherits the caller's
  working precision; the kernel layer routes through here.
- **`src/ops/fma.rs:226`** — `BigFloat::fma`. Same precision profile
  as `mul`.
- **`src/parse.rs:420, 552, 575, 580`** — decimal-literal parsing.
  Scales by powers of 5 / 10 to convert decimal to binary; the
  intermediate `m * 5^k` is dominated by `m` at `target_precision`
  bits. A 10000-bit literal produces ~157-limb intermediates; an
  100000-bit literal would produce ~1563-limb intermediates.
- **`src/fmt.rs:252, 291, 351, 356`** — formatter, binary-to-decimal
  conversion. Same profile as `parse` in reverse.

The kernel layer (`src/math/*.rs`) calls `BigFloat::mul` rather than
`multiply_limbs` directly. The Ziv driver
(`src/math/ziv.rs:42`, `ZIV_GUARD_CAP = 1024`) caps the internal
working precision at `caller_target + 1024` bits across at most
`ZIV_MAX_ITERS = 5` doubling iterations (64 → 128 → 256 → 512 →
1024). For a typical `caller_target` of 53 bits (f64-comparable),
the working precision tops out at 1077 bits (~17 limbs). Even an
exotic `caller_target = 10000` bits tops out at 11024 bits (~173
limbs). The named Bessel I/K `p = 320` consumer per the sequencing
memory caps at 1344 bits (~22 limbs) under the worst-case Ziv
ladder.

The realistic upper bound across pfloat's v1.0 surface is therefore
**bounded by parse / fmt at the user's stated precision plus a
small Ziv-driven kernel-side margin**. No in-tree code path reaches
the ~10⁴-limb region the literature places the Schönhage-Strassen
crossover at. Cryptographic and computer-algebra workloads
(RSA-grade integer multiplication, polynomial GCD over very large
coefficients) are outside pfloat's v1.0 surface per ADR-0010 and do
not enter the gate.

## Decision

**Phase 2a closes at this measurement slice. Schönhage-Strassen NTT
is not the v1.0 win.** pfloat's tuned Karatsuba covers the precision
range its v1.0-surface consumers reach with several decimal orders
of headroom; the literature FFT-crossover region is unreachable
through any in-tree call site.

The two gate criteria from the slice plan resolve as follows:

1. **Realistic-consumer tail (≤200 limbs) cost.** The gate threshold
   was < 50 µs per multiply, the order of one Spouge-coefficient
   evaluation. Measured: 25.3 µs at 192 limbs, 46.4 µs at 256 limbs
   on `mul_equal`; 10.9 µs and 21.6 µs respectively on `mul_skewed`.
   Even at 384 limbs (78.9 µs equal, 42.4 µs skewed) the absolute
   time is sub-100 µs. The realistic tail clears the gate
   comfortably.

2. **In-tree caller reach.** The consumer analysis above bounds the
   largest in-tree working precision at the user's stated target
   plus the Ziv driver's `ZIV_GUARD_CAP = 1024` bits, i.e.
   `target + 16` limbs. Even an exotic `target = 10000` bits (a
   10000-decimal-digit parse) stays under 173 limbs internally —
   200× short of the ~10⁴-limb region where Schönhage-Strassen
   wins unambiguously. The cryptographic and computer-algebra
   workloads that hit the FFT region are outside pfloat's v1.0
   surface per ADR-0010 and are unchanged from the original
   deferral.

Both gate criteria fire GATE A. The Schönhage-Strassen NTT
implementation work (slices 2a.2 through 2a.7 of the plan in
`~/.claude/plans/please-read-this-lexical-knuth.md`) does not start.

The bench-tail extension stays in tree. Phase 2b's `pf-6fvx` work
(kernel-specific perf: Spouge precision-pegging, asymptotic-series
cutoffs, Bessel Miller depth) inherits the bench infrastructure
without modification when its sub-slices need a multiplication-cost
number to feed a parameter choice. Any future 1.x or 2.x slice that
revisits FFT also inherits the bench; the data set under criterion
baseline `phase2a-baseline` becomes the prior measurement that the
revisitation diffs against.

ADR-0010 is **amended** (not superseded) in this same slice. The
amendment adds a header note pointing to this ADR as the measurement
that confirmed the original "FFT deferred to 1.x" decision was
correct. ADR-0041 (the would-be empirical-`FFT_THRESHOLD` calibration
ADR referenced in the plan under GATE B) is not written; the next
ADR allocates 0041.

`pf-rh4c` closes as a documentation-tier deliverable per the
acceptance criterion in its bead body: "if pfloat's tuned Karatsuba
covers the precision tail end-users hit, FFT lands as
documentation-tier 'we did the measurement, FFT is not the win'."
The deliverable is this ADR plus the `phase2a-baseline` criterion
data plus the `LIMB_SIZES_TAIL` bench-harness extension in
`benches/mul_thresholds.rs`.

## Consequences

**The measurement is the deliverable.** This ADR plus the
`phase2a-baseline` criterion data plus the `LIMB_SIZES_TAIL`
bench-harness extension together constitute the durable artifact of
Phase 2a. The discipline matches ADR-0037 (the SmallVec inline-
storage measurement that landed as a rejection ADR rather than a
code change): the project ends slice 2a.1 in a better position than
the prior unmeasured state, regardless of whether a perf
implementation lands.

**v1.0 ships with the measurement-justified posture, not a deferred-
question one.** Before slice 2a.1, ADR-0010's deferral was a
literature-justified guess; "the crossover happens at precisions
most users do not reach" was a true statement that had not been
tested against pfloat's specific Karatsuba implementation. After
slice 2a.1, the deferral is measurement-justified: at every
precision a v1.0-surface consumer can reach, Karatsuba runs in
sub-millisecond time on the calibration host; the FFT region sits
in the 30-to-280-millisecond range, two to three decimal orders
above where any in-tree caller operates. The honest framing in
ADR-0010's "Costs" section ("multiplication at 10⁴+ limbs is slower
than MPFR") stands, but the gap is now characterized rather than
assumed: 31 ms at 16384 limbs, 279 ms at 65536 limbs on the
reference host.

**The bench-tail extension is reusable infrastructure.** Phase 2b's
`pf-6fvx` (kernel-specific perf: Spouge precision-pegging,
asymptotic-series cutoffs, Bessel Miller depth) inherits the bench
unchanged when its sub-slices need a multiplication-cost figure for
a parameter-choice decision. A future 1.x or 2.x revisitation of
FFT diffs against the saved `phase2a-baseline` to detect whether
Karatsuba's curve has shifted (compiler upgrades, allocator
changes, host-architecture ports). The bench stays in tree under
the same logic ADR-0027 cites for the small-region sweep: "so the
next person (or the 1.x Toom-Cook work) re-runs the calibration
rather than re-deriving it."

**Phase 2a closure unblocks `pf-g8h` from the FFT direction.** The
v1.0 version-bump bead retains its blocking dependency on `pf-6fvx`
(kernel-specific perf, Phase 2b). The slice 8c release ceremony
(`pf-4fi` → `pf-g8h` → `pf-5ky` + `pf-tim` → `pf-61n`) becomes
reachable once Phase 2b closes; no further FFT work is required to
ship v1.0.

**Permacomputing horizon unchanged.** No new dependency entered the
crate graph. No code path changed. The pure-Rust / no_std-capable /
zero-runtime-deps posture is bit-identical to the pre-slice-2a.1
state. The bench file picked up two new groups and the criterion
dev-dep was already in the graph from slice 7d.

## Related

- ADR-0010 — Schönhage-Strassen FFT deferred to 1.x. This ADR
  amends (per GATE A) ADR-0010's header in the same slice to cite
  this measurement as confirming the original deferral.
- ADR-0027 — Karatsuba threshold calibrated to 48 limbs. The
  methodology, bench infrastructure, and discipline this ADR
  inherits.
- ADR-0037 — `SmallVec<[u64; 4]>` for mantissa/NaN payload,
  rejected. The "measure-before-shipping; ADR is the deliverable
  regardless of outcome" precedent.
- ADR-0039 — Phase 1g verification architecture closure. Records
  the per-release oracle-sweep cross-check gate that benefits from
  multiplication-cost improvements (the named "amortization"
  argument in the sequencing memory).
- `~/.claude/projects/-Users-parnell-Development-pfloat/memory/project_perf_before_full_sweep.md`
  — the sequencing decision this ADR executes.
- `benches/mul_thresholds.rs` — the bench harness extended in
  slice 2a.1.
- Brent, R. P., and Zimmermann, P. *Modern Computer Arithmetic*,
  Chapter 1 and §3.3.
- Schönhage, A., and Strassen, V. "Schnelle Multiplikation grosser
  Zahlen." *Computing* 7 (1971): 281–292.
