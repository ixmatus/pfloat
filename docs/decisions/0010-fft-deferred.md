# ADR-0010: Schönhage-Strassen FFT multiplication deferred to 1.x

- **Status**: accepted; amended at Phase 2a slice 2a.1 (2026-05-27)
  with a measurement-based reconfirmation per ADR-0040, and again on
  2026-06-04 when Toom-Cook 3-way landed (ADR-0061). The FFT deferral
  stands; the Toom-Cook 3-way deferral does not.
- **Date**: 2026-05-10 (amended 2026-05-27, 2026-06-04)

## Context

Multi-precision multiplication has a hierarchy of algorithms keyed
to operand size (Brent and Zimmermann, *Modern Computer Arithmetic*,
§1.3 and §3.3):

| algorithm | regime (limbs) | complexity |
|---|---|---|
| schoolbook | small (≤ ~30) | O(n²) |
| Karatsuba | medium (~30 to ~1000) | O(n^{1.585}) |
| Toom-Cook 3-way | large (~10² to ~10³) | O(n^{1.465}) |
| Toom-Cook k-way (k ≥ 4) | larger | O(n^{log_k(2k − 1)}) |
| Schönhage-Strassen FFT | very large (≳ 10⁴) | O(n log n log log n) |

MPFR ships all of these. `astro-float` ships schoolbook plus
Karatsuba.

The crossover from Karatsuba to Toom-Cook to FFT happens at
precisions most users do not reach. Cryptographic and
computer-algebra workloads hit FFT regularly; numerical and
financial workloads almost never do.

FFT multiplication is also where a multi-precision library is
hardest to implement correctly. Schönhage-Strassen requires modular
arithmetic in `Z/(2^N + 1)`, careful root-of-unity selection,
recursive applications, and rigorous bounding of round-off through
the negacyclic transform. The literature has worked out every
detail; the engineering effort to land it is still weeks to
months, with a heavy testing burden because the algorithm only
fires above the crossover.

## Decision

1.0 ships schoolbook plus Karatsuba multiplication.

Toom-Cook 3-way landed post-1.0 (ADR-0061, 2026-06-04): the measured
A/B inverted the going-in prior, because pfloat's allocation-bound
Karatsuba makes Toom-3's shallower split-by-three recursion allocate
less, winning ~40% at the consumer tail. Schönhage-Strassen FFT
remains deferred (no v1.0-surface call site reaches the ~10^4-limb
region where it would pay; ADR-0040). The thresholds for Karatsuba
and Toom-3 are tuned empirically against the bench harness; the upper
bound at which "no faster algorithm available" kicks in is documented
honestly.

A user who needs multiplication faster than Karatsuba at 10⁴+ limbs
in 1.0 has two options: use MPFR via `gmp-mpfr-sys` for that
specific path, or wait for the 1.x release that lands FFT.
Documented in the README's "Known limitations" section.

## Consequences

**Wins:**

- Phase 1 ships in weeks, not months. Schoolbook + Karatsuba is a
  well-understood pattern with literature support and reference
  implementations to differential-test against.
- The 1.0 surface is correct, just not asymptotically optimal at
  the precision tail. Honesty about the gap (in the README and in
  this ADR) preserves the conviviality property: a user reading the
  source learns where to pick up if they want to land FFT.
- The 1.x roadmap has a concrete next item beyond bug fixes. The
  cadence stays healthy.

**Costs:**

- Multiplication at 10⁴+ limbs is slower than MPFR. The differential
  lane will surface this not as a correctness gap but as a
  performance gap; CI does not fail, but bench reports show it.
  Documented.
- Users coming from MPFR who happen to use it at very high precision
  (cryptographic key generation, computer-algebra coefficient
  arithmetic) will not switch to pfloat 1.0 for those workloads.
  This is fine; the 1.0 selection criterion is "every user the 1.0
  surface serves is served correctly," not "every user MPFR serves."
- Phase 7's perf calibration has to land empirically-tuned
  thresholds for the schoolbook-Karatsuba crossover. Reasonable
  cost; the bench harness is built up alongside the kernels.

## Future work

The FFT path is large enough to merit its own ADR when it lands. A
full implementation involves:

- A negacyclic transform implementation in `Z/(2^N + 1)` for chosen
  `N`.
- Root-of-unity selection and the recursive structure of the
  transform.
- The modular reduction routine and its differential-test coverage.
- Threshold tuning for the schoolbook → Karatsuba → FFT crossovers
  on the platforms pfloat targets.
- Kani-discharged bounding of round-off through the transform, if
  tractable. Under the Schönhage-Strassen NTT variant the
  round-off bound is moot: integer arithmetic in `Z/(2^N + 1)` is
  exact, so the verification story reduces to modular-arithmetic
  correctness plus convolution-to-product carry propagation.

Tracked under: GitHub issue (to be filed when 1.0 is tagged) and
the `2026-mm-dd-fft.md` plan that will produce the corresponding
ADR.

## Amendment (2026-05-27, Phase 2a slice 2a.1)

The sequencing decision recorded in the project MEMORY entry
`project_perf_before_full_sweep` (2026-05-26) folded Phase 2 perf
work into pre-v1.0 scope, reopening the FFT question with a
measurement-first discipline. Slice 2a.1 of Phase 2a extended
`benches/mul_thresholds.rs` with a `LIMB_SIZES_TAIL` sweep from 768
to 65536 limbs on the calibration host (`aarch64-apple-darwin`) and
landed ADR-0040 with the measurement results.

The measurement reconfirms this ADR's original deferral: pfloat's
tuned Karatsuba covers the precision range its v1.0-surface
consumers reach with several decimal orders of headroom (25.3 µs at
192 limbs equal-size; the FFT-region cost sits in the 30 to 280 ms
range at 16384 to 65536 limbs, two to three decimal orders above
where any in-tree caller operates). The original "crossover happens
at precisions most users do not reach" claim now has a measured
quantification: the Ziv driver's `ZIV_GUARD_CAP = 1024` (16 limbs)
bounds the internal working precision at `caller_target + 16
limbs`, so even an exotic 10000-bit user request stays under 173
limbs internally, 200× short of the literature crossover.

The "Costs" section's "multiplication at 10⁴+ limbs is slower than
MPFR" statement stands and is now quantified: 31 ms at 16384 limbs,
279 ms at 65536 limbs on the reference host.

Phase 2a closes at slice 2a.1; the bead `pf-rh4c` closes as a
documentation-tier deliverable. No code path changes. The
bench-tail extension stays in tree for Phase 2b (`pf-6fvx`,
kernel-specific perf) and for any future 1.x or 2.x revisitation.

## Related

- DESIGN.md, "Multiplication" subsection.
- Brent, R. P., and Zimmermann, P. *Modern Computer Arithmetic*,
  Chapter 1.
- Schönhage, A., and Strassen, V. "Schnelle Multiplikation grosser
  Zahlen." *Computing* 7 (1971): 281–292.
