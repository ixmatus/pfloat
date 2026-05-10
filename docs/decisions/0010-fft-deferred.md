# ADR-0010: Schönhage-Strassen FFT multiplication deferred to 1.x

- **Status**: proposed
- **Date**: 2026-05-10

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

Toom-Cook 3-way and Schönhage-Strassen FFT are tracked as a 1.x
issue. The threshold for Karatsuba is tuned empirically in Phase 7
against the bench harness; the upper bound at which "no faster
algorithm available" kicks in is documented honestly.

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
  tractable.

Tracked under: GitHub issue (to be filed when 1.0 is tagged) and
the `2026-mm-dd-fft.md` plan that will produce the corresponding
ADR.

## Related

- DESIGN.md, "Multiplication" subsection.
- Brent, R. P., and Zimmermann, P. *Modern Computer Arithmetic*,
  Chapter 1.
- Schönhage, A., and Strassen, V. "Schnelle Multiplikation grosser
  Zahlen." *Computing* 7 (1971): 281–292.
