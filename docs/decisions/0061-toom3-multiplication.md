# ADR-0061: Toom-Cook 3-way multiplication, the rung above Karatsuba

Status: Accepted (2026-06-04)

## Context

pfloat multiplies limb magnitudes by schoolbook below
`KARATSUBA_THRESHOLD` (48 limbs) and Karatsuba above it (ADR-0027). The
next algorithmic rung is Toom-Cook 3-way, which splits each operand into
three parts and reduces the nine schoolbook part-products to five at the
cost of an interpolation step. Slice 2a.1 (ADR-0040) deferred the
Schönhage-Strassen FFT to v2.0 after measuring that Karatsuba covers the
v1.0 consumer reach; Toom-3 is the cheaper intermediate rung that the same
measurement left open.

The going-in prior was that Toom-3 would **lose** for pfloat. The
in-tree Karatsuba is allocation-bound (each recursion node heap-allocates
several `Vec`s; the 30 to 48 threshold raise in ADR-0027 was a response to
exactly that), and Toom-3 does more work and allocates more temporaries
per node, so the GMP-literature crossover near 20 to 30 limbs (which
assumes an allocation-free Karatsuba) could not be expected to transfer.
The honest expectation was a measure-and-revert, with the negative result
as the deliverable. The slice was run as full-implementation-then-measure
to give Toom-3 its best case before the decision.

## Decision

**Land Toom-3.** The prior was wrong, and instructively so.

### The measurement

A same-session A/B on `benches/mul_thresholds.rs` (aarch64-apple-darwin)
compared `BigFloat::mul` with Toom-3 active against a Karatsuba baseline
(the dispatch threshold raised out of range) at identical sizes. The
small region (≤ 256 limbs) sits inside the run-to-run measurement drift
(~2 to 3%, visible as a uniform shift on sizes that use Karatsuba in both
arms). Above it the win is unambiguous and grows with size:

| limbs | Toom-3 vs Karatsuba |
|------:|--------------------:|
|   384 | ~ −7%               |
|   512 | ~ −5%               |
|   768 | −12%                |
|  1024 | −13%                |
|  1536 | −43%                |
|  2048 | −41%                |
|  3072 | −47%                |
|  4096 | −45%                |

The step change near 1536 is the recursion compounding: once an operand
exceeds three times the threshold its Toom-3 sub-products are themselves
large enough to take the Toom-3 path.

### Why the prior inverted

The allocation-bound nature that was expected to sink Toom-3 is exactly
what floats it. Karatsuba splits by two, so reaching a base case of `b`
limbs from `n` limbs takes `log2(n/b)` levels and a recursion tree of
about `3^{log2(n/b)} = (n/b)^{1.585}` nodes, each allocating. Toom-3
splits by three: `log3(n/b)` levels and about `(n/b)^{1.465}` nodes. For
the same operand size Toom-3 visits **fewer total allocating nodes**, so
the dominant cost here, allocation, falls even though each node does
more arithmetic. In an allocation-free implementation the comparison
would be the textbook one (a marginal asymptotic edge); under pfloat's
allocating multiplier the gap is large and widens with size.

### Threshold

`TOOM3_THRESHOLD = 176`, calibrated to the measured crossover. Below it
the small/noisy band and the small in-tree consumers (a Ziv-lifted Bessel
call lands near 173 limbs) stay on proven Karatsuba; above it the large
operands (a 100,000-bit decimal parse reaches ~1563 limbs, where the win
is ~−43%) take Toom-3. Sub-products recurse through `multiply_limbs`, so
they bottom out in Karatsuba below 176. A lower threshold was measured to
regress: at 48 the sub-products recurse into small-operand Toom-3, which
loses by 12 to 49%.

### Implementation

`multiply_limbs_toom3` in `src/ops/limbs.rs` evaluates the product
polynomial at `{0, 1, −1, 2, ∞}` and interpolates the five coefficients.
The interpolation was derived from the Vandermonde inverse at those
points, cross-checked against Brent & Zimmermann, *Modern Computer
Arithmetic* §1.3.3 and Bodrato & Zanoni (ISSAC 2007); it was not
transcribed from GMP's `mpn_toom33_mul`. Two supporting primitives:

- `divexact_by3`, exact division by three by the modular-inverse method
  (Jebelean, 1993), differentially checked against `divmod_limbs` by `[3]`.
- A `Signed` magnitude type carrying an explicit sign field, because the
  evaluation at `−1` and the interpolation differences go negative even
  though the operands and the final coefficients do not. The sign is a
  represented field, not a two's-complement encoding over a dynamic width.

Correctness rests on a differential test against schoolbook across all
three dispatcher regimes and a dominant-middle-third case that forces the
signed path, plus the existing `differential_mul` MPFR lane (which now
exercises Toom-3 wherever a product exceeds 176 limbs) and the full oracle
suite re-certifying unchanged.

## Consequences

- Large multiplications get materially faster (the consumer tail by ~40%),
  with no value change: the kernels and the oracle status are untouched.
- The win is a property of the allocating multiplier, not of Toom-3 in the
  abstract. The same insight names the real lever for any further gain:
  eliminating Karatsuba's per-node allocation (a scratch-buffer or arena
  multiplier) would shift every rung, and is the prerequisite for an FFT
  rung paying off. That rework is a separate, larger effort.
- ADR-0010's "Toom-Cook 3-way deferred to 1.x" is superseded by this slice.

## Related

- pf-8fda (this work); ADR-0027 (Karatsuba threshold), ADR-0040 (FFT
  deferred, the measurement this continues from), ADR-0052 (the
  Burnikel-Ziegler divider, the other sub-quadratic limb routine).
- Brent & Zimmermann, *Modern Computer Arithmetic*, CUP 2010, §1.3.3.
- Bodrato & Zanoni, "Integer and Polynomial Multiplication: Towards
  Optimal Toom-Cook Matrices", ISSAC 2007, ACM. (Originally cited here
  as WAIFI 2007, conflating this paper with Bodrato's solo WAIFI 2007
  paper; corrected per `docs/references/bodrato-zanoni-2007.md`.)
- Jebelean, "An algorithm for exact division", J. Symbolic Computation 15
  (1993).
