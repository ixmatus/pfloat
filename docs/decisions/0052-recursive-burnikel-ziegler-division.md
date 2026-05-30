# ADR-0052: Recursive (Burnikel-Ziegler) integer division

- **Status**: accepted
- **Date**: 2026-05-30

## Context

`divmod_limbs` (`src/ops/limbs.rs`) was Knuth Algorithm D throughout,
`O(qlen · dlen)` in the quotient and divisor limb counts. That is the
right constant for small divisors but quadratic for large ones. Two
callers reach large divisors: the decimal base conversion in `fmt`
(dividing by powers of ten near the operand's magnitude) and `div` /
`parse` at high precision. The quadratic division was the asymptotic
behind the formatter's `O(digits^2)` `Display` (ADR-0051, review finding
12) and a latent cost on every high-precision divide; `div.rs` had long
carried a note that a Newton or recursive divider was future work.

A sub-quadratic decimal base conversion (ADR-0051) is only sub-quadratic
if the division it rests on is too: divide-and-conquer over a quadratic
divider is still quadratic. So the conversion needed a sub-quadratic
divider underneath it.

## Decision

Add a Burnikel-Ziegler recursive divider and dispatch to it from
`divmod_limbs` by divisor size. It reduces division to multiplication,
`O(M(n)·log n)` over the existing Karatsuba `multiply_limbs`, and bottoms
out in the existing Knuth Algorithm D core for small sizes.

- `divmod_limbs` is now a size dispatcher. Divisors below
  `RECURSIVE_DIV_DISPATCH` (64 limbs) take the Knuth core (renamed
  `divmod_knuth`), which keeps the smaller constant there; larger
  divisors take `divmod_recursive`.
- `divmod_recursive` normalizes the divisor to a power-of-two limb count
  with its top bit set, via a combined whole-limb and sub-limb left
  shift. The shift leaves the quotient unchanged and only scales the
  remainder, which a final right shift undoes. A power-of-two limb count
  makes every recursion halve to an even size down to the Knuth base
  case, so no odd-size special case is needed.
- The mutually recursive `div_2n_1n` and `div_3n_2n` are the standard
  "two digits by one" and "three halves by two" steps. They use only the
  existing add, subtract, compare, shift, and multiply primitives, so the
  new code is the recursion structure rather than new arithmetic.

The implementation was derived from the published algorithm (the
Burnikel-Ziegler report and *Modern Computer Arithmetic* §1.4.3), not
adapted from a particular implementation; the decomposition, names, and
file layout are chosen fresh.

## Consequences

- The decimal base conversion and high-precision `div` / `parse` are now
  sub-quadratic in operand size. Normal-precision arithmetic is
  unchanged: its divisors stay below the dispatch threshold on the Knuth
  path.
- The recursive divider is verified differentially against the
  independent bit-at-a-time oracle (`divmod_limbs_bitwise`), the same
  oracle the Knuth core is checked against, plus a structural
  reconstruction sweep and adversarial divisors (all-ones, single-bit,
  alternating) that exercise the correction branches. The recursion's
  invariants are guarded by debug assertions: at most two corrections,
  remainder below divisor, no add or subtract overflow.
- Blast radius is contained by the threshold: only large-divisor
  divisions change path, and those are exactly the ones the new
  differential tests cover. The existing `div` and `parse` suites and the
  MPFR differential lanes pass unchanged.
- The normalization to a power-of-two limb count costs up to a factor of
  two in padding for a divisor just over a power of two. That constant is
  acceptable against the asymptotic win and keeps the recursion uniform.

## Related

- Plan: `plans/buzzing-stargazing-elephant.md` (pf-vbm2).
- Commit: `f2f79b0`.
- Code: `src/ops/limbs.rs` (`divmod_limbs`, `divmod_knuth`,
  `divmod_recursive`, `div_2n_1n`, `div_3n_2n`).
- References: Burnikel, C., and Ziegler, J. "Fast Recursive Division."
  MPI-I-98-1-022, 1998. Brent, R. P., and Zimmermann, P. *Modern Computer
  Arithmetic*, §1.4.3.
- Issues: `pf-vbm2`.
- Other ADRs: ADR-0051 (the formatter conversion that depends on this);
  ADR-0010 (the defer-invasive-work posture, the precedent for landing a
  general primitive when a v1.0 caller needs it).
