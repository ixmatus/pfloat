# ADR-0030: Beta function on the negative real domain (sign and pole convention)

- **Status**: accepted
- **Date**: 2026-05-18

## Context

Roadmap slice 8a closes the inline TODO in `src/math/beta.rs`. Slice
4c restricted `beta(a, b)` to `a, b > 0`: the kernel computes the
magnitude through `lgamma` as `exp(lgamma(a) + lgamma(b) -
lgamma(a + b))`, and `lgamma` returns only `ln|Gamma|`, so any input
that could make a `Gamma` factor negative or hit a pole was coerced
to `qNaN + INVALID`. The TODO scoped the follow up as "negative non
integer inputs that produce a well defined result", deferred for
want of explicit sign tracking. v1.0 is a completeness boundary, so
the extension also covers the integer pole cancellation values
(decided below), giving `beta` the full real domain rather than a
narrow slice of it.

The sign and pole structure is recalled math until derived from the
primary source. This is the sixth instance of the pfloat derive do
not recall discipline (LN2, 2/sqrt(pi), the Airy `u_k` recurrence,
the Bessel `K` recurrence sign, the zeta algorithm assumption, now
the beta sign and pole rule). It was derived from DLMF, not recalled,
and every case was pinned against mpmath before any kernel code was
written.

### Primary sources (DLMF, fetched 2026-05-18)

- **5.12.1**: `B(a,b) = integral_0^1 t^(a-1) (1-t)^(b-1) dt =
  Gamma(a) Gamma(b) / Gamma(a+b)`. The integral needs `Re a > 0`,
  `Re b > 0`; the `Gamma` quotient is the analytic continuation that
  defines `B` elsewhere.
- **5.5.1**: `Gamma(z+1) = z Gamma(z)`.
- **5.5.3**: `Gamma(z) Gamma(1-z) = pi / sin(pi z)`, for
  `z != 0, +-1, ...`.
- **5.2**: `Gamma` is meromorphic with **no zeros**, and simple
  poles of residue `(-1)^n / n!` at `z = -n` for `n = 0, 1, 2, ...`.
  Equivalently `1 / Gamma` is entire with simple zeros at the non
  positive integers.

### Derived facts

**Sign of `Gamma` on the reals.** For `x > 0` the Euler integral
(5.2.1) has a positive integrand, so `Gamma(x) > 0`. For `x < 0`
non integer, 5.5.3 gives `Gamma(x) = pi / (sin(pi x) Gamma(1-x))`;
here `1 - x > 1`, so `Gamma(1-x) > 0`, and `pi > 0`, hence
`sign Gamma(x) = sign sin(pi x)`. This is exactly the existing
private `gamma_sign_of` in `src/math/gamma.rs` (positive for `x > 0`,
`sign sin(pi x)` otherwise); the derivation confirms that helper
rather than reproducing it.

**Sign of `B`.** Wherever `B` is finite and non zero,
`sign B(a,b) = sign Gamma(a) . sign Gamma(b) . sign Gamma(a+b)`,
a product of three `gamma_sign_of` values. Verified against mpmath
for the finite cases below (e.g. `B(-0.5, 0.75) < 0`,
`B(-2.5, -1.25) < 0`, `B(-0.5, 0.25) > 0`).

**Pole and zero structure**, from 5.2 (poles only at non positive
integers, no zeros):

- `a + b` a non positive integer with `a, b` not poles: `Gamma(a+b)`
  is a denominator pole, `1/Gamma(a+b) = 0`, so `B = 0` exactly
  (mpmath: `B(-0.5, 0.5) = B(-0.5, -0.5) = B(-1.5, 0.5) = 0`).
- `a` (or `b`) a non positive integer, no compensating `a+b` pole:
  numerator `Gamma` pole, `B` diverges. The two sided limit changes
  sign across the integer, so no single signed infinity is correct.
- `a = -n` (`n >= 0` integer) and `b = m` (`m >= 1` integer) with
  `1 <= m <= n` so `a + b = m - n` is also a non positive integer:
  numerator and denominator poles cancel and `B` is finite. Using
  the 5.2 residues, `Gamma(a) ~ (-1)^n / (n! (a+n))` and
  `Gamma(a+b) ~ (-1)^(n-m) / ((n-m)! (a+n))` as `a -> -n`, the
  `(a+n)` factors cancel and

  ```text
  B(-n, m) = (-1)^m (m-1)! (n-m)! / n!,   integers n >= 0, 1 <= m <= n
  ```

  Symmetric in its arguments. Verified exactly against mpmath for
  `(n,m) in {(3,2),(3,1),(5,3),(2,2),(4,4),(6,1),(1,1),(5,5),
  (10,3),(7,4)}` and `B(b,a) = B(a,b)`; `m > n` correctly diverges.

## Decision

`beta` accepts the full real domain. The kernel classifies inputs
and the public behavior is, with `Zle = {0, -1, -2, ...}`:

| # | Condition | Result | Status |
|---|-----------|--------|--------|
| 1 | `a > 0`, `b > 0` | `exp(lgamma sum)`, positive | `OK` |
| 2 | `a not in Zle`, `b not in Zle`, `a+b not in Zle` | `sign . exp(lgamma(a)+lgamma(b)-lgamma(a+b))`, `sign = gamma_sign_of(a).gamma_sign_of(b).gamma_sign_of(a+b)` | `OK` |
| 4 | `a` or `b` in `Zle`, the other a positive integer, `a+b in Zle` (pole cancellation) | `(-1)^m (m-1)! (n-m)! / n!` via the closed form (order independent by symmetry) | `OK` |
| 5 | `a+b in Zle`, `a not in Zle`, `b not in Zle` | `+0` (denominator pole) | `OK` |
| 3 | `a` or `b` a negative integer, not a case 4 cancellation | `qNaN + INVALID` (pole, two sided sign ambiguous) | `INVALID` |
| 0 | `a = +-0` or `b = +-0`, not a case 4 cancellation | `+-inf + DIV_BY_ZERO` (signed, mirrors `gamma(+-0)`) | `DIV_BY_ZERO` |
| 6 | `a in Zle` and `b in Zle` | `qNaN + INVALID` (net pole) | `INVALID` |

Rationale for the pole rows: they mirror the established gamma family
convention already shipped (`src/math/gamma.rs`:
`gamma(negative integer) = qNaN + INVALID` because the two sided
limit is sign ambiguous; `gamma(+-0) = +-inf + DIV_BY_ZERO` because a
signed zero carries a side). `beta` inherits that convention so the
two functions agree at their shared poles. Case 4 takes precedence
over rows 3 and 0 when its condition holds (the cancellation makes
the value finite).

Reuse, no duplication: `gamma_sign_of` is widened from `fn` to
`pub(super) fn` (both modules are children of `src/math/`); the
reflection derivation lives in one place. Integer detection reuses
`super::lgamma::is_integer_test` (already used by the gamma kernel),
and the case-4 `(-1)^m` sign reuses `super::pow::integer_parity`
(also widened to `pub(super)`, the same precedent). The magnitude
path for cases 1, 2 is the existing `lgamma` composition, unchanged.
Case 4 must not pass a non-positive integer to `lgamma`
(`lgamma(-n)` is `+inf`), but it *does* route through `lgamma` of
the three *positive*-integer factorials: with
`B(-n,m) = (-1)^m (m-1)!(n-m)!/n!`, the magnitude is
`exp(lgamma(m) + lgamma(n-m+1) - lgamma(n+1))` (arguments `m`,
`n-m+1`, `n+1`, all `>= 1`, never a `Gamma` pole) and the sign is
`(-1)^m` from the parity of `m`. This is `O(1)`.

Correction (robustness fix, recorded not silently applied). This
ADR originally specified case 4 as the reciprocal product
`(-1)^m (1/m) prod_{i=0}^{m-1} (i+1)/(n-i)`, "computed in
reciprocal-product form to avoid forming the huge binomial ...
for the rare very large `m` this trades speed for exactness, the
correct call here". That was exact but ran `m` iterations, and `m`
is a caller-supplied integer: `B(-2e18, 1e18)` is a legal call that
spins ~`1e18` BigFloat iterations and does not terminate. The
"trades speed for exactness" framing understated a caller-reachable
unbounded-resource-consumption defect (CLAUDE.md security posture).
The `lgamma`-of-factorials form above is the same value, still
total, and `O(1)`; it replaces the loop. The accuracy of case 4
changes from exact rational arithmetic to the same `~p`-bit
`lgamma`-composition accuracy as the negative-domain magnitude path
(cases 1, 2) — the right trade-off, since the alternative is
non-termination. This is the same recalled-design failure mode the
derive-do-not-recall discipline targets: a plan/ADR algorithm
choice (here "use the reciprocal product") is recalled math too;
its cost bound (`O(m)`, `m` unbounded) was not derived against the
input domain until the fuzz slice forced it. Pinned by
`beta_case4_large_m_terminates` (`m = 1e12`, returns) and
`beta_case4_factorial_exact_rational` (`B(-10,4) = 1/840`,
hand-derived).

The `NaN` / signaling `NaN` / `+-inf` handling already in
`beta_kernel` is unchanged. `B` with an infinite operand and a
negative finite operand is not part of this extension and keeps the
existing conservative behavior (only the documented
`beta(+inf, finite positive) = +0` is asserted); this is recorded as
a known narrow edge in DESIGN.md rather than silently widened.

## Consequences

- `beta` is correctly signed and total over the reals at v1.0; no
  user has to route around a positive only restriction. Completeness
  is a frugality property, so the integer pole cancellation values
  ship now rather than as a deferred 1.x slice.
- The differential oracle must track sign. MPFR has no native `beta`;
  the existing lane composes `ln_gamma_ref`, which loses sign. The
  oracle is extended to carry the reflection sign (or MPFR's
  `lgamma` sign output) so the negative domain is checked, not just
  the positive arm.
- `lgamma` is unchanged. The sign lives entirely in `beta`; pushing
  a sign out of `lgamma` would perturb every `lgamma` consumer for
  no benefit.
- A second 8a concern, the parser exponent cap
  (`MAX_DECIMAL_EXPONENT`), is a separate behavior change in the same
  slice; it closes an inline TODO and is documented in DESIGN.md, it
  does not get its own ADR. Recorded here so the slice's two behavior
  changes are discoverable from one place.
- Case 4's implementation approach was corrected post-ship from an
  `O(m)` reciprocal-product loop (unbounded on a caller-supplied
  `m`) to the `O(1)` `lgamma`-of-factorials form, as a labelled
  robustness fix. See the "Correction" note in the Decision section;
  it is also the recalled-algorithm-cost variant of the derive do
  not recall lesson (the cost bound of a recalled algorithm choice
  must be derived against the input domain, not assumed).
- This is the sixth derive do not recall instance and the second
  (after zeta) where the resolution required reproducing the source's
  own structure (here the 5.2 residues) and pinning every case
  against mpmath before writing code. The lesson is recorded in the
  derive do not recall memory, not only here.

## Related

- Closes the `src/math/beta.rs` slice 4c TODO.
- Reuses `src/math/gamma.rs` `gamma_sign_of` (visibility widened)
  and `src/math/lgamma.rs` `is_integer_test`.
- Convention precedent: ADR-0019 (integral specials pole
  convention), ADR-0021 (Airy), ADR-0025 (Bessel `I`/`K` sign),
  ADR-0026 (zeta domain convention).
- Plan: `plans/abundant-yawning-badger.md` (slice 8a).
- Other ADRs: none superseded.
