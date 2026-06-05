# ADR-0064: INEXACT flag fidelity for the proven-transcendence special functions

Status: Accepted (2026-06-05)

## Context

ADR-0060 corrected the IEEE 754-2019 §7.6 `INEXACT` flag for exp/log and
sin/cos; ADR-0063 extended the correction across the rest of the
elementary transcendental surface and `erf`/`erfc`. The pattern is: a
kernel dispatches its decidable exact-input set to `Status::OK` before the
Ziv loop, then forces `INEXACT` on the finite-normal fall-through. The
force is sound because, outside the dispatched set, the result is
irrational (Lindemann-Weierstrass / Gelfond-Schneider give it
unconditionally for the elementary functions), so the flag is correct even
where the working-precision evaluation collapses onto a grid value and the
rounding goes unobserved.

ADR-0063 deferred the heavier special functions; their exact-output sets
are subtler and a mistabulation would clear `INEXACT` wrongly. pf-cjnk
takes them up. It splits them by the strength of the irrationality
guarantee available, because that strength is no longer the clean
unconditional theorem the elementary surface enjoyed. This ADR covers the
families for which **no named open problem** obstructs the guarantee: the
integrals `Si`/`Ci`/`Ei`/`li`, the Bessel family `J`/`Y`/`I`/`K`, and the
Airy family `Ai`/`Bi`/`Ai′`/`Bi′`. The families whose guarantee rests on
an unsolved problem (`gamma`, `lgamma`, `digamma`, `beta`, `zeta`) are
deferred to a separate bead and its own residual-gap ADR (see Scope).

## Decision

Apply the ADR-0060 pattern to each kernel in scope: force `INEXACT` on the
finite-normal Ziv fall-through, guarded on a `Class::Normal` result so the
domain `qNaN` / `INVALID`, the poles' `±∞` / `DIV_BY_ZERO`, and the exact
non-finite limits keep their status. No new exact-input dispatch is needed
for any of these families; every dyadic-output input is already returned
by an existing special-case arm before the Ziv loop (see the audit below).

### The dispatched exact sets (the load-bearing audit)

The force is sound only if every input whose true value is exactly
representable at the target precision is dispatched to `OK` first. For
these families the dyadic-output inputs are few and all already
special-cased before the Ziv loop:

- Integrals: `Si(±0) = ±0` and `li(0) = 0`. (`Si(±∞) = ±π/2` is returned
  via `pi_over_2_at_round`, which already sets `INEXACT`; the poles
  `Ci(+0)`, `Ei(±0)`, `li(1)` and the exact limits `Ei(±∞)`, `Ci(+∞)`,
  `li(+∞)` are non-finite and so not `Class::Normal`.)
- Bessel: `J₀(0) = 1`, `Jₙ(0) = 0` (n ≠ 0), `I₀(0) = 1`, `Iₙ(0) = 0`
  (n ≠ 0); the `±∞` envelope limits; the `Y`/`K` poles at `+0`. The
  `J₀(0) = 1` and `I₀(0) = 1` values are themselves `Class::Normal`, but
  they are returned on the zero-input arm before the Ziv loop, so the
  force never sees them.
- Airy: none. Every Airy value is transcendental, including at `x = 0`
  (the boundary constants are rational multiples of `1/Γ(1/3)`,
  `1/Γ(2/3)`). The `x = 0` arm rounds the boosted-precision constant and
  is given the same force so a finite-normal result there cannot regress
  to a cleared flag.

### Edit shapes

(i) Most kernels bind `(result, status) = ziv_round(...)` then call
`auto_raise(status)`; the force is one guarded OR inserted between them
(`si`/`ci`/`ei`/`li`; `bessel_y`/`bessel_i`/`bessel_k`; the Airy main
path). (ii) `bessel_j` returns `ziv_round(...)` directly; it is
restructured to bind, force, raise, and return, which is the shape the
other kernels already use. (iii) The Airy `x = 0` arm is an early return
yielding a `Class::Normal` transcendental constant; the same guarded OR
is applied there before its `auto_raise`.

### Per-function classification

| Family | Fall-through value | Dispatched exact outputs |
|---|---|---|
| Si, Ci, Ei, li | transcendental at nonzero algebraic argument | Si(±0)=±0, li(0)=0 |
| J, I (integer order) | transcendental at nonzero algebraic argument | J₀(0)=1, Jₙ(0)=0, I₀(0)=1, Iₙ(0)=0 |
| Y, K (integer order) | transcendental at nonzero algebraic argument | (poles / limits only) |
| Ai, Bi, Ai′, Bi′ | transcendental at every algebraic argument incl. 0 | (none) |

### Soundness basis

The strongest rigorous statement is the one the force needs: for every
family in scope, no dyadic output is known at any non-dispatched dyadic
input, and no irrationality question for any of their structural values is
an open problem. This is the line separating this bead from the deferred
families, where a named open problem (the irrationality of Euler's
`γ = −ψ(1)`, of `ζ(5)`, and the transcendence of `Γ` at dyadic
denominators ≥ 8) sits exactly at a structural input.

For two of the families the guarantee is the classical transcendence
theorem, not merely the absence of a counterexample:

- **Bessel `Jₙ`, `Iₙ` (integer order)** are E-functions in the sense of
  Siegel: entire, of exponential type, with algebraic Taylor
  coefficients. The Siegel-Shidlovsky theorem gives `Jₙ(α)` and `Iₙ(α)`
  transcendental for every nonzero algebraic `α`. (Siegel proved the `J₀`
  case in 1929.)
- **Airy `Ai`, `Bi`** reduce to `₀F₁`, an E-function, at argument
  `x³ / 9`, so `Ai(α)`, `Bi(α)`, and their derivatives are transcendental
  for nonzero algebraic `α` by the same theory. At `x = 0` the boundary
  constants are rational multiples of `1/Γ(1/3)`, `1/Γ(2/3)`, with
  `Γ(1/3)` transcendental (Chudnovsky); they are therefore irrational.

The remaining families (`Y`/`K`, the second-kind Bessel solutions with a
logarithmic term at the origin, and the integrals `Si`/`Ci`/`Ei`/`li`,
which carry a logarithmic singularity) are not E-functions, so the clean
one-line argument does not apply. Their transcendence at nonzero algebraic
arguments is established through the transcendence theory of the
associated linear differential systems rather than the E-function theorem,
and, decisively for soundness, none of their structural values poses an
open irrationality question. Because the force only changes behavior on
the collapse cases (where the working-precision evaluation lands on a grid
value with no observed rounding), a false positive there would require an
actual dyadic value at a non-dispatched dyadic input; no such value is
known for any function in scope.

### Scope

Deferred to a follow-up bead (`pf-umlm`) and its own residual-gap ADR:
`gamma`, `lgamma`, `digamma`, `beta`, `zeta`. For these the force rests on
believed-but-unproven irrationality: `ψ(1) = −γ` (digamma), `ζ(5)`
(zeta), and `Γ` at dyadic denominators ≥ 8 (gamma family) are named open
problems. `beta` additionally needs a new construct-and-check exact
dispatch, because its positive-integer values `β(1, 2ᵏ) = 2⁻ᵏ` and
`β(1, 1) = 1` are genuinely dyadic and must be dispatched, not forced.
Isolating these keeps the open-problem disclosure out of the proven work.

### libm gate

None of the families in scope are in the v0.1 `LibmFnId` surface, so the
`pfloat-libm` `inexact_is_gated` widening from ADR-0060/0063 does not
extend to them. Validation is pfloat-side (`tests/inexact_fidelity.rs`
plus the `differential_*` MPFR lanes).

## Consequences

- `INEXACT` is reliable across the integrals, Bessel, and Airy families.
  No value changes; the fix is metadata, so the oracle `worst_ulp` is
  untouched and the `differential_*` lanes (which ignore status) confirm
  value preservation.
- The three families land as three signed-merge slices
  (integrals → Bessel → Airy) sharing this one ADR. The clean-first order
  establishes the edit and test pattern before the deferred families that
  need new dispatch logic.
- The weak-soundness special-function tables remain the open follow-up.

## Related

- pf-cjnk (this work); pf-umlm (the deferred weak-soundness families);
  ADR-0060 and ADR-0063 (the pattern and the elementary-surface
  discharge); ADR-0019 (integrals), ADR-0021 (Airy), ADR-0023 (Bessel) for
  the kernel designs and the `±∞` conventions.
- C. L. Siegel, *Über einige Anwendungen diophantischer Approximationen*,
  Abh. Preuss. Akad. Wiss., 1929 (E-functions; the `J₀` transcendence).
- A. B. Shidlovskii, *Transcendental Numbers*, de Gruyter, 1989 (the
  Siegel-Shidlovsky theorem).
- G. V. Chudnovsky, *Algebraic independence of the values of elliptic
  functions at algebraic points*, 1980 (`Γ(1/3)`, `Γ(1/4)` transcendental).
- A. Baker, *Transcendental Number Theory*, Cambridge University Press,
  1975.
