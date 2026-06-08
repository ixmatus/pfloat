# pfloat-complex rounding status

The per-operation rounding posture of pfloat-complex at 1.0. The complex analog
of `docs/rounding-status.md`: where the scalar table records correct rounding per
function across the five IEEE 754-2019 modes, this table records the
*componentwise* correct-rounding claim, the C99/C11 Annex G branch-cut
convention, and how each row is verified.

## The claim

Each operation rounds the real and imaginary parts **each correctly under their
own real rounding mode** (the model MPC uses; the only coherent strong rounding
claim for a type that carries no total order). A single `RoundingMode` argument
is applied to both components. The result `Status` is the OR-merge of the two
component statuses.

Branch selection and signed-zero discrimination are a documented **Annex G
convention layered on top of rounding**, not a rounding guarantee. The
load-bearing failure mode is a wrong-branch result when a caller supplies an
unsigned zero where the sign of zero was the only distinguishing information
(`csqrt(-4 + 0i) = +2i` but `csqrt(-4 - 0i) = -2i`).

## How rounding is achieved

- **A (delegation):** the component is exactly one already-verified
  correctly-rounded scalar kernel, so it is correctly rounded for free, with no
  Ziv loop.
- **B (directed-pair enclosure):** the component is a composition of
  transcendentals, enclosed by a `TowardNegative` / `TowardPositive` directed
  pair at a growing working precision (`GUARDS = [64, 128, 256, 512, 1024]`, cap
  five) and rounded once when both ends agree in value and sign (ADR-0091,
  shared in `src/enclosure.rs`, the same pattern as `div`). `INEXACT` is computed
  from whether the bracket collapsed, never forced, so exact algebraic outputs
  report `OK`.
- **F (fused single rounding):** one fused two-product (`mul_add_mul` /
  `mul_sub_mul`, ADR-0088), correctly rounded with a single rounding, no loop.

| Operation | Components | How | Annex G | Feature |
| --- | --- | --- | --- | --- |
| `add` / `sub` | re, im | scalar `add`/`sub` per part (A) | — | `big` |
| `neg` / `conj` | re, im | exact sign-bit flip (no rounding, no `Status`) | — | `big` |
| `norm_sqr` | real `T` | `re² + im²`, one fused rounding (F) | — | `big` |
| `mul` | re, im | `ac − bd`, `ad + bc`, each one fused rounding (F) | §G.5.1 infinity recovery | `big` |
| `div` | re, im | directed-pair enclosure Ziv loop (B) | §G.5.1 infinity recovery | `big` |
| `abs` | real `T` | `hypot(re, im)` (A) | §9.2.1 inf-dominates-NaN | `exp-log` |
| `sqrt` | re, im | Kahan robust enclosure (B) + axis-exact zeros | §G.6.4.2 branch cut | `exp-log` |
| `arg` | real `T` | `atan2(im, re)` (A) | §9.2.1 signed-zero / cut | `trig` |
| `to_polar` | `(abs, arg)` | A + A | as `abs`, `arg` | `trig` |
| `exp` | re, im | sign-aware product enclosure (B) + axis-exact zeros | §G.6.3.1 entire | `trig` |
| `log` | re | `ln(hypot)` enclosure (B); im = `atan2` (A) | §G.6.3.2 branch cut | `trig` |

## Verification

Every row is verified by the five-lane posture of ADR-0092:

- **Enumerated Annex G tables** (`tests/annex_g_special_values.rs`): every
  special-value row through the public API, across precisions and all five
  modes; the primary branch-cut and signed-zero guard.
- **Dispatch totality** (`tests/dispatch_totality.rs`): exhaustive over the
  finite IEEE class grid (a complete no-gap proof), plus the Annex G conjugation
  symmetry.
- **Algebraic identities** (`tests/identities.rs`): `csqrt(z)² = z`,
  `cexp(clog z) = z`, `clog(cexp z) = z`, `cexp(z+w) = cexp(z)cexp(w)`,
  oracle-free.
- **acb componentwise certified-rounding differential**
  (`tests/differential_acb.rs`, per release): the independent numeric pin
  against python-flint's rigorous complex Arb, bit-for-bit per component.
- **Kani** (`src/kani_harness.rs`, advisory): the componentwise `Status`-merge
  monoid; the BigFloat kernels are CBMC-hostile (ADR-0062) and rest on the lanes
  above.

## The Ziv-cap caveat

The (B) enclosures cap the working-precision schedule at five steps. A
hard-to-round near-zero component (the `clog` near `|z| = 1` and `cexp` near
`y = kπ/2` regimes) could in principle exhaust the cap and return a best-effort
value, the documented measure-zero MPFR caveat shared with `div`. The acb
differential probe in those exact regimes measured **zero** Ziv-cap residuals at
the 1.0 cut: every probed case certified and matched. This is a measured
quantity, not an assumed one.

## Deferred to v1.x (additive)

`sin` / `cos` / `tan`, the hyperbolics, and inverse trig with their Annex G
cuts; `pow` / `cis` / `from_polar`; the `clog` `log1p` reformulation tightening
the `|z| ≈ 1` band; per-component `Status` and per-component rounding modes
(ADR-0093).
