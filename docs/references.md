# References

This is the reference catalog for the pfloat family: the standards, papers, and
function definitions the implementations derive from. It gathers into one place
the citations the source already carries in module doc comments and the
architecture decision records, so a maintainer can check a kernel against its
primary source without hunting through the tree.

How to read the tables. A citation in the "cited in" or "module" column is a
place the source actually names the reference. The special functions table
distinguishes a DLMF pointer quoted from the module from a chapter pointer
curated here for a family that the module leaves implicit; the footnotes say
which is which. Where a function cites no DLMF section, the table says so rather
than inventing one: pfloat anchors a function's domain and special case
behavior to IEEE 754-2019, and reaches for a numerical handbook only for the
method.

## Standards

pfloat implements IEEE 754-2019
([registry entry](references/ieee-754-2019.md): pointer only, paywalled,
with the free proxy literature named). The clauses below are the ones
the source cites; each governs a load bearing behavior.

| Clause | Governs | Cited in |
|---|---|---|
| 4.3, 4.3.1 | The five rounding direction attributes and their rules | `src/rounding.rs`, ADR-0057 |
| 5.3.1 | `remainder` (the integer factor rounds to nearest even; the result is exact) | `src/ops/remainder.rs`, ADR-0069 |
| 5.4.1 | Correctly rounded square root | `src/ops/sqrt.rs` |
| 5.11 | Partial comparison (no order when an operand is NaN) | `pfloat-ball/src/scalar.rs` |
| 6.2, 7.2, 7.3 | Special value and exception dispatch (NaN, signed zero, poles) | `sqrt.rs`, `cbrt.rs`, `log1p.rs` |
| 6.5 | Rounding behavior, the funnel every kernel discharges | `src/rounding.rs` |
| 7.6 | The INEXACT exception, raised when the delivered value differs from the exact one | ADR-0057, ADR-0060 |
| 9.2, 9.2.1 | Recommended function special cases; `hypot`, `rootn`, `cbrt` | the elementary modules, `pfloat-libm` |
| 9.4, 9.4.3 | `gamma` special cases; `erf` and `erfc` special cases | `gamma.rs`, `erf.rs`, `erfc.rs` |

IEEE 1788-2015 ([registry entry](references/ieee-1788-2015.md)) is the
interval arithmetic standard. pfloat-ball v1.0 does not implement its
unbounded inf sup interval face; 1788 is named only as the boundary
marker for that future work (ADR-0076), and no 1788 clause numbers are
in tree. The complex surface anchors its branch cut and special value
behavior to C11 Annex G ([registry entry](references/iso-c-annex-g.md));
the special function methods trace to the DLMF
([registry entry](references/dlmf.md)).

## Algorithms and papers

Each named external work has a registry entry under `docs/references/`
carrying its full provenance: the verified citation, canonical and
archived URLs, license, rot risk class, and the code paths that consume
it. This list is the crosswalk; the entry is the record.

- Ziv's correct rounding strategy, the driver every transcendental lane
  runs through: [ziv-1991](references/ziv-1991.md). `math/ziv.rs`.
- The gamma family approximation: [spouge-1994](references/spouge-1994.md),
  with the error analysis in [pugh-2004](references/pugh-2004.md) and the
  evaluation pattern from [toth-gamma-2005](references/toth-gamma-2005.md).
  `gamma_stirling.rs`.
- Zeta acceleration:
  [cohen-rodriguez-villegas-zagier-2000](references/cohen-rodriguez-villegas-zagier-2000.md)
  and [borwein-zeta-1995](references/borwein-zeta-1995.md). `zeta.rs`, ADR-0026.
- Multiplication and division kernels:
  [brent-zimmermann-mca-2010](references/brent-zimmermann-mca-2010.md),
  [bodrato-zanoni-2007](references/bodrato-zanoni-2007.md) (the registry entry
  corrects a venue conflation the tree carried),
  [jebelean-1993](references/jebelean-1993.md),
  [knuth-taocp-v2](references/knuth-taocp-v2.md), and
  [burnikel-ziegler-1998](references/burnikel-ziegler-1998.md).
  `ops/limbs.rs`, ADR-0052, ADR-0061.
- Trig argument reduction: [payne-hanek-1983](references/payne-hanek-1983.md).
  `math/trig_reduce.rs`.
- AGM constants: [brent-1976](references/brent-1976.md) (with Salamin's
  companion result) and [brent-mcmillan-1980](references/brent-mcmillan-1980.md).
  `math/agm.rs`, `math/agm_constants.rs`, ADR-0015, ADR-0017, ADR-0018.
- The shortest decimal formatter:
  [steele-white-1990](references/steele-white-1990.md) and
  [burger-dybvig-1996](references/burger-dybvig-1996.md). `fmt.rs`, ADR-0071.
- The hard to round corpus and its lineage:
  [core-math-wc-corpus](references/core-math-wc-corpus.md) (the data pfloat
  samples), descending from
  [lefevre-muller-2001](references/lefevre-muller-2001.md) through
  [vinc17-testlibm](references/vinc17-testlibm.md) and
  [stehle-lefevre-zimmermann-2005](references/stehle-lefevre-zimmermann-2005.md),
  surveyed in
  [sibidanov-zimmermann-core-math-2023](references/sibidanov-zimmermann-core-math-2023.md)
  and consolidated in print in
  [muller-handbook-fp-2018](references/muller-handbook-fp-2018.md).
  `docs/lefevre-muller-corpus-provenance.md`, ADR-0083.
- The INEXACT discipline's transcendence anchor:
  [baker-transcendental-1975](references/baker-transcendental-1975.md).
  ADR-0060, ADR-0063, ADR-0064.
- The verification oracles: [mpfr](references/mpfr.md),
  [arb](references/arb.md), [mpmath](references/mpmath.md),
  [maxima](references/maxima.md); the complex design reference
  [mpc](references/mpc.md); and the deliberately absent conformance suite
  [berkeley-testfloat](references/berkeley-testfloat.md), whose entry records
  why no third party vector set applies at arbitrary precision.

## Special functions

The canonical reference and the implemented method, per function. The DLMF
column quotes the chapter or section the module cites; `none cited` means the
module names no DLMF section and the IEEE clause in parentheses is the reference
for the function's special case behavior. See the footnotes for the gamma family
and Airy.

| Function | DLMF | Method | Module | ADR |
|---|---|---|---|---|
| `exp` | none cited (IEEE 9.2) | ln(2) range reduction, Maclaurin, recompose by `2^k` | `src/math/exp.rs` | 0022 |
| `ln` | none cited (IEEE 9) | exponent reduction, `atanh` series `ln(m) = 2 atanh((m-1)/(m+1))` | `src/math/ln.rs` | 0022 |
| `log2`, `log10` | none cited (via `ln`) | `ln(x)` divided by `ln(2)` or `ln(10)` at target plus 64 bits | `src/math/log2.rs`, `log10.rs` | 0022 |
| `expm1`, `log1p` | none cited (IEEE 9.2, 7.3) | cancellation boosted core with a tiny argument short circuit | `src/math/expm1.rs`, `log1p.rs` | 0038 |
| `pow` | none cited (IEEE 9.2.1) | square and multiply for integer exponent, `exp(y ln x)` otherwise | `src/math/pow.rs` | 0022 |
| `sin`, `cos`, `tan` | none cited (IEEE 9.2) | Payne Hanek reduction, quadrant dispatch, Maclaurin | `src/math/sin.rs`, `cos.rs`, `tan.rs` | 0038 |
| `cot`, `sec`, `csc` | 4.14 | direct kernels, shared Payne Hanek reduction, single rounding | `src/math/trig_reciprocal.rs` | 0032, 0056 |
| `asin`, `acos`, `atan` | none cited (IEEE 9.2) | half angle `atan` identities, reduction, Maclaurin | `src/math/asin.rs`, `acos.rs`, `atan.rs` | 0038 |
| `sinh`, `cosh`, `tanh` | none cited (IEEE 9.2) | stable `expm1` and `exp` compositions | `src/math/sinh.rs`, `cosh.rs`, `tanh.rs` | 0022, 0050 |
| `asinh`, `acosh`, `atanh` | none cited (IEEE 9.2) | `log1p` reformulations | `src/math/asinh.rs`, `acosh.rs`, `atanh.rs` | 0038, 0059 |
| `sqrt` | none cited (IEEE 5.4.1) | integer square root of the shifted mantissa, single rounding | `src/ops/sqrt.rs` | (none) |
| `cbrt` | none cited (IEEE 9.2.1) | integer cube root, single rounding | `src/ops/cbrt.rs` | 0056 |
| `gamma`, `lgamma`, `digamma` | Chapter 5 (curated) | reflection, recurrence shift, Stirling or Spouge | `src/math/gamma.rs`, `lgamma.rs`, `digamma.rs` | 0022, 0038, 0066 |
| `beta` | 5.12.1 (with 5.2, 5.5.3) | `lgamma` composition with a negative domain dispatch | `src/math/beta.rs` | 0030, 0038 |
| `erf`, `erfc` | none cited (IEEE 9.4.3) | Maclaurin for small `x`, divergent asymptotic for large `x` | `src/math/erf.rs`, `erfc.rs` | 0022, 0063 |
| `J`, `Y`, `I`, `K` (Bessel) | Chapter 10 | Maclaurin or log series, Miller or upward recurrence, Hankel asymptotic | `src/math/bessel_{j,y,i,k}.rs` | 0023, 0024, 0025 |
| `Ai`, `Bi`, `Ai'`, `Bi'` (Airy) | Chapter 9 | Maclaurin for small `x`, exponential and oscillatory asymptotics | `src/math/airy.rs` | 0021, 0048 |
| `Ei`, `Si`, `Ci` | 6.2 | convergent series, divergent or auxiliary function asymptotic | `src/math/ei.rs`, `si.rs`, `ci.rs` | 0038, 0064 |
| `zeta` | Chapter 25 (25.2.3, 25.4.2) | Borwein and Cohen Rodriguez-Villegas Zagier acceleration, functional equation | `src/math/zeta.rs` | 0026 |

The arithmetic and format kernels carry their own references: `remainder`
(IEEE 5.3.1, ADR-0069), `convert` (IEEE 4.3, ADR-0055), `rounding` (IEEE 4.3,
4.3.1, 6.5), the limb multiplier (Toom-3, Karatsuba, Knuth Algorithm D,
ADR-0061), and the decimal formatter (Steele White, Burger Dybvig, ADR-0071).

Footnotes.

- Gamma family chapter pointer is curated, not quoted. `gamma.rs`, `lgamma.rs`,
  and `digamma.rs` cite no DLMF section in tree; the family is unambiguously
  DLMF Chapter 5, so the chapter is listed as a reading aid. Only `beta.rs`
  carries quoted section numbers (5.2, 5.5.3, 5.12.1).
- The Airy asymptotic `u_k` recurrence has a known transcription pitfall. The
  defining product is DLMF 9.7.2; the correct ratio carries a `(2k - 1)` factor
  in the denominator, which the code implements. One doc comment in `airy.rs`
  states the bare `216 k` form without that factor; trust the code and DLMF
  9.7.2, not that one comment.

## Interval arithmetic (pfloat-ball)

The ball references are theorems and a design oracle, not equation numbers.

- The Fundamental Theorem of Interval Arithmetic is the soundness law: an
  interval extension of a function contains the function's image over the input
  interval. pfloat-ball states its consequences as five laws in
  `pfloat-ball/src/spec.rs` (ADR-0076, ADR-0077). There is no numbered clause to
  cite; the in tree spec is the canonical statement.
- Arb is the principal midpoint and radius design reference (Fredrik Johansson,
  "Arb: Efficient Arbitrary-Precision Midpoint-Radius Interval Arithmetic", IEEE
  Transactions on Computers, 2017;
  [registry entry](references/johansson-arb-2017.md)). The source names
  Johansson's work by its behavior and identifiers (`mag_t`, relative accuracy
  bits, `printn`) rather than by a transcribed title; the bibliographic anchor
  is supplied here.
  pfloat-ball diverges from Arb deliberately: a 64 bit radius significand where
  Arb uses 30 bits, for tightness and Kani verifiability (ADR-0074).
- The transcendence theorems (Lindemann Weierstrass, Gelfond Schneider; Baker
  1975) justify the INEXACT discipline a priori: the result of a transcendental
  kernel on a non trivial algebraic input is irrational, so it cannot land
  exactly on the target grid. See ADR-0060, ADR-0063, ADR-0064.

## See also

- `docs/references/INDEX.md` for the per source registry this file links
  into (one entry per external source, with archived URLs, licenses, rot
  risk, and consumers; schema in `docs/references/README.md`, decision in
  ADR-0094).
- `docs/references/coverage-gaps.md` for the reconciliation of the README
  disclosures' named failure modes against the recorded vector coverage.
- `docs/algorithms.md` for the narrative reading guide into these algorithms.
- `docs/decisions/` for the architecture decision records cited above.
- `docs/lefevre-muller-corpus-provenance.md` for the hard to round corpus
  provenance and its MIT license posture.
