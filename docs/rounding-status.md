# pfloat rounding status

Per-function correct-rounding status across the five IEEE 754-2019
rounding modes (NE NearestEven, NA NearestAway, TZ TowardZero, TP
TowardPositive, TN TowardNegative), for the 63 functions the
verification oracle tracks. Generated from the status records under
`tests/oracle/status/` by `scripts/rounding-status-table.sh` and
checked in CI; the records are the single source of truth.

Across all rows the worst observed error is 0 ULP, with 0
mismatches and 0 panics over the sampled input grids.

Legend:

- `CR`: correctly rounded, certified across the exhaustive binary32
  input grid. Every binary32 value is computed at high working
  precision and rounded to the target, then compared bit for bit
  against the oracle.
- `CR(d)`: correctly rounded, certified by the five-mode differential
  lanes against MPFR (with the synthesized NearestAway oracle, which
  MPFR lacks a primitive for) and reconfirmed by the per-release
  cross-check sweep (ADR-0049). The exhaustive binary32 oracle's
  bf-to-f32 bridge carries NearestEven only, so it does not sweep this
  directed mode; the guarantee for these cells rests on the lanes.

The oracle column names the primary rigorous backend: MPFR where it has
a primitive for the function, Arb otherwise.

| Function | Oracle | NE | NA | TZ | TP | TN |
| --- | --- | --- | --- | --- | --- | --- |
| acos | MPFR | CR | CR | CR | CR | CR |
| acosh | MPFR | CR | CR | CR | CR | CR |
| Ai | MPFR | CR | CR | CR | CR | CR |
| Ai_prime | Arb | CR | CR | CR | CR | CR |
| asin | MPFR | CR | CR | CR | CR | CR |
| asinh | MPFR | CR | CR | CR | CR | CR |
| atan | MPFR | CR | CR | CR | CR | CR |
| atanh | MPFR | CR | CR | CR | CR | CR |
| Bi | Arb | CR | CR | CR | CR | CR |
| Bi_prime | Arb | CR | CR | CR | CR | CR |
| Ci | Arb | CR | CR | CR | CR | CR |
| cos | MPFR | CR | CR | CR | CR | CR |
| cosh | MPFR | CR | CR | CR | CR | CR |
| digamma | MPFR | CR | CR | CR | CR | CR |
| Ei | MPFR | CR | CR | CR | CR | CR |
| erf | MPFR | CR | CR(d) | CR(d) | CR(d) | CR(d) |
| erfc | MPFR | CR | CR | CR | CR | CR |
| exp | MPFR | CR | CR(d) | CR(d) | CR(d) | CR(d) |
| exp10 | MPFR | CR | CR | CR | CR | CR |
| exp2 | MPFR | CR | CR | CR | CR | CR |
| expm1 | MPFR | CR | CR | CR | CR | CR |
| gamma | MPFR | CR | CR | CR | CR | CR |
| I0 | Arb | CR | CR | CR | CR | CR |
| I1 | Arb | CR | CR | CR | CR | CR |
| In (n=2) | Arb | CR | CR | CR | CR | CR |
| In (n=5) | Arb | CR | CR | CR | CR | CR |
| In (n=10) | Arb | CR | CR | CR | CR | CR |
| In (n=25) | Arb | CR | CR | CR | CR | CR |
| In (n=100) | Arb | CR | CR | CR | CR | CR |
| J0 | MPFR | CR | CR(d) | CR(d) | CR(d) | CR(d) |
| J1 | MPFR | CR | CR(d) | CR(d) | CR(d) | CR(d) |
| Jn (n=2) | MPFR | CR | CR | CR | CR | CR |
| Jn (n=5) | MPFR | CR | CR(d) | CR(d) | CR(d) | CR(d) |
| Jn (n=10) | MPFR | CR | CR | CR | CR | CR |
| Jn (n=25) | MPFR | CR | CR | CR | CR | CR |
| Jn (n=100) | MPFR | CR | CR | CR | CR | CR |
| K0 | Arb | CR | CR | CR | CR | CR |
| K1 | Arb | CR | CR | CR | CR | CR |
| Kn (n=2) | Arb | CR | CR | CR | CR | CR |
| Kn (n=5) | Arb | CR | CR | CR | CR | CR |
| Kn (n=10) | Arb | CR | CR | CR | CR | CR |
| Kn (n=25) | Arb | CR | CR | CR | CR | CR |
| Kn (n=100) | Arb | CR | CR | CR | CR | CR |
| lgamma | MPFR | CR | CR(d) | CR(d) | CR(d) | CR(d) |
| li | Arb | CR | CR | CR | CR | CR |
| ln | MPFR | CR | CR(d) | CR(d) | CR(d) | CR(d) |
| log10 | MPFR | CR | CR(d) | CR(d) | CR(d) | CR(d) |
| log1p | MPFR | CR | CR | CR | CR | CR |
| log2 | MPFR | CR | CR(d) | CR(d) | CR(d) | CR(d) |
| Si | Arb | CR | CR | CR | CR | CR |
| sin | MPFR | CR | CR | CR | CR | CR |
| sinh | MPFR | CR | CR | CR | CR | CR |
| sqrt | MPFR | CR | CR(d) | CR(d) | CR(d) | CR(d) |
| tan | MPFR | CR | CR | CR | CR | CR |
| tanh | MPFR | CR | CR(d) | CR(d) | CR(d) | CR(d) |
| Y0 | MPFR | CR | CR | CR | CR | CR |
| Y1 | MPFR | CR | CR | CR | CR | CR |
| Yn (n=2) | MPFR | CR | CR | CR | CR | CR |
| Yn (n=5) | MPFR | CR | CR | CR | CR | CR |
| Yn (n=10) | MPFR | CR | CR | CR | CR | CR |
| Yn (n=25) | MPFR | CR | CR | CR | CR | CR |
| Yn (n=100) | MPFR | CR | CR | CR | CR | CR |
| zeta | MPFR | CR | CR(d) | CR(d) | CR(d) | CR(d) |
