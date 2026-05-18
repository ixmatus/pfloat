# ADR-0024: Bessel functions of the second kind Y0, Y1, Yn

- **Status**: accepted
- **Date**: 2026-05-17

## Context

Roadmap slice 6p is the companion to slice 6o (Bessel `J`,
ADR-0023). It adds `Y0`, `Y1`, and `Yn` for integer order and real
argument under the existing `bessel = ["specials"]` cluster feature,
closing the ordinary-Bessel pair (6q adds modified `I`/`K`, 6r adds
zeta). ADR-0023 is the structural template: the same public-API
shape, binary-exponent regime dispatch, smallest-term asymptotic
truncation, decaying-envelope limit convention, and bit-exact MPFR
tiered oracle. 6p also re-enables the J/Y Wronskian cross-tie that
6o could not exercise because `Y` did not yet exist.

Unlike `J`, which is entire on the real line, `Y` has a logarithmic
branch point at the origin and is complex for `x < 0`, so the
domain is `x > 0` only. `Y` is also the *dominant* solution of the
Bessel three-term recurrence, which inverts the recurrence strategy
relative to `J`.

All formulas were derived from the DLMF primary source
(user-authorized `WebFetch` of DLMF §10.8, §10.17, §10.5 at plan
time), not recalled. The 6n Airy `(2k−1)`-divisor defect is the
precedent that makes this discipline load-bearing.

## Decision

**Two-regime base pair plus upward recurrence**, dispatched on the
binary exponent of `x` (the `airy`/`si`/`bessel_j` integer-exponent
selector idiom):

- *Below the asymptotic cut*: the DLMF 10.8.1 logarithmic series for
  the base pair `Y0`/`Y1`.
- *At or above the cut*: the DLMF 10.17.4 Hankel asymptotic.
- *Order `n ≥ 2`*: the DLMF 10.6.1 upward recurrence from
  `(Y0, Y1)`.

**Upward recurrence, chosen over Miller backward descent** (this
resolves the 6p analog of 6o open-decision #4). `Y` is the dominant
solution of `𝒞_{ν−1}+𝒞_{ν+1} = (2ν/z)𝒞_ν`, so the forward form
`Y_{k+1}(x) = (2k/x)·Y_k(x) − Y_{k−1}(x)` is numerically stable: the
`(2k/x)·Y_k` term dominates the `−Y_{k−1}` subtraction, which
therefore never cancels catastrophically. `J`'s Miller backward
descent and its sum-rule normalization do not apply here: there is
no recessive solution to renormalize. The climb is a rolling pair,
`O(m)` time and `O(1)` space, no `Vec`. As a consequence `Y` has no
cheap middle regime: where `J` uses Miller for moderate `|x|`, `Y`
runs the DLMF 10.8.1 series with an `≈ x·log₂e` cancellation guard
(the `ci.rs`/`si.rs` capped form) across the entire sub-asymptotic
range, which is costlier than `J`'s moderate path. This is the
unavoidable price of the dominant solution lacking a normalization
trick.

**The DLMF 10.8.1 logarithmic series, with the digamma terms
reduced to elementary sums.** For integer `n ≥ 0`,

```text
Y_n(x) = −(x/2)^{−n}/π · Σ_{k=0}^{n−1} (n−k−1)!/k! (x²/4)^k
       + (2/π) ln(x/2) · J_n(x)
       − (x/2)^n/π · Σ_{k≥0} (ψ(k+1)+ψ(n+k+1)) (−x²/4)^k/(k!(n+k)!)
```

The finite head sum is empty for `n = 0` and is the source of the
`n ≥ 1` pole. The digamma terms reduce to running harmonic sums:
`ψ(k+1) = −γ + H_k` and `ψ(n+k+1) = −γ + H_{n+k}`, so
`ψ(k+1)+ψ(n+k+1) = −2γ + H_k + H_{n+k}` (`H_j = Σ_{i=1}^{j} 1/i`,
`H_0 = 0`). No digamma kernel is needed: only harmonic partial sums
and the in-tree Euler–Mascheroni `γ` (slice 6m0, also used by
`ei`/`ci`/`li`). The `J_n(x)` piece is 6o's `bessel_j` kernel. The
reduction is cross-pinned by the DLMF 10.8.1 ↔ 10.8.2 algebraic
identity: substituting `−2γ + 2H_k` (the `n = 0` case) into the
tail and using `Σ (−x²/4)^k/(k!)² = J_0(x)` folds the `−2γ J_0`
into a `+γ` and (since `H_0 = 0`) kills the `k = 0` tail term,
exactly reproducing DLMF 10.8.2 `Y_0 = (2/π)(ln(x/2)+γ)J_0 +
(2/π)[(x²/4)/(1!)² − …]`. The two DLMF forms agreeing is the
worked check.

**The DLMF 10.17.4 asymptotic reuses 6o's coefficients verbatim.**
The `a_k(n)` recurrence (`a_0=1`,
`a_k = a_{k−1}(4n²−(2k−1)²)/(8k)`) and the phase
`ω = x − nπ/2 − π/4` are identical to `J`'s (DLMF 10.17.1/10.17.2)
and were already derived and Pochhammer-cross-pinned at `k=1,2` in
ADR-0023; 6p references that pin rather than re-deriving. `Y`
differs from `J` only in the trig combination: DLMF 10.17.4 is
`√(2/πx)·[sinω·ΣP + cosω·ΣQ]` where `J`'s 10.17.3 is
`√(2/πx)·[cosω·ΣP − sinω·ΣQ]` (`ΣP = Σ(−1)^k a_{2k}/x^{2k}`,
`ΣQ = Σ(−1)^k a_{2k+1}/x^{2k+1}`). Folding the explicit `(−1)^k`
into the trig assignment yields `Y`'s period-4 factor
`[+sinω, +cosω, −sinω, −cosω]` on `a_j(n)/x^j` (vs `J`'s
`[+cosω, −sinω, −cosω, +sinω]`), summed to the smallest term. The
asymptotic threshold is shared with `J` (`bessel_j_threshold`,
already `pub(super)`): the DLMF 10.17.4 error has the same
`e^{−2|x|}` order as DLMF 10.17.3, so the same conservative cut is
strictly more than enough.

**Domain and limit conventions, the `Ci`/`li` precedent.** `Y` is
real only for `x > 0`. The table, Kani-pinned:

| input | result | status |
|-------|--------|--------|
| `qNaN` | `qNaN` (payload propagated) | `OK` |
| `sNaN` | `qNaN` | `INVALID` |
| `Y_n(+0)` | `−∞` | `DIV_BY_ZERO` (a pole, DLMF 10.8.1) |
| `x < 0`, `−0`, `−∞` | `qNaN` | `INVALID` |
| `Y_n(+∞)`, any `n` | `+0` | `OK` |

`Y_n(+0) = −∞` raising `DIV_BY_ZERO` is the pole convention shared
with `Ci(+0)` and `li(1)`. `x < 0` (and `−0`, `−∞`) giving `qNaN` +
`INVALID` is the "complex in the reals" convention shared with
`Ci`/`li`. `Y_n(+∞) = +0` is the decaying-envelope convention
(ADR-0021/0023, the Airy and `J` precedent): the true behaviour at
`+∞` is a bounded decaying oscillation with no limit, so the
conservative total result keeps the function total. Negative order
reduces before regime dispatch: `Y_{−n}(x) = (−1)^n Y_n(x)`
(DLMF 10.4.1), so the kernel evaluates `Y_m(x)` for `m = |n|` and
applies one parity sign. There is no argument-parity reduction (a
negative argument is `INVALID`, not folded to `|x|`).

**NearestEven-only differential tier (no Ziv).** Like 6o, the `Y`
kernel uses the fixed working-precision guard, not the 7c `pow_ziv`
interval-test driver. The differential lane asserts bit-exact
(`assert_eq!`) under `NEAREST_EVEN_ROUNDING_MODES`.

**Tiered oracle.** rug 1.30 exposes MPFR's `mpfr_y0`/`mpfr_y1`/
`mpfr_yn` (`y0_ref`/`y1_ref`/`yn_ref`), so all three functions get
a true bit-exact MPFR differential lane under NearestEven over the
full `TRANSCENDENTAL_PRECISIONS` including `p = 1024` (the user
decision for 6p: mirror `differential_jn`, accept the deliberately
slow CI tier). The cross-order property is the **J/Y Wronskian**
`J_{n+1}(x)·Y_n(x) − J_n(x)·Y_{n+1}(x) = 2/(πx)` (DLMF 10.5.2),
binding the 6o `J` kernel to the new `Y` kernel; the property test
checks it in the π-free invariance form (the product is constant in
`n` and `x`), and the in-module unit test pins the actual
`2/(πx)` constant. This is the deliverable 6o explicitly could not
produce: it had only the `J`-recurrence property because `Y` did
not exist.

## Consequences

- `Y0`/`Y1`/`Yn` are correct for `x > 0` and integer order, pinned
  bit-exact against MPFR under NearestEven over a bounded sweep
  including `p = 1024`, with the Wronskian, parity, pole/domain,
  regime-continuity, and self-consistency properties green.
- The ordinary-Bessel pair `J`/`Y` is complete; the J/Y Wronskian
  is now a live cross-oracle for both. The `a_k(n)` machinery and
  the asymptotic threshold are shared, so 6q (`I`/`K`) inherits
  them.
- Cost: `Y` composes the 6o `J` Miller kernel plus `γ`/`ln`, and
  `n ≥ 2` adds the recurrence, so the differential lane is markedly
  heavier than `differential_jn` and `p = 1024` is the cost driver;
  mitigated by running it as a deliberately slow CI tier with the
  kernel correctness already pinned cheaply by the in-module
  `p = 1024` and Wronskian unit tests, and by the conservative
  asymptotic taking over for large `x`.
- The lack of a recessive-normalization analog means `Y`'s
  sub-asymptotic path is the log series with the `x·log₂e` guard
  across the whole range, costlier than `J`'s moderate Miller path.
  This is inherent to the dominant solution, not a deficiency to
  revisit.
- The NE-only tier remains a documented limitation: `pow` is still
  the only transcendental with five-mode correct rounding.
  Extending the `pow_ziv` driver to Bessel is a clean future slice.

## Related

- Commits: `fda562d` (skeleton) … `8f0bec2` (fuzz), the 6p.1–6p.9
  per-concern commits on `slice-6p-bessel-y`.
- Other ADRs: ADR-0023 (Bessel `J`; the regime-dispatch,
  smallest-term-truncation, decaying-envelope, tiered-oracle shape,
  and the `a_k`/`ω` derivation 6p reuses), ADR-0021 (Airy; the
  decaying-envelope convention), ADR-0019 (`Ci`/`li`; the
  pole-and-complex-domain convention), ADR-0022 (the Ziv driver
  this kernel deliberately does not yet use).
- Primary source: DLMF Chapter 10 (§10.4, §10.5, §10.6, §10.8,
  §10.17).
