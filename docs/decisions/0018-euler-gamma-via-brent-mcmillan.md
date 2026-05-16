# ADR-0018: Euler–Mascheroni constant via Brent–McMillan

Status: accepted (slice 6m0)

## Context

Slice 6m adds the integral special functions `Ei`, `Si`, `Ci`, and
`li`. Three of the four (`Ei`, `Ci`, and `li` transitively through
`Ei`) need the Euler–Mascheroni constant γ in their small-argument
series. pfloat had no γ: a search of the tree found only a short
literal inside a `digamma` test.

Slices 7b, 7b1, and 7b2 established and then defended the principle
that pfloat's mathematical constants are correct at any requested
precision. 7b2 specifically removed the last precision cap. A γ that
only reached a fixed width would reintroduce exactly the regression
those slices closed, so γ has to meet the same bar: correct at any
precision, with a fast hardcoded path for the common range.

The constants already in `agm_constants` (π via Brent–Salamin, ln 2
via the atanh series) each have an arithmetic–geometric-mean or
series form with a clean error bound. γ has no AGM form. Its
canonical high-precision algorithm is Brent and McMillan's, which
expresses γ through Bessel-type sums with an exponentially small
truncation error.

## Decision

1. **Add `euler_gamma_via_bm(prec)` to `src/math/agm_constants.rs`,
   derived from Brent–McMillan algorithm B1.** With
   `t_k = (nᵏ / k!)²` and `H_k` the kth harmonic number (`H_0 = 0`),

   ```text
   I(n) = Σ_{k≥0} t_k
   S(n) = Σ_{k≥0} t_k · H_k
   γ    = S(n)/I(n) − ln(n) + O(π·e^{−4n})
   ```

   The implementation derives from the published identity, not from
   any reference implementation's source. `t_k` is carried by the
   recurrence `t_k = t_{k−1}·n²/k²` and `H_k` by `H_{k−1} + 1/k`, so
   no factorial or power is materialized. `n` is chosen so the
   `e^{−4n}` truncation error sits below the working precision; the
   integer ratio `7/40 > (ln 2)/4` keeps the choice no_std-clean and
   conservative. The inner sums rise to a peak near `k = n` and then
   decay super-geometrically; the loop runs until a term is
   negligible relative to the running `I(n)` and past the peak, with
   a hard iteration cap so a pathological input cannot spin.

2. **Memoize γ like the other constants.** `Kind::EulerGamma` joins
   the slice-7b1 `(kind, precision)` thread-local cache. The `std`
   and `no_std` arms are already generic over `Kind`; no cache
   surgery is needed.

3. **Hardcoded table plus dispatch in `src/math/mod.rs`, mirroring
   `ln_2_at`.** `EULER_GAMMA_LIMBS_1024` holds the correctly-rounded
   1024-bit value; `euler_gamma_at(prec)` returns the rounded table
   for `prec ≤ 1024` and `euler_gamma_via_bm(prec)` above it. There
   is no precision cap: the table is a fast path, not a ceiling.

4. **Pin the table to an authoritative decimal and the independent
   computation.** The limbs are the mantissa of the bit-exact parse
   of an authoritative γ decimal (OEIS A001620, treated as a
   mathematical fact, the `LN2_REFERENCE` pattern). A regression
   test asserts the table equals the parsed reference and equals
   `euler_gamma_via_bm` computed with generous headroom and rounded
   to 1024 bits, all three bit-for-bit. Because Brent–McMillan is an
   independent code path anchored at low precision against the
   universally-established short value of γ, a transcription error
   in the reference cannot pass this test.

## Consequences

- `Ei`, `Ci`, and `li` get a γ that is correct at any precision,
  preserving the no-cap property restored in 7b2.
- γ shares the memoization, so a high-precision integral sweep pays
  the Brent–McMillan cost once per distinct working precision rather
  than once per call, the same amortization 7b1 gave π and ln 2.
- Unlike ln 2 (cross-checked by two independent in-repo algorithms,
  atanh and AGM), γ has only one practical high-precision algorithm.
  The independent check is therefore the published-algorithm
  derivation plus the authoritative decimal, not a second in-repo
  method. This is the strongest check γ admits and is stronger than
  the self-consistency fallback the roadmap reserved for constants
  without an oracle.
- Algorithm B1 is used rather than the refined B3 (which halves `n`
  for the same accuracy). B1 is the simpler, more obviously correct
  form; the n-selection over-computes by a small constant factor.
  Tightening to B3 is a candidate performance slice, gated on a
  measured win per the performance-patch discipline.

## References

- Brent, R. P. and McMillan, E. M. *Some new algorithms for
  high-precision computation of Euler's constant.* Mathematics of
  Computation 34 (1980), pp. 305–312. The B1 identity and its
  `e^{−4n}` error bound.
- DLMF §5.2 — definition and basic properties of γ.
- ADR-0017 — the AGM constant kit and the table-versus-on-the-fly
  dispatch pattern this constant follows.
- OEIS A001620 — the authoritative decimal digit reference for γ.
- Plan: `let-s-review-the-backlog-vast-harbor.md` — slice 6m;
  γ scoped as the preceding sub-slice 6m0.
