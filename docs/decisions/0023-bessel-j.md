# ADR-0023: Bessel functions of the first kind J0, J1, Jn

- **Status**: accepted
- **Date**: 2026-05-17

## Context

Roadmap slice 6o opens the Bessel family, the largest remaining
tier-2 surface on the road to v1.0 (6p adds `Y`, 6q adds `I`/`K`).
It adds `J0`, `J1`, and `Jn` for integer order and real argument
behind a new `bessel = ["specials"]` cluster feature. It is the
first kernel to carry a genuine recurrence rather than a per-point
series or asymptotic alone, so it builds the recurrence and
normalization machinery the later Bessel slices reuse. Slice 7c
(pow Ziv) shipped first specifically so an integer-power capability
was available for recurrence work. The just-landed slice 6n (Airy,
ADR-0021) is the structural template: the same public-API shape,
binary-exponent regime dispatch, smallest-term asymptotic
truncation, decaying-envelope limit convention, and tiered oracle.

`J0`/`J1`/`Jn` are all entire on the real line (no poles, no domain
restriction), which keeps the special-value handling simple. `Y` is
out of scope here and deferred to 6p, so the J/Y Wronskian cannot
yet be the cross-tie; the three-term recurrence plays that role
instead.

## Decision

**Three-regime dispatch on the binary exponent of `|x|`** (the
`airy`/`si` integer-exponent selector idiom):

- *Tiny `|x| < 1`*: the convergent Maclaurin series, DLMF 10.2.2,
  `J_m(x) = (x/2)^m Σ_{k≥0} (−1)^k (x/2)^{2k}/(k!·(m+k)!)` (integer
  order ⇒ `Γ(m+k+1) = (m+k)!`). This is the small-argument backstop
  and the reason the `2k/x` recurrence never sees `x → 0`.
- *Moderate `|x|`*: Miller backward recurrence.
- *Large `|x|`*: the Hankel asymptotic, DLMF 10.17.3.

**Miller backward recurrence with sum-rule normalization, chosen
over direct per-order series** (this resolves roadmap open-decision
#4, "Bessel `Jn` recurrence direction"). The DLMF 10.6.1 three-term
relation `𝒞_{ν−1}+𝒞_{ν+1} = (2ν/z)𝒞_ν` rearranged downward is
`f_{k−1} = (2k/x)·f_k − f_{k+1}`. Started at a high seed index `M`
with `f_{M+1}=0`, `f_M=1`, the descent converges to a fixed multiple
`c·J_k(x)` of the recessive solution. The DLMF 10.12.4 identity
`1 = J_0(x) + 2J_2(x) + 2J_4(x) + ⋯` (re-derived independently by
setting `t = 1` in the 10.12.1 generating function and using
`J_{−n}=(−1)^n J_n`) gives `S = f_0 + 2(f_2+f_4+⋯) = c·1 = c`, so
`J_m(x) = f_m / S`, every order from one descent and one normalizing
constant. Forward recurrence is simpler but loses accuracy for the
recessive `J`; direct per-order series would not share the machinery
6p/6q need. The descent is carried with a rolling triple
(`f_{k+1}`, `f_k`, `f_{k−1}`) accumulating `S` and capturing `f_m`
on the fly, so it is `O(M)` time and `O(1)` memory, no `Vec`.

**Seed index `M` derived, not guessed.** DLMF 10.19.1 gives
`J_M(x) ∼ (1/√(2πM))·(eX/(2M))^M` as `M→∞`; the prefactor only
shrinks the bound, so requiring `(eX/(2M))^M < 2^{−P}`,
`P = target+64`, is conservative. In natural logs that is
`M·(1 + ln(x/(2M))) < −P·ln2`, solved by an exponential search at
low working precision plus a small fixed step guard. Overshoot
`≤ 2×` the optimal `M` is a deliberate robustness/cost trade; the
crossover is not perf-tuned without a bench.

**Asymptotic coefficient provenance discipline (the 6n lesson).**
The DLMF 10.17.1 coefficients are derived from the spec, never
recalled: `a_0(m)=1`, `a_k(m) = a_{k−1}(m)·(4m²−(2k−1)²)/(8k)`. The
`8k` divisor was cross-checked two independent ways against the
primary source (user-authorized `WebFetch` of DLMF §10.17): the
closed-form ratio `[(4m²−1²)…(4m²−(2k−1)²)]/(k!·8^k)` and the
Pochhammer form `(½−m)_k(½+m)_k/((−2)^k k!)`, which agree at
`k=1,2` (`a_1=(4m²−1)/8`, `a_2=(4m²−1)(4m²−9)/128`). DLMF 10.17.3
is `J_m(x) ∼ √(2/(πx))·[cosω·Σ(−1)^k a_{2k}/x^{2k} −
sinω·Σ(−1)^k a_{2k+1}/x^{2k+1}]`, `ω = x − mπ/2 − π/4`; folding the
explicit `(−1)^k` into the trig assignment yields the period-4
factor `[+cosω, −sinω, −cosω, +sinω]` on `a_j(m)/x^j`, summed to the
smallest term. The 6n Airy `(2k−1)`-divisor defect, invisible at
`k=1`, is the precedent: a coefficient recurrence must be derived
and cross-pinned, not recalled.

**Conservative asymptotic threshold.** The optimally-truncated
Bessel-J asymptotic has error of order `e^{−2|x|}`, so reaching
`target+64` bits needs `|x| ≳ 0.347·(target+64)`. Requiring
`2^{e_x} ≥ target+64` is strictly more than enough; Miller (always
correct, if slower) carries everything below. Deliberately
conservative, not tuned.

**Cancellation guard.** The Miller descent and the sum-rule share
cancellation, so working precision is boosted by `≈ |x|·log₂e` bits
(the `si.rs` rational `23/16`, capped form) plus the `+64` base. The
kernel rounds once at the end via `round_to_precision` and raises
flags via `auto_raise`.

**NearestEven-only differential tier (no Ziv).** The J kernel uses
the fixed working-precision guard like `ei`/`si`/`erf`, not the 7c
`pow_ziv` interval-test driver. The differential lane asserts
bit-exact (`assert_eq!`) under `NEAREST_EVEN_ROUNDING_MODES`. A
later slice could extend the `pow_ziv` driver to Bessel for full
five-mode correct rounding; that is out of scope here (no perf or
correctness machinery beyond the slice without a measurement).

**Limit and parity conventions.** `J_n` is entire, so there are no
poles. The table, Kani-pinned:

| input | result | status |
|-------|--------|--------|
| `qNaN` | `qNaN` (payload propagated) | `OK` |
| `sNaN` | `qNaN` | `INVALID` |
| `J_0(±0)` | `1` (exact, DLMF 10.2.2) | `OK` |
| `J_n(±0)`, `n ≠ 0` | `+0` (exact) | `OK` |
| `J_n(±∞)`, any `n` | `+0` | `OK` |

`J_n(±∞) = +0` is the decaying-envelope convention (ADR-0021, the
Airy precedent): the true behaviour at `±∞` is a bounded decaying
oscillation with no limit, so the conservative total result keeps
the function total. Negative argument and negative order reduce
before regime dispatch: `J_n(−x) = (−1)^n J_n(x)` (DLMF 10.11.1) and
`J_{−n}(x) = (−1)^n J_n(x)` (DLMF 10.4.1), so the kernel evaluates
`J_m(|x|)` for `m = |n|` and applies one parity sign, negating
exactly when `m` is odd and exactly one of `{n<0, x<0}` holds.

**Tiered oracle.** rug 1.30 exposes MPFR's `mpfr_j0`/`mpfr_j1`/
`mpfr_jn` for all real arguments, so all three functions get a true
**bit-exact** MPFR differential lane under NearestEven (the
`differential_ei` idiom), unlike Airy where only `Ai` had an MPFR
oracle. The cross-order property is the recurrence
`J_{n−1}+J_{n+1} = (2n/x)·J_n`, binding three independently
descended orders (the J/Y Wronskian is re-enabled in 6p). The
differential sweep is a small bounded representative set, not a wide
random sweep: the integer arguments land in Miller, whose seed index
grows with precision, so `p = 1024` is the cost driver and the lane
is a deliberately slow CI tier (the `differential_si`/`_ci`
posture).

## Consequences

- `J0`/`J1`/`Jn` are correct across the whole real line for integer
  order, pinned bit-exact against MPFR under NearestEven over a
  bounded sweep including `p = 1024`, with the recurrence, parity,
  boundary, and self-consistency properties green.
- The recurrence and sum-rule normalization machinery is in place
  for slices 6p (`Y`) and 6q (`I`/`K`), and roadmap open-decision #4
  is resolved (Miller backward, normalized).
- Cost: the Miller seed index `M` grows with target precision, so
  `p = 1024` is the cost driver; mitigated by the conservative
  asymptotic taking over for genuinely large `|x|`, the tiny regime
  capping the `2k/x` blow-up near zero, and a bounded differential
  sweep validated as a fast subset locally with the full lane in CI.
- The NE-only tier is a documented limitation: `pow` is still the
  only transcendental with five-mode correct rounding. Extending the
  `pow_ziv` driver to Bessel is a clean future slice.

## Related

- Commits: `13956bc` (skeleton) … `5cf7ca7` (fuzz), the 6o.1–6o.8
  per-concern commits on `slice-6o-bessel-j`.
- Other ADRs: related to ADR-0021 (Airy; the same regime-dispatch,
  smallest-term-truncation, decaying-envelope, and tiered-oracle
  shape) and ADR-0022 (the Ziv driver this kernel deliberately does
  not yet use).
- Primary source: DLMF Chapter 10 (§10.2, §10.4, §10.6, §10.11,
  §10.12, §10.17, §10.19).
