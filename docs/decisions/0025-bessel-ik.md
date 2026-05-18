# ADR-0025: Modified Bessel functions I0, I1, In, K0, K1, Kn

- **Status**: accepted
- **Date**: 2026-05-17

## Context

Roadmap slice 6q closes the modified-Bessel pair, the companion to
slice 6o (Bessel `J`, ADR-0023) and slice 6p (Bessel `Y`,
ADR-0024). It adds `I0`/`I1`/`In` and `K0`/`K1`/`Kn` for integer
order and real argument under the existing `bessel = ["specials"]`
cluster feature (no new feature; the 6p precedent, recorded here as
the resolved fork). ADR-0023/0024 are the structural templates: `I`
reuses the `bessel_j` Miller recessive-solution template, `K` reuses
the `bessel_y` dominant-solution upward-recurrence-plus-log-series
template. The shared `a_k(ν)` asymptotic coefficients (DLMF 10.17.1,
Pochhammer-cross-pinned in ADR-0023) are reused a third time.

The biggest structural delta from 6o/6p: **rug 1.30 / MPFR expose no
modified-Bessel `I`/`K` primitive** (only `j*`/`y*`). There is no
bit-exact MPFR oracle, so 6q uses the Airy-style tiered oracle (the
6n precedent) rather than 6o/6p's bit-exact `differential_*` lane.

All formulas were derived from the DLMF primary source
(user-authorized `WebFetch` of DLMF §10.25, §10.27, §10.28, §10.29,
§10.30, §10.31, §10.35, §10.40 at plan time), not recalled. The 6n
Airy `(2k−1)`-divisor defect is the precedent that makes this
discipline load-bearing; this slice caught a second instance of it
(the `K` recurrence sign, below).

## Decision

**`I` (recessive in order) reuses the `J` template; `K` (dominant
in order) reuses the `Y` template.** DLMF 10.30.1
`Iν ∼ (½z)ν/Γ(ν+1) → 0` and DLMF 10.30.2 `Kν ∼ ½Γ(ν)(2/z)ν → ∞` as
`ν → ∞` at fixed `z` classify `I` as the recessive and `K` as the
dominant solution in order, exactly the roles `J` and `Y` play for
the ordinary recurrence.

- *`I`, three-regime dispatch on the binary exponent of `|x|`:* the
  DLMF 10.25.2 convergent Maclaurin series for tiny `|x|`
  (all-positive, monotone — no cancellation, so no `x·log₂e` boost,
  unlike `J`'s alternating 10.2.2); Miller backward recurrence for
  moderate `|x|`; the DLMF 10.40.1 asymptotic for large `|x|`.
- *`K`, base pair plus upward recurrence:* the DLMF 10.31.1
  logarithmic series below the asymptotic cut and the DLMF 10.40.2
  asymptotic above it for `K0`/`K1`; `Kn (n ≥ 2)` climbs by upward
  recurrence. `K` has no recessive-normalisation analog, so (like
  `Y`) there is no cheap middle regime.

**The three-term recurrence sign — a derive-don't-recall catch.**
DLMF 10.29.1 is the unified `𝒵_{ν−1} − 𝒵_{ν+1} = (2ν/z)·𝒵_ν` where,
by the §10.25(ii) standard-solution convention,
`𝒵_ν = I_ν(z)` **or `e^{νπi} K_ν(z)`**. Applied to `I` directly it
gives `I_{ν−1} = (2ν/z)I_ν + I_{ν+1}` — the **`+ I_{ν+1}`** backward
form (Miller, normalised by the DLMF 10.35.5 sum rule
`eˣ = I_0 + 2·Σ_{k≥1} I_k`, *every* order, the modified generating
function DLMF 10.35.1 at `t = 1`; not `J`'s even-only DLMF 10.12.4).
The `e^{νπi}` factor on the `K` form injects a `(−1)` that **flips
`K`'s sign**: `e^{−πi}K_{ν−1} − e^{πi}K_{ν+1} = (2ν/z)K_ν` collapses
to `K_{ν+1} = (2ν/z)K_ν + K_{ν−1}` — the **`+ K_{ν−1}`** upward
form, the *opposite* sign to a naive reading of the unified relation
and to `I`. The plan and the bead text had recalled the `K`
recurrence as `K_{k+1} = K_{k−1} − (2k/x)K_k`; a numeric spot-check
against textbook values (`K₂ − K₀ = (2/x)K₁` at `x = 1`) caught the
wrong sign before it shipped, and it is now pinned by
`recurrence_spot_check`. This is the 6n `(2k−1)` defect pattern
recurring at the recurrence-sign level; recorded as a named lesson,
not a silent fix.

**Small-argument forms.** `I_n` is the all-positive DLMF 10.25.2
Maclaurin. `K_n` is the DLMF 10.31.1 logarithmic series, the
modified analog of `Y`'s 10.8.1 with three sign differences derived
from §10.31: the finite head carries `+½` and an **alternating**
`(−x²/4)^k` (vs `Y`'s `−1/π`, plain `(x²/4)^k`); the log term is
`(−1)^{n+1}` with **no `2/π`** (vs `Y`'s `+2/π`); the tail carries
`(−1)^n ½` and a **positive** `(x²/4)^k` (vs `Y`'s `−1/π`,
alternating). The digamma terms reduce to harmonic sums exactly as
for `Y` (`ψ(k+1)+ψ(n+k+1) = −2γ + H_k + H_{n+k}`, no digamma
kernel, the in-tree `γ`); the `I_n(x)` piece is the new `bessel_i`
kernel. The reduction is cross-pinned by the worked DLMF
10.31.1 ↔ 10.31.2 identity (`n = 0` collapses to
`K_0 = −(ln(x/2)+γ)I_0 + Σ_{k≥1} H_k(x²/4)^k/(k!)²`), runnable as
`k0_dlmf_10_31_2_crosscheck`.

**Asymptotics reuse 6o's coefficients; no trig.** DLMF 10.40.1
`I_n ∼ eˣ/√(2πx)·Σ (−1)^k a_k(n)/x^k` (alternating); DLMF 10.40.2
`K_n ∼ √(π/2x)·e^{−x}·Σ a_k(n)/x^k` (**all positive**, no
`(−1)^k`). `a_k(n)` is the same sequence as ordinary Bessel
(DLMF 10.40 ≡ §10.17(i)), reused from the ADR-0023 pin rather than
re-derived. Neither has a trig factor — markedly simpler than
`J`/`Y`'s `sinω·ΣP + cosω·ΣQ`. DLMF 10.40.5's `e^{−x}` companion to
`I` is `O(e^{−2x})` relative on the positive real axis, below the
optimal-truncation error, so a single series suffices. The
asymptotic cut **reuses `bessel_j_threshold`**, and the reuse is
*derived*, not reflexive (the CLAUDE.md "derive the cut" reflex):
the accuracy-controlling quantity is the optimal-truncation
*relative* error of the shared `a_k(ν)` divergent series,
`O(e^{−2x})`, identical to the ordinary-Bessel 10.17 series; the
`e^{±x}` prefactor is computed exactly and does not enter the
relative error, so the conservative `2^{e_x} ≥ target+64` cut is
strictly more than enough.

**Parity and domain — both even in order, no sign.** DLMF 10.27.1
`I_{−n} = I_n` and DLMF 10.27.3 `K_{−n} = K_n`: **even in order
with no sign**, unlike `J`/`Y`'s `(−1)^n`. The table, Kani-pinned:

| input | `I` result | `K` result |
|-------|------------|------------|
| `qNaN` | `qNaN`, `OK` | `qNaN`, `OK` |
| `sNaN` | `qNaN`, `INVALID` | `qNaN`, `INVALID` |
| `±0` | `I_0 = 1`, `I_n = 0` (n≠0), `OK` | `+∞`, `DIV_BY_ZERO` (`+0`); `qNaN`, `INVALID` (`−0`) |
| `x < 0` | finite, `(−1)^n I_n(|x|)`, `OK` (entire) | `qNaN`, `INVALID` (complex) |
| `+∞` | `+∞`, `OK` | `+0`, `OK` |
| `−∞` | `(−1)^n·∞`, `OK` | `qNaN`, `INVALID` (complex) |

`I` is entire (the `bessel_j` domain shape): `I_n(−x) = (−1)^n
I_n(x)` argument parity, `I_n(±∞) = ±∞` a **genuine infinite
limit** (`Status::OK`, the `exp(+∞) = +∞` precedent — explicitly
**not** the decaying-envelope convention, which covers a bounded
non-converging oscillation). `K` is `x > 0` only (the
`Y`/`Ci`/`li` shape): `K_n(+0) = +∞` raising `DIV_BY_ZERO` is a
pole — **positive**, the opposite of `Y_n(+0) = −∞`. `K_n(+∞) = +0`
is a **genuine exponential-decay limit** (DLMF 10.40.2
`√(π/2x)·e^{−x} → 0`), `Status::OK` — explicitly distinguished from
the decaying-envelope *convention* of `J`/`Y`/Airy: there the
function oscillates with a shrinking but non-converging envelope and
`+0` is a conservative total-keeping choice; here `K` actually
converges to `0`, so `+0` is the true mathematical limit.

**Tiered oracle (the biggest delta from 6o/6p).** With no MPFR
`I`/`K` primitive, `tests/differential_ik.rs` is the Airy-style
tiered oracle ([`differential_bi`] precedent,
`feedback_differential_lane_cost`): (1) a checked-in `mpmath`
`besseli`/`besselk` reference table at 86 digits, capped `p ≤ 256`,
`close_within(p−2)`; (2) the **DLMF 10.28.2 I/K cross-tie**
`I_ν K_{ν+1} + I_{ν+1} K_ν = 1/x` (a **plus** and `1/x`, the
modified-Bessel analog of 6p's J/Y Wronskian — π-free, so the
property form needs no constant); (3) precision self-consistency on
**dyadic** arguments (the pf-ok9 lesson, applied from the start).
`p = 1024` is pinned cheaply by the in-module `high_precision_pin`
and `ik_wronskian_10_28_2` unit tests. NearestEven-only, no Ziv
(the 6o/6p posture).

## Consequences

- The modified-Bessel pair `I`/`K` is complete for integer order
  and real argument; with `J`/`Y` (6o/6p) the entire ordinary-and-
  modified Bessel surface for v1.0 is shipped. The DLMF 10.28.2 I/K
  cross-tie is a live cross-oracle binding both new kernels.
- The `a_k(ν)` machinery and `bessel_j_threshold` are now shared
  four ways (`J`, `Y`, `I`, `K`); the derived (not reflexive)
  threshold-reuse argument is recorded so the next reuse does not
  re-litigate it.
- Cost: `K` composes `bessel_i` plus `γ`/`ln` and (for `n ≥ 2`) an
  upward recurrence, so the tiered lane is heavier than
  `differential_jn`; mitigated by the `p ≤ 256` cap with the
  `p = 1024` path pinned by in-module unit tests, and by the
  conservative asymptotic taking over for large `x`. The lane is a
  deliberately slow CI tier.
- The recurrence-sign catch is the second instance of the
  recalled-coefficient failure mode (after 6n Airy); the lesson is
  promoted to the derive-don't-recall memory so the reflex covers
  recurrence *signs*, not only coefficient *recurrences*.
- The NE-only tier remains a documented limitation: `pow` is still
  the only transcendental with five-mode correct rounding.
  Extending the `pow_ziv` driver to Bessel is a clean future slice.

## Related

- Commits: the 6q.1–6q.11 per-concern commits on
  `slice-6q-bessel-ik`.
- Other ADRs: ADR-0023 (Bessel `J`; the Miller/recessive template
  `I` reuses, the `a_k`/threshold derivation reused a third time),
  ADR-0024 (Bessel `Y`; the dominant/upward-plus-log-series template
  `K` reuses, the J/Y Wronskian whose I/K analog is DLMF 10.28.2),
  ADR-0021 (Airy; the decaying-envelope convention `K_n(+∞)`
  explicitly is *not*, and the no-MPFR-primitive tiered-oracle
  precedent), ADR-0019 (`Ci`/`li`; the pole-and-complex-domain
  convention `K` shares), ADR-0022 (the Ziv driver this kernel
  deliberately does not yet use).
- Primary source: DLMF Chapter 10 (§10.25, §10.27, §10.28, §10.29,
  §10.30, §10.31, §10.35, §10.40).
