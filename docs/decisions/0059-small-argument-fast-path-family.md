# ADR-0059: small-argument fast-path family for the odd elementary kernels

- **Status**: accepted
- **Date**: 2026-06-03
- **Disposition**: accepted for `atanh`, `asinh`, `sinh`, `tanh`;
  measured and rejected for `atan`, `sin`, `tan`.

## Context

The pf-lm3 exhaustive `f32` verification sweep found `atanh` on
tiny/subnormal inputs (`|x| < 2^-95`) costing ~0.6-1.5 ms per input,
20-40x the normal-magnitude cost, and a contiguous band of these
dominated the sweep wall-clock and spend. The `atanh` kernel drives the
cancellation-resistant identity `atanh(x) = (log1p(x) − log1p(−x))/2`
through the Ziv loop for every input. For tiny `x` the true value is
`x` plus a cubic tail `x³/3 + …` that sits far below a ULP, so the
entire composition computes a value that is `x` to within rounding, at
high working precision, for nothing.

All seven odd elementary kernels share this shape: `f(x) = x + c·x³ + …`
for small `x` (the even-power terms vanish). The same waste is latent in
each. The fix already has a proven precedent in the codebase: `log1p`
(`src/math/log1p.rs`) and `expm1` short-circuit their own tiny-x regime
by rounding `x` perturbed by a known-sign infinitesimal through
`round_with_infinitesimal` (`src/rounding.rs`), a single mode-aware
round that is already certified across the `f32` exhaustive sweep.

The trap to avoid is the one ADR-0050 documents: a bare "return `|x|`
rounded under mode" short-circuit is correct under round-to-nearest but
**wrong under the directed modes**, because it drops the sign of the
correction term. `tanh` carried exactly such a short-circuit; ADR-0050
removed it (returning the grid point `|x|` also defeated the Ziv
interval test's convergence). `round_with_infinitesimal` is the correct
mechanism: it carries the dropped direction, so the directed modes round
to the right neighbour.

## Decision

Add a tiny-x short-circuit to the odd-kernel family, gated on the same
threshold `log1p`/`expm1` use, `x.exponent ≤ −(target_precision + 2)`,
returning `round_with_infinitesimal(x, x.sign(), subtracts_magnitude,
target, mode)`. The `subtracts_magnitude` flag is a per-function
compile-time constant fixed by the sign of the cubic correction relative
to `x`:

| fn | series near 0 | magnitude | `subtracts_magnitude` |
|---|---|---|---|
| `atanh` | `x + x³/3 + …` | grows | `false` |
| `asinh` | `x − x³/6 + …` | shrinks | `true` |
| `sinh` | `x + x³/6 + …` | grows | `false` |
| `tanh` | `x − x³/3 + …` | shrinks | `true` |

Because these kernels are odd, the correction term **follows** `x`'s
sign, so the magnitude direction is unconditional (unlike `log1p`, whose
`−x²/2` term is sign-independent and so needs a sign-dependent flag).

**Threshold soundness.** The cubic tail is `O(x³)`, smaller than
`log1p`'s `O(x²)`, so the shared threshold is conservative. Round-to-
nearest needs `|tail| < half-ULP(x) = 2^(e−target)`; with
`|tail| ≈ 2^(3e)/3` and `e ≤ −(target+2)` the headroom is `~2^(target+4)`.
The directed modes need only a nonzero, known-sign tail, which holds.
Sound for all five IEEE modes; this is the same argument the `f32` sweep
already certified for `log1p`/`expm1`.

**Strict revert (ADR-0041 precedent).** A new `benches/small_arg.rs`
measured all seven against a normal-magnitude control before any kernel
edit. The four functions above grind in the moderate-tiny band because
they compose through another Ziv kernel (`atanh`/`asinh` through
`log1p`, `sinh`/`tanh` through `expm1`); they land the fast-path. The
remaining three (`atan`, `sin`, `tan`) evaluate a direct Taylor series
that terminates in ~one term for tiny `x`, already cost ~1-2 µs, and are
**rejected as measured-neutral**: a short-circuit would save a sub-µs
constant for the verification and ADR cost of a fourth path, which is
not frugal.

## Consequences

- **The atanh sweep hotspot closes.** In-band cost drops to a constant
  ~170-280 ns floor (the `round_with_infinitesimal` cost), from
  per-function baselines of 59-130 µs (`atanh`), 4-11 µs (`asinh`),
  10-19 µs (`sinh`), and 5-10 µs (`tanh`): a 95-99.8% reduction, deepest
  on `atanh`. The bench confirms the floor is precision-independent (the
  helper is `O(1)`).
- **No correctness change.** Every function stays correctly rounded in
  all five modes; the returned values are unchanged, only the path that
  produces them. The `tests/oracle/status/*.toml` rows and
  `docs/rounding-status.md` are byte-identical, so those gates pass
  unchanged.
- **Out-of-band inputs are untouched.** The threshold is conservative,
  so inputs just above the edge stay on the Ziv path. The added guard is
  a single untaken branch (extract the exponent, one integer compare)
  ahead of the unchanged composition, so the out-of-band code path and
  its result are identical by construction; this rests on the code
  structure, not on the bench (whose out-of-band cells drift with
  thermal state across a long run).
- **tanh distinguished from ADR-0050.** The removed short-circuit
  returned the grid point `|x|` (directed-mode-wrong, Ziv-divergent).
  The restored one rounds through `round_with_infinitesimal`, which is
  mode-aware and never returns a bare grid point, so neither ADR-0050
  failure mode recurs. The `expm1` stable form from ADR-0050 remains the
  out-of-band path.
- **A dedicated verification lane exists** because the integer-input
  differential lanes cannot reach the tiny-x band; see Verification.

## Verification

- `benches/small_arg.rs`: per-function tiny-vs-control ratios, the
  measure-then-decide evidence (the rejected three show no win).
- `tests/differential_small_arg.rs`: a new MPFR differential lane
  sweeping exact dyadic `±2^exp` inputs straddling each function's
  activation edge across precisions {24, 53, 113, 200} and all five
  modes. The integer-input lanes only sweep `|x| ≥ 1` and never reach
  the band; this lane adds the `bigfloat_pow2`/`rug_pow2` generators.
- Per-function hermetic unit tests: `atanh_tiny_input_directed_modes`
  pins all five modes against the exact away-neighbour (the grow
  signature of `subtracts_magnitude = false`); the existing
  `tanh_tiny_input_round_to_nearest_returns_input` (ADR-0050) still
  passes under the restored short-circuit.

## Related

- ADR-0050 — removed `tanh`'s bare-`|x|` short-circuit; this ADR
  restores a directed-mode-correct one and contrasts the mechanisms.
- ADR-0041 — strict-revert "measured, rejected" precedent for the three
  neutral siblings.
- ADR-0049 / ADR-0038 — the Ziv driver and cross-check the kernels run
  under.
- `src/math/log1p.rs`, `src/math/expm1.rs`,
  `src/rounding.rs::round_with_infinitesimal` — the precedent mechanism.
- pf-lm3 — the exhaustive `f32` sweep finding that motivated the work.
