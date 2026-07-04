# ADR-0115: complex status and sign fidelity — qNaN/sNaN INVALID, exact-quotient OK, zero-dividend sign, csqrt directed imaginary

- **Status**: accepted
- **Date**: 2026-07-03

## Context

The epic pf-8iji review of pfloat-complex 1.0.0 found four behaviour
defects where the returned value or its IEEE `Status` disagreed with C11
Annex G and IEEE 754. The complex 1.0.0 public API is frozen (ADR-0093);
every fix here is behaviour-only (a 1.0.1 patch), touching no signature,
no public type, and neither the merged `Status` nor the single
`RoundingMode` surface.

The four defects share a theme: a signed zero, a signaling flag, or an
exactness flag was decided from the wrong source (the rounding mode, a
boxed operand, or a non-collapsing bracket) instead of from the operands
and the spec.

- **pf-pz9r** (`div.rs`). `(−0 − 0i)/(3 + 0i)` under NearestEven returned
  `re = +0`. Componentwise IEEE division gives `−0`: the real numerator
  `ac + bd = (−0)·3 + (−0)·0 = −0` divided by the non-negative denominator
  `c² + d²` is `−0` (the sign flows from the inputs, ADR-0091's
  cross-function signed-zero rule). The exact-zero short-circuit stamped
  the sign from the rounding mode, which is right only for the z/z
  cancellation it was written for.
- **pf-yprp** (`csqrt.rs`). `csqrt(−2 − 0i)` under TowardPositive returned
  `im = −1.41421356237309504887…`, more negative than `−√2`, on the wrong
  side for a mode that must not overshoot toward `−∞`. The negative-axis
  imaginary part is `copysign(√|x|, y)`; the magnitude was rounded under
  `mode` and then negated, which lands a directed result on the wrong side.
- **pf-hdq1** (`specials.rs`). `(qNaN + 1i)/(2 + 0i)` raised INVALID, but a
  quiet NaN must propagate silently (IEEE 754 §7.2 reserves INVALID for a
  signaling NaN or a genuine invalid operation). Conversely
  `mul((∞ + sNaN·i), (2 + 3i))` returned `Status::OK`: the §G.5.1 recovery
  boxes the signaling NaN to a signed zero and discarded the INVALID it
  must raise.
- **pf-bv2i** (`div.rs`, the hardest, where exactness meets the enclosure
  contract).
  - `z/z` for `z = (1 + 2^-199) + i` at `p = 200` returned `(1, +0)` with
    INEXACT, contradicting ADR-0090 (an exact quotient is OK). The real
    numerator `N = ac + bd` equals the denominator `D = c² + d²`
    structurally, so `N/D = 1` exactly; but at the working precision `N`
    and `D` are each separately inexact, the directed quotient bracket
    straddles `1` without collapsing, and the collapse-only exactness test
    under-reported.
  - `(1 + 1i)/(c + ci)` with `c = 2^-(2^62 + 100)` returned a wrong,
    saturating value with `Status::OK`. Forming `c² + d²` drives the
    exponent below `i64::MIN`; pfloat has no emin, so the exponent
    saturates and the directed pair no longer brackets the true (tiny)
    denominator. The enclosure is unsound there, yet the OVERFLOW/UNDERFLOW
    the formation raised was dropped by the `.0` on the bracket helpers.

## Decision

### pf-pz9r: the zero-dividend sign follows the operands, not the mode

`div.rs::resolve`'s exact-zero short-circuit now distinguishes two sources
for the sign of a zero numerator over a positive denominator, by inspecting
the directed numerator pair:

- The pair AGREES in sign (`[−0, −0]` or `[+0, +0]`): the numerator is an
  exactly-representable signed zero fixed by the inputs. Because
  `D = c² + d² ≥ 0`, IEEE `(±0)/(+D) = ±0` takes the numerator's sign.
- The pair STRADDLES (`[−0, +0]`): the numerator is a genuine cancelling
  difference (the imaginary part of `z/z`, `bc − ad = 0`), whose zero sign
  is mode-determined (`−0` under TowardNegative, `+0` otherwise). This is
  the pre-existing behaviour, now scoped to the case it is correct for.

The agree/straddle test separates the two cleanly: a numerator that is
exactly `±0` rounds to that same signed zero in both directed modes, while
a cancellation `x − x` rounds to `−0` toward `−∞` and `+0` toward `+∞`.

### pf-yprp: negate after the mirrored directed round

The negative-real-axis imaginary root is `copysign(√|x|, y)`. When `y` is
negative the value is `−√|x|`; rounding a negated value under `mode` means
rounding the magnitude under the *mirrored* mode and then negating, since
`round(−v, mode) = −round(v, mirror(mode))`. `csqrt` now rounds `√|x|`
under `mirror_mode_for_negation(mode)` before the `copysign`. The mirror
swaps TowardNegative and TowardPositive and passes the nearest and
toward-zero modes through (they are symmetric under negation). This is the
crate-local twin of pfloat's scalar `signed_constant_at_round` idiom
(ADR-0101). Only the exact real-axis dispatch needed it; the interior
Kahan enclosure already rounds both bracket ends under `mode` and is
correct by convergence.

### pf-hdq1: a signaling-NaN INVALID floor over the §G.5.1 recovery

Both `complex_div_special` and `recover_mul` now wrap a core that pins the
Annex G *values* and OR a `signaling_invalid` floor onto every recovered
row: INVALID iff any of the four operand parts is a signaling NaN, else OK.
This reinstates the flag the boxing (`sNaN → signed 0`) would otherwise
drop, mirroring the `csqrt`/`cexp`/`clog` pattern (a signaling NaN raises
INVALID even where an infinity overrides the value). The divide's
NaN-operand fall-through now returns a base `Status::OK`, so a quiet NaN
propagates silently and only the floor adds INVALID for a signaling one.
The genuine `∞/∞` and `∞·0` rows keep their INVALID (an invalid operation,
not a NaN-propagation).

### pf-bv2i (a): a sound residual exactness certificate

When the quotient bracket converges but does not collapse, `resolve` falls
back to the residual certificate `N − r·D = 0`, with `N` the component
numerator (`ac + bd` or `bc − ad`), `D = c² + d²`, and `r` the rounded
result, all recomputed exactly from the operands. A cheap screen first
brackets the residual at the working precision from the pairs already in
hand and bails unless it can contain zero (the common inexact case, where
`|N − r·D| ≈ ulp(r)·D` sits far from zero). Only then is the exact residual
formed at a generous bounded precision; if any step cannot be represented
exactly (INEXACT or exponent saturation), the certificate conservatively
returns false. The check therefore never reports a false OK: a `true`
result means every step was exact and `N − r·D = 0`, hence `N/D = r`
exactly.

### pf-bv2i (b): saturation carries the status, not a clean flag

The `mul_add_mul` / `mul_sub_mul` bracket helpers now return their status,
and the loop OR the OVERFLOW/UNDERFLOW of the numerator and denominator
formation into the result (the INEXACT of a directed bracket rounding is
filtered out, as it is expected). When forming `c² + d²` saturates the i64
exponent, the enclosure contract's assumption ("a finite product of finite
operands never saturates the exponent before the round", ADR-0091's M-OVF
exemption) fails for a *squared tiny* operand on the underflow side, so the
result now honestly reports not-OK (UNDERFLOW here) rather than a wrong
value with a clean flag.

The value in that saturating case stays best-effort. Making it correct
would require a scaling reformulation of the whole divide (Smith-style
pre-scaling to keep `c² + d²` in range), a new algorithm outside a
status-fidelity patch and outside the frozen API's intent. The honest
signal (a raised flag) is the fix; the exact value on a divisor whose
square underflows `i64::MIN` is a documented, flagged limitation, the dual
of ADR-0091's overflow exemption.

## Consequences

- Four value/flag disagreements with Annex G and IEEE 754 are closed under
  a behaviour-only patch; the frozen 1.0 API is untouched.
- The residual exactness certificate costs, on the converged iteration of
  an inexact quotient, one cheap residual bracket at the working precision;
  the expensive exact residual (a few fused ops at `2p + 4096` bits) runs
  only when the cheap screen cannot exclude zero, which happens for exact
  or hard-to-round quotients, a measure-zero set. No unbounded cost: the
  exact residual is capped and bails to a spurious INEXACT on a wild
  exponent gap, the same caveat ADR-0090 already accepts for the value.
- The saturation flag is `OK` on every normal-magnitude division (no
  spurious OVERFLOW/UNDERFLOW), so ordinary results are unchanged; it fires
  only when the i64 exponent genuinely saturates while forming a bracket.
- Verification: a new `tests/regression_pf8iji.rs` pins all four defects
  and their sign/mode variants; the existing complex unit, identity,
  Annex G enumeration, dispatch-totality, and the independent `acb`
  componentwise certified-rounding differential (cmul/cdiv 840, csqrt 480,
  cexp/clog 480 each, 0 mismatches) all stay green, confirming the divide
  and csqrt values remain correctly rounded.

### Inversion (refuted alternatives)

- **Stamp the zero-dividend sign from the mode always** (the status quo).
  Refuted: it returns `+0` for `(−0 − 0i)/(3 + 0i)` under NearestEven,
  where componentwise IEEE gives `−0`. The mode is the right source only
  for a cancelling difference, told apart by the straddle test.
- **Round the csqrt magnitude under `mode` and negate.** Refuted: it lands
  `−√2` on the `−∞` side under TowardPositive, overshooting the value a
  directed mode must not pass. The mirror is not optional for a negated
  directed round.
- **Keep `nan_pair` (INVALID) for every NaN operand.** Refuted: it raises
  INVALID for a quiet NaN, which IEEE 754 requires to propagate silently.
  The floor must key on *signaling*, not on *NaN*.
- **Trust the naive recovery status for multiply.** Refuted: the boxing
  turns a signaling NaN into a signed zero, so the naive status is OK where
  IEEE mandates INVALID. The floor reads the original operands, before
  boxing.
- **Decide exactness only by bracket collapse** (the status quo).
  Refuted: for `z/z` the bracket straddles the exact `1` without collapsing
  (N and D separately inexact), so an exact quotient reports INEXACT. The
  residual certificate reads exactness from `N − r·D`, not from the
  bracket's width.
- **Grow the working precision until the quotient bracket collapses to
  decide exactness.** Refuted: a genuinely inexact quotient never collapses,
  so the loop would run to the cap on every inexact divide and still
  over-report. Exactness is a property of the residual, tested directly.
- **Compute the exact numerator `N = ac + bd` unconditionally to divide
  once.** Refuted (already by ADR-0090): the exact sum's bit-length tracks
  the exponent gap and can exceed any representable precision. The residual
  is only formed when the cheap screen admits it and bails when it cannot be
  formed exactly.
- **Make the saturating `c² + d²` divide return the correct value.**
  Refuted for this patch: it needs a Smith-style scaling reformulation of
  the divide, a new algorithm, not a flag fix. Raising the honest
  OVERFLOW/UNDERFLOW is the in-scope, in-API remedy; the value stays a
  flagged best-effort.

## Related

- Beads: pf-pz9r, pf-yprp, pf-hdq1, pf-bv2i (epic pf-8iji, review
  remediation R4)
- Other ADRs: ADR-0090 (mul/div and the directed-pair enclosure; the
  exact-quotient-is-OK contract), ADR-0091 (Annex G branch cuts, the
  §G.5.1 recovery, the cross-function signed-zero rule, and the M-OVF
  exponent exemption this extends to the underflow side), ADR-0093 (the
  frozen 1.0 API these fixes respect), ADR-0101 (the scalar
  negate-after-mirrored-round idiom `csqrt` borrows)
- Code: `pfloat-complex/src/div.rs`, `pfloat-complex/src/csqrt.rs`,
  `pfloat-complex/src/specials.rs`
- Tests: `pfloat-complex/tests/regression_pf8iji.rs`
