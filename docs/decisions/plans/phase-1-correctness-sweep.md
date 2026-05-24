# Phase 1: pfloat correctness sweep

- **Status**: accepted plan, runs to completion before the v1.0 tag (re-sequenced 2026-05-22 per ADR-0033; was originally queued behind v1.0).
- **Date**: 2026-05-22.

## Why this phase

pfloat's README headline claim is "correctly rounded
transcendentals." That claim is per function, not surface wide,
and today it is not true for most of the surface. The basic
arithmetic operations are correctly rounded. `pow` is correctly
rounded under all five IEEE rounding modes (ADR-0022, the Ziv
interval test). Nearly every special function is `NearestEven`
only, with the differential lane using a `close_within(p − 2)`
tolerance (faithful, not correctly rounded), and the distinction
is not surfaced anywhere a user can read it. A user evaluating
"can I depend on this" cannot today get a per function answer
without reading the source.

A sharper concern. "Matches MPFR" does not prove "correctly
rounded." MPFR's correct rounding guarantee covers the basic ops
and a documented subset of functions. For specials where MPFR is
itself the implementation (`mpfr_ai`, `mpfr_zeta`, `mpfr_eint`)
MPFR uses Ziv internally and is believed correctly rounded in
practice, but it does not ship the formal guarantee. For
modified Bessel I and K MPFR has no primitive at all; the
existing in tree lane uses an mpmath table plus a DLMF identity
for exactly this reason. Differential bit exact agreement with
MPFR is therefore the strongest statement the current
verification can make for half the special surface, and that
statement is not equivalent to "correctly rounded."

A third concern. There is no external conformance vector suite,
only differential and property tests. ferrodec ships decTest
(8721 external vectors authored by the standards body) and that
authority is the load bearing part of ferrodec's credibility.
pfloat's analog is the Lefèvre and Muller worst case rounding
tables, cited in the README as a provenance source but not yet
wired as integration vectors. The differential sweep is bounded
random within representable inputs; the worst case rounding
inputs are precisely the adversarial cases least likely to be
hit by random or quasi random sampling. The highest value test
layer is the one that is not there.

Phase 1 closes all three at once. A single harness exhausts the
unary surface over `f32` against a certified oracle, fixes what
is wrong, and publishes a per function correct vs faithful
status table that users can read directly. The same harness is
the verifier the Phase 2 libm reuses verbatim, which is the
second reason this is the one phase that must come before any
libm spinoff.

## One line discipline

Freeze the surface, sweep it, fix it, ship the table. Do not add
algorithmic functions mid phase; a moving surface means the
table never stabilizes.

## Scope

In. Every unary (single real argument, `f32 → f32` shaped)
elementary and special function currently in pfloat. These have
a 2³² input space and so can be verified by enumeration, not
sampling. The boundary cases (signed zero, subnormals, infinity,
NaN, domain edges like `asin` at ±1, poles like `gamma` at non
positive integers) are covered automatically; that is the real
prize of exhaustion, not just the happy path.

Out, deferred:

- Multi argument functions (`pow`, `atan2`, `beta`, `fma`,
  `agm`). 2⁶⁴ plus input space; cannot be exhausted over `f32`.
  Stay on differential plus worst case vectors; their rounding
  rigor is a later track.
- New algorithmic functions (incomplete gamma and beta,
  Lambert W). These are Phase 3 additions; each gets swept as it
  lands using this same harness. Building them now widens an
  unverified surface on top of an unverified surface.
- The libm shell itself. Phase 2 work, after this harness is
  green and the status table is published.
- Pre freeze "trivial aliases" (cot, sec, csc, cbrt, hypot,
  rootn). The original plan proposed adding these as a few line
  compositions over the existing surface (`cot = 1 / tan`,
  `hypot = sqrt(a * a + b * b)`, etc.) so the frozen surface
  looks complete to a user scanning against libm or MPFR. That
  step is removed. A correctly rounded `tan` followed by a
  correctly rounded reciprocal is not a correctly rounded `cot`
  for hard to round inputs; the two roundings compose into up to
  1 ULP error in the reciprocal direction. The frozen surface
  for Phase 1 lists these six as `absent, deferred to the libm
  phase`, and Phase 2 ships them as direct primary kernels with
  their own derivations. See ADR-0032 for the full reasoning.

## The frozen unary surface

Sweep targets (confirm against `src/math/` at freeze time):

```
exp expm1 exp2 exp10
ln log1p log2 log10
sin cos tan
asin acos atan
sinh cosh tanh
asinh acosh atanh
sqrt
erf erfc
gamma lgamma digamma
zeta Ei Si Ci li
Ai Bi Ai_prime Bi_prime
J0 J1 Jn   Y0 Y1 Yn
I0 I1 In   K0 K1 Kn
```

Bessel note. Unary per fixed integer order. Exhaust `J0`, `J1`,
`Y0`, `Y1`, `I0`, `I1`, `K0`, `K1` over all `f32`; sweep a
chosen set of higher orders (`n = 2..N`) and record `N` in the
status table. The order parameter is discrete, so "exhaustive
in `x` for each swept order" is the honest claim.

## Oracle layer: two backends behind one trait

The verifier never compares pfloat to "the oracle's rounded
answer." It compares pfloat to the unique `f32` inside the
oracle's certified enclosure. The trait returns a proven
bracket, not a value. This makes MPFR and Arb interchangeable
behind one interface and makes the "is the rounding determined"
question uniform.

```rust
/// A proven bracket of the true value: lo <= true <= hi.
/// Endpoints are exact at the oracle's working precision; the
/// bracket tightens as working precision rises.
pub struct Enclosure { pub lo: OracleReal, pub hi: OracleReal }

pub trait OracleBackend {
    /// Enclose f(input) at `working_prec` bits. Result is a
    /// proven bracket, never a rounded scalar.
    fn enclose(&self, f: FnId, input: F32Bits, working_prec: u32)
        -> Enclosure;
    fn name(&self) -> &'static str;
}
```

Every finite `f32` converts exactly into the oracle's high
precision type, so the input is never itself a source of
rounding.

MPFR backend (in process, via `rug`, behind the existing
`differential-mpfr` feature lane). Produce the bracket by
evaluating the function twice, once with directed rounding down
and once up, at `working_prec`: `[f_RNDD(input), f_RNDU(input)]`.
Those two values provably bracket the true result; the bracket
tightens as `working_prec` grows. Fast, in process, ideal for
the bulk 2³² enumeration. Covers every function MPFR has a
primitive for.

Arb backend (subprocess, via `python-flint`). Arb is ball
arithmetic; it returns `[midpoint − radius, midpoint + radius]`
natively, a proven enclosure by construction, exactly the
`Enclosure` shape. Use it for functions with no MPFR primitive
(I and K Bessel, Si, Ci, li, Airy, future Lambert W and
incomplete gamma or beta), where Arb is the primary enclosure
oracle, and to spot check a subset of the MPFR covered
functions.

Do not spawn a process per input. Run a long lived `python-flint`
worker that reads input batches over a pipe and streams back
`(lo, hi)` pairs. The Rust harness owns enumeration and the
boundary check; the worker only encloses.

Binding caveat to verify before committing. Confirm the chosen
binding exposes the ball's radius (or both endpoints), not just
a rounded midpoint. The whole enclosure straddle check depends
on reading the certified radius out. `python-flint`'s `arb`
type carries midpoint plus radius and surfaces both; confirm in
the version pinned.

Licensing. FLINT and Arb are LGPL; the oracle is a CI lane,
test time only, never linked into the shipped crate. Same
posture as the existing MPFR lane. Behind a feature flag, out
of the published dependency tree.

## Verification core

```rust
/// The unique f32 the true value must round to under `mode`,
/// or None if the enclosure is too loose to decide.
fn certified_round_f32(enc: &Enclosure, mode: RoundingMode)
    -> Option<f32>
{
    let lo = round_to_f32(&enc.lo, mode);
    let hi = round_to_f32(&enc.hi, mode);
    if lo == hi { Some(lo) } else { None }   // bracket straddles
}

/// Verify one input under one rounding mode, with Ziv at oracle.
fn verify_input(
    o: &dyn OracleBackend, f: FnId, x: F32Bits, mode: RoundingMode,
) -> Verdict {
    let mut p = START_PREC;          // comfortably above f32's 24 bits
    loop {
        let enc = o.enclose(f, x, p);
        match certified_round_f32(&enc, mode) {
            Some(expected) => {
                let got = pfloat_round_f32(f, x, mode);
                return if got == expected {
                    Verdict::Ok
                } else {
                    Verdict::Mismatch { x, mode, expected, got }
                };
            }
            None => {
                if p >= MAX_PREC {
                    return Verdict::OracleInconclusive { x, mode };
                }
                p *= 2;
            }
        }
    }
}
```

`OracleInconclusive` is rare and is not a pfloat failure. It
means the oracle could not certify the rounding at `MAX_PREC`.
Those inputs are captured separately; they are candidates for
the worst case vector set and for raising `MAX_PREC` for that
function in a future sweep.

## Per function oracle routing

Each function declares its route. Three classes:

```rust
enum OracleRoute {
    /// MPFR has a primitive: MPFR is primary, Arb spot checks.
    MpfrPrimary,
    /// No MPFR primitive: Arb is primary, plus an independent
    /// identity cross check pfloat already has.
    ArbPrimary { identity: IdentityCheck },
    /// Hardest specials: Arb plus identity plus the checked in
    /// mpmath reference table as a third independent signal.
    ArbPlusTable { identity: IdentityCheck, table: TableId },
}
```

| Function group | Route | Independent cross check |
|---|---|---|
| exp / ln / log family, sin / cos / tan plus inverses, sinh / cosh / tanh plus inverses, sqrt, erf / erfc, gamma / lgamma / digamma, zeta, J and Y Bessel, Ei | `MpfrPrimary` | Arb spot check |
| I and K Bessel | `ArbPrimary` | I and K cross tie `I_{ν+1} K_ν + I_ν K_{ν+1} = 1 / z` |
| Si, Ci, li | `ArbPrimary` | existing reference table + dyadic self consistency |
| Airy Ai, Bi, Ai′, Bi′ | `ArbPrimary` | Wronskian `Ai · Bi′ − Ai′ · Bi = 1 / π` |
| (Phase 3) Lambert W, incomplete gamma / beta | `ArbPlusTable` | defining identity (`W · e^W = x`, etc.) |

The cross check on the Arb only rows is load bearing. For those
functions Arb is the primary oracle, not a second opinion, so it
must not be treated as infallible ground truth. The identity
already computed in pfloat's own tests is the cheap independent
signal that catches an oracle error.

## Exhaustive `f32` harness

- Enumeration. Iterate all 2³² bit patterns per function.
  Include every special and boundary input explicitly. Exhaustion
  means the special case and exception handling is verified too,
  not just the happy path; that is the real prize of an
  exhaustive sweep.
- Sharding. Partition the 2³² space across cores; embarrassingly
  parallel. A shard coordinator detects "this shard has been
  running an order of magnitude longer than its peers, the
  function's hard region is over represented here" and rebalances
  rather than running serially through a function's worst inputs.
- Exhaustive vs sampled, recorded honestly. Cheap functions get
  literal 2³² enumeration. For the most expensive Arb only
  specials, where IPC plus Arb evaluation over 2³² is
  prohibitive, a function may instead run dense sample plus all
  boundary and special inputs; the status table records
  `exhaustive` vs `sampled(N)` per function. Never claim
  exhaustive where sampled.
- Cadence. Per release CI lane, not per commit (it is CPU hours).
  Per commit keeps the existing fast differential plus property
  suite, plus the regression corpus (see below).

## Status table output schema

The sweep emits one row per (function, rounding mode set):

```
function          : sin
order             : -            # for Bessel: the swept order(s)
kernel_kind       : primary      # primary | derived_alias
domain_coverage   : exhaustive   # exhaustive | sampled(2^28)
oracle            : MPFR-primary # MPFR / Arb / Arb+table
oracle_independence: independent # independent | shared_algorithm_class
rounding_modes    : RNE,RNA,RZ,RP,RM  # the modes claimed
rounding_status   : correctly-rounded # correctly-rounded | faithful | has-errors
worst_ulp         : 0.0          # max observed error in ULP
mismatch_count    : 0
inconclusive_count: 0            # OracleInconclusive inputs
panic_count       : 0            # pfloat-side panics seen
vectors           : tests/vectors/sin_regression.bin
```

This artifact is published in the README and docs as the per
function table (the credibility document; the thing that lets a
user adopt the crate per function) and is the machine readable
input to CI gating.

`rounding_status` is the verdict that matters:

- `correctly-rounded`: zero mismatches across the swept space.
- `faithful`: mismatches all `≤ 1 ULP`; honestly downgraded in
  the table rather than over claimed.
- `has-errors`: anything worse; a Phase 1 blocker until fixed.

`correctly-rounded` is structurally unavailable for
`derived_alias` rows; the highest verdict an alias can earn is
`faithful` (ADR-0032).

## Failing input capture, regression corpus

Every `Mismatch` and every `faithful` not `correct` input has
its exact `f32` bits (and mode) appended to a per function
checked in vector file. Fixes are then regression tested against
the captured inputs, and the expensive exhaustive sweep does not
have to re discover known hard cases on every run.
`OracleInconclusive` inputs go to a separate file (worst case
candidates), not the failure corpus. `Panic` inputs go to a
third file (panic regression) and run on every CI push, not only
at the next per release sweep.

## Sequencing within Phase 1

1. Surface freeze. Confirm the surface list against `src/math/`;
   freeze it. No trivial aliases added (ADR-0032).
2. Oracle trait plus MPFR backend. `Enclosure`, `OracleBackend`,
   `enclose` via the directed rounding bracket through `rug`.
3. Verification core. `certified_round_f32`, `verify_input`
   with Ziv at oracle; harness scaffold (sharding, per function
   driver, status row emitter).
4. Run on MPFR covered functions first. Fast feedback; validates
   the harness end to end before the Arb backend exists. The
   standard elementary functions should come back clean.
5. Layer the Lefèvre and Muller worst case vectors onto the
   standard elementary functions during step 4. They target
   precisely the hard to round inputs a random or enumerated
   sweep is least likely to stress, and they let those functions
   claim provable (not merely exhaustive against an oracle)
   correct rounding.
6. Arb backend. Long lived `python-flint` worker plus streaming
   protocol; wire the `ArbPrimary` and `ArbPlusTable` routes.
7. Run the Arb only specials. I and K, Si, Ci, li, Airy. These
   are the functions most likely to surface `faithful` not
   `correct` or a Ziv ladder bug, because they are the ones
   with no MPFR oracle today.
8. Triage, fix pfloat, re run. Fix `has-errors` to
   `correctly-rounded`, or knowingly document as `faithful`.
   Re run the sweep as the regression gate.
9. Publish the status table. README and docs; CI gates on it.

## Exit criteria

Phase 1 is done when:

- Every frozen unary function has a status table row with a
  definitive `rounding_status`.
- No function is `has-errors`. Each is either correctly rounded
  or documented `faithful` with rationale.
- The regression corpus (mismatch, inconclusive, panic) is
  checked in and CI gated.
- The per function status table is published.
- The oracle trait plus harness are reusable as is for Phase 2
  (the libm shell re points the same harness at the shell output
  to catch `BigFloat → f32` double rounding).

Only then does Phase 2 (libm spinoff) begin.

## Slice p1.1 + p1.2 closure (2026-05-23)

Slice p1.1 extended the slice 8b L-M tier from nine functions to
twenty (added `exp10`, `expm1`, `log1p`, `sinh`, `cosh`, `asinh`,
`acosh`, `atanh`, `erf`, `erfc`, `gamma` clean; deferred `log2`,
`log10`, `tanh`, `lgamma` behind their fix slice).

Slice p1.2 closed the five-finding `has-errors` class on the v1.0
surface:

- Lifted `pow`'s Ziv interval-test driver (ADR-0022) into
  `src/math/ziv.rs` as a shared `pub(crate) fn ziv_round`. Pure
  refactor.
- Wired `exp`, `ln`, `tanh`, `lgamma` through `ziv_round`. `log2`
  and `log10` inherit via the existing `ln_round` composition.
- Reinstated CORE-MATH's `exp.wc` leading underflow block as the
  regression guard for the slice 8b documented `exp` underflow
  defect. The L-M corpus generator's marker-scan filter was
  removed entirely; the `domain_ok` filter handles per-function
  edge cases (±0 dropped; `log2`/`log10` reject `x ≤ 0`; `lgamma`
  rejects negative integers).
- Test driver compares at the f64 bit-pattern level rather than
  BigFloat-level (required for binary64 subnormal expected outputs
  whose p=53 BigFloat representation carries more precision than
  the binary64 source). Inputs constructed bit-exact via
  `bf53_of_bits` (integer mantissa plus chained 2^30 mul/div
  scalings).

The L-M corpus now covers twenty-four functions at 1200 bit-exact
cases. The five known `has-errors` findings on the v1.0 surface
(documented exp underflow plus log2/log10/tanh/lgamma corpus
deferrals) are closed. The exhaustive `f32` sweep (planned p1.6 -
p1.8) will surface any remaining latent trips; the corpus tier
catches the obvious ones.

## Slice p1.3 closure (2026-05-23)

Slice p1.3 implemented the Phase 1 plan's steps 1-4 (surface
freeze, Oracle trait + MPFR backend, verification core, run on
MPFR-covered functions). The slice landed eleven unsigned branch
commits closing the harness scaffold:

- **Surface freeze + ADR-0034**: `docs/v1.0-surface.md` enumerates
  the 47 frozen v1.0 unary surface entries (21 elementary + 10
  specials + 4 Airy + 12 Bessel fixed-order); ADR-0034 records the
  Oracle layer architecture (`Enclosure` bracket type,
  `OracleBackend` trait, MPFR backend's RNDD / RNDU directed-
  rounding bracket, `FnId` enum dispatch, Ziv-at-oracle
  precision doubling, status table schema, regression corpus
  capture, Arb-backend posture for the next slice, LGPL
  isolation).
- **Types and dispatch**: `tests/oracle/` ships `Enclosure`,
  `OracleBackend`, `FnId` (47 entries), `Verdict`,
  `MpfrOracle::enclose` with 35 MPFR-primary dispatch entries via
  a `bracket!` macro plus an `lgamma_bracket` helper, and the
  parallel pfloat-side `pfloat_kernel` dispatch routing all 47
  FnIds.
- **Verification core**: `convert.rs` adds `bf24_of_bits` (bit-
  exact f32-to-BigFloat construction; mirrors slice-p1.2's
  `bf53_of_bits` sized for binary32), `bf_to_f32_bits`, and
  `round_f32` (rug Float to f32 under any of the five IEEE
  rounding modes, with NearestAway synthesized to compensate for
  MPFR's absent roundTiesToAway primitive). `verify.rs` adds
  `START_PREC = 64`, `MAX_PREC = 1024`, `certified_round_f32`
  (with NaN-aware (NaN, NaN) handling that returns
  `Some(f32::NAN)`), and the `verify_input` Ziv-at-oracle loop.
- **Driver and status emitter**: `driver.rs` adds the per-
  function `run_function` runner with `std::panic::catch_unwind`
  capture, `outcome_to_status_row` builder, and
  `write_mismatch_corpus` binary serializer. `status.rs` ships
  the `StatusRow` schema and hand-written TOML emitter matching
  ADR-0034's schema verbatim.
- **Smoke gate**: `tests/oracle_smoke_gate.rs` runs each of the
  33 MPFR-primary `FnId` variants at 64 representative inputs
  under NE on every `differential-mpfr` CI push; ~20 seconds in
  debug, ~5 seconds in release. All 2112 verdicts return Ok at
  slice close.
- **Standalone runner**: `examples/oracle_sweep.rs` runs the same
  harness at a larger budget for per-release sweeps;
  configurable via `--function NAME`, `--exhaustive` /
  `--sample N`, `--mode MODES`, output paths.
- **First sweep in-tree**: 33 per-function TOML status rows in
  `tests/oracle/status/` plus three binary regression corpora
  in `tests/vectors/` from a 65536-input subnormal-range sweep
  under NE. Output:
    - 30 of 33 MPFR-primary FnIds: `correctly-rounded`.
    - tanh: `has-errors` (65535 mismatches), pf-7d7. Root cause:
      the `1 - exp(-2x)` cancellation underflows to exact 0 at
      every working precision for subnormal inputs; the slice-
      p1.2 `ziv_round` envelope cannot recover because
      `half_width(0)` is 0 and the interval test certifies the
      (wrong) 0.
    - erf: `has-errors` (111 mismatches), pf-z0f. Likely the
      slice-3a fixed-guard convention surviving in erf; the
      slice-p1.2 Ziv upgrade did not extend to erf.
    - J1: `has-errors` (16388 mismatches, all 1-ULP), pf-n5d.
      Likely faithful-not-correctly-rounded in the small-
      argument formula.
- Three defect beads (pf-7d7, pf-z0f, pf-n5d) and one
  deferred-feature bead (pf-r2b9: L-M corpus inputs as
  adversarial seeds in the runner) are filed for slice p1.4+.

Slice p1.3 ships the harness; p1.4 closes the three findings and
extends the runner with L-M adversarial seeds. Phase 1's Arb
backend (`ArbPrimary` and `ArbPlusTable` routes per
ADR-0034) ships in a follow-up slice; until then the 12
Arb-primary FnIds (`Si`, `Ci`, `Li`, `Bi`, `Ai_prime`,
`Bi_prime`, modified Bessel `I` and `K` families) verify only
through the existing differential lane plus identity cross
checks.

## Slice p1.4 closure (2026-05-23)

Slice p1.4 closed all three slice p1.3 has-errors findings
(pf-7d7, pf-z0f, pf-n5d), extended the runner with L-M
adversarial seeds (pf-r2b9), and wired three more elementary
kernels onto the canonical Ziv envelope (erf, the Bessel J
family). The slice landed seven unsigned branch commits:

- **tanh tiny-input short circuit (slice p1.4.1, closes pf-7d7).**
  The composition `(1 - exp(-2|x|)) / (1 + exp(-2|x|))` collapses
  to exactly zero for tiny inputs (when `2|x| < 2^-w` at working
  precision `w`, `exp(-2|x|)` rounds to one and the numerator
  becomes zero); the slice p1.2 Ziv envelope cannot recover
  because `half_width(0) = 0` makes the interval test certify
  zero as the answer. `src/math/tanh.rs::tanh_at_w` now short
  circuits to `|x|` for any input whose binary exponent is below
  `-ceil((working_prec - 22) / 2)`, the threshold derived from
  Taylor's theorem (`tanh(|x|) = |x| - |x|^3/3 + O(|x|^5)`) so
  the truncation error fits the Ziv driver's error guard
  `|y| * 2^-(working_prec - ZIV_ERROR_GUARD)`. Closed the 65535
  / 65536 f32 subnormal mismatches on the slice p1.3 sweep.
- **Harness bf -> f32 17-digit Display + initial p = 53 bump (slice
  p1.4.2).** Diagnosed pf-z0f as a harness conversion bug rather
  than a kernel defect: pfloat's `erf` at `p = 24` is correctly
  rounded at that precision (verified by probing the kernel at
  `p = 24, 53, 64, 113, 200`; all give identical decimal), but
  the bf -> f32 conversion through Display + parse loses
  information on f32-subnormal-grid midpoints because the
  Display digit count (9 for `p = 24`) carries enough precision
  for f32 normals but not for the subnormal grid spacing. Slice
  p1.4.2 bumped the digit count to 17 (= f64 round-trip) and
  raised the verification precision to `p = 53` globally; the
  second move turned out to introduce a directed-mode
  regression (see slice p1.4.6) and was scoped back to per-
  function in the same slice.
- **erf Ziv envelope (slice p1.4.3, closes pf-0qp9).**
  Architectural cleanup. The slice p1.2 pattern (exp / ln / tanh
  / lgamma) wraps the working-precision body in `ziv_round` to
  inherit the canonical correctness criterion of ADR-0022;
  `src/math/erf.rs` now follows. Regime decision is fixed at
  kernel entry; both `erf_maclaurin` and
  `super::erfc::erfc_asymptotic` flow under the Ziv envelope.
- **Bessel J Ziv envelope (slice p1.4.4, closes pf-ydna).**
  Same pattern around `bessel_j_eval_normal`: the regime
  dispatcher (tiny / Miller / asymptotic) runs under one Ziv
  envelope so J0, J1, Jn share the slice p1.2 correctness
  scaffolding.
- **Per-function verification precision (slice p1.4.5, closes
  pf-n5d).** J1's residual mismatch was a different harness
  bug: the Maclaurin's first correction term sits at relative
  `~2^-298` for the smallest f32 subnormal exponent, far below
  the slice p1.4.2 verification precision; the kernel's final
  round to `p = 53` stripped the correction and left the value
  on the exact f32-subnormal-grid midpoint where the conversion
  ties to even instead of tracking the true sub-midpoint
  position. The fix is per-function verification precision in
  the harness: `BesselJ1` and `BesselJn` route through
  `BESSEL_TINY_VERIFICATION_PRECISION = 320` (sized so the
  worst-case cubic correction at relative `2^-298` survives the
  kernel's final round with a 22-bit headroom).
- **Tighten verification precision (slice p1.4.6, closes
  pf-kg12).** Slice p1.4.2's global `p = 53` default introduced
  a directed-mode regression: the bf -> f32 bridge always uses
  NE rounding (Rust's f32 parser is NE-only), so at `p = 53`
  the kernel's directed-mode rounding (`TowardPositive` /
  `TowardNegative` / `TowardZero` / `NearestAway`) is silently
  overridden by the bridge's NE re-encode. At `p = 24` the
  kernel returns a value that lands exactly on the f32 grid, so
  the bridge is lossless under every mode. Slice p1.4.6 set
  `DEFAULT_VERIFICATION_PRECISION = 24` and reified the erf
  bump as `ERF_VERIFICATION_PRECISION = 53`; the Bessel J
  family keeps its `p = 320` bump from slice p1.4.5. Both
  bumped paths run under NE only in the f32 sweep, so the
  NE-only bridge correctly rounds them.
- **L-M adversarial seeds (slice p1.4.7, closes pf-r2b9).**
  `examples/oracle_sweep.rs` includes
  `tests/differential/lefevre_muller_data.rs` via a relative
  `#[path]` and prepends per-function L-M inputs (cast from f64
  bit patterns to f32) to the linear sweep iterator. The
  status row schema gained an `lm_seeds_run` field so the
  verification posture records the adversarial-seed count
  per row.

Slice p1.4 closure leaves the 33 in-tree status rows all reading
`correctly-rounded` and no regression corpora under
`tests/vectors/`. The per-push smoke gate stays under one
minute (~20 seconds debug for the default path, plus the J1
and Jn bumped runs at `p = 320`). Next: pf-cvs (1.x smallvec
inline BigFloat storage; do not pick up under the Phase 1
posture) remains the only open work bead aside from the
slice 8c.* parked items.

## Design notes from the 2026-05-22 critique pass

The following refinements landed during a critique of the
original work breakdown, drawing on the slice 8a evidence (the
beta `case 4` hang, the mul / div / fma exponent panic
fuzz found via Airy `bi_prime`, the parser cap reframe). They
are part of the accepted plan.

**Panic as a finding category, not a missing path.** The slice
8a fuzz find (a pre existing panic in `src/ops/mul.rs:189`
reachable via Airy `bi_prime` with extreme exponent inputs) is
the exact class of defect an exhaustive `bi_prime` sweep finds
in the first kilo inputs, and it is not a "rounding mismatch."
The harness must treat a pfloat side panic as `Verdict::Panic`,
emit a `panic_count` column, and capture the input to the
regression corpus. A function that panics on its seventeenth
swept input is wrong before any oracle is consulted; that is
the cheapest finding class and the harness must surface it.

**`oracle_independence` is a status table column.** MPFR for
`pow` is genuinely independent (different algorithm). MPFR for
`gamma` runs Stirling alongside pfloat running Stirling;
agreement is suggestive, not proof. The slice 8a differential
beta oracle used MPFR's `ln_abs_gamma` precisely to get MPFR's
own independent sign tracking; that reflex is right and the
status table must surface it per function so readers know what
the green check actually proves. The Arb only rows already
carry an independent identity cross check; the column makes
the same epistemic question explicit for the MPFR primary rows.

**`kernel_kind` is a status table column.** `correctly-rounded`
is structurally unavailable for `derived_alias` rows. ADR-0032
locks the libm reciprocal and root kernels as direct primary
kernels; the column is the type system level guard against a
future contributor sneaking an alias in under the assumption
that the sweep would catch the problem.

**`rounding_modes` is a status table column.** `pow` is
correctly rounded under all five IEEE rounding modes
(ADR-0022). Most pfloat specials are `NearestEven` only. The
column records which modes each function claims. Without it
`correctly-rounded sin` is ambiguous (under which modes?), and
the table cannot honestly mix functions of different rigor.

**The `Enclosure` should carry a `ProvablyExact` outcome.** A
bracket whose both endpoints are the same oracle value at the
working precision means the true value is exactly representable
at the target precision: every rounding mode produces the same
result, and the input warrants a stronger per row claim. The
slice 7c Ziv interval test (ADR-0022) hit the close cousin: the
original "recompute and compare" composition false converged on
hard to round inputs, replaced by the interval test. Making
`Verdict::ProvablyExact` an explicit third outcome lets the
status table emit a stronger per input claim where the input
warrants it.

**`OracleInconclusive` inputs feed per commit CI, not only the
next per release sweep.** Once the corpus of hard inputs is
bounded (the count is small in practice), running it on every
CI push is cheap. They also feed back into raising `MAX_PREC`
for that specific function in the next sweep. The function
specific Ziv ladder learns from its own hard cases; the manual
session level pattern (slice 8a's parser boundary test that
took 49 minutes at exponent `3 * 10^6`, shrunk to `1.1 * 10^6`)
is the kind of feedback that wants to be automatic.

**The CPU budget is realistic only at order of magnitude
granularity.** Today's experience (parser boundary differential
at `p = 113` against MPFR ran 5 minutes for two strings; an
earlier draft at exponent `3 * 10^6` ran 49 minutes) confirms
MPFR per call cost is bursty at high working precision. With
Ziv at oracle on top, tail latency on adversarial inputs at
high `working_prec` is seconds per input. The honest budget:

- Cheap functions (`sqrt`, basic exp / ln family): hours.
- Trig: a day or two, dominated by range reduction near
  multiples of π.
- Specials with Stirling or series tails (gamma family, zeta,
  erfc near `+∞`): days per function.
- Arb only specials with subprocess overhead: dense sample plus
  boundary inputs is the only feasible mode; the status table
  column making that explicit is doing real work.

This does not kill the plan. It makes the per release cadence
load bearing (which the plan already states) and makes shard
coordinator quality matter more than the original framing
suggested.

**Coordinate L-M with in repo Phase 8b.** The in repo Phase 8b
slice (the v1.0 conformance preparation) wires the Lefèvre and
Muller corpus as a checked in resource with a hard licensing
stop ask gate, ahead of the v1.0 tag. Phase 1 reuses that
corpus as the adversarial seed inside the exhaustive harness
rather than re vetting and re importing it. The corpus is the
durable artifact; the basic differential tier in 8b plus the
worst case adversarial seed in Phase 1 is the program.

**Sequencing relative to the v1.0 tag (revised 2026-05-22 per
ADR-0033).** Phase 1 runs to completion before the v1.0 tag, not
after. The original sequencing (Phase 1 as v1.0 → v1.1 work)
was reversed after slice 8b surfaced two findings that re-weighted
the credibility-vs-timeline trade-off:

- pfloat's `exp` kernel mis-rounds inputs at the binary64
  underflow boundary (the slice-8b L-M corpus excluded CORE-MATH's
  underflow block to keep the suite green; `src/math/exp.rs` lines
  23-30 carry the durable note). A library claiming
  correctly-rounded transcendentals cannot honestly ship 1.0 with
  a known wrong-rounding case in its headline kernel.
- The tier-2 specials' differential verification posture
  (mpmath table at `p ≤ 256` plus identity cross-ties) is
  structurally weaker than the README's correct-rounding claim
  implies. The Arb backend specified above closes that gap, but
  only once it is actually integrated and the sweep is run.

The published 1.0 on crates.io is immutable (yank-only, no
replace). It is the version users will quote, link to in their
own changelogs, and judge the project's permacomputing-horizon
discipline against permanently. The credibility cost of
"shipped a wrong-rounding case in the headline kernel" is
permanent in a way the timeline cost of the audit is not.

Revised reconciliation:

- v1.0 ships when every frozen unary function has a definitive
  `rounding_status` (correctly-rounded or honestly-downgraded
  faithful); no function is `has-errors`. The per-function status
  table replaces the count-only `## Conformance evidence` section
  at the v1.0 cut.
- Slice 8c (the v1.0 tag + `cargo publish` slice) is parked
  indefinitely behind Phase 1. The 8c beads stay on file; the
  disclosure-correction diff artifact at
  `docs/disclosure-correction-v1.0.diff` stays in tree (the two
  factual corrections it makes are still required at the eventual
  tag, just later).
- `has-errors` findings in the sweep are v1.0 blockers, not v1.1
  fixes. The exp-underflow defect is the first known one; the
  CORE-MATH corpus expansion across all pfloat-supported binary64
  functions will surface more.

The accepted timeline cost (six to twelve months per the CPU
budget above) is borne entirely by the project; pfloat has no
published 0.x consumers paying a migration cost. ADR-0033 is the
durable record of this re-sequencing.

## Appendix A: Path to performance, post correctness

Performance is not part of Phase 1 and not a near term concern.
This appendix exists for two reasons: so the correct kernels
are written in a way that will not have to be thrown away (see
A.4, actionable during Phase 1), and so the eventual
performance phase has a stated shape and exit criterion.

### A.1 Framing

Correctness first is what unlocks aggressive performance work.
After Phase 1 the project holds the one asset performance
optimization requires and most projects lack: an exhaustive
`f32` oracle plus a per function regression corpus. That lets
any kernel be rewritten as violently as needed with proof that
behavior did not change.

Target (per DESIGN.md): "documented gap to MPFR, never absurd."
Not "beat MPFR" — MPFR rides GMP's decades of per architecture
assembly, and the project's principles forbid FFI to close the
last percentages. The honest goal is same algorithmic
complexity class as MPFR, a constant factor behind, with the
gap measured and published per function per precision.

### A.2 The ordered path: complexity before constants

The wins are wildly unequal; do them in this order.

1. Algorithmic complexity (the order of magnitude wins).
   - **FFT or NTT multiplication** (ADR-0010, currently deferred):
     the single largest performance item in the library. At
     high precision everything bottoms out in multiply: every
     transcendental, every Newton step, every AGM iteration.
     Schoolbook `O(n²)` → Karatsuba `O(n^1.58)` → FFT
     `O(n log n)` is the difference between "usable at
     thousands of digits" and "falls off a cliff." Implement
     after the sweep exists, with the oracle watching; it is
     the scariest kernel to get right.
   - **Newton based division and square root** computed via
     the reciprocal, so they ride fast multiplication instead
     of naive long division.
   - **Karatsuba thresholds** (ADR-0027) re calibrated once
     FFT lands (three way crossover: schoolbook to Karatsuba
     to FFT).

2. Transcendental algorithm choices (complexity again).
   - **Binary splitting** evaluation of series, rather than
     naive term by term summation; turns an N term rational
     series into something that rides fast multiply and is
     far better asymptotically.
   - **AGM** logarithm and constants (ADR 0015 and 0017) are
     already the right asymptotic choice; leave them.
   - **Tighten argument reduction** to cut the Ziv retry rate.

3. Constant factors (least impactful; do last).
   - Limb kernels: tight carry propagation, auto vectorizable
     inner loops.
   - **Allocation and inline storage** (ADR-0028): small
     precision `BigFloat`s should not heap allocate.
   - SIMD last, if at all; architecture specific, fights the
     pure Rust portability story, and is a constant factor on
     top of work that should already be complexity optimal.

### A.3 The instrument: a standing MPFR benchmark lane

Turn `benches/` into a tracked comparison: per function, across
a precision sweep (53, 113, 256, 1024, 4096 bits), pfloat vs
`rug` timing, tracked over time so regressions are visible and
the "documented gap" is a published number per function per
precision. The performance analog of the correctness status
table.

### A.4 Phase 1 guardrails (actionable now, cost nothing)

To avoid a rewrite when the performance phase comes, honor
these while writing the correct kernels:

- Multiplication behind a clean interface so an FFT backend
  can slot in without touching callers. (Slice 8a's mul, div,
  fma exponent saturation fix accidentally validated this:
  three kernels patched with the same shape because callers
  did not know the algorithm underneath.)
- No schoolbook or `O(n²)` assumption baked into callers (no
  loops that assume a particular multiply cost shape).
- Division and square root callers agnostic to whether the
  implementation is long division or Newton via reciprocal.
- The limb kernel API does not leak a fixed algorithm
  assumption upward.

Interface hygiene, not optimization. Free during correctness
work, saves a rewrite later.

### A.5 The hard discipline

Never optimize a kernel while correcting it. A kernel being
optimized and corrected at once is one where the engineer
cannot tell which change the oracle just caught. Correct first
(oracle green); then make it fast (oracle still green as proof
of no regression).

### A.6 Exit criterion (performance phase)

- Gap to MPFR measured and documented per function per
  precision.
- No function asymptotically worse than MPFR's algorithm class.
  In particular FFT multiply shipped, Newton division and
  square root, binary split series.
- The benchmark lane standing and regression gated.

The performance phase does not begin until correctness ships.

## Related

- `docs/ROADMAP.md` — Phase 1 in the wider ecosystem sequence.
- `docs/decisions/0032-libm-reciprocal-and-root-kernels-direct.md`
  — the one Phase 1 / Phase 2 discrete decision already discharged.
- `~/.claude/plans/abundant-yawning-badger.md` — in repo Phase 8
  per slice detail (where the L-M corpus first lands as a
  checked in resource).
- ADR-0014 (sweep size convention), ADR-0022 (Ziv interval test
  for `pow`), ADR-0029 (Dragon4 deferred to 1.x, picked up in
  Phase 3 of the roadmap).
