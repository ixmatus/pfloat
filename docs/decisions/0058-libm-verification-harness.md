# ADR-0058: The pfloat-libm verification harness (MPFR only, value hard, range sharded)

Status: Accepted (2026-06-01)

## Context

`pfloat-libm` (ADR-0056, ADR-0057) is a thin shell that rounds pfloat's
correctly rounded kernels to a hardware float through a directed pair outer
Ziv loop. The shell's output is therefore already an `f32` or `f64` paired
with a `Status`. Nothing yet proves that output correct against an oracle
independent of pfloat. This ADR records the verification harness that
supplies the proof: an exhaustive sweep of the 2^32 `binary32` inputs to
each unary function, a structured differential for `binary64`, and the
per function status records that publish the result.

Three forces shape the design. The dependency direction
(`pfloat-libm` depends on pfloat) forbids importing pfloat's own test
oracle, which lives in pfloat's `tests/` and is not part of pfloat's
public API. The exhaustive `f32` claim needs a sweep far larger than any
prior pfloat sweep, so it must shard across machines. And the shell
returns a `Status`, so the harness can check the IEEE flags, not only the
value.

## Decision

### 1. The oracle is MPFR only, reimplemented inside pfloat-libm

The harness lives entirely in `pfloat-libm` (a `tests/harness/` shared
module, the `tests/libm_*.rs` consumers, and `examples/libm_sweep.rs`)
and reimplements the thin oracle bridge there over `rug` as a unix only
dev dependency. It does not import pfloat's `tests/oracle/`. The libm
bridge is simpler than pfloat's: the shell already produced the hardware
float, so the oracle only encloses `f(x)` in two `rug::Float` endpoints
(`f(x)` rounded toward minus infinity and toward plus infinity at a
working precision) and rounds both ends to the hardware width under the
requested mode. Agreement of the two ends certifies the expected bits;
disagreement grows the working precision (the Ziv at oracle loop,
`START_PREC = 64`, doubling, `MAX_PREC = 1024`). There is no `BigFloat`
to hardware bridge, no `OracleBackend` trait, and no Arb backend.

Finding: pfloat's own oracle routes cot, sec, and csc through the Arb
subprocess on the stated belief that "MPFR has no cot/sec/csc primitive"
(`tests/oracle/meta.rs`, `tests/oracle/mpfr.rs`, ADR-0056 notes). That
belief is mistaken. MPFR ships `mpfr_cot`, `mpfr_sec`, and `mpfr_csc`,
and `rug` 1.30 exposes them as `cot_ref`, `sec_ref`, and `csc_ref`. The
libm harness verifies all three directly against MPFR, which makes its
oracle for those functions more independent than pfloat's (it shares
neither pfloat's kernel nor pfloat's Arb path). Every v0.1 function is
an MPFR primitive in `rug`, including cbrt, hypot, and the IEEE
754-2019 n-th root (`root_ref` for positive order, `root_i_ref` for
signed order), so the harness is Python free with no Arb venv. A
pfloat side cleanup to drop the unnecessary Arb routing for cot, sec,
csc is filed discovered-from pf-lm3.

### 2. Value is the hard gate; INVALID and DIV_BY_ZERO are also hard

The bit pattern is the headline claim, so value bit exactness (NaN
aware, with signed zero distinguished by the bit pattern) is the hard
gate: zero mismatches under nearest even, and under the directed modes
wherever the enclosure certifies a unique float. A measure zero hard to
round input that still straddles at `MAX_PREC` is recorded as
inconclusive, not a failure, the same caveat pfloat's own sweep carries.

Two status flags are gated hard alongside the value, and both
expectations are derived independently of the shell rather than from a
hand written domain table:

- INVALID is read straight from the enclosure: a non NaN input whose
  true value is NaN is a domain error (`ln` of a negative, `asin` out of
  range, a trigonometric function of an infinity), and the enclosure
  being NaN is the witness. A NaN input producing NaN is propagation,
  not INVALID.
- DIV_BY_ZERO is read from a small exact pole set: the flag is expected
  only at an exactly representable pole (`ln`, `log2`, `log10` at zero;
  `log1p` at minus one; `atanh` at plus or minus one; `cot`, `csc` at
  zero; `rootn` of zero with negative order). Irrational poles (`cot`,
  `sec` at multiples of pi or pi over two) are never a representable
  float, so the value there is finite. Deriving the flag from "the
  enclosure is infinite" instead would misfire on the overflow regime,
  where MPFR's own exponent range overflows `exp(huge)` to infinity
  though the true value is finite; the exact input list avoids that.

The remaining flags (INEXACT, OVERFLOW, UNDERFLOW) are recorded but not
gated. INEXACT specifically cannot be gated: the directed pair shell
conservatively over reports it on composed exact results such as
`log10(1000) = 3` (pf-njs5, ADR-0057), so an exact INEXACT gate would
fail a correctly rounded value. The smoke battery enumerates every
domain and pole input and confirms the derived flags agree with the
shell, so a future disagreement is surfaced as a real signal rather than
absorbed.

### 3. Exhaustive is unary only; binary is differential

The 25 unary functions are swept exhaustively over the 2^32 `binary32`
grid. The two binary functions (`hypot`, `rootn`) cannot be exhausted
over either axis and rest on a structured differential plus worst case
inputs at both widths, alongside the Lefevre and Muller hard to round
corpus for the 20 elementary functions the corpus covers.

### 4. Range sharding, NearestEven exhaustive plus directed sampled

`examples/libm_sweep.rs` splits one function's `[0, 2^32)` input space
into contiguous shards with `--shard-index` and `--shard-count`, the
capability pfloat's pf-hcz4 runner lacked (it sharded only one function
per instance and never split a function's space). The shard arithmetic
is `u64` throughout so the final shard reaches exactly 2^32 without
wrapping a `u32`; a unit test asserts the shards partition `[0, 2^32)`
with no gap or overlap.

Each shard verifies NearestEven over its full range (the headline "every
`f32` input correctly rounded under nearest even" claim) and the four
directed modes over a strided subsample of the range. This keeps the
exhaustive cost near the nearest even cost rather than five times it,
while still certifying the directed modes against the oracle across the
range. It matches pfloat's own posture and the kernel list's claim
("nearest even, and the directed modes wherever the enclosure determines
them"). The runner emits a status TOML row mirroring pfloat's schema and
a JSON sidecar the aggregator merges across shards.

### 5. EC2 ceremony adapted from pf-hcz4, Arb venv dropped

The fan out, poll, relaunch, and aggregate scripts adapt pfloat's
pf-hcz4 ceremony (ADR-0049). The Arb venv bootstrap is dropped (the
harness is MPFR only), so cloud init installs only the GMP and MPFR
build prerequisites. Every Noble AMI mitigation is kept: `m4` for the
`gmp-mpfr-sys` vendored build, `HOME` set under `set -u`, rustup install
with retries, the self terminating shutdown trap, and on demand
instances over spot. The launch is gated behind an explicit pay
confirmation and a `--dry-run` that prints the fan out and cost without
launching; the aggregator validates that the merged shard ranges cover
`[0, 2^32)` with no hole before writing a function's exhaustive row.

## Consequences

- The harness is generic over the width through an `Hw` trait (the
  analogue of the shell's `Shell` trait): one Ziv at oracle loop, one
  driver, and one sweep runner serve both `f32` and `f64`, so the two
  lanes cannot drift apart.
- The oracle is genuinely independent of pfloat for every v0.1 function,
  including the reciprocal trigonometric functions, because MPFR
  computes them with its own algorithm. The exhaustive `f32` sweep is
  the strongest correctness claim a pure Rust libm can make.
- The flag gate makes the harness sensitive to a pfloat conformance gap
  on INVALID or DIV_BY_ZERO, not only a value error. The cost is a small
  risk of a false flag failure if the derived expectation and the shell
  disagree on an edge; the smoke battery is where that is caught before
  any paid sweep.
- The exhaustive run is a paid EC2 campaign sized by a local micro
  benchmark and gated behind an explicit confirmation, not an automatic
  cost.

## References

- ADR-0057: the directed pair outer Ziv loop the harness verifies.
- ADR-0056: the six direct libm kernels.
- ADR-0049: the pf-hcz4 EC2 cross check ceremony this adapts.
- ADR-0034: pfloat's enclosure based oracle harness the bridge mirrors.
- pf-njs5: the pfloat INEXACT flag conservatism that makes INEXACT a
  recorded, not gated, flag.
- IEEE 754-2019 §4 (rounding), §7 (exceptions), §9.2 (recommended
  functions and the n-th root sign rules).
