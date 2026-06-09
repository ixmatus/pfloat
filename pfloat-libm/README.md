# pfloat-libm

Pure Rust correctly-rounded `libm`.

## Overview

pfloat-libm is a correctly-rounded implementation of the elementary math
functions (the `libm` surface) for `f32` and `f64`. It is a thin shell over
[pfloat](https://github.com/ixmatus/pfloat): a call widens the hardware float to
an arbitrary-precision `BigFloat`, evaluates the function correctly-rounded
through pfloat's kernel, and rounds the result back to the hardware width under
an outer Ziv loop that commits a value only once an enclosure proves it. The
second rounding is therefore never blind, which is what closes the
`BigFloat` to float double-rounding gap a naive widen, compute, round shell
leaves open.

The niche pfloat-libm fills sits between the existing pure-Rust
[`libm`](https://crates.io/crates/libm) crate (a faithful, roughly one ULP port
of fdlibm and musl, not correctly rounded) and the C-backed correctly-rounded
stacks (CRlibm, CORE-MATH, and MPFR via [`rug`](https://crates.io/crates/rug),
all of which need a C toolchain). pfloat-libm is pure Rust, `no_std` with `alloc`,
permissively licensed, and its unary `f32` surface is verified
exhaustively: every one of the 2^32 `binary32` inputs is checked against
an independent oracle under all five rounding modes, a claim no competing
pure-Rust libm makes. The `f64` surface rests on differential testing
against a reference oracle plus the published worst-case hard-to-round
vectors, since the 2^64 input space cannot be enumerated.

pfloat-libm is `no_std` + `alloc`. There is no alloc-free profile: correct
rounding grows the working precision at runtime past any compile-time width, so
the computation allocates.

## Alternatives

The honest comparison. pfloat-libm does not claim to be the fastest or the most
complete correctly-rounded libm; CORE-MATH and CRlibm are mature, faster, and
broader. What pfloat-libm offers is a combination none of the others does: pure
Rust with no C toolchain, a permissive license, correct rounding, and an
exhaustive `f32` verification.

| Library | Pure Rust | Correctly rounded | `no_std` | No C toolchain | Exhaustive `f32` check | License |
| --- | --- | --- | --- | --- | --- | --- |
| **pfloat-libm** | yes | yes | yes (`alloc`) | yes | yes (unary) | MIT OR Apache-2.0 |
| [`libm`](https://crates.io/crates/libm) | yes | no (~1 ULP) | yes | yes | no | MIT OR Apache-2.0 |
| CORE-MATH | no (C) | yes | — | no | yes (its own) | MIT-style |
| CRlibm | no (C) | yes | — | no | proven, not enumerated | LGPL |
| MPFR via [`rug`](https://crates.io/crates/rug) | no (C / FFI) | yes | — | no | n/a (it is the oracle) | LGPL |

The differentiator is the row, not any single column. The only other pure-Rust
member, the `libm` crate, is a faithful `~1` ULP port and not correctly rounded;
every correctly-rounded alternative is C-backed and pulls a toolchain (and, for
CRlibm and MPFR, an LGPL dependency). pfloat-libm is the pure-Rust member that is
correctly rounded, and among pure-Rust libms the only one whose claim is checked
by enumerating the whole `binary32` grid rather than by sampling. The price is
performance: every call routes through arbitrary-precision `BigFloat`, so this is
a correctness-first libm, not a speed-first one.

## Usage

```rust
use pfloat_libm::{f32 as lm, RoundingMode};

// Correctly rounded to nearest even: the last bit is always right.
let r = lm::exp(1.5_f32);

// The directed form picks the rounding mode and returns the sticky flags.
// exp(1.5) is transcendental, so its upward and downward roundings differ by
// exactly one ULP.
let (up, status) = lm::exp_round(1.5, RoundingMode::TowardPositive);
let (down, _) = lm::exp_round(1.5, RoundingMode::TowardNegative);
assert_eq!(up.to_bits() - down.to_bits(), 1);
assert!(status.inexact());
```

Every unary function has a bare form (`lm::exp(x) -> f32`, rounded to nearest
even) and a directed form (`lm::exp_round(x, mode) -> (f32, Status)`). The
`f64` surface mirrors it under `pfloat_libm::f64`. A runnable walk through,
including the saturation fast path, lives in
`examples/correctly_rounded_exp.rs`.

## Verified surface

The unary `f32` surface is verified exhaustively: every one of the 2^32
`binary32` inputs is evaluated and compared against an independent MPFR oracle
(via `rug`, kept out of the shipped link graph) under all five IEEE 754 rounding
modes. Across all 25 unary kernels the result is uniform: **0 mismatches, 0
panics, worst error 0 ULP**, correctly rounded on every input in every mode. The
per-function status records are checked in under `tests/harness/status/`.

| Family | Functions | `f32` (2^32 x 5 modes) | `f64` |
| --- | --- | --- | --- |
| exponential | `exp` `exp2` `exp10` `expm1` | correctly rounded | differential |
| logarithm | `ln` `log2` `log10` `log1p` | correctly rounded | differential |
| root | `sqrt` `cbrt` | correctly rounded | differential |
| forward trig | `sin` `cos` `tan` | correctly rounded | differential |
| reciprocal trig | `cot` `sec` `csc` | correctly rounded | differential |
| inverse trig | `asin` `acos` `atan` | correctly rounded | differential |
| hyperbolic | `sinh` `cosh` `tanh` | correctly rounded | differential |
| inverse hyperbolic | `asinh` `acosh` `atanh` | correctly rounded | differential |

The 25 unary kernels are the 21 elementary functions pfloat shipped correctly
rounded at 1.0 plus the four net-new direct kernels `cbrt`, `cot`, `sec`, `csc`
(ADR-0032: a correctly-rounded reciprocal or root is a direct kernel, never an
alias of a composed one). The two multi-argument kernels `hypot` and `rootn` are
binary, so they cannot be enumerated over either width and rest on differential
testing plus the published worst-case vectors. The `f64` surface rests on the
same MPFR differential lane over a structured sample plus the Lefevre-Muller
hard-to-round vectors, since the 2^64 input space cannot be enumerated.
`docs/kernel-list.md` records the per-kernel implementation site and tier.

For the development process that produced this code, read the disclosure
immediately below before deciding whether to adopt.

## How pfloat-libm is developed

This is an open disclosure of the development process so users can judge for themselves whether the resulting code meets their bar.

**Authorship and collaboration.** Parnell Springmeyer is the author of record. pfloat-libm is developed in collaboration with Claude, an AI coding agent from Anthropic. Parnell owns architecture, acceptance criteria, test and verification strategy, and release boundaries. Claude drafts the implementation, writes and runs tests and verification harnesses, and produces analysis under that direction. **Parnell does not review the generated code line by line.** Human oversight operates at the level of design, strategy, and outcomes: does the architecture make sense, are the right invariants being checked, does the verification strategy cover the risk surface, do the tests and proofs pass. Merges to main are GPG signed by Parnell to attest to that level of review, not to an audit of every line.

**Provenance.** Implementations derive from primary sources: IEEE 754-2019 for floating point semantics and rounding modes, the DLMF for the reciprocal and root function definitions and their range reduction, and the published Lefèvre and Muller worst case rounding tables for the hard to round vectors. pfloat-libm is mostly a shell and a verification harness over pfloat; the arbitrary precision mathematics, including the direct reciprocal and root kernels, lives in pfloat and carries its provenance in pfloat's own architecture decision records. The agent is instructed to cite recalled sources rather than reproduce verbatim, to surface provenance uncertainty rather than hide it, and to choose surface forms (identifiers, helper decomposition, file layout) fresh for idiomatic Rust rather than copying from existing C and C++ reference implementations (CRlibm and fdlibm are behavioral oracles, not templates).

These are instructions to the agent, not guarantees about every line of output. A verbatim reproduction or an unflagged derivation could slip through. The project's defense against that is the instruction discipline above plus the human reviewer's ability to notice architectural smells that suggest a problem upstream, not a clean room audit. If you spot a passage that reads like a copy from a source it should not be copied from, please open an issue.

**Verification.** The verification places the central claim in an exhaustive enumeration of the `f32` input space: every `binary32` value of the unary surface is evaluated and compared against an independent reference oracle under all five rounding modes, so a wrong rounding on any single input is a test failure rather than a sampling gap; the inputs the oracle's own exponent range cannot certify are recorded inconclusive rather than passed. The `f64` surface rests on differential testing against the same reference oracle over a structured sample plus the published worst case rounding vectors as adversarial seeds. The shell commits a hardware float only when an enclosure of the true value determines it, so that the `BigFloat` to float rounding step cannot double round silently. CI runs the usual lints and the test suite on every change; the full exhaustive sweep runs out of band on a release cadence. Significant decisions are recorded as architecture decision records. The unary `f32` sweep has run; the crate remains under active construction, so this is not a claim that every planned function and surface is complete today.

**Scope.** pfloat-libm is a personal project. The intended consumer is the broader Rust scientific and embedded ecosystem: anyone who needs a correctly rounded `libm` without a C toolchain dependency. Durability and quality are goals, but this is not a funded library with a maintenance team behind it. The crate is at v0.1 and under active construction; the surface grows toward completeness in stated increments, and the API will break without notice before a 1.0 tag. pfloat-libm is not yet published to crates.io, and depends on pfloat, which is also not yet published. The repository remains public for users who want to read or follow the work.

**What this does not promise.** AI collaboration does not transfer responsibility. The author is accountable for what ships under his name. The disciplines above narrow the failure surface; they do not eliminate it. In particular, this process is most exposed to subtle bugs that a careful human reading of the code would catch but tests, types, and verification would not. For a correctly rounded `libm` that specifically includes roundings in the wrong direction at the `BigFloat` to hardware float step on pathological inputs the sweep did not reach, cancellation near the poles of the reciprocal kernels, and boundary cases in the `f32` subnormal range where the available mantissa shrinks below the working precision. Issues are welcome and will be triaged as time allows; no SLA is offered. This README describes the project's development process and is not a warranty; see the LICENSE file for the legal terms governing use.

## Status

v0.1, under active construction. The shell exposes the elementary
surface over pfloat's kernels; the unary `f32` surface is verified
exhaustively against an independent MPFR oracle under all five rounding
modes, with the per-function tiers recorded in `docs/kernel-list.md`.
The public API is unstable and will break without notice before a 1.0
tag. pfloat-libm is not published to crates.io.

## License

Dual-licensed under either:

- Apache License, Version 2.0 (`LICENSE-APACHE`)
- MIT License (`LICENSE-MIT`)

at your option.
