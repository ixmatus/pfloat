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
permissively licensed, and its unary `f32` surface is designed to be
exhaustively verified: every one of the 2^32 `binary32` inputs checked against
an independent oracle, a claim no competing pure-Rust libm makes. The `f64`
surface rests on differential testing against a reference oracle plus the
published worst-case hard-to-round vectors, since the 2^64 input space cannot be
enumerated.

pfloat-libm is `no_std` + `alloc`. There is no alloc-free profile: correct
rounding grows the working precision at runtime past any compile-time width, so
the computation allocates.

For the development process that produced this code, read the disclosure
immediately below before deciding whether to adopt.

## How pfloat-libm is developed

This is an open disclosure of the development process so users can judge for themselves whether the resulting code meets their bar.

**Authorship and collaboration.** Parnell Springmeyer is the author of record. pfloat-libm is developed in collaboration with Claude, an AI coding agent from Anthropic. Parnell owns architecture, acceptance criteria, test and verification strategy, and release boundaries. Claude drafts the implementation, writes and runs tests and verification harnesses, and produces analysis under that direction. **Parnell does not review the generated code line by line.** Human oversight operates at the level of design, strategy, and outcomes: does the architecture make sense, are the right invariants being checked, does the verification strategy cover the risk surface, do the tests and proofs pass. Merges to main are GPG signed by Parnell to attest to that level of review, not to an audit of every line.

**Provenance.** Implementations derive from primary sources: IEEE 754-2019 for floating point semantics and rounding modes, the DLMF for the reciprocal and root function definitions and their range reduction, and the published Lefèvre and Muller worst case rounding tables for the hard to round vectors. pfloat-libm is mostly a shell and a verification harness over pfloat; the arbitrary precision mathematics, including the direct reciprocal and root kernels, lives in pfloat and carries its provenance in pfloat's own architecture decision records. The agent is instructed to cite recalled sources rather than reproduce verbatim, to surface provenance uncertainty rather than hide it, and to choose surface forms (identifiers, helper decomposition, file layout) fresh for idiomatic Rust rather than copying from existing C and C++ reference implementations (CRlibm and fdlibm are behavioral oracles, not templates).

These are instructions to the agent, not guarantees about every line of output. A verbatim reproduction or an unflagged derivation could slip through. The project's defense against that is the instruction discipline above plus the human reviewer's ability to notice architectural smells that suggest a problem upstream, not a clean room audit. If you spot a passage that reads like a copy from a source it should not be copied from, please open an issue.

**Verification.** The verification design places the central claim in an exhaustive enumeration of the `f32` input space: every `binary32` value is intended to be evaluated and compared against an independent reference oracle, so a wrong rounding on any single input is a test failure rather than a sampling gap. The `f64` surface is designed to rest on differential testing against the same reference oracle over a structured sample plus the published worst case rounding vectors as adversarial seeds. The shell is designed to commit a hardware float only when an enclosure of the true value determines it, so that the `BigFloat` to float rounding step cannot double round silently. CI is intended to run the usual lints and the test suite on every change, with the full exhaustive sweep on a release cadence. Significant decisions are recorded as architecture decision records. This is a design statement for a crate under active construction, not a claim that the full verification has run on every function today.

**Scope.** pfloat-libm is a personal project. The intended consumer is the broader Rust scientific and embedded ecosystem: anyone who needs a correctly rounded `libm` without a C toolchain dependency. Durability and quality are goals, but this is not a funded library with a maintenance team behind it. The crate is at v0.1 and under active construction; the surface grows toward completeness in stated increments, and the API will break without notice before a 1.0 tag. pfloat-libm is not yet published to crates.io, and depends on pfloat, which is also not yet published. The repository remains public for users who want to read or follow the work.

**What this does not promise.** AI collaboration does not transfer responsibility. The author is accountable for what ships under his name. The disciplines above narrow the failure surface; they do not eliminate it. In particular, this process is most exposed to subtle bugs that a careful human reading of the code would catch but tests, types, and verification would not. For a correctly rounded `libm` that specifically includes roundings in the wrong direction at the `BigFloat` to hardware float step on pathological inputs the sweep did not reach, cancellation near the poles of the reciprocal kernels, and boundary cases in the `f32` subnormal range where the available mantissa shrinks below the working precision. Issues are welcome and will be triaged as time allows; no SLA is offered. This README describes the project's development process and is not a warranty; see the LICENSE file for the legal terms governing use.

## Status

v0.1, under active construction. This is an early scaffold: the repository
carries the kernel list (`docs/kernel-list.md`) recording the planned surface
and its verification tiers, the dual license, and the shell crate skeleton. The
public API is unstable and will break without notice before a 1.0 tag. pfloat-libm
is not published to crates.io.

## License

Dual-licensed under either:

- Apache License, Version 2.0 (`LICENSE-APACHE`)
- MIT License (`LICENSE-MIT`)

at your option.
