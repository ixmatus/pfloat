# pfloat-complex

Componentwise correctly-rounded complex arithmetic in pure Rust, built on
[pfloat](https://github.com/ixmatus/pfloat). A `Complex<T>` is a pair `(re, im)`
of pfloat scalars ([`BigFloat`] or [`FixedFloat<PREC>`]); each operation rounds
the real and imaginary parts each correctly under their own real rounding mode,
the model MPC uses and the only coherent strong rounding claim for complex
numbers (which carry no total order).

This is the MPC analog in the pfloat family: arbitrary-precision complex
arithmetic with no C toolchain, where `num-complex` is a bare container and
`rug`/MPC require FFI.

## Componentwise correct rounding

The component type is constrained by a sealed `RealScalar` trait, implemented
only for `BigFloat` and `FixedFloat<PREC>`, so a `Complex` here is always built
over a verified, correctly rounded pfloat scalar; the crate's own surface cannot
be made to construct one over an unverified type.

On top of that componentwise rounding, branch selection and signed-zero
discrimination follow the C99/C11 Annex G convention. The branch cuts are a
semantic layer above rounding: `csqrt(-4 + 0i) = +2i` while `csqrt(-4 - 0i) =
-2i`, the sign of the input zero choosing the half of the cut. The §G.5.1
complex-infinity rules are honored too, so `(1 + 0i) / (0 + 0i)` is a directed
complex infinity rather than a NaN.

## Features

- `big` (default): `Complex<BigFloat>`, the headline dynamic-precision type.
- `fixed`: `Complex<FixedFloat<PREC>>`, compile-time precision.
- `std` (default): pfloat's thread-local sticky flags.
- `exp-log`: the magnitude `abs` and the square root `sqrt`.
- `trig` (implies `exp-log`): the phase `arg`, `to_polar`, and the complex
  `exp` and `log`.

A `Complex` is a pair of pfloat scalars, and the scalar engine needs `alloc`, so
there is no public item below `big`: a bare `--no-default-features` build is
empty.

## Usage

```rust
use pfloat::{BigFloat, RoundingMode};
use pfloat_complex::Complex;

let re = BigFloat::try_from_i64_exact(-7, 128).unwrap();
let im = BigFloat::try_from_i64_exact(24, 128).unwrap();
let z = Complex::new(re, im);

// Principal square root with the C99 Annex G branch cut, componentwise
// correctly rounded. csqrt(-7 + 24i) = 3 + 4i, exact (an OK status).
let (w, status) = z.sqrt(RoundingMode::NearestEven);
assert!(status.is_ok());
println!("sqrt = {} + {}i", w.re, w.im);
```

`Complex::sqrt` needs the `exp-log` feature, and `exp`/`log`/`arg` need `trig`;
the arithmetic surface (`add`, `sub`, `mul`, `div`, `neg`, `conj`, `norm_sqr`)
is always available under `big`.

## How pfloat-complex is developed

This is an open disclosure of the development process so users can judge for themselves whether the resulting code meets their bar.

**Authorship and collaboration.** Parnell Springmeyer is the author of record. pfloat-complex is developed in collaboration with Claude, an AI coding agent from Anthropic. Parnell owns architecture, acceptance criteria, test and verification strategy, and release boundaries. Claude drafts the implementation, writes and runs tests and verification harnesses, and produces analysis under that direction. **Parnell does not review the generated code line by line.** Human oversight operates at the level of design, strategy, and outcomes: does the architecture make sense, are the right invariants being checked, does the verification strategy cover the risk surface, do the tests and proofs pass. Merges to main are GPG signed by Parnell to attest to that level of review, not to an audit of every line.

**Provenance.** Implementations derive from primary sources: the C99/C11 Annex G specification for the complex branch cuts and the §G.5.1 infinity recovery, IEEE 754-2019 §9.2.1 for the signed-zero `atan2` and `hypot` tables the magnitude and phase ride, the componentwise correct rounding model (MPC is the design reference for the model, not a template for code), and Kahan's cancellation robust reformulation for the square root. The complex arithmetic is built on pfloat's correctly rounded arbitrary precision scalars rather than on an FFI binding to MPC. The agent is instructed to cite recalled sources rather than reproduce verbatim, to surface provenance uncertainty rather than hide it, and to choose surface forms (identifiers, helper decomposition, file layout) fresh for idiomatic Rust rather than copying from existing C and C++ reference implementations.

These are instructions to the agent, not guarantees about every line of output. A verbatim reproduction or an unflagged derivation could slip through. The project's defense against that is the instruction discipline above plus the human reviewer's ability to notice architectural smells that suggest a problem upstream, not a clean room audit. If you spot a passage that reads like a copy from a source it should not be copied from, please open an issue.

**Verification.** The verification design places what it can in the type system: a sealed scalar trait keeps each component a correctly rounded pfloat value the crate's own surface cannot break, so a `Complex` here is always built over a verified scalar. Above the types, the C99 Annex G special-value rows are enumerated through the public API across precisions and all five rounding modes, the primary branch-cut and signed-zero guard; the special-value dispatch is checked exhaustively over the finite class grid; algebraic identities cross-tie the elementary kernels to one another; and a formal proof harness (Kani) discharges the componentwise status-merge contract. An independent lane checks each result bit for bit against an established arbitrary precision complex interval library, reached out of process so no complex library enters the link graph, for the stronger claim that the rounding is correct rather than merely self-consistent; it runs per release. The crate forbids `unsafe`. CI runs the usual lints and the test and verification suite; specific harness counts change as the project evolves. Significant decisions are recorded as ADRs in the repo.

**Scope.** pfloat-complex is a personal project and the complex face of a rigorous numeric tower over pfloat. The intended consumers are the broader Rust scientific and embedded ecosystem (anyone who needs arbitrary precision complex arithmetic without a C toolchain) and the pfloat family itself. Durability and quality are goals, but this is not a funded library with a maintenance team behind it. The crate is at 1.0: componentwise correctly rounded complex arithmetic and the elementary core (`sqrt`, `exp`, `log`) with their C99 Annex G branch cuts are in place and the public API is meant to be settled from here, while the trigonometric, hyperbolic, and inverse functions, `pow` / `cis` / `from_polar`, and the `ComplexBall` join with pfloat-ball are planned as later work. It is not yet published to crates.io; the signed 1.0 tag is the milestone. The repository remains public for users who want to read or follow the work.

**What this does not promise.** AI collaboration does not transfer responsibility. The author is accountable for what ships under his name. The disciplines above narrow the failure surface; they do not eliminate it. In particular, this process is most exposed to subtle bugs that a careful human reading of the code would catch but tests, types, and formal verification would not. For componentwise complex arithmetic that specifically includes the wrong half of a branch cut returned for an unsigned zero on inputs no signed-zero test happened to cover (the sign of a result zero is input determined and load bearing: `csqrt(-4 + 0i)` and `csqrt(-4 - 0i)` differ only in that sign), catastrophic cancellation in the `ac - bd` and `ad + bc` cross products on inputs no random sweep lands on, and a near-zero component of `log` near the unit circle that a bounded enclosure rounds a hair wrong. Issues are welcome and will be triaged as time allows; no SLA is offered. This README describes the project's development process and is not a warranty; see the LICENSE file for the legal terms governing use.

## Status

Version `1.0`: the public API is settled, and the aim is to keep `1.x` changes
additive (semver); as a personal project this is an intent rather than a
contractual guarantee (ADR-0093). Part of the pfloat workspace; built and
verified alongside it.

## License

MIT OR Apache-2.0, at your option.
