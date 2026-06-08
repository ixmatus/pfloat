# pfloat-ball

Rigorous arbitrary-precision real **ball** (midpoint-radius) arithmetic
in pure Rust, built on [pfloat](https://github.com/ixmatus/pfloat).

A ball `[m ± r]` denotes the closed real interval `[m − r, m + r]`. The
midpoint `m` is a full-precision correctly-rounded pfloat scalar; the
radius `r` is a small upward-rounded magnitude. Every operation computes
the midpoint with pfloat's verified kernel and bounds the radius by the
rounding error that kernel already produces, so the result is a *sound*
enclosure: the true mathematical result, applied to any point of the
input ball, lies inside the output ball.

This is the first cut of pfloat's rigorous-enclosure tower. v1.0 scopes
to real-ball arithmetic and the elementary functions over `BigFloat`;
special functions and the complex / IEEE-1788 faces are separate later
work.

## Soundness by construction

The radius type [`Mag`] makes an unsound radius unrepresentable: it has
no sign (no negative radius) and no NaN, and every `Mag` operation rounds
toward `+∞` (no inward-rounded radius). The remaining soundness
obligations live in the in-tree enclosure spec and the verification
lanes.

## Features

- `big` (default): `Ball<BigFloat>`, the headline dynamic-precision type.
- `fixed`: `Ball<FixedFloat<PREC>>`, compile-time precision.
- `std` (default): pfloat's thread-local sticky flags.
- `exp-log`, `trig`: the matching ball elementary functions.
- `serde`: `Serialize`/`Deserialize` for `Mag` and `Ball`.

A bare `--no-default-features` build exposes only `Mag` (no_std,
alloc-free).

## Usage

```rust
use pfloat::BigFloat;
use pfloat_ball::{Ball, Mag};

// A tight input interval: 1.0 plus or minus 2^-20.
let mid = BigFloat::try_from_i64_exact(1, 128).unwrap();
let x = Ball::new(mid, Mag::from_pow2(-20)).unwrap();

// Enclose sin over the whole interval. The result is a sound superset:
// sin(t) for every t in the input is provably inside [lower, upper].
let (s, _status) = x.sin();
println!("sin is enclosed by [{}, {}]", s.lower(), s.upper());
println!("certified accuracy: {} bits", s.rel_accuracy_bits());
```

`Ball::sin` needs the `trig` feature; the arithmetic surface (`add`, `sub`,
`mul`, `div`, `sqrt`, `cbrt`) is always available. Two runnable programs sit
under `examples/`: `sin_over_interval` (the enclosure contract) and
`ball_oracle` (a ball as a rigorous self-oracle for a `pfloat` scalar).

## How pfloat-ball is developed

This is an open disclosure of the development process so users can judge for themselves whether the resulting code meets their bar.

**Authorship and collaboration.** Parnell Springmeyer is the author of record. pfloat-ball is developed in collaboration with Claude, an AI coding agent from Anthropic. Parnell owns architecture, acceptance criteria, test and verification strategy, and release boundaries. Claude drafts the implementation, writes and runs tests and verification harnesses, and produces analysis under that direction. **Parnell does not review the generated code line by line.** Human oversight operates at the level of design, strategy, and outcomes: does the architecture make sense, are the right invariants being checked, does the verification strategy cover the risk surface, do the tests and proofs pass. Merges to main are GPG signed by Parnell to attest to that level of review, not to an audit of every line.

**Provenance.** Implementations derive from primary sources: the midpoint-radius interval arithmetic literature (Johansson's work on Arb is the principal design reference), IEEE 754-2019 for the scalar rounding semantics the radius is built on, and the Fundamental Theorem of Interval Arithmetic for the enclosure laws. The radius engine is built on pfloat's correctly rounded arbitrary precision scalars rather than on an FFI binding to a C interval library. The agent is instructed to cite recalled sources rather than reproduce verbatim, to surface provenance uncertainty rather than hide it, and to choose surface forms (identifiers, helper decomposition, file layout) fresh for idiomatic Rust rather than copying from existing C and C++ reference implementations.

These are instructions to the agent, not guarantees about every line of output. A verbatim reproduction or an unflagged derivation could slip through. The project's defense against that is the instruction discipline above plus the human reviewer's ability to notice architectural smells that suggest a problem upstream, not a clean room audit. If you spot a passage that reads like a copy from a source it should not be copied from, please open an issue.

**Verification.** The verification design places soundness in the type system where it can: the radius type makes a negative or not-a-number radius unrepresentable and rounds outward by construction, and a sealed scalar trait keeps the midpoint a correctly rounded value the crate's own surface cannot break. The enclosure laws are written down as an in-tree specification rather than left implicit. Above the types, a self-consistency property lane samples inside each ball and re-evaluates at higher precision to check that the radius covers the kernel's own residual, formal proof harnesses (Kani) discharge the allocation-free radius invariants, and the decimal parser entry point is fuzzed as an adversarial boundary. An independent containment lane samples ball results against an established arbitrary precision interval library, reached out of process so no interval library enters the link graph, for the stronger claim that an enclosure is sound rather than merely self-consistent; it runs per release. The crate forbids `unsafe`. CI runs the usual lints and the test and verification suite; specific harness counts change as the project evolves. Significant decisions are recorded as ADRs in the repo.

**Scope.** pfloat-ball is a personal project and the first cut of a rigorous enclosure tower over pfloat. The intended consumers are the broader Rust scientific and embedded ecosystem (anyone who needs a sound real enclosure at arbitrary precision without a C interval library) and pfloat itself, which can use the ball arithmetic as a rigorous oracle for its own scalar kernels. Durability and quality are goals, but this is not a funded library with a maintenance team behind it. The crate is at an early version: real ball arithmetic and the elementary functions are in place, while special functions and the complex and IEEE 1788 interval faces are planned as later work, and it is not yet published to crates.io. The repository remains public for users who want to read or follow the work.

**What this does not promise.** AI collaboration does not transfer responsibility. The author is accountable for what ships under his name. The disciplines above narrow the failure surface; they do not eliminate it. In particular, this process is most exposed to subtle bugs that a careful human reading of the code would catch but tests, types, and formal verification would not. For a rigorous enclosure library that specifically includes a radius that under-estimates the true error on inputs no test happened to cover: a directed rounding at a domain or saturation boundary that lands a hair too small, a propagation bound that holds across the sampled inputs but under-covers at an unsampled corner of the input box, or a special-case value that drifts outside its enclosure while the bulk of the surface stays sound. Issues are welcome and will be triaged as time allows; no SLA is offered. This README describes the project's development process and is not a warranty; see the LICENSE file for the legal terms governing use.

## Status

Pre-release (`0.x`); the API will change until `1.0`. Part of the pfloat
workspace; built and verified alongside it.

## License

MIT OR Apache-2.0.
