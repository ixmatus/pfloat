# pfloat

Pure Rust correctly-rounded arbitrary-precision floats.

## How pfloat is developed

This is an open disclosure of the development process so users can judge for themselves whether the resulting code meets their bar.

**Authorship and collaboration.** Parnell Springmeyer is the author of record. pfloat is developed in collaboration with Claude, an AI coding agent from Anthropic. Parnell owns architecture,
acceptance criteria, test and verification strategy, and release boundaries. Claude drafts the implementation, writes and runs tests and verification harnesses, and produces analysis under that
direction. **Parnell does not review the generated code line by line.** Human oversight operates at the level of design, strategy, and outcomes: does the architecture make sense, are the right
invariants being checked, does the verification strategy cover the risk surface, do the tests and proofs pass. Merges to main are GPG signed by Parnell to attest to that level of review, not to
an audit of every line.

**Provenance.** Implementations derive from primary sources: IEEE 754-2019 for floating point semantics and rounding modes, the published Lefèvre and Muller worst case rounding tables and the
table maker's dilemma literature for correctly rounded transcendentals, and published open algorithm work with stated licenses for the special functions (gamma, erf, Bessel, zeta, and related).
The agent is instructed to cite recalled sources rather than reproduce verbatim, to surface provenance uncertainty rather than hide it, and to choose surface forms (identifiers, helper
decomposition, file layout) fresh for idiomatic Rust rather than copying from existing C and C++ reference implementations (MPFR is the principal behavioral oracle).

These are instructions to the agent, not guarantees about every line of output. A verbatim reproduction or an unflagged derivation could slip through. The project's defense against that is the instruction discipline above plus the human reviewer's ability to notice architectural smells that suggest a problem upstream, not a clean room audit. If you spot a passage that reads like a copy from a source it should not be copied from, please open an issue.

**Verification.** The verification design places correctness in the type system where it can (precision and rounding mode as types, compile-time precision parameterization that lets `no_std`
consumers escape `alloc`, IEEE 754 exception semantics surfaced as observable state), in formal proof harnesses (Kani) where the cost is justified, in IEEE 754-2019 conformance vectors and the
worst case rounding tables run as integration tests, in differential tests against a trusted reference oracle on a separate CI lane, and in fuzz coverage of parser entry points. CI runs the usual lints and the full test and verification suite; specific harness counts and conformance counts change as the project evolves. Significant decisions are recorded as ADRs in the repo. `unsafe` blocks carry a written justification at the call site.

**Scope.** pfloat is a personal project. The intended consumer is the broader Rust scientific and embedded ecosystem (anyone who needs more than `f64` with correctly rounded results and no C
toolchain dependency); durability and quality are goals, but this is not a funded library with a maintenance team behind it. The crate has not reached 1.0 and is unpublished; the design, the
architecture decision records, and the CI scaffolding are in place, but the algorithmic kernels are not yet implemented. The repository remains public for users who want to read or follow the
work.

**What this does not promise.** AI collaboration does not transfer responsibility. The author is accountable for what ships under his name. The disciplines above narrow the failure surface; they
do not eliminate it. In particular, this process is most exposed to subtle bugs that a careful human reading of the code would catch but tests, types, and formal verification would not. For
numerical code that specifically includes roundings in the wrong direction on pathological inputs, mantissa shifts that lose a bit, and special function values that drift outside the claimed
error bound on inputs no test happened to cover. Issues are welcome and will be triaged as time allows; no SLA is offered. This README describes the project's development process and is not a
warranty; see the LICENSE file for the legal terms governing use.

## Status

Pre-1.0. The repository carries the design (`DESIGN.md`), the
architecture decision records (`docs/decisions/`), the CI
scaffolding, and the algorithmic kernels (arithmetic, the elementary
transcendental and special-function surface listed below, both
precision profiles). The Phase 1 correctness sweep is complete per
ADR-0033 (exhaustive f32 audit closed by slices p1.1 through p1.11
and the follow-ups pf-jn1y, pf-cvs, pf-06sw). The public API is
unstable and will break without notice until 1.0; slice 8c is the
v1.0 tag ceremony.

## Quickstart

```rust
use pfloat::{BigFloat, FixedFloat, RoundingMode};

// BigFloat: runtime precision, heap-allocated mantissa.
// Square root of two at 200-bit precision, correctly rounded to
// nearest even. The call returns the result and a Status carrying
// any IEEE 754-2019 sticky exception flags raised by the operation.
let two = BigFloat::try_from_i64_exact(2, 200).unwrap();
let (sqrt2, _status) = two.sqrt(RoundingMode::NearestEven);

// FixedFloat: precision fixed at compile time via a const generic.
// Stack-allocated; works without `alloc`. `FixedFloat<113>` is the
// binary128 mantissa width.
type F128 = FixedFloat<113>;
let two_f128 = F128::try_from_i64_exact(2).unwrap();
let (sqrt2_f128, _) = two_f128.sqrt(RoundingMode::NearestEven);
```

Five IEEE 754-2019 rounding modes are available
(`RoundingMode::NearestEven`, `NearestAway`, `TowardZero`,
`TowardPositive`, `TowardNegative`); every kernel returns a
`(value, Status)` pair so callers can inspect or accumulate sticky
flags without thread-local state. Under the `std` feature, the same
flags also accumulate into a thread-local set accessible via
`pfloat::flags`.

## Scope target

v1.0 covers the MPFR-equivalent surface:

- IEEE 754-2019 arithmetic with all five rounding modes (RNE, RNA, RZ, RP, RM) and sticky exception flags.
- Correctly-rounded elementary transcendentals: `exp`, `log` family, trig and inverses, hyperbolic and inverses, `pow`.
- Special functions: `gamma`, `lgamma`, `digamma`, `beta`, `erf`, `erfc`, Bessel `J/Y/I/K`, `zeta`, `Ei`, `Si`, `Ci`, Airy `Ai/Bi/Ai′/Bi′`, AGM.
- Two precision profiles in one crate: `BigFloat` (runtime precision, needs `alloc`) and `FixedFloat<const PREC: u32>` (compile-time precision, stack-allocated, runs without `alloc`).
- `no_std`-first, embedded-friendly. CI cross-compiles to `thumbv6m-none-eabi`.

## Installation

pfloat is pre-1.0 and unpublished; depend on it as a git dependency
until the v1.0 tag and crates.io release land:

```toml
[dependencies]
pfloat = { git = "https://github.com/ixmatus/pfloat" }
```

The crate requires a nightly Rust toolchain for
`feature(generic_const_exprs)` (the bit-level const-generic mantissa
storage that lets `FixedFloat<const PREC: u32>` avoid `alloc` per
ADR-0011). The pfloat repository pins its toolchain channel in
`rust-toolchain.toml`; downstream consumers need a matching nightly
in their own workspace.

For embedded targets, build with `--no-default-features
--features=fixed` (or `--features=alloc,big` if the runtime
precision profile is wanted). CI exercises the
`thumbv6m-none-eabi` and `thumbv8m.main-none-eabi` cross targets.

## Feature flags

The crate ships in a small number of feature clusters so embedded
and `alloc`-free consumers can pick only what they need.

| Cluster | Feature | Implies | Adds |
|---|---|---|---|
| Defaults | `std`, `fmt`, `big` | `alloc` | `std::error::Error` impls, thread-local sticky flags, `core::fmt::Write` parse and format, `BigFloat` |
| Storage | `alloc` | | Heap-allocated multi-limb mantissa buffers |
| | `fixed` | `big` | `FixedFloat<const PREC: u32>` (stack-allocated) |
| Operator sugar | `ops` | `big` | `core::ops::Add`/`Sub`/`Mul`/`Div`/`Neg` overloads (default `NearestEven`, discards `Status`) |
| Elementary | `exp-log` | `big`, `agm` | `exp`, `expm1`, `exp2`, `exp10`, `ln`, `log1p`, `log2`, `log10` |
| | `trig` | `exp-log` | `sin`, `cos`, `tan`, `asin`, `acos`, `atan`, `atan2` |
| | `agm` | `big` | Arithmetic-geometric mean |
| Tier-1 specials | `specials` | `trig` | `erf`, `erfc`, `gamma`, `lgamma`, `digamma`, `beta` |
| Tier-2 specials | `integrals` | `specials` | `Ei`, `Si`, `Ci`, `li` |
| | `airy` | `specials` | `Ai`, `Bi`, `Ai′`, `Bi′` |
| | `bessel` | `specials` | `J₀`, `J₁`, `Jₙ` (and `Y/I/K` siblings) |
| | `zeta` | `specials` | Riemann ζ on the real axis |
| Verification (dev) | `kani` | | Compile Kani proof harnesses (off in normal builds) |
| | `differential-mpfr` | full feature union | Differential testing against MPFR via `rug` (Unix-only) |
| | `differential-arb` | `differential-mpfr` | Arb oracle backend via Python subprocess |

## Why

`rug` and `gmp-mpfr-sys` force a C toolchain on every Rust project that needs more than `f64` with correct rounding. `astro-float` is the closest pure-Rust alternative and covers basic arithmetic plus elementary transcendentals well; the special-function surface (gamma, erf, Bessel, zeta, etc.) and shipped formal-verification artifacts are gaps that pfloat fills directly. The companion goal is to displace the GMP/MPFR build dependency for scientific, financial, and symbolic-computation crates that want WebAssembly or embedded targets.

## Verification posture

- IEEE 754-2019 conformance vectors and Lefèvre–Muller worst-case-rounding tables run as integration tests.
- Kani harnesses discharge no-panic, rounding-direction, and sign-of-zero properties on the arithmetic core.
- `gmp-mpfr-sys` runs as a feature-gated dev-dependency on a separate Linux CI lane for primary differential testing against MPFR. The default lane stays pure Rust.
- Phase 1's exhaustive `f32` audit (ADR-0035) cross-checks three independent oracles. Arb via `python-flint` (LGPL, subprocess-only) is the primary oracle for the FnIds MPFR cannot cover (`Si`, `Ci`, `li`, `Bi`, `Ai′`, `Bi′`, the Bessel `I`/`K` family). mpmath (BSD, pure Python) cross-checks Arb's certified outputs in the same Python venv. Maxima (GPL, invoked through `nix-shell`) supplies a sampling third opinion on a small pinned corpus of worker outputs. The per-push CI gate stays MPFR-only and Python-free; the per-slice gate exercises the Arb worker plus a per-push diff against the pinned-corpus snapshot; the per-release sweep exercises Arb + mpmath agreement at full f32 coverage plus three-way agreement on the pinned corpus.
- `cargo-fuzz` covers the parser entry points.

## Conformance evidence

Per-bucket counts of pfloat's verification harnesses. Each bucket
is asserted independently by `scripts/conformance-counts.sh` and
gated in CI; one bucket shrinking while another grows cannot hide
under an aggregate floor.

- **Kani proof harnesses:** 354 `#[kani::proof]`
  attributes across 61 files in `src/verify/`.
- **Differential lanes:** 34 `tests/differential_*.rs`
  files. CI sweep 10⁴ inputs per (op × precision × rounding
  mode); `PFLOAT_DEEP=1` escalates to 10⁶ (ADR-0014).
- **Fuzz targets:** 7 targets under `fuzz/fuzz_targets/`.
- **Property tests:** 33 `tests/property_*.rs` files
  (addsub agm ai bi ci classify digamma_beta div ei erf exp exp_family fma fmt gamma hyperbolic ik jn li ln log_family mul parse partial_cmp pow rounding si sqrt total_cmp trig trig_inverse yn zeta).

## Documentation

The in-tree documents are the primary reference; the published
rustdoc is a thin layer on top of them.

- `DESIGN.md` — full architectural design (numeric representation, arithmetic algorithms, transcendental and special-function strategy, verification stack, feature gating, phase plan).
- `docs/ROADMAP.md` — current direction, Phase 1/2 sequencing, and the post-ADR-0035 follow-up narrative.
- `docs/decisions/` — architecture decision records (37 ADRs covering load-bearing choices, with provenance and consequences).
- `tools/alloc-profile/` — standalone allocation profiler crate (path-dep on pfloat) used for ADR-0028 and ADR-0037 measurements.
- `benches/mul_thresholds.rs` — criterion bench backing the Karatsuba threshold calibration in ADR-0027.

To generate the rustdoc locally:

```sh
cargo +nightly doc \
  --features=std,fmt,big,fixed,ops,exp-log,trig,specials,agm,integrals,airy,bessel,zeta \
  --no-deps --open
```

## Issues

Issue reports are welcome. See the disclosure at the top of this
README for the project's maintenance posture (personal project, no
SLA, triage as time allows). A useful issue carries the smallest
input that exhibits the problem, the rounding mode, the precision,
the expected and observed value (or status), and the pfloat commit
the reproducer was run against. Reproducers that fit in a single
property test or a `cargo test` invocation are the easiest to act
on.

If a passage in the source looks like a verbatim copy from a
license-encumbered reference implementation (the provenance
discipline in the disclosure above lists what the agent is
instructed to avoid; verbatim reproduction is the failure mode it
guards against), please flag it in an issue. The honest framing is
"this looks like it might be lifted from X" rather than waiting for
certainty; the project would rather investigate a false positive
than miss a real one.

## License

Dual-licensed under either:

- Apache License, Version 2.0 (`LICENSE-APACHE`)
- MIT License (`LICENSE-MIT`)

at your option.
