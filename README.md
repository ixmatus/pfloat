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
decomposition, file layout) fresh for idiomatic Rust rather than copying from existing C and C++ reference implementations (MPFR is the principal behavioral oracle) that serve as oracles for
behavior.

These are instructions to the agent, not guarantees about every line of output. A verbatim reproduction or an unflagged derivation could slip through. The project's defense against that is the
instruction discipline above plus the human reviewer's ability to notice architectural smells that suggest a problem upstream, not a clean room audit. If you spot a passage that reads like a copy   from a source it should not be copied from, please open an issue.

**Verification.** Correctness lives in the type system where it can (precision and rounding mode as types, the `FixedFloat<const PREC: u32>` parameterization that lets `no_std` consumers escape
`alloc`, sticky exception flags as observable state), in formal proof harnesses (Kani) where the cost is justified, in IEEE 754-2019 conformance vectors and the worst case rounding tables run as
integration tests, in differential tests against MPFR on a separate CI lane behind a feature flag, and in fuzz coverage of every parser entry. CI runs the usual lints and the full test and
verification suite; specific harness counts and conformance counts change as the project evolves. Significant decisions are recorded in ADRs under `docs/decisions/`. `unsafe` blocks carry a
written justification at the call site.

**Scope.** pfloat is a personal project. The intended consumer is the broader Rust scientific and embedded ecosystem (anyone who needs more than `f64` with correctly rounded results and no C
toolchain dependency); durability and quality are goals, but this is not a funded library with a maintenance team behind it. The crate has not reached 1.0 and is unpublished; the design
(`DESIGN.md`), the architecture decision records under `docs/decisions/`, and the CI scaffolding are in place, but the algorithmic kernels are not yet implemented. The repository remains public
for users who want to read or follow the work.

**What this does not promise.** AI collaboration does not transfer responsibility. The author is accountable for what ships under his name. The disciplines above narrow the failure surface; they
do not eliminate it. In particular, this process is most exposed to subtle bugs that a careful human reading of the code would catch but tests, types, and formal verification would not. For
numerical code that specifically includes roundings in the wrong direction on pathological inputs, mantissa shifts that lose a bit, and special function values that drift outside the claimed
error bound on inputs no test happened to cover. Issues are welcome and will be triaged as time allows; no SLA is offered.


## Status

Pre-1.0. The repository carries the design (`DESIGN.md`), the
architecture decision records (`docs/decisions/`), and the CI
scaffolding. The algorithmic kernels are not yet implemented. The
public API is unstable and will break without notice until 1.0.

## Scope target

v1.0 covers the MPFR-equivalent surface:

- IEEE 754-2019 arithmetic with all five rounding modes (RNE, RNA, RZ, RP, RM) and sticky exception flags.
- Correctly-rounded elementary transcendentals: `exp`, `log` family, trig and inverses, hyperbolic and inverses, `pow`.
- Special functions: `gamma`, `lgamma`, `digamma`, `beta`, `erf`, `erfc`, Bessel `J/Y/I/K`, `zeta`, `Ei`, `Si`, `Ci`, Airy, AGM.
- Two precision profiles in one crate: `BigFloat` (runtime precision, needs `alloc`) and `FixedFloat<const PREC: u32>` (compile-time precision, stack-allocated, runs without `alloc`).
- `no_std`-first, embedded-friendly. CI cross-compiles to `thumbv6m-none-eabi`.

## Why

`rug` and `gmp-mpfr-sys` force a C toolchain on every Rust project that needs more than `f64` with correct rounding. `astro-float` is the closest pure-Rust alternative and covers basic arithmetic plus elementary transcendentals well; the special-function surface (gamma, erf, Bessel, zeta, etc.) and shipped formal-verification artifacts are gaps that pfloat fills directly. The companion goal is to displace the GMP/MPFR build dependency for scientific, financial, and symbolic-computation crates that want WebAssembly or embedded targets.

## Verification posture

- IEEE 754-2019 conformance vectors and Lefèvre–Muller worst-case-rounding tables run as integration tests.
- Kani harnesses discharge no-panic, rounding-direction, and sign-of-zero properties on the arithmetic core.
- `gmp-mpfr-sys` runs as a feature-gated dev-dependency on a separate Linux CI lane for differential testing. The default lane stays pure Rust.
- `cargo-fuzz` covers every parser entry.

## License

Dual-licensed under either:

- Apache License, Version 2.0 (`LICENSE-APACHE`)
- MIT License (`LICENSE-MIT`)

at your option.
