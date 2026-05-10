# pfloat

Pure Rust correctly-rounded arbitrary-precision floats.

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
