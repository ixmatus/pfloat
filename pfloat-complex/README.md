# pfloat-complex

Componentwise correctly-rounded complex arithmetic in pure Rust, built on
[`pfloat`](..). A `Complex<T>` is a pair `(re, im)` of pfloat scalars
([`BigFloat`] or [`FixedFloat<PREC>`]); each operation rounds the real and
imaginary parts each correctly under their own real rounding mode, the model
MPC uses and the only coherent strong rounding claim for complex numbers
(which carry no total order).

This is the MPC analog in the pfloat family: arbitrary-precision complex
arithmetic with no C toolchain, where `num-complex` is a bare container and
`rug`/MPC require FFI. The component type is constrained by a sealed
`RealScalar` trait, so a `Complex` here is always built over a verified,
correctly-rounded pfloat scalar.

## Status

Pre-release (`0.x`); the API will change until `1.0`. This cut ships the
`Complex<T>` type and componentwise additive arithmetic (`add`, `sub`,
`neg`, `conj`). Multiplication and division (the cancellation-safe fused
two-product forms), magnitude and phase, and the elementary functions with
their C99 Annex G branch cuts are later slices. Part of the pfloat
workspace; built and verified alongside it. The full development-process
disclosure lands with the 1.0 cut, as it did for `pfloat` and `pfloat-ball`.

## License

MIT OR Apache-2.0, at your option.
