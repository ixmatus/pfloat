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

## Status

Pre-release (`0.x`); the API will change until `1.0`. Part of the pfloat
workspace; built and verified alongside it.

## License

MIT OR Apache-2.0.
