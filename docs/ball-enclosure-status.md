# pfloat-ball enclosure-accuracy posture

Per-operation enclosure shape and accuracy posture for the pfloat-ball
real-ball surface. This is the ball analog of `docs/rounding-status.md`,
but it declares a different kind of property. The scalar table reports a
*correct-rounding* verdict per function (a point claim). A ball returns an
*enclosure* of the true result over the whole input ball, so the property
that matters is two-layered:

- **Soundness (the hard guarantee).** Every operation returns a ball that
  contains the true result of the real function applied to every point of
  the input ball (the Fundamental Theorem, Law 1 in `src/spec.rs`). The
  radius may over-estimate; it must never under-estimate. This is not on a
  spectrum: it holds for every operation in the table, backed by the
  blocking FTIA self-consistency property lane and the per-release
  independent Arb containment plus BRACKETI range-soundness lane
  (ADR-0078, ADR-0082).
- **Accuracy (the posture below, a quality property).** *How tight* the
  sound enclosure is. Tightness is **measured, not asserted**: the
  `differential_arb` tightness lane logs the ball width against Arb's
  rigorous image per `(function, precision, magnitude)` bucket
  (`tests/arb_tightness_expected.txt`), regression-guarded per bucket
  rather than on an aggregate floor. The posture column records the
  *expected* tightness class that follows from each operation's enclosure
  shape; the measured lane is the source of truth that a regression would
  trip. ADR-0087.

Legend (enclosure shapes, `src/elem.rs` / `src/arith.rs`):

- **directed-pair.** The midpoint is the correctly-rounded scalar kernel
  result; the radius bounds half the directed `(TowardNegative,
  TowardPositive)` spread (the kernel's own residual), rounded up under
  `Mag`, plus the propagated input width (Law 2).
- **monotone endpoints.** For a function monotone on the input interval,
  the enclosure is the directed-kernel-correct endpoints `[f(alo), f(ahi)]`
  (or reversed for a decreasing function). The tightest representable
  interval for the monotone range.
- **1-Lipschitz.** For `|f'| <= 1` (sin, cos, and hypot per argument) the
  radius is the midpoint kernel residual plus the input radius, since the
  function's variation across the input ball is at most the input radius.
- **composed.** Built from sub-operations through ball arithmetic (tan =
  sin/cos via ball division; atan2 via the gradient-magnitude bound; cosh
  via a magnitude interval through the monotone route).

Legend (accuracy posture):

- **tightest.** The enclosure is the narrowest representable sound ball for
  the operation's own rounding: a correctly-rounded midpoint (directed
  pair) or directed-kernel-correct monotone endpoints. Any residual width
  beyond this is the inherent dependency/input-propagation term, not a
  rounding loss.
- **accurate.** Sound and tight in the near-linear or monotone region, but
  the shape's variation bound conservatively over-covers near a local
  extremum (where `|f'| -> 0` but the radius still carries the full input
  radius). The exact straddles where Arb's image is tightest are exactly
  where the Lipschitz/composed radius is loosest; this is the
  `|f'| -> 0` effect ADR-0082 records, and the reason tightness is measured.

| Operations | Enclosure shape | Posture |
| --- | --- | --- |
| `add` `sub` `mul` `div` `sqrt` `cbrt` | directed-pair | tightest |
| `exp` `expm1` `exp2` `exp10` `ln` `log2` `log10` `log1p` | monotone endpoints | tightest |
| `sinh` `tanh` `asinh` `acosh` `atanh` | monotone endpoints | tightest |
| `asin` `acos` `atan` | monotone endpoints | tightest |
| `sin` `cos` | 1-Lipschitz | accurate |
| `cosh` | composed (magnitude, then monotone) | accurate |
| `hypot` | 1-Lipschitz per argument | accurate |
| `tan` `atan2` | composed | accurate |

The completeness of this table against the public surface is gated in CI
by `scripts/ball-enclosure-status.sh --check`: every public ball operation
(`pub fn` in `src/arith.rs` and `src/elem.rs`) must appear here, so a new
operation cannot ship without a declared posture. The shape and posture
cells themselves are authored, not generated: they are a design
declaration that the measured tightness lane substantiates and would
regress against.

Special functions (gamma, zeta, Bessel, Airy, the integrals) are not part
of this surface; their ball lift is later work (the oscillatory and
pole-bearing families need per-family variation bounds the four shapes
above do not provide), and they will extend this table when they land.
