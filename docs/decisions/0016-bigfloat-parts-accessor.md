# ADR-0016: Public `BigFloat::parts` accessor and bit-exact MPFR converter

Status: accepted (slice 7a shipped)

## Context

ADR-0014's slice 6h status update flagged a 1-ULP divergence
between pfloat and MPFR under non-NearestEven rounding for a
handful of operations (`div(-966132233652331, 1233101814760529)`
at `p=53, NearestAway`, `sqrt(2473446)` at the same cell, and
several `fma` cases). The slice-6h diagnosis attributed all of
these to the differential lane's converter: `bigfloat_to_rug`
routed through pfloat's `core::fmt::Display` output and
`rug::Float::parse`, both of which round under NearestEven, so
values produced under non-NE rounding lost the 1 ULP separation
through the round trip. The proposed fix was a public `raw_parts`
accessor on `BigFloat` that lets the converter read the limbs
directly.

Slice 7a lands that fix and surfaces a second finding that the
slice-6h diagnosis did not separate cleanly from the converter
loss: pfloat's kernels for `div`, `sqrt`, `fma`, and `parse_str`
are correctly rounded under NearestEven but not (necessarily)
under the other four IEEE modes. With the bit-exact converter in
place, the differential lane runs `add`, `sub`, and `mul` under
all five rounding modes cleanly; `div`, `sqrt`, `fma`, and
`parse_str` continue to fail by 1 ULP under non-NE rounding on
the same input cells slice 6h captured. The kernel gap is a real
correctness gap and is tracked separately for slice 7c (or a
dedicated successor slice).

## Decision

1. **Public `BigFloat::parts` accessor and `pfloat::Parts<'a>` enum.**
   Mirror of the internal `Class` tagged union, with borrowed slice
   storage. Each variant carries the IEEE-shaped fields for its
   value kind; `Normal` includes the precision so downstream tools
   can build representations without a separate `precision()` call.
   No converse constructor is exposed: pfloat's construction paths
   enforce top-bit-set mantissa normalization, the
   `limbs_for(precision)` storage shape, and the precision-bound
   payload length as compile-time and construction-time invariants,
   and a raw-parts constructor would bypass those.

2. **Bit-exact differential converter on top of `parts()`.** The
   helper `bigfloat_to_rug` in `tests/differential/mod.rs` builds a
   `rug::Float` from a `BigFloat` by constructing a `rug::Integer`
   from the little-endian mantissa limbs, setting a `Float` at the
   originating precision (exact because the integer's
   precision-bit-wide top-bit-set value fits), and applying the
   exponent shift via `mul_2si` (the `Float << i32` operator on
   rug). NaN payload is not preserved because MPFR does not expose
   payload bits via the public API; differential tests do not
   compare NaN values for bit-equality (NaN != NaN under IEEE) so
   the loss is intentional and matches the differential lane's
   pre-existing semantics. Specials route through MPFR's
   `Special::Infinity` and `Special::Nan` constructors. The
   `Display + parse` path is removed.

3. **Per-op rounding-mode coverage.** The constants
   `BIT_EXACT_ROUNDING_MODES` (all five) and
   `NEAREST_EVEN_ROUNDING_MODES` (NE only) replace the old
   blanket `ALL_ROUNDING_MODES` constant; an alias of the same
   name pointing at the NE-only list preserves the existing
   imports for tests still gated to NE. Each
   `tests/differential_<op>.rs` picks the list that matches its
   kernel's guarantees:

   | Op | Coverage | Why |
   | --- | --- | --- |
   | `add`, `sub`, `mul` | All five | Bit-exact correctly rounded under any IEEE mode; the rounding pipeline produces the unique correctly-rounded result for each mode. |
   | `div`, `sqrt`, `fma`, `parse_str` | NE only | Bit-exact correctly rounded under NearestEven; the non-NE rounding decision can disagree with MPFR by 1 ULP on tie-adjacent inputs. Tracked as a Phase 7 follow-up alongside slice 7c (the `pow` Ziv retry). |
   | All transcendentals, tier-1 specials, AGM | NE only | Fixed working-precision guard (`target + 64`) is not Ziv-strategy retry; non-NE modes can diverge by 1 ULP under tie cases. Slice 7c addresses `pow`; broader transcendental correctness lands in subsequent slices. |

4. **Converter sanity tests live with the helper.** A nested
   `mod converter_tests` in `tests/differential/mod.rs` covers the
   round-trip for arithmetic results under NearestEven (the case
   the converter must hit exactly), the slice-6h problem case
   under NearestEven (regression guard against my own converter
   logic), and integer round-trips at the four CI precisions
   (53, 113, 256, 1024).

## Consequences

- The 1-ULP false negatives slice 6h documented for `div`, `sqrt`,
  and `fma` under non-NE are now correctly attributed: they are
  real kernel correctness gaps in those ops, not converter
  artifacts. The fix is per-kernel work that lifts each op from
  "correctly rounded under NE" to "correctly rounded under all
  five modes." That is tracked as a Phase 7 task that mirrors
  slice 7c's `pow` work.
- The differential lane now exercises five-mode coverage on three
  ops (`add`, `sub`, `mul`), which is the largest single
  expansion of MPFR-against-pfloat coverage since slice 6a. Total
  per-cell coverage on those ops scales by ×5; the CI cost
  scaling is modest because `add`/`sub`/`mul` on i64 operands are
  individually cheap.
- The `pfloat::Parts` API is public and stable. Future downstream
  tools (serializers, alternative formatters, MPFR adapters) can
  build on it without taking a dependency on test-only helpers.
- `rug`'s `integer` feature joins the dev-dependency set
  (alongside `float`). The differential lane needs
  `rug::Integer::from_digits`. No change to the published-crate
  dependency surface; `rug` remains a Unix-only dev dependency.
- ADR-0014's slice 6h status update overstated what the converter
  fix would close. ADR-0016 narrows the claim: the converter
  closed the converter-side false negatives; the kernel-side
  non-NE gaps for `div`, `sqrt`, `fma`, and `parse_str` remain
  open.

## References

- Plan: `let-s-review-the-backlog-vast-harbor.md` — slice 7a.
- ADR-0014 — MPFR differential gating. The slice-6h status update
  flagged the limitation this ADR partially closes; the remaining
  kernel-correctness piece is the open work.
- Slice 7c (planned) — `pow` Ziv-strategy retry. The pattern
  `pow` adopts will guide the kernel-correctness work on `div`,
  `sqrt`, `fma`, and `parse_str` for non-NE rounding.
