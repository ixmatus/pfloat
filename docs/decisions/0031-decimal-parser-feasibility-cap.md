# ADR-0031: Decimal parser feasibility cap and the intrinsic pow5 cost

- **Status**: accepted
- **Date**: 2026-05-21

## Context

Roadmap slice 8a closes the inline TODO in `src/parse.rs:349`.
`BigFloat::parse_str` rounds `digits * 10^exponent` correctly by
building the *exact* multi precision integer
`m * 5^exponent * 2^exponent` and rounding once through the
universal rounding pipeline with an exact sticky bit
(`src/parse.rs:385` onward). For huge `|exponent|`,
`pow5(|exponent|)` is enormous: at `|e| = 10^6` it is already on the
order of a few megabytes; at `|e| = 10^9` it is on the order of a
gigabyte and a multi second multiply.

To keep that under control, slice 4 capped `|exponent|` at the
magic constant `MAX_DECIMAL_EXPONENT = 1_000_000` and short
circuited any `|e|` past that to `±inf + OVERFLOW + INEXACT` or
`±0 + UNDERFLOW + INEXACT` without allocating. The doc on the
constant called this a "computational feasibility" bound and the
inline TODO proposed a Phase 7 follow up to "replace the explicit
`pow5` build with a logarithm based dispatch so this cap can move
from computational feasibility to exponent fits the binary exponent
type only".

That TODO framing rests on two assumptions that do not hold under
derivation:

1. The parsed decimal exponent is `i32`. The lexer rejects exponents
   that do not fit `i32` with `ParseError::ExponentOutOfRange`
   (`src/parse.rs:250`). Therefore the largest exponent the kernel
   can see is `i32::MAX`, and `i32::MAX * log2(10) ~ 7.1e9` bits is
   far below the `i64::MAX ~ 9.2e18` binary exponent range pfloat
   carries (ADR-0006). No valid parsed input ever overflows pfloat's
   binary exponent. The "exponent fits the binary exponent type"
   cap is therefore unreachable for any valid input; adopting it as
   the operative cap would never short circuit and every large
   exponent would fall through to the full `pow5` allocation. That
   reintroduces the very resource exhaustion the 1e6 cap exists to
   prevent (the doc at `src/parse.rs:81` already states the cap is
   the only defence pfloat offers against an oversized exponent).
2. The big `pow5` is intrinsic to correctly rounded decimal to
   binary conversion. The hard to round cases sit arbitrarily close
   to the halfway point between two adjacent binary floats, and
   resolving them requires the exact rational value, hence the
   exact `5^|e|`. This is the result of the Clinger / Steele and
   White / Gay line of work (Clinger 1990, "How to Read
   Floating Point Numbers Accurately", PLDI; Gay 1990 / `dtoa`).
   Modern fast parsers keep a fast 128 bit path for typical inputs
   but explicitly fall back to an exact bignum on near tie inputs;
   Lemire and Lemire's `fast_float` is the canonical example. There
   is no algorithm that is simultaneously correctly rounded and
   bounded allocation for adversarial huge exponent inputs.

Together, those two facts say the cap is intrinsically a *resource
budget*, not a deferrable algorithm bug. The 1e6 number itself was
a recalled magic constant; the *existence* of a cap is structural.

## Decision

Cap `|decimal exponent|` at a bound *derived in code* from an
explicit `pow5` storage budget; beyond, short circuit to
`±inf + OVERFLOW + INEXACT` or `±0 + UNDERFLOW + INEXACT` as before.

Budget: **16 MiB** for the `pow5(|e|)` intermediate. Derivation in
code:

```text
bits(pow5(e)) = ceil(e * log2(5))
log2(5) > 232 / 100 = 2.32          (since log2(5) = 2.321928...)
e * 2.32 <= e * log2(5) <= budget_bits
=> e <= budget_bits / 2.32 = budget_bits * 100 / 232
```

The rational lower bound on `log2(5)` makes the integer division
conservative (rounds the cap down, so the worst case `pow5` stays
inside the budget). At 16 MiB:

```text
budget_bits = 16 * 1024 * 1024 * 8 = 134_217_728
MAX_DECIMAL_EXPONENT = floor(134_217_728 * 100 / 232) = 57_852_468
```

That is ~57 times the recalled 1e6 magic, while staying inside a
documented memory budget the reader can audit and raise. The
existing doc remark on the cap (real workloads do not approach
`10^9`) still applies; the new cap simply replaces a recalled
number with a derived one and lets a much wider band of inputs
parse to correct finite values.

The short circuit branches in `finite_to_bigfloat` are unchanged in
shape; only the constant `MAX_DECIMAL_EXPONENT` changes, and its
doc records the derivation and the intrinsic pow5 fact.

## Consequences

- The cap is principled and derived. The comment in `parse.rs`
  carries the derivation and a pointer to this ADR; a future reader
  does not re encounter the "magic number" objection or the
  unsound "exponent fits the binary exponent type" framing.
- Behavior change at the input boundary: strings with `|exponent|`
  in `(10^6, ~5.785 * 10^7)`, previously short circuited to `±inf /
  ±0`, now parse to correct finite values. Strings with `|exponent|`
  past the new cap still short circuit; strings beyond `i32` still
  return `ParseError::ExponentOutOfRange` from the lexer.
- The DoS protection character is preserved: an attacker supplied
  oversized exponent still cannot drive an unbounded allocation;
  the bound is the explicit 16 MiB storage budget rather than the
  arbitrary 1e6. Within the cap, allocation is at most a few times
  the budget; beyond it, allocation is zero.
- A Lemire / `fast_float` style fast path (128 bit common path with
  an exact slow fallback) is a Phase 7 style performance follow up.
  It does not change the cap or the correctness guarantee: the
  slow path bignum stays, for the same intrinsic reason, and the
  same storage budget applies. The TODO is therefore considered
  closed; the algorithmic follow up is performance, not
  correctness.
- The 8a slice's parser coverage (8a.7) pins the new boundary on
  both sides and confirms parse to format round trip still holds
  through the conversion.

## Related

- `src/parse.rs:81` (untrusted input doc) and `src/parse.rs:336` to
  `353` (the constant and its derivation), now updated to point at
  this ADR.
- ADR-0006 (i64 exponent, the i64 binary range that the parser
  always fits within for valid i32 decimal exponents).
- ADR-0005 (tagged enum, the reason `i64::MAX` is a saturating
  sentinel not infinity; see also the DESIGN.md "Exponent" note
  updated alongside the 8a.5b `mul`/`div`/`fma` saturation fix).
- Primary references for the intrinsic pow5 result: Clinger 1990
  "How to Read Floating Point Numbers Accurately" (PLDI); Gay 1990
  "Correctly Rounded Binary Decimal and Decimal Binary Conversions"
  (Bell Labs NA report, the `dtoa` lineage); Lemire and Lemire
  2021 "Number parsing at a gigabyte per second" (`fast_float`).
  Cited as the body of work establishing the property; pfloat's
  current implementation is its own derivation from the exact
  rational, not adapted from any of these.
- Plan: `plans/abundant-yawning-badger.md` (slice 8a.6); refined
  from the plan's literal "derive from i64::MAX / log2(10)" framing
  because that derivation does not survive the i32 exponent fact
  and would reintroduce the documented resource exhaustion. Fork
  resolved with the user on 2026-05-21.
