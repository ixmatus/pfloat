# ADR-0128: the decimal cap governs value magnitude, not the point-baked exponent

- **Status**: accepted
- **Date**: 2026-07-08
- **Amends**: ADR-0031 / ADR-0051 (the `10^6` resource cap stands; what it is measured against changes)

## Context

The libFuzzer `parse` target asserts that any string `parse` accepts as a
finite value round-trips: `parse(Display(parse(s))) == parse(s)` under
nearest-even. On the R5.2+R5.3 merge (`7584632`) it failed with a finite
value comparing GREATER than its own re-parsed Display output. The failure
was pre-existing (the fuzz smoke is non-deterministic; earlier merges'
seeds missed it) and unrelated to the merged math kernels.

Root cause. `parse` builds `value = m · 10^exponent` from the integer
mantissa `m` and bakes the decimal-point position into `exponent`
(`exponent = e_part − frac_digits`), then caps `|exponent|` at
`MAX_DECIMAL_EXPONENT = 10^6` (a pow5 resource budget, ADR-0031). But
`Display` renders `round_trip_digit_count(precision)` significant digits
(36 for `p = 113`), so re-parsing its output uses a baked exponent
`value_magnitude − (digits − 1)`, which sits up to `digits − 1` BELOW the
value's decimal magnitude. For a value within `~digits` of the cap, that
re-parse exponent exceeds the cap and `parse` saturates its own faithful
output:

- `9e-999990` → `Display` `8.9999…987e-999990` (36 digits, correct) →
  re-parse baked exponent `−1000025 < −10^6` → underflow to `0` →
  `parsed > reparsed`. An independent `Fraction` oracle confirms 36 digits
  round-trips this value mathematically, so `parse` is wrong to decline it.
- `99e1000000` = `9.9e1000001` stayed finite (its own `exponent = 10^6` is
  within the cap) yet its value magnitude `1000001` exceeds the FORMAT cap,
  so it rendered as the approximate `1e{D}` token (ADR-0120) and re-parsed
  to `inf`: a finite value that could not render exactly.
- `1e-1000000` (magnitude exactly `−10^6`): the `30103/100000` `log10`
  estimate is off by ~`0.014` at `log2 ≈ 3.3·10^6`, so its floor landed at
  `−1000001`, one past the FORMAT cap, wrongly routing it to the
  approximate token.

The failures form a narrow band within `round_trip_digit_count` orders of
`±10^6` — astronomically extreme magnitudes, but a genuine
`parse`/`Display` inconsistency the fuzz correctly caught.

The resource cost is `pow5(|exponent|)`, and `|exponent|` differs from the
value magnitude by only the digit count; at these magnitudes the pow5 for
`9e-999990` and its 36-digit form differ by ~35 in the exponent, a
negligible cost difference. So capping on magnitude rather than the baked
exponent is a CONSISTENCY fix, not a change to the resource budget.

## Decision

Separate the two things the cap governs.

1. **Representability is capped on the value MAGNITUDE**
   `L = exponent + (digits − 1)` (the leading significant digit's place).
   A value with `|L| > MAX_DECIMAL_EXPONENT` saturates by `L`'s sign
   (overflow / underflow), keeping the mode-aware saturation of pf-mw6u
   (ADR unchanged). A large positive `exponent` implies `L ≥ exponent`, so
   overflow needs only this magnitude test. `|L| ≤ MAX` is NECESSARY for a
   finite result but not sufficient: the pow5 guard of item 2 can still
   saturate a small-magnitude input carrying an extreme digit count. So
   representability is magnitude-driven (precision-independent at the
   overflow boundary), while the underflow boundary additionally depends
   on the input's digit count through `round_trip_digit_count`.

2. **The pow5 COST budget is widened to `MAX + round_trip_digit_count`.**
   A value's own `Display` output re-parses with `|exponent| ≤ |L| + rt −
   1`, so this budget admits every renderable value while still bounding a
   pathological small-magnitude, many-fractional-digit direct input (whose
   `|L|` stays small but whose `|exponent|` blows up) — that case
   saturates on the widened `|exponent|` guard, as before.

3. **`fmt`'s exact-render cap is widened by the same `round_trip_digit_count`.**
   Every value `parse` produces (`|L| ≤ MAX`) must render EXACTLY so it
   round-trips; the `log10` estimate is off by up to ~1 at the cap, and a
   bare `MAX` cutoff wrongly routed magnitude-`MAX` values to the
   approximate token. The wider bound is cost-safe (the exact render's
   pow5 stays `O(pow5(MAX + rt + digits))`) and still emits the ADR-0120
   approximate token for genuinely-huge arithmetic-built magnitudes.

`parse` carries a local `round_trip_digit_count` copy (the formula
`ceil(p·log10 2) + 1`) so it does not depend on the `fmt` feature; it is
kept in sync with `fmt::BigFloat::round_trip_digit_count`.

## Consequences

- Every value `parse` produces round-trips through its own `Display` under
  nearest-even, at any magnitude within the cap and any input precision.
  The fuzz `parse` target passes; the regression
  (`tests/regression_review_2026_07_08_bpaq.rs`) pins the tiny-near-cap
  family, the magnitude-past-cap saturation, and precision-independence.
- **Behaviour change:** a value whose magnitude exceeds the cap because of
  its mantissa digits (e.g. `99e1000000` = `9.9e1000001`) now saturates to
  `inf` (overflow) rather than remaining finite-but-approximately-rendered.
  This makes `parse` and `fmt` agree on representability. Values `2×` past
  the cap (`1e2000000`) saturate exactly as before (pf-mw6u regression
  green).
- **Inversion: the round-trip invariant is nearest-even only.** `Display`
  targets nearest-even (its digit count is round-trip-safe under NE, not
  the directed modes), so `parse(Display(v), directed)` is NOT expected to
  equal `v`. The regression sweep and the fuzz both round in NE for both
  legs; a first-draft directed-mode sweep failed on exactly this and was a
  test error, not a defect.
- **Inversion: `parsed > reparsed` was a PARSE fault, not a `Display`
  one.** `Display` rendered the correct 36 digits; `parse` saturated them.
  The oracle (mathematical 36-digit round-trip) located the fault, against
  the intuition that a rendering mismatch is a formatter bug.
- **Named boundary:** a direct input at magnitude exactly `−cap` carrying
  `round_trip_digit_count + 2` or more significant digits saturates on the
  pow5 guard: its baked exponent `−cap − (digits − 1)` drops below the
  strict `−(cap + rt)` bound (an `rt + 1`-digit input at `−cap` sits
  exactly on the bound and stays finite). A value's own `Display` output,
  carrying exactly `rt` digits, never reaches this. This is the documented
  resource limit at an astronomically extreme magnitude, essentially
  unchanged from the pre-fix cap. The `<` (not `≤`) on the guard is
  load-bearing: a value that rounds down across a power of ten renders at
  `decimal_exp = −(cap + 1)`, whose `rt`-digit re-parse bottoms out at
  `−(cap + rt)`, which the strict bound admits with zero margin.

## References

- pf-bpaq (the bug), ADR-0031 (the pow5 resource cap), ADR-0051 / ADR-0120
  (the format cap and its approximate-magnitude token), pf-mw6u (the
  mode-aware saturation).
- `tests/regression_review_2026_07_08_bpaq.rs`; oracle
  `scratchpad/rt_oracle.py`.
