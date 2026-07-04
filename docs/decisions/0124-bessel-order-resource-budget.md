# ADR-0124: a resource budget on Bessel integer order (the pf-ap01 DoS)

- **Status**: accepted
- **Date**: 2026-07-03

## Context

pf-ap01 (code-read, not executed). The integer-order Bessel kernels
`Y_n`/`J_n`/`I_n` evaluate a normal argument by an O(|n|)-step recurrence
(the `Y` upward recurrence; the `J`/`I` Miller backward descent) at a
working precision that grows with the order (`bessel_y`'s recurrence uses
`~32 + 4·|n|` extra bits). The order `n` is a caller-supplied `i32`, so an
attacker-controlled order up to `i32::MAX ≈ 2.1e9` is an **unbounded
resource DoS**: `~4·|n|` bits is a multi-hundred-MB mantissa and O(|n|)
steps of arithmetic at that width is hours. The current code is *correct*
— it would eventually return the right value — but the cost is not
bounded by anything the caller controls that is checked.

## Decision

Impose a shared **resource budget** `MAX_BESSEL_ORDER = 2^14 = 16384` on
`|n|` for the recurrence kernels. In the normal-argument path (after the
O(1) special cases — poles, `x = 0` exact values, ±∞ — which are *not*
capped), refuse an order with `|n| > MAX_BESSEL_ORDER` and return
`NaN + INVALID`.

`2^14` bounds the worst case to ~2 s and ~32 KB of mantissa while staying
far above any order that arises in practice (physical/engineering Bessel
orders are `< ~10^3`). A larger order is *representable but not computable
within budget* — the same posture as parse/format past their cost caps
(ADR-0031/0051). Refusal (not saturation) is chosen because the value
does **not** overflow pfloat's `i64` exponent range for any `i32` order at
`x ≥ 1` (`|n|·log₂|n| ≤ 2.1e9·31 ≈ 6.5e10 ≪ i64::MAX ≈ 9.2e18`), so there
is no honest `±∞`-saturation to return; `NaN + INVALID` signals "outside
the supported computable range." (Only an extreme-tiny-`x`, large-`n`
corner genuinely overflows; distinguishing it would need its own
derivation and is not worth a special case for an already-refused input.)

Applied at `bessel_y_kernel`, `bessel_j_kernel`, and `bessel_i_kernel`;
the constant lives in `bessel_j.rs` (`pub(super)`) and the `Y`/`I` kernels
reference it.

## Consequences

- The order-driven DoS is closed: cost is bounded by a checked, caller-
  independent constant. In-cap orders (every order that arises in
  practice) are unchanged — the differential lanes (`differential_jn`,
  `differential_yn`, `differential_ik`) stay green, and the exact `x = 0`
  cases are unaffected (they short-circuit before the cap).
- Refusal is a behaviour change for `|n| > 16384`: a mathematically finite
  (and, up to ~`10^17`, exponent-representable) value now returns
  `NaN + INVALID` instead of grinding. This is the deliberate
  security/frugality-over-convenience trade (the value hierarchy), matched
  to the cost-cap precedent.
- The value at `|n| > cap` remains recoverable in principle by a future
  large-order uniform (Debye) asymptotic that is O(1) in the order rather
  than O(|n|); that is an additive enhancement (a new kernel), out of
  scope here and not required to close the DoS.

### Inversion (failure paragraphs considered)

- *"Cap only the working precision, not the order."* Insufficient: capping
  the `4·|n|` boost bounds memory but the recurrence still runs O(|n|)
  steps, so the time DoS survives. The order itself is the lever.
- *"Saturate to ±∞ instead of refusing."* Rejected: the value does not
  overflow the exponent range for an `i32` order at `x ≥ 1`, so `±∞` would
  be a *wrong* value with a clean-looking flag — the exact certified-wrong
  posture the whole review remediation is closing. `NaN + INVALID` is the
  honest "unsupported" signal.
- *"16384 is too low — it refuses computable orders."* Orders in
  `(16384, ~50000]` are computable but slow (seconds to minutes), and
  above that they are genuine DoS territory; a cap in the middle of that
  band is the least-surprising bound. Real Bessel orders are three orders
  of magnitude below the cap, so the refusal band is exotic.

## References

- pf-ap01 (epic pf-8iji); ADR-0031/0051 (the parse/format cost-cap
  precedent), the bignum DoS-budget posture, the beta case-4 "derive the
  cost bound against the input domain" lesson.
- `src/math/bessel_j.rs` (`MAX_BESSEL_ORDER`, `bessel_j_kernel`),
  `src/math/bessel_y.rs` (`bessel_y_kernel`), `src/math/bessel_i.rs`
  (`bessel_i_kernel`); test
  `bessel_y::tests::bessel_order_resource_budget_refuses_exotic_orders`.
