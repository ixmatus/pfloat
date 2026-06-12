# ADR-0108: u32 arithmetic at the precision ceiling

- **Status**: accepted
- **Date**: 2026-06-12

## Context

pf-9wb2 (verified): `parse_str("0.5", u32::MAX, NearestEven)`
**panicked** at the `expect("non-zero quotient")` — at precisions
near the documented `u32::MAX` ceiling (ADR-0002) the buffer-width
add `m_bits + l` wrapped, the under-sized buffer truncated the
mantissa's top bits, and the quotient came out zero. The wrap made
the failure cheap to hit: a three-character literal and a precision
parameter. The review's code-read found six sibling expressions in
the ops (`mul`, `div`, `sqrt`, `cbrt`, `addsub`, `fma`) that wrap or
panic in the same regime.

## Decision

**parse**: clamp `l` itself (`min(u32::MAX − m_bits)`) so every
downstream use — the buffer size, the shift, the exponent
bookkeeping — stays self-consistent; the first draft saturated only
`total_bits`, which desynchronized the buffer from the shift and
silently produced a WRONG VALUE for the dyadic `"0.5"` (caught by
the timed verification run, not by the panic's absence — saving then
reloading is the hedge). At the ceiling the clamp trims only guard
tail. The two `limbs · 64 − precision` placement shapes compute in
u64. Verified by run: `parse_str("0.5", u32::MAX)` now returns
exactly 0.5 with Status OK in ~0.6 s release (the honest ~512 MB
cost of the request); release-gated lane row.

**div/sqrt/cbrt**: `l` is clamped at its SOURCE
(`min(u32::MAX − p_a)`, after cbrt's mod-3 bump and sqrt's parity
bump) so the buffer, the shift, and the exponent bookkeeping stay
consistent. The first draft saturated only `total_bits` — repeating,
in the same slice, the exact desync recorded below as parse's
lesson; the slice's adversarial verifier refuted it by run
(`div_round(3@u32::MAX, 7, 53)` and `cbrt_round(2@u32::MAX, 53)`
still panicked in both profiles). sqrt's 1-bit parity bump happened
to fit `limbs_for`'s rounding slack; it gets the clamp anyway so the
invariant does not rest on a coincidence.

**mul/fma**: the product of two near-ceiling operands spans more
bits than a u32 can name; the raw `(top_bit + 1) as u32` wrapped
silently. Saturating conversion plus a debug assertion documenting
the envelope: reaching it requires two ~2^31-bit operands (half a
gigabyte of mantissa each), and full support would need
u64-precision plumbing through the rounding pipeline — recorded here
as the ceiling's edge rather than silently wrong. **addsub**: the
aligned-span `expect` could panic for two ceiling-precision operands
with a small exponent gap (legal input); now saturating with the
same debug assertion, dropping at most the bottom few bits below the
sticky horizon.

## Consequences

- Operations stay total and panic-free across the documented
  precision domain in both profiles; at the very edge
  (operand-precision sums past `u32::MAX`) behavior is
  self-consistent saturation with a debug-build assertion, recorded
  as the supported envelope.
- The `limbs · 64 − precision` placement shape was a FAMILY, not two
  sites: the slice verifier found the same expression in
  `round_finite_to_precision` (a debug overflow panic across the
  whole near-ceiling parse band — release was self-correcting
  modular arithmetic, so the release-gated lane row could not see
  it) and in `mantissa::storage_shift` (a debug panic at mere
  construction of a ceiling-precision value); the post-fix probe
  then hit `ops::limbs::extract_as_integer`, and the full sweep
  found nine more across the ops (div/mul/fma/sqrt/cbrt/addsub/
  remainder/adjacent dst-placement shapes). All thirteen now compute
  in u64; the grep-able invariant (no u32 `· 64 −` arithmetic) holds
  crate-wide, probe-verified in debug across the ceiling band.
- Failure modes considered (inverted):
  1. **Saturation could silently change mid-domain results** — every
     clamp engages only within ~64 bits of the ceiling; the full lib
     and lane suites run unchanged.
  2. **The parse clamp could trim accuracy a caller relies on** — it
     trims guard bits only at requests within `m_bits` of
     `u32::MAX`, where the guard's headroom still dominates the
     trim for any literal short enough to construct in memory.
  3. **The first draft's lesson, learned twice in one slice**: a
     saturating fix not applied at the SOURCE of the quantity
     desynchronizes its consumers. Parse's first draft hit it (the
     wrong-value regression, caught by reloading the parsed value
     rather than trusting the panic's absence) — and the div/cbrt
     edits in the SAME slice repeated it, caught only by the
     adversarial verifier's run. The pattern is now uniform: clamp
     `l` where it is born, never its derived widths.
  4. **A "sweep" that enumerates sites from one defect's
     neighborhood misses the family**: the placement-shape fix
     covered parse's two expressions; the verifier found the same
     shape in the rounding funnel and the mantissa helper. The
     grep-able invariant (`· 64 −` arithmetic on u32) now holds
     nowhere in the crate.

## Related

- Issues: pf-9wb2 (closed by this ADR), epic pf-8iji.
- Other ADRs: ADR-0002 (the bit-level precision ceiling this
  hardens), ADR-0107 (the bottom-rim sibling, same slice), the
  parse-OOM budget posture (`feedback_bignum_dos_budget`).
