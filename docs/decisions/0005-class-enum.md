# ADR-0005: Special-value encoding via tagged `Class` enum

- **Status**: accepted
- **Date**: 2026-05-10

## Context

IEEE 754 values include `±0`, `±∞`, signaling and quiet NaNs (with
payload), and finite normals. A representation has to encode all of
these unambiguously.

MPFR encodes special values by reserving exponent values: a sentinel
exponent flags zero, another flags infinity, another flags NaN. Sign
rides as a separate field. The mantissa contents for special values
are unspecified. The wins are storage compactness (no discriminant
byte) and a single struct shape.

The cost is fragility. Every code path that touches the value has to
remember which exponent values are sentinels. A bug that lets a
finite value slip into the sentinel range (or vice versa) becomes a
correctness incident that the type system did nothing to catch. MPFR
mitigates this with disciplined accessor functions; pfloat is not
constrained to mirror MPFR's layout choices.

The Rust alternative is a tagged enum. The discriminant byte costs
storage; the upside is that the compiler refuses programs that mix
the variants incorrectly. Pattern matches are exhaustive; the
classifier is "what variant is this" rather than "what range does
the exponent fall into."

## Decision

Use tagged enums for both `BigFloat` and `FixedFloat<PREC>`.

```rust
pub enum Class {
    Zero { sign: Sign },
    Infinity { sign: Sign },
    Nan { quiet: bool, sign: Sign, payload: NanPayload },
    Normal { sign: Sign, exponent: i64, mantissa: Vec<u64> },
}

pub enum ClassFixed<const PREC: u32> {
    Zero { sign: Sign },
    Infinity { sign: Sign },
    Nan { quiet: bool, sign: Sign, payload: NanPayload },
    Normal { sign: Sign, exponent: i64, mantissa: [u64; limbs_for(PREC)] },
}
```

Sign rides inside each variant rather than at the outer level. This
keeps `+0` and `−0` distinct (IEEE-required), lets `+NaN` and `−NaN`
remain distinct (IEEE-required for sign-bit propagation through
operations like `copysign`), and matches the spec's mental model
where the sign is part of the value's identity.

The `quiet` flag on `Nan` is required by IEEE 754-2019 §6.2.1: a
signaling NaN traps on most operations and is converted to a quiet
NaN under specified circumstances. The `payload` carries the
diagnostic bits per §6.2.2.

## Consequences

**Wins:**

- Pattern matches on the value are exhaustive. The compiler refuses
  to let a code path forget that NaN exists.
- `Normal` is the only variant that carries a mantissa, so the
  arithmetic kernels destructure once and work on the relevant
  fields. Special-value dispatch is a `match` at the top of each
  kernel, not a series of magic-number comparisons.
- Sign rides with the variant, so the IEEE-required distinctions
  (signed zero, signed NaN) are representable by construction.
- The enum is `Copy`-free for `BigFloat` (because of `Vec`) but
  `Copy` for `FixedFloat<PREC>` at any reasonable PREC; that
  asymmetry mirrors the storage decision in ADR-0004.

**Costs:**

- One discriminant byte per value, plus padding to the alignment
  boundary. For `BigFloat` (which already carries a `Vec`'s 24-byte
  footprint) the cost is invisible; for `FixedFloat<53>` (single
  limb mantissa) the discriminant is a measurable fraction of the
  total. Acceptable; correctness wins over packing density at this
  scale.
- The Kani harness in ADR-0009 has to enumerate variants explicitly.
  The cost is one extra `match` arm per harness; the win (catching
  variant-mishandling at proof time) is the point.
- Differential testing against MPFR has to translate at the boundary:
  pfloat's `Class::Nan { quiet: true }` corresponds to MPFR's
  `mpfr_nan_p` returning true. The translation is mechanical and
  lives in the `differential-mpfr` test module.

## Related

- ADR-0001 (limb representation, where the sign-magnitude rule is
  set)
- ADR-0004 (mantissa storage)
- ADR-0009 (verification scaffolding)
- DESIGN.md, "Special values" subsection.

## Update (2026-05-10)

The `payload: NanPayload` field shown in the variant declarations is replaced with **raw limb storage matching the mantissa shape**: `payload: alloc::vec::Vec<u64>` inside `Class::Nan` and `payload: [u64; limbs_for(PREC)]` inside `ClassFixed::<PREC>::Nan`. The `NanPayload(u64)` newtype originally implied by the ADR body is dropped from 1.0; the variable-width per-precision representation is faithful to IEEE 754 §6.2.2's "diagnostic information in the trailing significand" framing and avoids forcing a u64 ceiling on the payload of a 4096-bit mantissa.

The rest of the ADR (tagged enum over reserved-exponent encoding, sign rides inside each variant, quiet/signaling distinction in the `Nan` variant) is unchanged. The accessor APIs on `BigFloat` / `FixedFloat<PREC>` expose the payload as `&[u64]` for read and as `&[u64]` (with length validation) for construction.
