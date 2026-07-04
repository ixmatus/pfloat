# Architecture Decision Records

This directory holds the record of *why* pfloat is the way it is.
Each significant choice (numeric representation, API shape, feature
gating, verification posture, performance tradeoffs) gets one
Architecture Decision Record. Together they form the audit log a
future reviewer would otherwise have to reconstruct from commit
messages and release notes.

The format is borrowed from ferrodec, which borrowed it from the
broader ADR community.

## Conventions

- **Filenames**: `NNNN-short-slug.md`, four-digit zero-padded sequence
  number, lowercase slug. Numbers are never re-used; superseded ADRs
  keep their slot and link forward.
- **Format**: see `template.md`. Each ADR is short. A single page is
  the target; the form matters more than the length.
- **Status lifecycle**:
  - `proposed` — drafted, not yet acted on. Avoid this for
    retroactive ADRs.
  - `accepted` — the decision is in effect.
  - `superseded by ADR-NNNN` — replaced; keep the file as a
    historical record, link forward.
  - `rejected` — considered and decided against. Document for the
    next person who wonders the same thing.
- **Plans**: planning artifacts archive under `plans/` with a date
  prefix (`YYYY-MM-DD-slug.md`). They are snapshots of the state at
  decision time, not living documents. ADRs reference the plan that
  produced them when applicable.

## Writing a new ADR

1. Pick the next available number.
2. Copy `template.md` to `NNNN-your-slug.md`.
3. Fill in: status, date, context, decision, consequences, related
   references.
4. If the decision supersedes a prior one, edit the prior ADR's
   status line to `superseded by ADR-NNNN`.

Decisions that are reversible or local in scope do not need an ADR.
These are for choices that matter to future contributors deciding
whether to revisit a path.

## Index

This list is generated from the ADR files by `scripts/adr-index.sh`
and checked in CI (`scripts/adr-index.sh --check`); each ADR file is
the source of truth for its own title and status. Edit the script (the
preamble and this note live in its heredocs), not this file by hand.

- [ADR-0001: `u64` limb representation, sign-magnitude, top-bit-set normalization](0001-limb-representation.md) (accepted)
- [ADR-0002: Bit-level precision granularity](0002-bit-level-precision.md) (accepted)
- [ADR-0003: Dual API, `BigFloat` (dynamic) and `FixedFloat<const PREC: u32>` (const-generic)](0003-dual-api.md) (accepted)
- [ADR-0004: Mantissa storage, `Vec<u64>` for `BigFloat`, `[u64; N]` for `FixedFloat`](0004-mantissa-storage.md) (accepted)
- [ADR-0005: Special-value encoding via tagged `Class` enum](0005-class-enum.md) (accepted)
- [ADR-0006: `i64` exponent](0006-exponent-type.md) (accepted)
- [ADR-0007: Rounding mode and exception flags](0007-rounding-and-flags.md) (accepted)
- [ADR-0008: Differential testing oracle, `gmp-mpfr-sys` on a feature-gated CI lane](0008-differential-oracle.md) (accepted)
- [ADR-0009: Verification scaffolding, copy-paste from ferrodec, no shared crate](0009-verification-scaffolding.md) (accepted)
- [ADR-0010: Schönhage-Strassen FFT multiplication deferred to 1.x](0010-fft-deferred.md) (accepted)
- [ADR-0011: MSRV moves to nightly to use `generic_const_exprs`](0011-msrv-nightly-for-generic-const-exprs.md) (accepted)
- [ADR-0012: Kani harness architecture and CI gating](0012-kani-harness-architecture.md) (accepted)
- [ADR-0013: Fuzz harness architecture](0013-fuzz-harness-architecture.md) (accepted)
- [ADR-0014: MPFR differential CI gating and implementation choice](0014-mpfr-differential-ci-gating.md) (accepted)
- [ADR-0015: AGM kernel uses Gauss's iteration with an independent `agm` feature flag](0015-agm-formulation.md) (accepted)
- [ADR-0016: Public `BigFloat::parts` accessor and bit-exact MPFR converter](0016-bigfloat-parts-accessor.md) (accepted)
- [ADR-0017: Transcendental constants computed on the fly via AGM](0017-agm-based-transcendental-constants.md) (accepted)
- [ADR-0018: Euler–Mascheroni constant via Brent–McMillan](0018-euler-gamma-via-brent-mcmillan.md) (accepted)
- [ADR-0019: Integral special functions (Ei, Si, Ci, li)](0019-integral-special-functions.md) (accepted)
- [ADR-0020: Bit-audit and pin the 1024-bit specials constants (2/√π, ln 2π)](0020-audit-and-pin-1024-bit-specials-constants.md) (accepted)
- [ADR-0021: Airy functions Ai, Bi, Ai′, Bi′](0021-airy-functions.md) (accepted)
- [ADR-0022: `pow` Ziv interval-test retry and integer-exponent fast path](0022-pow-ziv-retry-integer-fast-path.md) (accepted)
- [ADR-0023: Bessel functions of the first kind J0, J1, Jn](0023-bessel-j.md) (accepted)
- [ADR-0024: Bessel functions of the second kind Y0, Y1, Yn](0024-bessel-y.md) (accepted)
- [ADR-0025: Modified Bessel functions I0, I1, In, K0, K1, Kn](0025-bessel-ik.md) (accepted)
- [ADR-0026: Riemann zeta function on the real line](0026-zeta.md) (accepted)
- [ADR-0027: Karatsuba threshold calibrated to 48 limbs](0027-karatsuba-threshold-calibration.md) (accepted)
- [ADR-0028: Allocation profiling, and `BigFloat` inline storage deferred to 1.x with data](0028-allocation-profiling-and-inline-storage.md) (accepted)
- [ADR-0029: Dragon4 / Steele-White shortest formatter deferred to 1.x](0029-dragon4-shortest-formatter-deferred.md) (superseded)
- [ADR-0030: Beta function on the negative real domain (sign and pole convention)](0030-beta-negative-domain.md) (accepted)
- [ADR-0031: Decimal parser feasibility cap and the intrinsic pow5 cost](0031-decimal-parser-feasibility-cap.md) (accepted)
- [ADR-0032: Libm reciprocal and root kernels (cot, sec, csc, cbrt, hypot, rootn) ship as direct primary kernels, not derived aliases](0032-libm-reciprocal-and-root-kernels-direct.md) (accepted)
- [ADR-0033: Phase 1 correctness sweep runs to completion before the v1.0 tag](0033-phase1-correctness-sweep-precedes-v1.0.md) (accepted)
- [ADR-0034: Oracle layer for the Phase 1 exhaustive `f32` sweep](0034-oracle-layer.md) (accepted)
- [ADR-0035: Oracle worker reports certified `f32` directly; three-way agreement architecture](0035-oracle-worker-protocol-and-three-way-agreement.md) (proposed)
- [ADR-0036: `property_jn::self_consistent` argument constrained to dyadic rationals (pf-ok9 lesson)](0036-property-self-consistent-dyadic-argument.md) (accepted)
- [ADR-0037: `SmallVec<[u64; 4]>` for `Class::Normal::mantissa` and `Class::Nan::payload`, rejected](0037-mantissa-payload-inline-storage-rejected.md) (rejected)
- [ADR-0038: Five-mode kernel completeness as the v1.0 strong-claim gate](0038-five-mode-kernel-completeness-as-v1.0-gate.md) (accepted)
- [ADR-0039: Phase 1g verification architecture closure (v1.0 blocker)](0039-phase-1g-verification-closure.md) (accepted)
- [ADR-0040: Schönhage-Strassen FFT multiplication — measured, not the v1.0 win](0040-fft-schoenhage-strassen-measurement.md) (accepted)
- [ADR-0041: Spouge precision-pegging — measured, rejected](0041-spouge-precision-pegging.md) (rejected)
- [ADR-0042: pf-1axr trig range-cap pre-check and bessel_y recurrence boost — root fix](0042-trig-range-cap-and-bessel-y-recurrence-boost.md) (accepted)
- [ADR-0043: Bessel per-kernel asymptotic threshold — accepted](0043-bessel-per-kernel-asymptotic-threshold.md) (accepted)
- [ADR-0044: Airy asymptotic threshold — already at the mathematical optimum](0044-airy-asymptotic-threshold-already-optimal.md) (accepted)
- [ADR-0045: Si/Ci asymptotic threshold — already at the mathematical optimum](0045-si-ci-asymptotic-threshold-already-optimal.md) (accepted)
- [ADR-0046: erf/erfc asymptotic threshold — already at the mathematical optimum](0046-erf-erfc-asymptotic-threshold-already-optimal.md) (accepted)
- [ADR-0047: Bessel Miller seed-index tightening — accepted, precision-gated](0047-bessel-miller-seed-index-tightening.md) (accepted)
- [ADR-0048: Airy asymptotic kernel working-precision boost reduction](0048-airy-asymptotic-working-precision-boost-reduction.md) (accepted)
- [ADR-0049: pf-hcz4 full pf-tqzz cross-check sweep — first execution and v1.0 baseline](0049-pf-hcz4-full-cross-check-sweep.md) (in-progress)
- [ADR-0050: tanh stable `expm1` form replaces the cancelling composition and its tiny-x short circuit](0050-tanh-expm1-stable-form.md) (accepted)
- [ADR-0051: Decimal formatter magnitude cap and sub-quadratic conversion](0051-formatter-magnitude-cap-and-sub-quadratic-conversion.md) (accepted)
- [ADR-0052: Recursive (Burnikel-Ziegler) integer division](0052-recursive-burnikel-ziegler-division.md) (accepted)
- [ADR-0053: CI gate widening for v1.0 (oracle smoke coverage and a feature-union drift guard)](0053-ci-gate-widening-v1.0.md) (accepted)
- [ADR-0054: v1.0 public API freeze](0054-v1.0-public-api-freeze.md) (accepted)
- [ADR-0055: Public f32/f64 conversion API](0055-public-f32-f64-conversion-api.md) (accepted)
- [ADR-0056: The six direct libm kernels (cot, sec, csc, cbrt, hypot, rootn)](0056-libm-direct-kernels-implementation.md) (accepted)
- [ADR-0057: The pfloat-libm outer Ziv loop (the enclosure determines the float)](0057-libm-outer-ziv-loop.md) (accepted)
- [ADR-0058: The pfloat-libm verification harness (MPFR only, value hard, range sharded)](0058-libm-verification-harness.md) (accepted)
- [ADR-0059: small-argument fast-path family for the odd elementary kernels](0059-small-argument-fast-path-family.md) (accepted)
- [ADR-0060: INEXACT flag fidelity for the transcendental exp/log and sin/cos kernels](0060-inexact-fidelity-transcendental-kernels.md) (accepted)
- [ADR-0061: Toom-Cook 3-way multiplication, the rung above Karatsuba](0061-toom3-multiplication.md) (accepted)
- [ADR-0062: The BoundedBigFloat Kani discharge of Ziv soundness is blocked at the Vec allocation level, not the arithmetic](0062-bounded-bigfloat-ziv-discharge-investigation.md) (accepted)
- [ADR-0063: INEXACT flag fidelity across the rest of the elementary transcendental surface](0063-inexact-fidelity-elementary-transcendentals.md) (accepted)
- [ADR-0064: INEXACT flag fidelity for the proven-transcendence special functions](0064-inexact-fidelity-special-functions-proven.md) (accepted)
- [ADR-0065: Beta exact-input construct-and-check (correctness bugfix)](0065-beta-exact-construct-and-check.md) (accepted)
- [ADR-0066: Defensive INEXACT guard on the gamma family and zeta](0066-defensive-inexact-guard-gamma-zeta.md) (accepted)
- [ADR-0067: Public constants-at-precision API](0067-public-constants-at-precision-api.md) (accepted)
- [ADR-0068: Optional serde serialization behind a feature](0068-serde-serialization-behind-a-feature.md) (accepted)
- [ADR-0069: IEEE 754-2019 remainder kernel](0069-ieee-remainder-kernel.md) (accepted)
- [ADR-0070: num-traits impls scoped to FixedFloat](0070-num-traits-impls-fixedfloat.md) (accepted)
- [ADR-0071: Dragon4 shortest-output formatter](0071-dragon4-shortest-formatter.md) (accepted)
- [ADR-0072: Exact power-of-two scaling (`scale_by_pow2`)](0072-scale-by-pow2.md) (accepted)
- [ADR-0073: Adjacent-representable primitives (`next_up`/`next_down`/`ulp`)](0073-adjacent-representable.md) (accepted)
- [ADR-0074: `pfloat-ball` crate and the `Mag` radius primitive](0074-pfloat-ball-crate-and-mag.md) (accepted)
- [ADR-0075: the sealed `RealScalar` trait](0075-realscalar-sealed-trait.md) (accepted)
- [ADR-0076: `Ball<T>`, the in-tree spec, and the conversion boundary](0076-ball-type-and-conversion-boundary.md) (accepted)
- [ADR-0077: ball arithmetic and the directed-pair radius-soundness law](0077-ball-arithmetic-radius-soundness.md) (accepted)
- [ADR-0078: `pfloat-ball` verification design](0078-pfloat-ball-verification.md) (accepted)
- [ADR-0079: directed-mode rounding verification posture](0079-directed-mode-verification-posture.md) (accepted)
- [ADR-0080: directed-mode rounding for the asymptote-saturating kernels](0080-asymptote-saturation-directed-rounding.md) (accepted)
- [ADR-0081: Ziv-certify the directed rounding of log2 and log10](0081-log2-log10-ziv-certification.md) (accepted)
- [ADR-0082: `pfloat-ball` interval-input Arb bracket: range soundness and tightness](0082-pfloat-ball-interval-bracket-range-soundness-tightness.md) (accepted)
- [ADR-0083: Lefèvre–Muller hardest-to-round seeding of the ball generators](0083-lefevre-muller-seeding-ball-generators.md) (accepted)
- [ADR-0084: Creusot containment-composition lemma feasibility spike](0084-creusot-containment-lemma-feasibility-spike.md) (accepted)
- [ADR-0085: directed f32 EC2 sweep for zeta and Bessel J — feasibility and go/no-go](0085-directed-f32-ec2-sweep-feasibility.md) (accepted)
- [ADR-0086: pfloat-ball 1.0 public API freeze](0086-pfloat-ball-public-api-freeze.md) (accepted)
- [ADR-0087: pfloat-ball enclosure-accuracy posture](0087-pfloat-ball-enclosure-accuracy-posture.md) (accepted)
- [ADR-0088: the fused two-product primitive and its single-rounding proof](0088-fused-two-product-primitive.md) (accepted)
- [ADR-0089: pfloat-complex and its independent sealed RealScalar trait](0089-pfloat-complex-sealed-realscalar.md) (accepted)
- [ADR-0090: complex multiply and divide, and the rounding-mode and Status API](0090-complex-multiply-divide.md) (accepted)
- [ADR-0091: complex magnitude, phase, and the elementary core with C99 Annex G branch cuts](0091-complex-magnitude-phase-elementary-annex-g.md) (accepted)
- [ADR-0092: the pfloat-complex verification posture](0092-complex-verification-posture.md) (accepted)
- [ADR-0093: pfloat-complex 1.0 public API freeze](0093-pfloat-complex-public-api-freeze.md) (accepted)
- [ADR-0094: Per-source reference registry at docs/references/](0094-reference-registry.md) (accepted)
- [ADR-0095: agm convergence floor is relative to the iterate magnitude](0095-agm-relative-convergence-floor.md) (accepted)
- [ADR-0096: exp at the exponent rim — certified dispatch and pinned reduction](0096-exp-exponent-ceiling.md) (accepted)
- [ADR-0097: input-encoded cancellation at ln/asin/lgamma/digamma — three fixes, one mechanism](0097-cancellation-boost-positive-roots.md) (accepted)
- [ADR-0098: input-structure-aware resolution for zeta, trig reduction, and beta](0098-input-structure-aware-dispatch.md) (accepted)
- [ADR-0099: Ball widening at the exponent rim and the parse-budget collapse](0099-ball-rim-saturation-soundness.md) (accepted)
- [ADR-0100: clog's enclosure guard derives its cap from the input structure](0100-clog-structure-derived-enclosure-cap.md) (accepted)
- [ADR-0101: the exponent-rim dispatch composed through the shell and the exp-family siblings](0101-rim-composition-through-shell-and-siblings.md) (accepted)
- [ADR-0102: input-encoded depth dispatches through the infinitesimal rounding — hypot, atan, atan2](0102-input-depth-infinitesimal-dispatch.md) (accepted)
- [ADR-0103: the Ziv driver takes a lazy input-derived certification depth](0103-ziv-depth-hint.md) (accepted)
- [ADR-0104: the tiny-x family gains input-precision arms and depth hints](0104-tiny-x-family-precision-arms.md) (accepted)
- [ADR-0105: acos's gap boost, and the INEXACT force covers the zero-result gap](0105-acos-gap-boost-zero-posture.md) (accepted)
- [ADR-0106: agm normalizes its operands toward exponent 0](0106-agm-exponent-normalization.md) (accepted)
- [ADR-0107: exponent arithmetic at the bottom rim — the lift, the Taylor floor, and the deep-rung ceiling](0107-bottom-rim-exponent-hardening.md) (accepted)
- [ADR-0108: u32 arithmetic at the precision ceiling](0108-precision-ceiling-arithmetic.md) (accepted)
- [ADR-0109: zeta_fe pre-boosts its working precision by the input's proximity to the negative integers](0109-zeta-fe-trivial-zero-conditioning.md) (accepted)
- [ADR-0110: interior-zero conditioning for li, Ei, and digamma — geometric cancellation growth, an Ei zero-window, and a Spouge digamma](0110-interior-zero-conditioning-li-ei-digamma.md) (accepted)
- [ADR-0111: route lgamma's Spouge regime through cancellation_boosted (the pf-rlrb certified-wrong-value sibling of ADR-0110)](0111-lgamma-spouge-cancellation-boost.md) (accepted)
- [ADR-0112: a near-zero dispatch for ζ(s), s → 0⁻ (the pf-qt7v 0·∞ form)](0112-zeta-near-zero-conditioning.md) (accepted)
- [ADR-0113: tiny-x dispatch for Si (the pf-31ql directed-mode 1-ulp fix)](0113-si-tiny-x-dispatch.md) (accepted)
- [ADR-0114: add/sub result fidelity — directed-mode sign mirroring and exponent-saturation flag parity](0114-addsub-flag-and-sign-fidelity.md) (accepted)
- [ADR-0115: complex status and sign fidelity — qNaN/sNaN INVALID, exact-quotient OK, zero-dividend sign, csqrt directed imaginary](0115-complex-status-sign-fidelity.md) (accepted)
- [ADR-0116: Ball kernel-status fidelity and validating serde deserialize](0116-ball-status-serde-fidelity.md) (accepted)
- [ADR-0117: conversion & parse fidelity — mode-aware parse-cap saturation, tininess-before-rounding underflow, and the sticky-flag lane](0117-conversion-saturation-and-flag-fidelity.md) (accepted)
- [ADR-0118: gamma/beta exact-dispatch precision — walk factorials at build precision, not at target](0118-gamma-beta-exact-dispatch-precision.md) (accepted)
- [ADR-0120: render a past-cap finite as an approximate magnitude `±1e{D}`, not `inf`/`0`](0120-fmt-approximate-magnitude-for-huge-finites.md) (accepted)
- [ADR-0121: tiny-x dispatch for the near-zero `x ± c·x³` trig/inverse family (sin, tan, asin)](0121-trig-tiny-x-dispatch-near-zero-family.md) (accepted)
- [ADR-0122: pow — mirror the mode under odd-parity negation, and dispatch the power-of-two exact set](0122-pow-mirror-mode-and-power-of-two-exact.md) (accepted)
