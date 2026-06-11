# Coverage gap reconciliation

This report reconciles the named failure mode in each crate README's
disclosure block against the coverage gaps recorded in the registry
entries. The disclosure text is read only for this program: every
clause is quoted from the committed READMEs, and any mismatch lands in
the owner review section rather than in a README edit.

Extraction date: 2026-06-11. README blobs at extraction (so future
drift is detectable): `README.md` `e69392690dd6`,
`pfloat-ball/README.md` `b8f7313bd6a5`, `pfloat-libm/README.md`
`b09e8a0a8a42`, `pfloat-complex/README.md` `efc65ed04fe4`.

The structural root, recorded in `berkeley-testfloat.md`: no third
party conformance vector set exists or can exist for the arbitrary
precision surface, so verification rests on differential oracles,
sampled corpora, one exhaustive fixed format enumeration, and property
plus Kani lanes. Every clause below is the honest residual of that
shape.

## pfloat (README.md, line 87)

| Clause (quoted) | Grounding gap entries | Verdict |
|---|---|---|
| "roundings in the wrong direction on pathological inputs" | `core-math-wc-corpus.md` (vectors are binary64, unary, sampled; directed mode coverage beyond f32 is differential, not exhaustive); `berkeley-testfloat.md` | reconciled; the class is demonstrably real (the shipped directed mode saturation defect fixed under ADR-0080 was exactly this clause) |
| "mantissa shifts that lose a bit" | `berkeley-testfloat.md` (no boundary vector set at arbitrary precision); property and Kani lanes cover the shift paths they state, no corpus targets shift boundaries per precision | reconciled |
| "special function values that drift outside the claimed error bound on inputs no test happened to cover" | `core-math-wc-corpus.md` (no corpus analog for Bessel, Airy, zeta, or the integrals); `arb.md` and `maxima.md` (sampled lanes and a finite pinned corpus, not sweeps) | reconciled |

## pfloat-libm (pfloat-libm/README.md, line 128)

| Clause (quoted) | Grounding gap entries | Verdict |
|---|---|---|
| "roundings in the wrong direction at the BigFloat to hardware float step on pathological inputs the sweep did not reach" | `core-math-wc-corpus.md` (f64 is corpus sampled at roughly fifty cases per function; the exhaustive enumeration is f32 unary only; multi argument kernels have no sampled vectors) | reconciled |
| "cancellation near the poles of the reciprocal kernels" | `core-math-wc-corpus.md` (the sampled twenty four functions exclude cot, sec, csc, so no hard to round vectors exercise the reciprocal poles at f64; the f32 surface is covered by the exhaustive sweep) | reconciled |
| "boundary cases in the f32 subnormal range" | none remaining for the verified surface: the exhaustive f32 sweep enumerates all 2^32 inputs, subnormals included, for the 25 unary functions | reconciled, with a note: this clause is now fully covered on the unary f32 surface it names; its honest residual is f64 subnormals (not enumerable) and the unverified multi argument tail |

## pfloat-ball (pfloat-ball/README.md, line 74)

| Clause (quoted) | Grounding gap entries | Verdict |
|---|---|---|
| "a directed rounding at a domain or saturation boundary that lands a hair too small" | `arb.md` (the containment lane is sampled; edge reconciliation rows exist but are finite); `berkeley-testfloat.md` | reconciled |
| "a propagation bound that holds across the sampled inputs but under-covers at an unsampled corner of the input box" | `arb.md` (point sampling was demonstrably blind to interior extrema once; the interval bracket lane closed that class but remains sampled, and interval input tightness is measured, not asserted) | reconciled |
| "a special-case value that drifts outside its enclosure" | `arb.md`, `maxima.md` (finite pinned corpus) | reconciled |

## pfloat-complex (pfloat-complex/README.md, line 76)

| Clause (quoted) | Grounding gap entries | Verdict |
|---|---|---|
| "the wrong half of a branch cut returned for an unsigned zero on inputs no signed-zero test happened to cover" | `iso-c-annex-g.md` (the Annex G rows are an enumerated finite table; the continuum approaching each cut is sampled only); `mpc.md` (the design contract is checked by a sampled acb lane) | reconciled |
| "catastrophic cancellation in the ac - bd and ad + bc cross products on inputs no random sweep lands on" | `arb.md` (the acb componentwise lane is a finite check set, 2940 checks at the 1.0 cut); no hard to round corpus exists for complex operations (`core-math-wc-corpus.md` is real and unary) | reconciled |
| "a near-zero component of log near the unit circle" | `arb.md` (sampled); the cancellation near the unit circle is exactly the class enumerated vectors cannot exhaust | reconciled |

## For owner review

No clause lacks a recorded gap, and no recorded gap contradicts a
disclosure claim. Two observations that do not rise to mismatches:

- The largest single named gap in the registry, the absence of multi
  argument hard to round vectors (atan2, pow; `core-math-wc-corpus.md`),
  is subsumed under the broad "inputs no test happened to cover"
  clauses rather than named. If a future authorized disclosure edit is
  ever undertaken, naming it would sharpen the libm failure mode.
- The libm subnormal clause predates the exhaustive f32 sweep and is
  now stronger than reality on its named surface (the sweep covers it).
  Harmless in the conservative direction; noted for the same future
  edit.
