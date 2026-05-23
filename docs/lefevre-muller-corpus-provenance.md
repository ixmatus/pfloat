# Hard-to-round-case corpus: provenance

This document records the provenance for pfloat's hard-to-round-case
test corpus, shipped as a differential tier under
`tests/differential_lefevre_muller.rs`. The file name honours the
historical Lefèvre-Muller lineage that the corpus descends from; the
direct upstream source is the modern CORE-MATH project, which
maintains the canonical correctly-rounded-binary64 test data and
itself cites the Lefèvre database for its origin.

## Direct source

- **Project.** CORE-MATH, an Inria-hosted open-source effort to
  provide correctly-rounded elementary functions for binary32,
  binary64, binary80, and binary128 IEEE 754 formats.
  - Homepage: <https://core-math.gitlabpages.inria.fr/>
  - Repository: <https://gitlab.inria.fr/core-math/core-math/>
  - Survey of the project's history, scope, and maintainers in
    Sibidanov and Zimmermann, *CORE-MATH: progress report*,
    FPBench 2023.
- **License.** MIT, declared per-source-file (CORE-MATH carries no
  top-level `LICENSE` file; every C source file in the repository
  begins with the MIT permission notice and a copyright line naming
  the contributor(s) responsible for that file). Copyright holders
  across the binary64 surface include Alexei Sibidanov, Paul
  Zimmermann, Tom Hubrecht, Stéphane Glondu, and Claude-Pierre
  Jeannerod, with Inria and CERN as institutional copyright holders
  on jointly authored files. Years span 2021 through 2026.
- **Data files.** The hard-to-round cases live in
  `src/binary64/<function>/<function>.wc` in the CORE-MATH repository.
  Each line is one binary64 input encoded as a C99 hex-float literal
  (`0x1.62e42fefa39efp+9`). Comment lines (`# ...`) annotate each
  block with the BaCSeL search parameters that produced it and the
  resulting minimum "identical bits after the round bit" guarantee
  for inputs in that block.

## Historical lineage

CORE-MATH's hard-to-round-case data is the modern endpoint of a
multi-decade research program. The relevant prior art that pfloat
acknowledges through CORE-MATH:

- **The foundational paper.** Vincent Lefèvre and Jean-Michel Muller,
  *Worst Cases for Correct Rounding of the Elementary Functions in
  Double Precision*, ARITH-15 (2001). Preprint as INRIA Research
  Report RR2000-35 at <https://inria.hal.science/inria-00072594/>.
  The paper introduces the search method used to find binary64
  worst cases and frames the resulting tabulated values as
  "good test cases for checking whether a library provides correct
  rounding or not" — the original library-testing intent that
  CORE-MATH's `.wc` files inherit and extend.
- **The 2000 snapshot.**
  <https://perso.ens-lyon.fr/jean-michel.muller/TMD.html>
  (last modified September 5, 2000) is the historical web page
  the ARITH-15 paper cites for "more examples." It is now a
  superseded snapshot; CORE-MATH covers a strictly larger surface
  with current data.
- **The maintained database.** Vincent Lefèvre maintains an
  updated worst-case archive at
  <https://www.vinc17.net/research/testlibm/> (the modern successor
  to the 2000 TMD page). CORE-MATH's per-function `.wc` files cite
  this database directly; for example, `src/binary64/sin/sin.wc`'s
  header reads `# worst cases from
  https://www.vinc17.net/research/testlibm/testlibm-data.xz (update
  from 2020-11-27), from 0 to pi, with 46 to 59 identical bits
  after the round bit`. CORE-MATH's role is to curate Lefèvre's
  database into per-function files, augment with CORE-MATH's own
  BaCSeL searches over additional input regions, and add
  CORE-MATH-derived non-regression cases.
- **The Stehlé-Lefèvre-Zimmermann algorithm.** Damien Stehlé,
  Vincent Lefèvre, and Paul Zimmermann's LLL-based search
  algorithm extends the original Lefèvre method to wider precisions.
  It is the basis for both the maintained vinc17 database and
  CORE-MATH's BaCSeL tool.
- **The textbook.** Jean-Michel Muller et al., *Handbook of
  Floating-Point Arithmetic*, 2nd edition (Birkhäuser, 2018) is
  the consolidated print reference for the binary64 worst cases.
  pfloat does not transcribe from the Handbook (CORE-MATH is the
  living source); the Handbook is cited so a future reader has
  the bibliographic anchor.

The lineage is L-M (2001 algorithm + paper) → vinc17 maintained
database → CORE-MATH per-function `.wc` files → pfloat subset.
pfloat takes from the last node directly and acknowledges every
prior node through this document and the file header below.

## What pfloat ships

The full CORE-MATH `.wc` corpus is on the order of 100 MB
uncompressed (e.g. `exp.wc` is 25 MB, ~1.1 million cases). Shipping
the full corpus verbatim is neither warranted as a test surface nor
proportionate to a Rust crate's size budget. pfloat ships a
representative subset:

- One Rust data block per covered function (approximately 50 cases
  per function, sampled from the CORE-MATH `.wc` blocks that
  exercise the hardest rounding decisions).
- The pfloat-supported binary64 surface that has a CORE-MATH analog
  and a passing pfloat kernel at `p = 53` `NearestEven`: `exp`, `ln`
  (CORE-MATH `log`), `sin`, `cos`, `tan`, `atan`, `asin`, `acos`,
  `exp2`, `exp10`, `expm1`, `log1p`, `sinh`, `cosh`, `asinh`,
  `acosh`, `atanh`, `erf`, `erfc`, `gamma` (CORE-MATH `tgamma`).
  Twenty functions as of slice p1.1; the original nine landed in
  slice 8b, eleven more in slice p1.1.
- Four functions whose CORE-MATH `.wc` data is available but whose
  pfloat kernel surfaces a `has-errors` finding at `p = 53` `NE` on
  at least one hard-to-round input are deferred to slice p1.2:
  `log2`, `log10`, `tanh`, `lgamma`. The pattern across all four is
  a 1-ULP miss at `p = 53` traceable to the elementary kernels'
  fixed-64-bit-guard convention (slice 3a, predating the Ziv retry
  the Phase 1 plan upgrades). The corpus addition for each of these
  is staged behind its kernel fix; the bead queue carries the
  pairing.
- Multi-argument functions in CORE-MATH (`atan2`, `pow`) and
  reciprocal/root primitives ADR-0032 reserves for the Phase 2 libm
  shell (`cbrt`, `hypot`, `rootn`) are intentionally absent from
  this corpus. They stay on differential + property tests as their
  v1.0 verification posture; the post-v1.0 rigor track is the
  multi-argument exhaustive sweep, separate from the Phase 1 unary
  surface.
- Inputs are transcribed from CORE-MATH `.wc` files. Outputs are
  NOT transcribed — pfloat computes each expected binary64 result
  independently by evaluating the function at `p = 200` bits via
  `mpmath` and rounding to binary64 under `NearestEven`. The
  differential test then asserts pfloat's own kernel rounds to the
  same binary64 value.
- The subset and its mpmath-computed expected outputs are committed
  as a static Rust constant in the test file; the
  CORE-MATH-original `.wc` files are not vendored, only sampled.

## What pfloat does not ship

- The full multi-MB `.wc` data files. pfloat does not redistribute
  CORE-MATH's data files verbatim; it derives a documented subset.
- The CORE-MATH C source code. pfloat is a pure-Rust crate; the
  CORE-MATH C implementation is the test oracle's algorithmic
  inspiration, not vendored code.
- Any CORE-MATH function's correctly-rounded output verbatim. The
  expected outputs in pfloat's test corpus are computed by pfloat's
  own oracle (mpmath at p=200, rounded to binary64), not transcribed
  from CORE-MATH's `--worst` mode output.

## License compliance

MIT requires the copyright notice and the permission notice to be
preserved in copies and substantial portions of the software. pfloat
satisfies this by:

1. Including the CORE-MATH MIT permission notice verbatim in the
   header of `tests/differential_lefevre_muller.rs` (the file that
   embeds the derived subset).
2. Naming the CORE-MATH project as the source and citing the
   repository URL so individual per-file copyrights remain accessible
   to readers of the derivative.
3. Listing the principal copyright holders by name (Alexei Sibidanov,
   Paul Zimmermann, Tom Hubrecht, Stéphane Glondu, Claude-Pierre
   Jeannerod) and the institutional holders (Inria, CERN) in the
   derivative's header. The exhaustive enumeration of which holder
   covers which function is delegated to the upstream repository,
   which is the canonical attribution surface.
4. Reproducing the MIT permission notice in this provenance document
   (above) so the license terms are visible without depending on the
   upstream remaining reachable.

## Attribution header for derived files

Every file derived from this corpus carries the following header
verbatim (the four-paragraph version is required at the head of
`tests/differential_lefevre_muller.rs`; an abbreviated single-line
citation suffices for cross-references in module documentation):

```text
// Hard-to-round cases for IEEE 754 binary64 elementary functions,
// transcribed from the CORE-MATH project at
// https://gitlab.inria.fr/core-math/core-math/ under the MIT
// license. Principal copyright holders for the source files this
// subset derives from: Alexei Sibidanov, Paul Zimmermann, Tom
// Hubrecht, Stéphane Glondu, and Claude-Pierre Jeannerod, with
// Inria and CERN as institutional copyright holders on jointly
// authored files; copyright years 2021 through 2026. See the
// upstream repository for per-function attribution and
// docs/lefevre-muller-corpus-provenance.md for the full provenance
// chain.
//
// Permission is hereby granted, free of charge, to any person
// obtaining a copy of this software and associated documentation
// files (the "Software"), to deal in the Software without
// restriction, including without limitation the rights to use,
// copy, modify, merge, publish, distribute, sublicense, and/or
// sell copies of the Software, and to permit persons to whom the
// Software is furnished to do so, subject to the following
// conditions:
//
// The above copyright notice and this permission notice shall be
// included in all copies or substantial portions of the Software.
//
// THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND,
// EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES
// OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND
// NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT
// HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY,
// WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING
// FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR
// OTHER DEALINGS IN THE SOFTWARE.
//
// Inputs are transcribed from CORE-MATH .wc files. Outputs are
// computed independently by pfloat's mpmath oracle at p=200 bits;
// they are not transcribed from CORE-MATH and represent pfloat's
// own verification rather than upstream's claim.
```

## If the upstream position changes

If a future maintainer learns CORE-MATH's license posture has
changed (a relicense, a more restrictive grant, a written objection
from a copyright holder), the appropriate response is to remove
the corpus, fall back on pfloat-derived hard-to-round cases (a
randomized or exhaustive-in-small-range search using pfloat's own
kernels and mpmath as oracle), and update this document with the
event date and reason. Nothing in pfloat's design depends
structurally on the CORE-MATH corpus being present; it is one
tier among several.
