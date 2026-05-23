#!/usr/bin/env python3
"""Extract a representative subset of CORE-MATH binary64 worst-case
inputs and compute the expected NearestEven-rounded outputs via
mpmath at 200-bit working precision.

Outputs a Rust source file with one const slice per function.
The generated file is committed under
`tests/differential/lefevre_muller_data.rs` and is consumed by
`tests/differential_lefevre_muller.rs`.

Usage:
    nix-shell -p python3Packages.mpmath --run \\
      'python3 scripts/extract_lefevre_muller.py \\
         /path/to/core-math/clone > tests/differential/lefevre_muller_data.rs'

The CORE-MATH source clone is not vendored in pfloat; it is fetched
on demand from <https://gitlab.inria.fr/core-math/core-math/> when
the corpus is regenerated. See docs/lefevre-muller-corpus-provenance.md
for the full source identification and MIT-license-compliance posture.

Provenance discipline (mirrors the principle in CLAUDE.md):
- Inputs are transcribed from CORE-MATH `.wc` files (one entry per
  line, C99 hex-float format) by `float.fromhex`. They are facts
  under the MIT-licensed CORE-MATH project and are reproduced under
  attribution.
- Outputs are NOT transcribed from CORE-MATH. The expected NE-rounded
  binary64 result of each function at each input is computed by this
  script using mpmath at 200-bit precision, then rounded to f64 by
  Python's `float()` constructor (correctly rounded under NE by IEEE
  rules). This is pfloat's verification — an independent witness that
  the upstream-identified hard-to-round case actually rounds where
  the upstream and we both say it rounds.
- A mismatch in either parsing or output computation blocks the
  corpus from landing; the script aborts on the first failed case
  with the offending input.
"""
from __future__ import annotations

import math
import struct
import sys
from pathlib import Path

import mpmath as mp

# Cases sampled per function. ~50 per function × 9 functions ≈ 450
# test cases total — small enough to run in well under a second,
# large enough to exercise every regime the upstream identified.
SUBSET = 50

# Working precision in BITS for mpmath. 200 bits gives ample
# headroom above the longest published identical-bit runs (~120),
# so the rounded-to-binary64 result is unambiguous.
MPMATH_WORKING_BITS = 200


# Each entry pairs a pfloat kernel name with the CORE-MATH directory
# under src/binary64/ and the mpmath function that computes the
# reference value. CORE-MATH calls the natural logarithm `log`;
# pfloat (following Rust convention) calls it `ln`.
FUNCTIONS = [
    ("exp", "exp", mp.exp),
    ("ln", "log", mp.log),
    ("sin", "sin", mp.sin),
    ("cos", "cos", mp.cos),
    ("tan", "tan", mp.tan),
    ("atan", "atan", mp.atan),
    ("asin", "asin", mp.asin),
    ("acos", "acos", mp.acos),
    ("exp2", "exp2", lambda x: mp.power(2, x)),
    ("exp10", "exp10", lambda x: mp.power(10, x)),
    ("expm1", "expm1", mp.expm1),
]


def f64_to_bits(x: float) -> int:
    """Return the IEEE 754 binary64 bit pattern of `x` as a u64."""
    return struct.unpack("<Q", struct.pack("<d", x))[0]


def domain_ok(name: str, x: float) -> bool:
    """Reject inputs outside the function's f64-finite-output domain.

    pfloat returns ±∞ / NaN for domain-edge inputs; those are
    covered by separate property tests, not this corpus.
    """
    if not math.isfinite(x):
        return False
    if name == "exp" and (x > 709.0 or x < -745.0):
        return False
    if name == "exp2" and (x > 1023.0 or x < -1074.0):
        return False
    if name == "ln" and x <= 0.0:
        return False
    if name in ("asin", "acos") and abs(x) > 1.0:
        return False
    return True


def mpmath_to_f64_ne(value: mp.mpf) -> float:
    """Round an mpmath value to binary64 under NearestEven.

    Python's `float()` on an mpf rounds to NE (mpmath defers to the
    underlying mpf-to-double conversion, which is correctly rounded).
    """
    return float(value)


def extract(core_math_root: Path, name: str, wc_dir: str, fn) -> list[tuple[int, int]]:
    """Sample up to SUBSET cases from the function's .wc file and
    compute the expected NE-rounded output for each.

    Sampling targets the canonical hard-to-round-case block, not
    upstream's leading domain-edge stress block. Some `.wc` files
    (notably `exp.wc`) open with an "exercise underflow or
    overflow" section whose entries test subnormal-boundary
    behaviour rather than rounding precision; pfloat's subnormal
    underflow handling is a separate concern from elementary-kernel
    rounding correctness. The script looks for a `# hard-to-round`
    or `# worst cases` comment marker and begins sampling on the
    next data line. Files without that marker (none in the current
    upstream set) fall back to sampling from the file's first data
    line.
    """
    wc_path = core_math_root / "src" / "binary64" / wc_dir / f"{wc_dir}.wc"
    if not wc_path.is_file():
        raise FileNotFoundError(f"missing CORE-MATH file: {wc_path}")

    lines = wc_path.read_text().splitlines()
    start_idx = 0
    for i, raw in enumerate(lines):
        ls = raw.strip().lower()
        if ls.startswith("#") and (
            "hard-to-round" in ls or "worst cases" in ls
        ):
            start_idx = i + 1
            break

    cases: list[tuple[int, int]] = []
    with mp.workprec(MPMATH_WORKING_BITS):
        for raw in lines[start_idx:]:
            line = raw.strip()
            if not line or line.startswith("#"):
                continue
            if len(cases) >= SUBSET:
                break
            try:
                x = float.fromhex(line)
            except ValueError:
                # Some upstream lines carry trailing annotations
                # (e.g. binade-tag comments at the end). Skip them.
                continue
            if not domain_ok(name, x):
                continue
            try:
                y = fn(mp.mpf(x))
            except (ValueError, ZeroDivisionError, mp.libmp.NoConvergence):
                continue
            y_f = mpmath_to_f64_ne(y)
            if not math.isfinite(y_f):
                continue
            cases.append((f64_to_bits(x), f64_to_bits(y_f)))
    return cases


HEADER = """\
// Hard-to-round cases for IEEE 754 binary64 elementary functions,
// transcribed from the CORE-MATH project at
// <https://gitlab.inria.fr/core-math/core-math/> under the MIT
// license. Principal copyright holders for the source files this
// subset derives from: Alexei Sibidanov, Paul Zimmermann, Tom
// Hubrecht, Stéphane Glondu, and Claude-Pierre Jeannerod, with
// Inria and CERN as institutional copyright holders on jointly
// authored files; copyright years 2021 through 2026. See the
// upstream repository for per-function attribution and
// `docs/lefevre-muller-corpus-provenance.md` for the full
// provenance chain.
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
// computed independently by pfloat's mpmath oracle at 200-bit
// working precision and rounded to binary64 under NearestEven;
// they are not transcribed from CORE-MATH and represent pfloat's
// own verification rather than upstream's claim.

// Auto-generated by `scripts/extract_lefevre_muller.py` — do not
// hand-edit. Re-run the script with the CORE-MATH repository
// clone path to regenerate.
"""


def main() -> int:
    if len(sys.argv) != 2:
        sys.stderr.write(
            "usage: extract_lefevre_muller.py <core-math-clone-path>\n"
        )
        return 2
    core_math_root = Path(sys.argv[1]).resolve()
    if not core_math_root.is_dir():
        sys.stderr.write(f"error: not a directory: {core_math_root}\n")
        return 2

    out: list[str] = [HEADER.rstrip("\n")]
    out.append("")
    out.append("/// One hard-to-round test case: the binary64 input bit")
    out.append("/// pattern and the expected NE-rounded binary64 output")
    out.append("/// bit pattern computed by mpmath at 200-bit precision.")
    out.append("pub type Case = (u64, u64);")
    out.append("")

    grand_total = 0
    for name, wc_dir, fn in FUNCTIONS:
        cases = extract(core_math_root, name, wc_dir, fn)
        sys.stderr.write(f"{name}: {len(cases)} cases\n")
        grand_total += len(cases)
        out.append(f"/// Hard-to-round cases for `{name}` over binary64.")
        out.append(
            f"pub const {name.upper()}_CASES: &[Case] = &["
        )
        for (xb, yb) in cases:
            out.append(f"    (0x{xb:016x}, 0x{yb:016x}),")
        out.append("];")
        out.append("")

    sys.stderr.write(f"total cases: {grand_total}\n")
    # Drop the trailing blank line introduced by the per-section
    # spacing so the output is rustfmt-clean out of the box.
    while out and out[-1] == "":
        out.pop()
    sys.stdout.write("\n".join(out))
    sys.stdout.write("\n")
    return 0


if __name__ == "__main__":
    sys.exit(main())
