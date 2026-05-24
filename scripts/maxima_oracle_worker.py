#!/usr/bin/env python3
"""Maxima oracle worker (ADR-0035 Tier 6 sampling oracle).

Inner script driven by the `maxima_oracle_worker.sh` nix-shell
wrapper. Reads requests from stdin, invokes Maxima per request
via subprocess, parses the bigfloat result, builds an exact
rational bracket, and certifies the f32 via the shared
``certified_round_f32`` routine. Same wire protocol as the Arb /
mpmath workers.

Per-request cost is ~500ms (Maxima startup); the Tier 6 use mode
is sampling (hand-derived corpus + tie-breakers + N-sample per
release), not full f32 sweep. The wrapper's docstring covers the
function-coverage caveats (li via ei composition, bessel_i
precision floor on the smallest subnormals).

Currently slice p1.10 ships this as scaffolding: the worker is
functional for the covered FnIds but is not yet wired into a
Rust-side cross-check test. The pinned corpus
(`tests/oracle/pinned/`) is the load-bearing artifact for slice
p1.10; the Maxima cross-check gets its own slice (p1.11+) once
the function-coverage gaps are fully characterized.
"""

import os
import subprocess
import sys
from fractions import Fraction
from typing import Optional

# Import the shared certified-rounding routine.
_SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(_SCRIPT_DIR, "oracle_workers"))

from certified_rounding import certified_round_f32  # noqa: E402


ZIV_START_PREC_DIGITS = 100
ZIV_MAX_PREC_DIGITS = 2000

_F32_POS_ZERO = 0x0000_0000
_F32_NEG_ZERO = 0x8000_0000
_F32_POS_INF = 0x7F80_0000
_F32_NEG_INF = 0xFF80_0000
_F32_QUIET_NAN = 0x7FC0_0000


def f32_to_maxima_bfloat(bits: int) -> str:
    """Build the Maxima expression for the exact f32 value at
    bit pattern ``bits``. Uses integer mantissa + power-of-two
    scale; both representable exactly in Maxima's bfloat."""
    sign = (bits >> 31) & 1
    exp_field = (bits >> 23) & 0xFF
    mant = bits & 0x7FFFFF
    if exp_field == 0xFF:
        if mant == 0:
            return "-inf" if sign else "inf"
        return "und"  # Maxima's undefined; the caller short-circuits NaN
    if exp_field == 0 and mant == 0:
        return "0"
    if exp_field == 0:
        int_mant = mant
        scale_exp = -149
    else:
        int_mant = mant | 0x800000
        scale_exp = exp_field - 127 - 23
    # exp can be positive or negative; bfloat(2)^N handles both.
    sign_prefix = "-" if sign else ""
    return f"{sign_prefix}bfloat({int_mant}) * bfloat(2)^{scale_exp}"


def maxima_call(fn_id: str, order: str, input_bits: int, prec_digits: int) -> str:
    """Build the Maxima batch command for the (fn_id, order,
    input_bits) request at the given decimal precision, returning
    the Maxima source string."""
    x_expr = f32_to_maxima_bfloat(input_bits)
    # Maxima's function names; some need composition.
    if fn_id == "si":
        # Sine integral.
        body = f"expintegral_si({x_expr})"
    elif fn_id == "ci":
        body = f"expintegral_ci({x_expr})"
    elif fn_id == "li":
        # Maxima has no direct logarithmic integral; compose via Ei.
        # li(x) = ei(ln x) for real x > 0 (and x != 1).
        body = f"expintegral_ei(log({x_expr}))"
    elif fn_id == "bi":
        body = f"airy_bi({x_expr})"
    elif fn_id == "ai_prime":
        body = f"airy_dai({x_expr})"
    elif fn_id == "bi_prime":
        body = f"airy_dbi({x_expr})"
    elif fn_id == "i":
        body = f"bessel_i({int(order)}, {x_expr})"
    elif fn_id == "k":
        body = f"bessel_k({int(order)}, {x_expr})"
    else:
        raise ValueError(f"unknown fn_id: {fn_id}")

    return (
        f"fpprec : {prec_digits}$ "
        f"display2d : false$ "
        f"y : bfloat({body})$ "
        f"disp(y)$ "
        f"quit()$"
    )


def parse_maxima_bfloat(line: str) -> Optional[Fraction]:
    """Parse a Maxima bigfloat output line (e.g.
    '7.006...b-46') into an exact Fraction.

    Returns None if the line is not parseable as a bigfloat
    (Maxima emitted an error message, a symbolic form, or a
    non-numeric result).
    """
    line = line.strip()
    if not line:
        return None
    # Bigfloats have the form `<mantissa>b<exp>` where mantissa is a
    # decimal like `7.0064...` and exp is a signed integer.
    if "b" not in line:
        return None
    parts = line.rsplit("b", 1)
    if len(parts) != 2:
        return None
    mantissa_str, exp_str = parts
    try:
        exp = int(exp_str)
    except ValueError:
        return None
    # Maxima bigfloat mantissa is decimal; convert to Fraction via
    # int-times-10^exp by walking the decimal point.
    if "." in mantissa_str:
        int_part, frac_part = mantissa_str.split(".", 1)
    else:
        int_part, frac_part = mantissa_str, ""
    sign = -1 if int_part.startswith("-") else 1
    int_part = int_part.lstrip("+-")
    try:
        mantissa_int = int(int_part + frac_part)
    except ValueError:
        return None
    decimal_shift = -len(frac_part)
    total_exp = exp + decimal_shift
    if total_exp >= 0:
        value = Fraction(sign * mantissa_int * (10**total_exp), 1)
    else:
        value = Fraction(sign * mantissa_int, 10 ** (-total_exp))
    return value


def run_maxima(source: str, timeout_sec: float = 30.0) -> str:
    """Invoke Maxima with the given source string; return stdout."""
    result = subprocess.run(
        ["maxima", "--very-quiet", "--batch-string=" + source],
        capture_output=True,
        text=True,
        timeout=timeout_sec,
    )
    if result.returncode != 0:
        raise RuntimeError(
            f"maxima exit {result.returncode}; stderr: {result.stderr.strip()}"
        )
    return result.stdout


def handle_request(line: str) -> str:
    line = line.strip()
    if line == "ready?":
        return "OK ready"
    parts = line.split()
    if len(parts) != 4:
        return f"ERR malformed request: expected 4 tokens, got {len(parts)}"
    fn_id, order, input_hex, mode = parts

    if mode not in ("NE", "RNA", "RZ", "RP", "RM"):
        return f"ERR malformed mode: {mode}"
    try:
        input_bits = int(input_hex, 16)
    except ValueError:
        return f"ERR malformed input_bits_hex: {input_hex}"

    # f32 +0 special cases match the Arb / mpmath workers.
    if input_bits == 0x0000_0000:
        if fn_id == "ci":
            return f"OK {_F32_NEG_INF:08x}"
        if fn_id == "k":
            return f"OK {_F32_POS_INF:08x}"

    sign = (input_bits >> 31) & 1
    exp_field = (input_bits >> 23) & 0xFF
    mant = input_bits & 0x7FFFFF
    if exp_field == 0xFF and mant != 0:
        return f"OK {_F32_QUIET_NAN:08x}"

    prec_digits = ZIV_START_PREC_DIGITS
    last_error = None
    while prec_digits <= ZIV_MAX_PREC_DIGITS:
        try:
            source = maxima_call(fn_id, order, input_bits, prec_digits)
            output = run_maxima(source)
        except Exception as e:
            return f"ERR {type(e).__name__}: {e}"

        # Maxima output may have multiple lines; the bfloat result
        # is the last non-empty parseable line.
        parsed = None
        for ln in reversed(output.splitlines()):
            parsed = parse_maxima_bfloat(ln)
            if parsed is not None:
                break
        if parsed is None:
            last_error = (
                f"unparseable Maxima output at fpprec={prec_digits}: "
                f"{output.strip()!r}"
            )
            prec_digits *= 2
            continue

        # Bracket: |y| * 10^-(prec_digits - 20) safety margin.
        # Maxima's bfloat accuracy is ~10^-prec_digits relative; the
        # 20-digit headroom matches the Arb worker's 64-bit headroom.
        if parsed == 0:
            margin = Fraction(0)
        else:
            abs_y = -parsed if parsed < 0 else parsed
            margin_digits = prec_digits - 20
            if margin_digits < 1:
                margin_digits = 1
            margin = abs_y / Fraction(10**margin_digits, 1)
        lo = parsed - margin
        hi = parsed + margin
        certified = certified_round_f32(lo, hi, mode)
        if certified is not None:
            return f"OK {certified:08x}"

        prec_digits *= 2

    if last_error is not None:
        return f"ERR Ziv exhausted ({ZIV_MAX_PREC_DIGITS} digits); last error: {last_error}"
    return "INC"


def main() -> None:
    for line in sys.stdin:
        response = handle_request(line)
        sys.stdout.write(response + "\n")
        sys.stdout.flush()


if __name__ == "__main__":
    main()
