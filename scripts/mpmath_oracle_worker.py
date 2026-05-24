#!/usr/bin/env python3
"""mpmath oracle worker for the pfloat Phase 1 Oracle harness.

Second independent oracle per ADR-0035, paired with the
`scripts/arb_oracle_worker.py` Arb backend. Same wire protocol as
the Arb worker (request `<fn_id> <order_or_dash> <input_bits_hex>
<mode>`; response `OK <f32_bits_hex>` | `INC` | `ERR <message>`)
so the Rust ``MpmathOracle`` struct mirrors ``ArbOracle`` modulo
the script path.

Why a second oracle
-------------------

The Arb backend's correctness rests on FLINT/Arb's ball arithmetic
and our shared certified-rounding routine. The slice p1.7
diagnostic showed that single-oracle protocols admit a silent
defect class (the slice p1.5 sweep certified the wrong f32 for
14030 BesselI1 inputs because two layers of the decimal-bridge
protocol miscarried sub-precision information). Even with the
ADR-0035 protocol fix on the Arb side, a future silent bug in
Arb itself would re-introduce the same risk. mpmath is a totally
independent multi-precision library (pure Python, no shared code
lineage with FLINT/Arb), and three-way agreement between Arb /
mpmath / Maxima (slice p1.10) is the defense against that class.

Approach
--------

mpmath's ``iv`` interval-arithmetic context covers some functions
but lacks several attributes the special functions rely on
(``si``, ``ci``, ``li``, ``airyai``, ``airybi``, ``besseli``,
``besselk`` all fail under ``iv``). The pragmatic alternative: use
``mpmath.mp`` (point arithmetic) at high precision and MANUALLY
bracket the result with a conservative relative error bound.
mpmath's relative error for these smooth functions at precision
``prec`` is empirically below ``2^-(prec - O(log(prec)))``; we use
``|y| * 2^-(prec - 64)`` as the bracket half-width, which gives
~64 bits of safety margin above mpmath's typical accuracy.

Worker Ziv loop: start at ``prec = 256``, double until the bracket
certifies a unique f32 via the shared ``certified_round_f32``
routine, capped at ``prec = 8192`` (same cap as the Arb worker).
For BesselI1 on f32 subnormal inputs the loop typically settles
at ``prec = 512`` (mpmath at ``prec = 200`` collapses the
sub-midpoint correction to zero; ``prec = 400`` captures it).

LGPL status
-----------

mpmath is BSD-licensed (pure Python, no compiled deps in the
shipped pip package). The subprocess pattern is identical to the
Arb worker's; the license is more permissive, so no isolation
discipline beyond the existing subprocess boundary is required.
"""

import os
import sys
from fractions import Fraction
from typing import Optional

# Import the shared certified-rounding routine.
_SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(_SCRIPT_DIR, "oracle_workers"))

from certified_rounding import certified_round_f32  # noqa: E402

import mpmath  # noqa: E402


ZIV_START_PREC = 256
ZIV_MAX_PREC = 8192

# Bracket safety margin: mpmath's relative error at precision p is
# below 2^-(p - O(log p)) for smooth special functions; 64 bits of
# headroom is conservative.
BRACKET_MARGIN_HEADROOM = 64


# f32 bit patterns for the special values the worker returns directly.
_F32_POS_ZERO = 0x0000_0000
_F32_NEG_ZERO = 0x8000_0000
_F32_POS_INF = 0x7F80_0000
_F32_NEG_INF = 0xFF80_0000
_F32_QUIET_NAN = 0x7FC0_0000


def mpf_from_f32_bits(bits: int) -> "mpmath.mpf":
    """Lift the exact f32 value with bit pattern ``bits`` to an
    mpmath ``mpf`` at the current ``mpmath.mp.prec``. Constructed
    via integer mantissa + power-of-two scale, exact in mpmath's
    binary representation.
    """
    sign = (bits >> 31) & 1
    exp_field = (bits >> 23) & 0xFF
    mant = bits & 0x7FFFFF
    if exp_field == 0xFF:
        if mant == 0:
            return mpmath.mpf("-inf") if sign else mpmath.mpf("inf")
        return mpmath.mpf("nan")
    if exp_field == 0 and mant == 0:
        return mpmath.mpf(0)
    if exp_field == 0:
        int_mant = mant
        scale_exp = -149
    else:
        int_mant = mant | 0x800000
        scale_exp = exp_field - 127 - 23
    # mpmath's mpf supports exact integer construction; the
    # power-of-two scale via `mpmath.power(2, k)` is exact for
    # integer k.
    if scale_exp >= 0:
        value = mpmath.mpf(int_mant) * mpmath.mpf(2) ** scale_exp
    else:
        value = mpmath.mpf(int_mant) / mpmath.mpf(2) ** (-scale_exp)
    return -value if sign else value


def dispatch(fn_id: str, order_or_dash: str, x: "mpmath.mpf") -> "mpmath.mpf":
    """Run the requested function on ``x`` and return the mpmath result."""
    if fn_id == "si":
        return mpmath.si(x)
    if fn_id == "ci":
        return mpmath.ci(x)
    if fn_id == "li":
        return mpmath.li(x)
    if fn_id == "bi":
        return mpmath.airybi(x)
    if fn_id == "ai_prime":
        return mpmath.airyai(x, derivative=1)
    if fn_id == "bi_prime":
        return mpmath.airybi(x, derivative=1)
    if fn_id == "i":
        return mpmath.besseli(int(order_or_dash), x)
    if fn_id == "k":
        return mpmath.besselk(int(order_or_dash), x)
    raise ValueError(f"unknown fn_id: {fn_id}")


def mpf_to_rational(y: "mpmath.mpf") -> Fraction:
    """Convert an exact mpmath ``mpf`` to a ``Fraction``. mpf stores
    ``(sign, man, exp, bc)`` where the value is exactly
    ``(-1)^sign * man * 2^exp``."""
    sign, man, exp, _bc = y._mpf_
    signed_man = -man if sign else man
    if exp >= 0:
        return Fraction(signed_man * (1 << exp), 1)
    return Fraction(signed_man, 1 << (-exp))


def special_case_at_zero(fn_id: str) -> Optional[int]:
    """Match the Arb worker's slice p1.5.6 limit-at-+0 special cases
    for `ci` and `k_n` so the two oracles return identical answers
    on these convention-divergence inputs.
    """
    if fn_id == "ci":
        return _F32_NEG_INF
    if fn_id == "k":
        return _F32_POS_INF
    return None


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
    if not 0 <= input_bits <= 0xFFFF_FFFF:
        return f"ERR input_bits out of u32 range: {input_hex}"

    # f32 +0 special cases for ci / k_n (limit-at-zero convention).
    if input_bits == 0x0000_0000:
        special = special_case_at_zero(fn_id)
        if special is not None:
            return f"OK {special:08x}"

    # NaN / inf input handling.
    sign = (input_bits >> 31) & 1
    exp_field = (input_bits >> 23) & 0xFF
    mant = input_bits & 0x7FFFFF
    if exp_field == 0xFF and mant != 0:
        return f"OK {_F32_QUIET_NAN:08x}"

    # Ziv-at-oracle loop.
    prec = ZIV_START_PREC
    last_error = None
    while prec <= ZIV_MAX_PREC:
        mpmath.mp.prec = prec
        try:
            x = mpf_from_f32_bits(input_bits)
            y = dispatch(fn_id, order, x)
        except Exception as e:
            return f"ERR {type(e).__name__}: {e}"

        # NaN / inf result handling.
        if mpmath.isnan(y):
            return f"OK {_F32_QUIET_NAN:08x}"
        if mpmath.isinf(y):
            return f"OK {(_F32_NEG_INF if y < 0 else _F32_POS_INF):08x}"

        # Build bracket: y +/- |y| * 2^-(prec - BRACKET_MARGIN_HEADROOM).
        try:
            y_rat = mpf_to_rational(y)
        except Exception as e:
            last_error = f"{type(e).__name__}: {e}"
            prec *= 2
            continue

        margin_exp = prec - BRACKET_MARGIN_HEADROOM
        if margin_exp < 1:
            margin_exp = 1
        # margin = |y_rat| / 2^margin_exp (conservative absolute
        # bound on mpmath's accumulated relative error)
        if y_rat == 0:
            # mpmath returned exact zero; bracket has no width (the
            # result is exactly representable). Certify directly.
            margin = Fraction(0)
        else:
            abs_y = -y_rat if y_rat < 0 else y_rat
            margin = abs_y / Fraction(1 << margin_exp, 1)
        lo = y_rat - margin
        hi = y_rat + margin

        certified = certified_round_f32(lo, hi, mode)
        if certified is not None:
            return f"OK {certified:08x}"

        prec *= 2

    if last_error is not None:
        return f"ERR Ziv exhausted ({ZIV_MAX_PREC} bits); last error: {last_error}"
    return "INC"


def main() -> None:
    for line in sys.stdin:
        response = handle_request(line)
        sys.stdout.write(response + "\n")
        sys.stdout.flush()


if __name__ == "__main__":
    main()
