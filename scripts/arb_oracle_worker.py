#!/usr/bin/env python3
"""Long-lived Arb oracle worker for the pfloat Phase 1 Oracle harness
under the ADR-0035 protocol.

The worker reads requests from stdin (one per line) and writes
responses to stdout (one per line). Each request specifies an
``FnId``, an optional order, an ``f32`` input bit pattern, and a
rounding mode; each response is the certified ``f32`` bit pattern
(or ``INC`` if the worker's Ziv loop cannot certify even at the
maximum Arb precision).

ADR-0035 replaces the slice p1.5 decimal-bridge protocol with the
worker-reports-certified-f32-directly architecture:

- The worker decodes the f32 input via exact integer mantissa and
  power-of-two scale (see ``arb_from_f32_bits``); no Python ``repr``
  truncation. This closes the input-side precision-loss class that
  caused the 14030 silent BesselI1 mismatches in slice p1.5.
- The worker runs the Ziv-at-oracle loop in-process: at increasing
  Arb working precision, compute the function's ball, extract the
  exact rational lower/upper bounds, hand off to the shared
  ``certified_round_f32`` routine. No decimal-bridge crosses the
  subprocess boundary. This closes the bracket-collapse class that
  caused the wrong f32 to silently certify at low Ziv precision.
- The shared ``certified_round_f32`` routine
  (``scripts/oracle_workers/certified_rounding.py``) is the
  load-bearing piece, library-agnostic and exhaustively
  property-tested in isolation.

Protocol
--------

Request (one per line)::

    <fn_id> <order_or_dash> <input_bits_hex> <mode>

where:

- ``fn_id`` is one of ``si``, ``ci``, ``li``, ``bi``,
  ``ai_prime``, ``bi_prime``, ``i`` (Bessel I, with order), or
  ``k`` (Bessel K, with order).
- ``order_or_dash`` is ``-`` for non-parametric variants and the
  integer order for Bessel ``i`` / ``k``.
- ``input_bits_hex`` is the binary32 input as an 8-character
  lowercase hex string of the f32 bit pattern (the natural
  human-readable big-endian representation: ``3f800000`` = 1.0).
- ``mode`` is one of ``NE``, ``RNA``, ``RZ``, ``RP``, ``RM``
  (IEEE 754 rounding modes).

The worker recognises one extra request, ``ready?``, which it
answers with ``OK ready`` so the Rust side can confirm the worker
has started.

Response (one per line)::

    OK <f32_bits_hex>

(the certified f32 bit pattern as 8 lowercase hex chars), or::

    INC

(the worker's Ziv loop reached its maximum precision without
certifying a unique f32 — the bracket straddles a rounding
boundary at every attempted precision), or::

    ERR <message>

(an error occurred while processing the request).

Ziv loop
--------

Internal precision sequence: ``64, 128, 256, 512, 1024, 2048,
4096, 8192``. At each precision the worker computes the function
in Arb's ball arithmetic, extracts exact rational bounds via
``arb.lower().fmpq()`` and ``arb.upper().fmpq()``, and calls
``certified_round_f32(lower, upper, mode)``. If the routine
returns ``Some(f32_bits)`` the worker emits ``OK <bits>``; if it
returns ``None`` the worker doubles precision and retries. The
8192-bit cap is well above the slice p1.5 Rust-side cap of 1024
because the worker pays only the in-process ball-arithmetic cost
(no decimal-bridge cost), so much higher precision is cheap.

LGPL isolation
--------------

FLINT and Arb are LGPL. The worker is an out-of-process Python
subprocess driven by the pfloat oracle harness; FLINT and Arb
never enter the shipped Rust crate's link graph (ADR-0034 +
ADR-0035). The venv that hosts ``python-flint`` lives outside the
pfloat repo (per ``scripts/setup_arb_oracle.sh``).
"""

import os
import sys
from fractions import Fraction
from typing import Optional

# Import the shared certified-rounding routine. The
# ``scripts/oracle_workers/`` package is sibling to this script
# (under ``scripts/``); we add the package directory to sys.path
# explicitly so the import works regardless of how the worker is
# invoked.
_SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(_SCRIPT_DIR, "oracle_workers"))

from certified_rounding import certified_round_f32  # noqa: E402

from flint import arb, ctx  # noqa: E402


# Ziv-at-oracle internal precision sequence and cap.
ZIV_START_PREC = 64
ZIV_MAX_PREC = 8192


def arb_from_f32_bits(bits: int) -> arb:
    """Lift the exact f32 value with bit pattern ``bits`` to an Arb
    point at the current ``ctx.prec``. The construction goes through
    integer mantissa and a power-of-two scale, both of which Arb
    represents exactly: ``value = signed_mantissa * 2^scale_exp``
    with no rounding.

    Handles all f32 special classes: ``+0`` / ``-0`` -> exact Arb
    zero (Arb has no signed-zero distinction; the sign is the f32
    side's concern, not the worker's); ``+inf`` / ``-inf`` -> Arb
    special infinities; quiet/signaling NaN -> Arb NaN.
    """
    sign = (bits >> 31) & 1
    exp_field = (bits >> 23) & 0xFF
    mant = bits & 0x7FFFFF
    if exp_field == 0xFF:
        if mant == 0:
            return arb("-inf") if sign else arb("inf")
        return arb("nan")
    if exp_field == 0 and mant == 0:
        return arb(0)
    if exp_field == 0:
        # Subnormal: value = (-1)^sign * mant * 2^-149
        int_mant = mant
        scale_exp = -149
    else:
        # Normal: value = (-1)^sign * (1.mant) * 2^(exp_field - 127)
        #               = (mant | 0x800000) * 2^(exp_field - 150)
        int_mant = mant | 0x800000
        scale_exp = exp_field - 127 - 23
    if scale_exp >= 0:
        value = arb(int_mant) * arb(2) ** scale_exp
    else:
        value = arb(int_mant) / arb(2) ** (-scale_exp)
    return -value if sign else value


def dispatch(fn_id: str, order_or_dash: str, x: arb) -> arb:
    """Run the requested function on ``x`` and return the Arb result.

    Bessel ``i`` / ``k`` take an integer order; the other six
    variants are non-parametric and ignore ``order_or_dash`` (the
    request format sends ``-`` for them so the line tokenisation
    stays uniform)."""
    if fn_id == "si":
        return x.si()
    if fn_id == "ci":
        return x.ci()
    if fn_id == "li":
        return x.li()
    if fn_id == "bi":
        return x.airy_bi()
    if fn_id == "ai_prime":
        return x.airy()[1]
    if fn_id == "bi_prime":
        return x.airy()[3]
    if fn_id == "i":
        return x.bessel_i(arb(int(order_or_dash)))
    if fn_id == "k":
        return x.bessel_k(arb(int(order_or_dash)))
    # Reciprocal circular functions (pfloat 1.1 libm kernels, ADR-0056).
    # Arb computes these natively in rigorous ball arithmetic, so the
    # enclosure stays independent of pfloat's own reduce+reciprocate path.
    if fn_id == "cot":
        return x.cot()
    if fn_id == "sec":
        return x.sec()
    if fn_id == "csc":
        return x.csc()
    raise ValueError(f"unknown fn_id: {fn_id}")


def arb_ball_to_rational_bounds(b: arb) -> tuple[Optional[Fraction], Optional[Fraction]]:
    """Extract exact rational lower / upper bounds from an Arb ball.

    Returns ``(None, None)`` for NaN balls (caller dispatches NaN
    out-of-band; this is not a rounding question).

    For finite balls, ``b.lower()`` and ``b.upper()`` return arbs
    that are conservative bounds (lower's upper-end is the original
    ball's lower-end; upper's lower-end is the original ball's
    upper-end). Both convert to exact rationals via ``.fmpq()``,
    which represents the bound as ``p / q`` with arbitrary-precision
    integers. We pull ``p`` and ``q`` into Python ints and wrap in
    ``Fraction`` for the shared routine.

    For infinite endpoints (the upper bound is +inf or the lower is
    -inf), ``.fmpq()`` raises; the caller handles the infinity case
    separately by examining the ball's structure."""
    if b.is_nan():
        return (None, None)
    lo_arb = b.lower()
    hi_arb = b.upper()
    # Both lower() and upper() return arbs whose magnitudes are
    # finite when the ball is finite. .fmpq() works on finite-valued
    # arbs; for the inf-bound case the caller short-circuits via
    # the is_finite check upstream.
    lo_q = lo_arb.fmpq()
    hi_q = hi_arb.fmpq()
    return (
        Fraction(int(lo_q.p), int(lo_q.q)),
        Fraction(int(hi_q.p), int(hi_q.q)),
    )


# f32 bit patterns for the special values the worker returns
# directly without going through the rounding routine.
_F32_POS_ZERO = 0x0000_0000
_F32_NEG_ZERO = 0x8000_0000
_F32_POS_INF = 0x7F80_0000
_F32_NEG_INF = 0xFF80_0000
# Canonical quiet NaN bit pattern (sign bit clear, exponent all 1s,
# quiet bit (MSB of mantissa) set; per IEEE 754 the payload of an
# oracle-emitted NaN is unspecified, so we pick a canonical value).
_F32_QUIET_NAN = 0x7FC0_0000


def special_case_at_zero(fn_id: str, sign: int) -> Optional[int]:
    """Handle the f32 ``+0`` / ``-0`` input cases for functions
    whose mathematical limit at zero is a definite infinity or NaN.

    Mirrors the slice p1.5.6 worker special case (Ci(+0) = -inf,
    K_n(+0) = +inf) which slipped in to align Arb's NaN-at-limit
    behavior with pfloat's IEEE-convention mathematical-limit
    behavior. The slice p1.7 reclassification confirms these are
    the correct convention; the worker preserves the alignment.

    Returns ``None`` for cases where the standard Ziv path should
    handle the input.
    """
    # ci / k limits at zero are sign-independent (the functions are
    # defined for x > 0); preserve the +0 behavior and let -0 fall
    # through to the standard path exactly as before this signed-pole
    # extension landed.
    if fn_id == "ci":
        return _F32_NEG_INF if sign == 0 else None  # ci(+0) = -inf
    if fn_id == "k":
        return _F32_POS_INF if sign == 0 else None  # K_n(+0) = +inf
    # cot / csc are odd with a pole at zero: cot(+0) = csc(+0) = +inf,
    # cot(-0) = csc(-0) = -inf (pfloat raises DIV_BY_ZERO there). Arb
    # collapses +/-0 to a single 0, so the worker must supply the sign.
    if fn_id in ("cot", "csc"):
        return _F32_NEG_INF if sign == 1 else _F32_POS_INF
    # Other functions at +/-0: sec(0) = 1, li(0) = 0, si(0) = 0, bi(0)
    # a specific finite value, etc. The standard path computes them.
    return None


def handle_request(line: str) -> str:
    """Process one request line. Returns the response line
    (without trailing newline)."""
    line = line.strip()
    if line == "ready?":
        return "OK ready"
    parts = line.split()
    # Dispatch on the first token if it is an explicit verb. The
    # pre-pf-tqzz protocol used implicit "CERTIFY" with the 4-token
    # form (fn_id, order, input_hex, mode); backward-compatible.
    if parts and parts[0] == "MIDPOINT":
        return handle_midpoint(parts[1:])
    if parts and parts[0] == "BRACKET":
        return handle_bracket(parts[1:])
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

    # f32 +0 special case for the ci / k limit-at-zero functions.
    # Slice p1.5.6 introduced this; ADR-0035 preserves it.
    if input_bits == 0x0000_0000:
        special = special_case_at_zero(fn_id, sign=0)
        if special is not None:
            return f"OK {special:08x}"
    # f32 -0: signed pole for cot / csc (ADR-0056). Arb has no signed
    # zero, so the generic infinite-ball path below would lose the sign;
    # supply it here. ci / k return None for sign=1, so their -0 still
    # falls through to the standard path (unchanged behavior).
    if input_bits == 0x8000_0000:
        special = special_case_at_zero(fn_id, sign=1)
        if special is not None:
            return f"OK {special:08x}"

    # NaN / inf input: dispatch directly to Arb at low precision
    # and emit the appropriate f32 special. Arb propagates NaN
    # through every function; Arb at inf gives function-specific
    # behavior.
    sign = (input_bits >> 31) & 1
    exp_field = (input_bits >> 23) & 0xFF
    mant = input_bits & 0x7FFFFF
    if exp_field == 0xFF:
        if mant != 0:
            # Input NaN: output NaN.
            return f"OK {_F32_QUIET_NAN:08x}"
        # Input +/-inf: let the Ziv path compute the function's
        # behavior at infinity (most functions saturate to +/-inf
        # or 0; the rounding routine handles both).

    # Run the Ziv-at-oracle loop.
    prec = ZIV_START_PREC
    last_error = None
    while prec <= ZIV_MAX_PREC:
        ctx.prec = prec
        try:
            x_arb = arb_from_f32_bits(input_bits)
            result_ball = dispatch(fn_id, order, x_arb)
        except Exception as e:
            return f"ERR {type(e).__name__}: {e}"

        # Handle NaN ball: output NaN as the certified answer.
        if result_ball.is_nan():
            return f"OK {_F32_QUIET_NAN:08x}"

        # Handle inf ball: the ball is rigorous; if both bounds
        # are +inf, certify +inf; both -inf, certify -inf;
        # otherwise the ball straddles inf and the rounding routine
        # would need a special handling. For pfloat's domain we
        # expect either both bounds inf with matching sign or
        # finite bounds.
        if not result_ball.is_finite():
            # mid is +inf or -inf; assume the entire ball is on
            # that side. (Arb compares ball to 0 via > / < which
            # return True only when the ball lies entirely on the
            # chosen side.)
            if result_ball > 0:
                return f"OK {_F32_POS_INF:08x}"
            if result_ball < 0:
                return f"OK {_F32_NEG_INF:08x}"
            # Sign indeterminate: try a higher precision (may
            # resolve once the ball tightens).
            prec *= 2
            continue

        # Finite ball: extract exact rational bounds and certify.
        try:
            lo, hi = arb_ball_to_rational_bounds(result_ball)
        except Exception as e:
            last_error = f"{type(e).__name__}: {e}"
            prec *= 2
            continue

        if lo is None or hi is None:
            # Unexpected NaN extraction; treat as needing more
            # precision.
            prec *= 2
            continue

        certified = certified_round_f32(lo, hi, mode)
        if certified is not None:
            return f"OK {certified:08x}"

        prec *= 2

    if last_error is not None:
        return f"ERR Ziv exhausted ({ZIV_MAX_PREC} bits); last error: {last_error}"
    return "INC"


def handle_midpoint(args: list[str]) -> str:
    """Process a MIDPOINT request.

    Request shape (verb already stripped)::

        <fn_id> <order_or_dash> <input_bits_hex> <oracle_prec>

    Computes the function at ``oracle_prec`` and returns the
    midpoint of the resulting Arb ball as a signed
    mantissa-and-binary-exponent triple. The midpoint is the centre
    of the ball; the Arb ball's radius bounds how far it can be
    from the true value, and at ``oracle_prec >= working_prec + 64``
    the midpoint is accurate to well within the pf-tqzz cross-check
    tolerance (``2^(error_guard - working_prec) * |midpoint|``).

    Response shape::

        OK <sign> <mantissa_hex> <exponent>

    where ``sign`` is ``+`` or ``-``, ``mantissa_hex`` is the
    absolute integer mantissa as a lowercase hex string (no ``0x``
    prefix), and ``exponent`` is a signed decimal integer such that
    ``value = sign * mantissa * 2^exponent``. The triple is the
    exact arf representation of the ball midpoint (``arf.man_exp``).

    Errors return ``ERR <message>``; ``INC`` is returned when the
    ball is non-finite (NaN or unbounded) and the midpoint has no
    finite representation.

    pf-tqzz (slice p1g.3, ADR-0039).
    """
    if len(args) != 4:
        return f"ERR MIDPOINT: expected 4 args, got {len(args)}"
    fn_id, order, input_hex, oracle_prec_str = args
    try:
        input_bits = int(input_hex, 16)
    except ValueError:
        return f"ERR MIDPOINT malformed input_hex: {input_hex}"
    if not 0 <= input_bits <= 0xFFFF_FFFF:
        return f"ERR MIDPOINT input_bits out of u32 range: {input_hex}"
    try:
        oracle_prec = int(oracle_prec_str)
    except ValueError:
        return f"ERR MIDPOINT malformed oracle_prec: {oracle_prec_str}"
    if not 1 <= oracle_prec <= ZIV_MAX_PREC:
        return f"ERR MIDPOINT oracle_prec out of range [1, {ZIV_MAX_PREC}]: {oracle_prec}"

    ctx.prec = oracle_prec
    try:
        x_arb = arb_from_f32_bits(input_bits)
        result_ball = dispatch(fn_id, order, x_arb)
    except Exception as e:
        return f"ERR MIDPOINT {type(e).__name__}: {e}"

    if result_ball.is_nan() or not result_ball.is_finite():
        return "INC"

    # mid() returns the arf centre of the ball. arf.man_exp() returns
    # an exact (mantissa, exponent) pair such that the arf value
    # equals mantissa * 2^exponent. The mantissa fits the arf's
    # storage precision (which is at least oracle_prec when the ball
    # was computed at ctx.prec = oracle_prec).
    try:
        mid_arf = result_ball.mid()
        man, exp = mid_arf.man_exp()
        man_int = int(man)
        exp_int = int(exp)
    except Exception as e:
        return f"ERR MIDPOINT midpoint extraction failed: {type(e).__name__}: {e}"

    if man_int == 0:
        return "OK + 0 0"
    sign_str = "-" if man_int < 0 else "+"
    abs_man_hex = format(abs(man_int), "x")
    return f"OK {sign_str} {abs_man_hex} {exp_int}"


# pfloat-ball elementary functions taking two operands. All others in
# `dispatch_elementary` are unary.
_BRACKET_BINARY = frozenset(("add", "sub", "mul", "div", "atan2", "hypot"))


def arb_from_dyadic(sign_str: str, man_hex: str, exp_str: str) -> arb:
    """Lift an exact dyadic ``sign * mantissa * 2^exp`` to an Arb point
    at the current ``ctx.prec``. Exact: an integer mantissa times a
    power-of-two scale, both of which Arb represents without rounding
    (the input side of the bracket bridge, the dyadic analogue of
    ``arb_from_f32_bits``)."""
    man = int(man_hex, 16)
    if sign_str == "-":
        man = -man
    exp = int(exp_str)
    if exp >= 0:
        return arb(man) * arb(2) ** exp
    return arb(man) / arb(2) ** (-exp)


def dispatch_elementary(fn_id: str, x: arb, y: Optional[arb]) -> arb:
    """Compute a pfloat-ball elementary function in Arb's rigorous ball
    arithmetic. The enclosure stays independent of pfloat's own series /
    argument-reduction path. Base-changed functions (``exp2`` / ``exp10``
    / ``log2`` / ``log10``) compose through Arb's own ``ln(2)`` / ``ln(10)``
    constants, so the result is still a rigorous ball (the constant is
    itself a ball)."""
    if fn_id == "exp":
        return x.exp()
    if fn_id == "expm1":
        return x.expm1()
    if fn_id == "ln":
        return x.log()
    if fn_id == "log1p":
        return x.log1p()
    if fn_id == "sqrt":
        return x.sqrt()
    if fn_id == "cbrt":
        # Real (odd) cube root over the whole real line. Arb's root(3) is
        # the PRINCIPAL root: NaN for negatives and NaN at 0. pfloat-ball's
        # cbrt is the real root (cbrt(-x) = -cbrt(x), cbrt(0) = 0), so the
        # bracket must extend by sign or the lane silently skips cbrt's
        # entire negative half-domain (a NaN bracket is dropped). BRACKET
        # inputs are exact dyadic points, so the sign test is exact.
        if x < 0:
            return -((-x).root(3))
        if x > 0:
            return x.root(3)
        return arb(0)
    if fn_id == "sin":
        return x.sin()
    if fn_id == "cos":
        return x.cos()
    if fn_id == "tan":
        return x.tan()
    if fn_id == "asin":
        return x.asin()
    if fn_id == "acos":
        return x.acos()
    if fn_id == "atan":
        return x.atan()
    if fn_id == "sinh":
        return x.sinh()
    if fn_id == "cosh":
        return x.cosh()
    if fn_id == "tanh":
        return x.tanh()
    if fn_id == "asinh":
        return x.asinh()
    if fn_id == "acosh":
        return x.acosh()
    if fn_id == "atanh":
        return x.atanh()
    if fn_id == "exp2":
        return (x * arb.const_log2()).exp()
    if fn_id == "exp10":
        return (x * arb.const_log10()).exp()
    if fn_id == "log2":
        return x.log() / arb.const_log2()
    if fn_id == "log10":
        return x.log() / arb.const_log10()
    # Binary. atan2 follows pfloat's (y, x) convention: atan2(self=y,
    # other=x). hypot composes in Arb so the radius stays rigorous.
    if fn_id == "add":
        return x + y
    if fn_id == "sub":
        return x - y
    if fn_id == "mul":
        return x * y
    if fn_id == "div":
        return x / y
    if fn_id == "atan2":
        return arb.atan2(x, y)
    if fn_id == "hypot":
        return (x * x + y * y).sqrt()
    raise ValueError(f"unknown bracket fn_id: {fn_id}")


def _arb_bound_to_dyadic(bound) -> tuple:
    """Exact ``(sign_str, abs_mantissa_hex, exp)`` of a finite Arb ball
    endpoint, via its ``man_exp`` (mantissa times a power of two). The
    value is reconstructed bit-exactly on the Rust side, so no decimal
    crosses the boundary in either direction."""
    man, exp = bound.man_exp()
    man_int = int(man)
    exp_int = int(exp)
    if man_int == 0:
        return "+", "0", 0
    sign = "-" if man_int < 0 else "+"
    return sign, format(abs(man_int), "x"), exp_int


def handle_bracket(args: list) -> str:
    """Process a BRACKET request: emit the rigorous rational enclosure
    ``[lo, hi]`` of a pfloat-ball elementary function over exact dyadic
    input(s), as dyadic triples. Unlike CERTIFY this does NOT collapse
    the bracket to a rounded f32 -- the *interval* is exactly what the
    ball-containment check needs (rounding f(x) to an f32 and asking
    whether the ball contains that f32 is a false backstop). pf-fe5f.2,
    ADR-0078 follow-up.

    Request shape (verb already stripped)::

        <fn_id> <oracle_prec> <s1> <m1_hex> <e1> [<s2> <m2_hex> <e2>]

    where each operand is the exact dyadic ``sign * mantissa * 2^exp``,
    and the second operand is present iff ``fn_id`` is binary
    (add/sub/mul/div/atan2/hypot).

    Response shape::

        OK <lo_s> <lo_m_hex> <lo_e> <hi_s> <hi_m_hex> <hi_e>

    (the exact dyadic lower and upper bounds with ``lo <= f(x) <= hi``),
    or ``NAN`` (result is NaN), ``POS_INF`` / ``NEG_INF`` (the ball lies
    entirely at +/-inf, e.g. at a pole), ``INC`` (sign indeterminate at
    this precision), or ``ERR <msg>``.
    """
    if not args:
        return "ERR BRACKET: missing fn_id"
    fn_id = args[0]
    binary = fn_id in _BRACKET_BINARY
    expected = 5 + (3 if binary else 0)
    if len(args) != expected:
        return f"ERR BRACKET {fn_id}: expected {expected} args, got {len(args)}"
    try:
        oracle_prec = int(args[1])
    except ValueError:
        return f"ERR BRACKET malformed oracle_prec: {args[1]}"
    if not 1 <= oracle_prec <= ZIV_MAX_PREC:
        return f"ERR BRACKET oracle_prec out of range [1, {ZIV_MAX_PREC}]: {oracle_prec}"

    ctx.prec = oracle_prec
    try:
        x = arb_from_dyadic(args[2], args[3], args[4])
        y = arb_from_dyadic(args[5], args[6], args[7]) if binary else None
        ball = dispatch_elementary(fn_id, x, y)
    except Exception as e:
        return f"ERR BRACKET {type(e).__name__}: {e}"

    if ball.is_nan():
        return "NAN"
    if not ball.is_finite():
        if ball > 0:
            return "POS_INF"
        if ball < 0:
            return "NEG_INF"
        return "INC"

    try:
        lo_s, lo_m, lo_e = _arb_bound_to_dyadic(ball.lower())
        hi_s, hi_m, hi_e = _arb_bound_to_dyadic(ball.upper())
    except Exception as e:
        return f"ERR BRACKET bound extraction: {type(e).__name__}: {e}"
    return f"OK {lo_s} {lo_m} {lo_e} {hi_s} {hi_m} {hi_e}"


def main() -> None:
    """Read requests from stdin until EOF, write responses to stdout.
    The worker flushes stdout after every response so the Rust side
    sees the line immediately."""
    for line in sys.stdin:
        response = handle_request(line)
        sys.stdout.write(response + "\n")
        sys.stdout.flush()


if __name__ == "__main__":
    main()
