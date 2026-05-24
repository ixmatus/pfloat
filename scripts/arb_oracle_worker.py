#!/usr/bin/env python3
"""Long-lived Arb oracle worker for the pfloat Phase 1 Oracle harness.

The worker reads requests from stdin (one per line) and writes
responses to stdout (one per line). Each request specifies an
``FnId``, an optional order, an ``f32`` input bit pattern, and a
working precision in bits; each response is an Arb enclosure
``[lo, hi]`` formatted as two decimals.

Protocol
--------

Request (one per line)::

    <fn_id> <order_or_dash> <input_bits_hex> <working_prec>

where:

- ``fn_id`` is one of ``si``, ``ci``, ``li``, ``bi``,
  ``ai_prime``, ``bi_prime``, ``i`` (Bessel I, with order), or
  ``k`` (Bessel K, with order).
- ``order_or_dash`` is ``-`` for non-parametric variants
  (``si``, ``ci``, ``li``, ``bi``, ``ai_prime``, ``bi_prime``)
  and the integer order for Bessel ``i`` and ``k`` (``0``,
  ``1``, or any other ``int``).
- ``input_bits_hex`` is the binary32 input as an 8-character
  lowercase hex string (little-endian byte order matches
  ``struct.unpack('<f', ...)``).
- ``working_prec`` is the requested Arb working precision in
  bits, set on ``flint.ctx.prec`` before the evaluation.

The worker recognises one extra request, ``ready?``, which it
answers with ``OK ready`` so the Rust side can confirm the
worker has started.

Response (one per line)::

    OK <lo_decimal> <hi_decimal>

or::

    ERR <message>

where each ``*_decimal`` is either:

- a scientific-notation decimal ``<mantissa>e<exponent>`` (the
  mantissa is an integer in decimal), or
- ``nan``, ``inf``, or ``-inf`` for non-finite endpoints.

The Rust side parses each via ``rug::Float::parse`` and assembles
the ``Enclosure`` for the verifier.

Bracket construction
--------------------

For each result the worker extracts the integer enclosure via
``arb.mid_rad_10exp(n)``, which returns ``(mid_mantissa,
rad_mantissa, exp)`` such that the ball is ``[mid_mantissa
+/- rad_mantissa] * 10^exp``. The endpoints then are::

    lo = (mid - rad - 1) * 10^exp
    hi = (mid + rad + 1) * 10^exp

The ``-1`` / ``+1`` absorb any sub-LSB rounding the Rust ``Float``
parser introduces converting decimal to binary, so the resulting
``rug::Float`` enclosure rigorously contains the true value at
the verifier's working precision. ``n`` is sized to carry the
requested binary precision (``n = working_prec * 0.31 + 5``)
so the decimal conversion does not itself widen the bracket
beyond Arb's native ball radius.

When ``rad == 0`` the Arb result is an exact decimal at the
chosen ``n``: the value equals ``mid * 10^exp`` exactly, and
mid_rad_10exp's documented post-condition guarantees rad accounts
for any decimal-rounding error in the conversion. The ``+/-1``
absorbed-rounding safety is purely decorative in that case and
artificially widens the bracket from a single point to
``[-10^exp, +10^exp]``, which the verifier reports as
``OracleInconclusive`` when both endpoints round to different
``f32`` values (the bracket ``[-1, +1]`` around an exact zero
straddles every ``f32`` boundary). Slice p1.6 skips the
``+/-1`` widening when ``rad == 0`` so exact results emit a
single-point bracket: ``li(0) = 0``, ``si(0) = 0``,
``i1(0) = 0`` all certify cleanly at the verifier.

LGPL isolation
--------------

FLINT and Arb are LGPL. The worker is an out-of-process Python
subprocess driven by the pfloat oracle harness; FLINT and Arb
never enter the shipped Rust crate's link graph (ADR-0034). The
venv that hosts python-flint lives outside the pfloat repo (per
``scripts/setup_arb_oracle.sh``).
"""

import struct
import sys

from flint import arb, ctx


def f32_from_hex(hex8: str) -> float:
    """Decode an 8-character hex string as a binary32 value (returned
    as a Python ``float``; the f32 value is exactly representable in
    f64 so no precision is lost in the cast). The hex string is the
    natural human-readable big-endian representation of the f32 bit
    pattern (`3f800000` = 1.0), matching what ``format!(\"{:08x}\",
    f32::to_bits(v))`` emits on the Rust side."""
    bits = int(hex8, 16)
    return struct.unpack(">f", struct.pack(">I", bits))[0]


def arb_from_f32(v: float) -> arb:
    """Lift a Python ``float`` (carrying an exact f32 value) to an
    Arb point. NaN and infinity get the Arb special values; finite
    values go through ``arb(str(v))`` so the decimal expansion (which
    round-trips through f64 and therefore through f32) is parsed
    exactly."""
    if v != v:
        return arb("nan")
    if v == float("inf"):
        return arb("inf")
    if v == float("-inf"):
        return arb("-inf")
    return arb(repr(v))


def dispatch(fn_id: str, order_or_dash: str, x: arb) -> arb:
    """Run the requested function on ``x`` and return the Arb result.

    Bessel I and K take an integer order; the other six variants are
    non-parametric and ignore ``order_or_dash`` (the request format
    sends ``-`` for them so the line tokenisation stays uniform)."""
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
    raise ValueError(f"unknown fn_id: {fn_id}")


def format_endpoint(value: arb, lower: bool) -> str:
    """Format one endpoint of an Arb result for the wire.

    ``lower=True`` emits the lower bound of ``value``'s ball;
    ``lower=False`` emits the upper bound. Returns either a
    ``<mantissa>e<exp>`` decimal or one of ``nan``, ``inf``,
    ``-inf``."""
    if value.is_nan():
        return "nan"
    if not value.is_finite():
        # Arb represents ±inf as a ball whose midpoint is ±inf;
        # at this branch we know the value is +inf or -inf.
        # `value > 0` is robust because Arb's comparison returns
        # True only when the bracket lies entirely on the chosen
        # side of zero.
        return "inf" if value > 0 else "-inf"

    # Heuristic: request enough decimal digits to carry the current
    # working precision plus a small headroom for the +/-1
    # safety adjustment.
    n_digits = max(20, int(ctx.prec * 0.31) + 5)
    mid, rad, exp = value.mid_rad_10exp(n_digits)
    # mid and rad come back as fmpz (FLINT big integers); convert
    # through ``int`` to Python integers so arithmetic and the
    # f-string formatter behave normally.
    mid_i = int(mid)
    rad_i = int(rad)
    # Exact result: rad == 0 means the value equals mid * 10^exp
    # exactly at the chosen n_digits, so no parser-rounding safety
    # widening is required. Emit a single-point bracket
    # ``[mid, mid] * 10^exp`` instead of ``[mid-1, mid+1] * 10^exp``
    # so the verifier sees the exact value rather than a 2-ULP
    # decimal ball that may straddle an f32 boundary (e.g. the
    # ``[-1, +1]`` ball around an exact zero straddles every f32
    # boundary). Slice p1.6 closes li(0), si(0), i1(0) inconclusive.
    if rad_i == 0:
        return f"{mid_i}e{exp}"
    if lower:
        mantissa = mid_i - rad_i - 1
    else:
        mantissa = mid_i + rad_i + 1
    return f"{mantissa}e{exp}"


def handle_request(line: str) -> str:
    """Process one request line. Returns the response line (without
    trailing newline)."""
    line = line.strip()
    if line == "ready?":
        return "OK ready"
    parts = line.split()
    if len(parts) != 4:
        return f"ERR malformed request: expected 4 tokens, got {len(parts)}"
    fn_id, order, input_hex, working_prec_s = parts
    try:
        working_prec = int(working_prec_s)
    except ValueError:
        return f"ERR malformed working_prec: {working_prec_s}"
    if working_prec < 1:
        return f"ERR working_prec must be >= 1, got {working_prec}"
    try:
        x_f = f32_from_hex(input_hex)
    except ValueError:
        return f"ERR malformed input_bits_hex: {input_hex}"

    # Special case: at f32 +0, Arb returns NaN for functions whose
    # IEEE-conventional value at +0 is the mathematical limit
    # (ci(+0) = -inf; K_n(+0) = +inf for any order). pfloat follows
    # the IEEE convention, so the Arb oracle's NaN-vs-limit
    # divergence here is an oracle-side bug, not a kernel-side bug;
    # special-casing here aligns the oracle with the IEEE
    # convention. Negative zero (0x80000000) is out of the f32 sweep
    # range that starts at 0; ci(-0) and K_n(-0) are domain errors
    # for both pfloat and Arb (both return NaN) so no special case
    # is needed there.
    if input_hex == "00000000":
        if fn_id == "ci":
            return "OK -inf -inf"
        if fn_id == "k":
            return "OK inf inf"

    ctx.prec = working_prec
    try:
        x_arb = arb_from_f32(x_f)
        result = dispatch(fn_id, order, x_arb)
        lo = format_endpoint(result, lower=True)
        hi = format_endpoint(result, lower=False)
    except Exception as e:
        return f"ERR {type(e).__name__}: {e}"
    return f"OK {lo} {hi}"


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
