#!/usr/bin/env python3
"""Long-lived complex-Arb (`acb`) oracle worker for the pfloat-complex C5
verification pass (ADR-0092).

The worker reads requests from stdin (one per line) and writes one response
line each. It computes a complex elementary function (or a complex arithmetic
op) in python-flint's rigorous `acb` ball arithmetic at a requested working
precision, then emits the EXACT dyadic enclosure of EACH component (real and
imaginary) separately, so the Rust side can run an independent componentwise
certified-rounding check against pfloat-complex's output.

This is the complex analogue of `scripts/arb_oracle_worker.py`'s `BRACKET`
verb. It exists as a separate worker because complex results are componentwise:
one request returns two component brackets, and Arb has no signed zero (so the
signed-zero / special-value branch rows are pinned by the Rust enumerated
tables, not here -- this worker certifies only the NUMERIC value of finite,
nonzero components).

Protocol
--------

Request (one per line)::

    CBRACKET <fn_id> <oracle_prec> <re_s> <re_m> <re_e> <im_s> <im_m> <im_e> \
        [<re2_s> <re2_m> <re2_e> <im2_s> <im2_m> <im2_e>]

where:

- ``fn_id`` is one of the unary ``csqrt`` / ``cexp`` / ``clog`` or the binary
  ``cadd`` / ``csub`` / ``cmul`` / ``cdiv`` (the second operand triples are
  present iff the op is binary).
- ``oracle_prec`` is the Arb working precision in bits.
- each operand component is the EXACT dyadic ``sign * mantissa * 2^exp``:
  ``sign`` is ``+`` / ``-``, ``mantissa`` is a lowercase hex integer, ``exp`` a
  signed decimal. The construction is exact (an integer times a power of two,
  both exact in Arb), so no decimal crosses the boundary.

The worker also answers ``ready?`` with ``OK ready`` for the startup handshake.

Response (one per line)::

    OK <re_component> <im_component>

where each ``<component>`` is one of:

- ``F <lo_s> <lo_m> <lo_e> <hi_s> <hi_m> <hi_e>`` -- a finite enclosure
  ``lo <= value <= hi`` as two exact dyadic triples;
- ``N`` -- the component is NaN;
- ``P`` / ``M`` -- the component is entirely +inf / -inf;
- ``Q`` -- the component's sign is indeterminate at this precision (the ball
  straddles a value the Rust side cannot certify).

or ``ERR <message>`` on a malformed request or evaluation error.

LGPL isolation
--------------

FLINT and Arb are LGPL. The worker is an out-of-process Python subprocess
driven by the pfloat-complex test harness; FLINT/Arb never enter the shipped
Rust crate's link graph (ADR-0034 + ADR-0035 posture, reused here). The venv
that hosts python-flint lives outside the repo (`scripts/setup_arb_oracle.sh`).
"""

import sys

from flint import acb, arb, ctx  # noqa: E402

_BINARY = frozenset(("cadd", "csub", "cmul", "cdiv"))


def arb_from_dyadic(sign_str: str, man_hex: str, exp_str: str) -> arb:
    """Lift an exact dyadic ``sign * mantissa * 2^exp`` to an Arb point.

    Exact: an integer mantissa times a power-of-two scale, both of which Arb
    represents without rounding (the input mantissa is at most the pfloat input
    precision, well under ``ctx.prec``)."""
    man = int(man_hex, 16)
    if sign_str == "-":
        man = -man
    exp = int(exp_str)
    if exp >= 0:
        return arb(man) * arb(2) ** exp
    return arb(man) / arb(2) ** (-exp)


def dispatch(fn_id: str, x: acb, y) -> acb:
    """Compute the requested complex function in rigorous acb ball arithmetic.

    csqrt / clog are the principal branches Arb implements, matching C99 Annex
    G off the cut; the harness keeps inputs off the negative-real cut so the
    branch choice is unambiguous."""
    if fn_id == "csqrt":
        return x.sqrt()
    if fn_id == "cexp":
        return x.exp()
    if fn_id == "clog":
        return x.log()
    if fn_id == "cadd":
        return x + y
    if fn_id == "csub":
        return x - y
    if fn_id == "cmul":
        return x * y
    if fn_id == "cdiv":
        return x / y
    raise ValueError(f"unknown fn_id: {fn_id}")


def _bound_to_dyadic(bound) -> str:
    """Exact ``<sign> <abs_mantissa_hex> <exp>`` of a finite Arb ball endpoint,
    via its ``man_exp`` (mantissa times a power of two). Reconstructed
    bit-exactly on the Rust side."""
    man, exp = bound.man_exp()
    man_int = int(man)
    exp_int = int(exp)
    if man_int == 0:
        return "+ 0 0"
    sign = "-" if man_int < 0 else "+"
    return f"{sign} {format(abs(man_int), 'x')} {exp_int}"


def component(b: arb) -> str:
    """Encode one Arb component (real or imaginary part) as a response token
    group: a finite dyadic enclosure, or a non-finite verdict."""
    if b.is_nan():
        return "N"
    if not b.is_finite():
        if b > 0:
            return "P"
        if b < 0:
            return "M"
        return "Q"
    lo = b.lower()
    hi = b.upper()
    return f"F {_bound_to_dyadic(lo)} {_bound_to_dyadic(hi)}"


def handle_cbracket(args: list) -> str:
    if not args:
        return "ERR CBRACKET: missing fn_id"
    fn_id = args[0]
    binary = fn_id in _BINARY
    expected = 2 + 6 + (6 if binary else 0)  # fn_id, prec, z (6), [w (6)]
    if len(args) != expected:
        return f"ERR CBRACKET {fn_id}: expected {expected} args, got {len(args)}"
    try:
        oracle_prec = int(args[1])
    except ValueError:
        return f"ERR CBRACKET malformed oracle_prec: {args[1]}"
    if not 1 <= oracle_prec <= 1 << 20:
        return f"ERR CBRACKET oracle_prec out of range: {oracle_prec}"

    ctx.prec = oracle_prec
    try:
        zre = arb_from_dyadic(args[2], args[3], args[4])
        zim = arb_from_dyadic(args[5], args[6], args[7])
        x = acb(zre, zim)
        y = None
        if binary:
            wre = arb_from_dyadic(args[8], args[9], args[10])
            wim = arb_from_dyadic(args[11], args[12], args[13])
            y = acb(wre, wim)
        result = dispatch(fn_id, x, y)
    except Exception as e:  # noqa: BLE001 -- report any Arb error to the caller
        return f"ERR CBRACKET {type(e).__name__}: {e}"

    try:
        re_tok = component(result.real)
        im_tok = component(result.imag)
    except Exception as e:  # noqa: BLE001
        return f"ERR CBRACKET component extraction: {type(e).__name__}: {e}"
    return f"OK {re_tok} {im_tok}"


def handle_request(line: str) -> str:
    line = line.strip()
    if line == "ready?":
        return "OK ready"
    parts = line.split()
    if parts and parts[0] == "CBRACKET":
        return handle_cbracket(parts[1:])
    return f"ERR unknown verb: {parts[0] if parts else '(empty)'}"


def main() -> None:
    """Read requests from stdin until EOF, flushing each response so the Rust
    side sees it immediately."""
    for line in sys.stdin:
        sys.stdout.write(handle_request(line) + "\n")
        sys.stdout.flush()


if __name__ == "__main__":
    main()
