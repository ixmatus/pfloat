#!/usr/bin/env python3
"""Verify the Arb worker's BRACKET verb (pf-fe5f.2, ADR-0078 follow-up).

The BRACKET verb emits the rigorous rational enclosure [lo, hi] of a
pfloat-ball elementary function over exact dyadic input(s), as dyadic
triples, WITHOUT collapsing to a rounded f32. The ball-containment lane
relies on two properties this test pins:

  1. CONTAINMENT: lo <= f(x) <= hi for the true mathematical value,
     cross-checked against a high-precision mpmath evaluation. If this
     ever fails the bracket is not rigorous and the whole backstop is
     unsound.
  2. NARROWING: the bracket width shrinks as oracle_prec rises, so the
     lane can make the enclosure far tighter than the ball under test
     (the tightness metric would be meaningless otherwise).

Run with the Arb venv python:
    ~/.cache/pfloat-arb-oracle/venv/bin/python scripts/tests/test_bracket_verb.py
"""

import os
import sys
from fractions import Fraction

import mpmath as mp

_HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.dirname(_HERE))  # scripts/

import arb_oracle_worker as w  # noqa: E402

mp.mp.prec = 800


def dyadic_tokens(man, exp):
    """(sign, abs_man_hex, exp) for the exact value man * 2^exp."""
    s = "-" if man < 0 else "+"
    return s, format(abs(man), "x"), str(exp)


def value(man, exp):
    """Exact Fraction for man * 2^exp."""
    return Fraction(man) * (Fraction(2) ** exp)


def parse_bracket(resp):
    """Parse an `OK lo_s lo_m lo_e hi_s hi_m hi_e` response into exact
    (lo, hi) Fractions, or return the sentinel string."""
    parts = resp.split()
    if parts[0] != "OK":
        return resp  # NAN / POS_INF / NEG_INF / INC / ERR ...
    assert len(parts) == 7, f"bad OK response: {resp}"
    def frac(s, m, e):
        v = Fraction(int(m, 16)) * (Fraction(2) ** int(e))
        return -v if s == "-" else v
    lo = frac(parts[1], parts[2], parts[3])
    hi = frac(parts[4], parts[5], parts[6])
    return lo, hi


def bracket(fn_id, prec, op1, op2=None):
    s1, m1, e1 = dyadic_tokens(*op1)
    line = f"BRACKET {fn_id} {prec} {s1} {m1} {e1}"
    if op2 is not None:
        s2, m2, e2 = dyadic_tokens(*op2)
        line += f" {s2} {m2} {e2}"
    return parse_bracket(w.handle_request(line))


# (fn_id, mpmath fn, op1=(man,exp), op2 or None). Inputs chosen in-domain
# and off-grid; op encodes the exact dyadic value man*2^exp.
UNARY = [
    ("exp",   mp.exp,                      (3, -1)),    # exp(1.5)
    ("expm1", mp.expm1,                   (1, -10)),   # tiny x
    ("ln",    mp.ln,                      (3, 0)),     # ln(3)
    ("log1p", lambda v: mp.log1p(v),     (1, -8)),
    ("sqrt",  mp.sqrt,                    (2, 0)),
    ("cbrt",  mp.cbrt,                    (5, 0)),
    # negative: the REAL odd root. mp.cbrt(-x) is the principal COMPLEX
    # root, so the oracle must take the real branch to match pfloat-ball.
    ("cbrt",  lambda v: mp.sign(v) * mp.cbrt(abs(v)), (-5, 0)),
    ("sin",   mp.sin,                     (1, 0)),
    ("cos",   mp.cos,                     (1, 0)),
    ("tan",   mp.tan,                     (1, 0)),
    ("asin",  mp.asin,                    (1, -1)),    # asin(0.5)
    ("acos",  mp.acos,                    (1, -1)),
    ("atan",  mp.atan,                    (7, -2)),
    ("sinh",  mp.sinh,                    (3, -1)),
    ("cosh",  mp.cosh,                    (3, -1)),
    ("tanh",  mp.tanh,                    (3, -1)),
    ("asinh", mp.asinh,                   (5, -1)),
    ("acosh", mp.acosh,                   (5, -1)),    # acosh(2.5)
    ("atanh", mp.atanh,                   (1, -2)),    # atanh(0.25)
    ("exp2",  lambda v: mp.mpf(2) ** v,  (5, -2)),
    ("exp10", lambda v: mp.mpf(10) ** v, (3, -2)),
    ("log2",  lambda v: mp.log(v, 2),    (5, 0)),
    ("log10", mp.log10,                   (7, 0)),
]

BINARY = [
    ("add",   lambda a, b: a + b,                 (3, -1), (5, -2)),
    ("sub",   lambda a, b: a - b,                 (3, -1), (5, -2)),
    ("mul",   lambda a, b: a * b,                 (3, -1), (5, -2)),
    ("div",   lambda a, b: a / b,                 (3, -1), (5, -2)),
    ("hypot", lambda a, b: mp.hypot(a, b),        (3, 0),  (4, 0)),
    ("atan2", lambda a, b: mp.atan2(a, b),        (1, 0),  (2, 0)),
]


def check_contains(fn_id, mpfn, op1, op2):
    v1 = mp.mpf(op1[0]) * mp.mpf(2) ** op1[1]
    if op2 is None:
        true = mpfn(v1)
        br = bracket(fn_id, 256, op1)
    else:
        v2 = mp.mpf(op2[0]) * mp.mpf(2) ** op2[1]
        true = mpfn(v1, v2)
        br = bracket(fn_id, 256, op1, op2)
    assert isinstance(br, tuple), f"{fn_id}: non-OK response {br!r}"
    lo, hi = br
    # Containment: lo <= true <= hi, compared at mpmath precision.
    lo_m, hi_m, true_m = mp.mpf(lo.numerator) / lo.denominator, mp.mpf(hi.numerator) / hi.denominator, true
    assert lo_m <= true_m <= hi_m, (
        f"{fn_id}: CONTAINMENT FAILED lo={float(lo_m):.17g} "
        f"true={float(true_m):.17g} hi={float(hi_m):.17g}")
    return float(hi - lo)


def check_narrows(fn_id, op1, op2):
    w64 = bracket(fn_id, 64, op1, op2)
    w256 = bracket(fn_id, 256, op1, op2)
    assert isinstance(w64, tuple) and isinstance(w256, tuple)
    width64 = w64[1] - w64[0]
    width256 = w256[1] - w256[0]
    assert width256 <= width64, f"{fn_id}: bracket did not narrow ({width256} > {width64})"


def main():
    n = 0
    for fn_id, mpfn, op1 in UNARY:
        check_contains(fn_id, mpfn, op1, None)
        check_narrows(fn_id, op1, None)
        n += 1
    for fn_id, mpfn, op1, op2 in BINARY:
        check_contains(fn_id, mpfn, op1, op2)
        check_narrows(fn_id, op1, op2)
        n += 1
    # At elementary-function poles / domain edges Arb returns a NaN ball
    # (it cannot bracket the value with a finite-or-signed-inf interval),
    # so BRACKET returns NAN. pfloat-ball uses a signed-inf / entire
    # convention there instead; reconciling the two is the S5 edge lane's
    # job (the lane skips containment on NAN and asserts the ball flags
    # INVALID / goes entire). The POS_INF / NEG_INF sentinels stay as a
    # defensive path but do not fire for the elementary surface.
    assert bracket("ln", 128, (0, 0)) == "NAN", "ln(0): Arb returns a NaN ball"
    assert bracket("atanh", 128, (1, 0)) == "NAN", "atanh(1): Arb returns a NaN ball"
    # cbrt over its full real domain: Arb's principal root(3) is NaN for
    # x <= 0, so the worker extends by the odd identity. Pin the exact-zero
    # and negative-real-root cases the unary positive sample never reaches;
    # without the extension the lane would silently skip cbrt's negative
    # half-domain (the bracket would come back NAN and be dropped).
    assert bracket("cbrt", 128, (0, 0)) == (Fraction(0), Fraction(0)), "cbrt(0) = exact 0"
    neg = bracket("cbrt", 256, (-27, 0))
    assert isinstance(neg, tuple), f"cbrt(-27): expected a finite bracket, got {neg!r}"
    nlo, nhi = neg
    assert nlo <= Fraction(-3) <= nhi, (
        f"cbrt(-27) must bracket -3, got [{float(nlo):.17g}, {float(nhi):.17g}]")
    print(f"BRACKET verb: {n} functions, containment + narrowing OK; pole NaN handling OK")


if __name__ == "__main__":
    main()
