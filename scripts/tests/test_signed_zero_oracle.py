#!/usr/bin/env python3
"""Pin the Arb worker's signed-zero handling at the origin (pf-0ja5).

Arb has no signed zero and collapses +/-0 to a single 0, so any function
that carries a signed zero through the origin must be handled by the
worker's ``special_case_at_zero`` rather than the generic Arb path, or the
lane records a false Mismatch against pfloat's (correct) signed-zero kernel.

The confirmed defect: ``si`` is odd with ``si(0) = 0``, so ``si(-0) = -0``;
the generic path returned +0, contradicting pfloat's si kernel. This test
pins ``si(+-0)`` and guards that the ci/k/cot/csc conventions are unchanged
and that non-odd zero-valued functions (li) still fall through.

Run with the Arb venv python (the module imports python-flint):
    ~/.cache/pfloat-arb-oracle/venv/bin/python scripts/tests/test_signed_zero_oracle.py
"""

import os
import sys

_HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.dirname(_HERE))  # scripts/

import arb_oracle_worker as w  # noqa: E402

_POS_ZERO = 0x0000_0000
_NEG_ZERO = 0x8000_0000
_NEG_INF = 0xFF80_0000
_POS_INF = 0x7F80_0000


def main() -> None:
    # si is odd through the origin: the input zero's sign is preserved.
    assert w.special_case_at_zero("si", 0) == _POS_ZERO, "si(+0) = +0"
    assert w.special_case_at_zero("si", 1) == _NEG_ZERO, "si(-0) = -0 (pf-0ja5)"
    # End-to-end through the request codec (the -0 input is 0x80000000).
    assert w.handle_request("si 0 80000000 NE") == "OK 80000000", "si(-0) request = -0"
    assert w.handle_request("si 0 00000000 NE") == "OK 00000000", "si(+0) request = +0"

    # Unchanged conventions: ci/k limit-at-zero (sign-independent, -0 falls
    # through), cot/csc signed pole.
    assert w.special_case_at_zero("ci", 0) == _NEG_INF, "ci(+0) = -inf"
    assert w.special_case_at_zero("ci", 1) is None, "ci(-0) falls through"
    assert w.special_case_at_zero("k", 0) == _POS_INF, "K_n(+0) = +inf"
    assert w.special_case_at_zero("cot", 1) == _NEG_INF, "cot(-0) = -inf"
    assert w.special_case_at_zero("csc", 0) == _POS_INF, "csc(+0) = +inf"

    # Non-odd zero-valued functions must NOT be given a signed zero: li(0)=0
    # but li is not odd, so its zero stays on the standard path.
    assert w.special_case_at_zero("li", 0) is None, "li(+0) falls through"
    assert w.special_case_at_zero("li", 1) is None, "li(-0) falls through"

    print("signed-zero oracle: si(+-0) signed, ci/k/cot/csc conventions intact, li pass-through OK")


if __name__ == "__main__":
    main()
