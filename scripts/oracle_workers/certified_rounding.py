"""Shared certified-rounding routine for the ADR-0035 Oracle workers.

Contract (load-bearing)
-----------------------

::

    certified_round_f32(lower, upper, mode) -> Optional[int]

Given exact rationals ``lower <= upper`` and a rounding mode ``m``,
return the ``f32`` bit pattern (as a Python ``int`` in ``[0, 2^32)``)
iff every value in ``[lower, upper]`` rounds to the same ``f32`` under
``m``; otherwise return ``None``.

The routine is library-agnostic. It does NOT compute the function
being verified; it only rounds a known bracket of the function's true
value. The function computation is the job of each individual oracle
worker (Arb, mpmath, Maxima), which extracts ``(lower, upper)`` as
exact rationals from its native ball/interval arithmetic and calls
this routine.

The routine works entirely on ``fractions.Fraction``. No floating
point arithmetic appears in the rounding logic; no decimal-bridge
conversions occur in either direction. This is the verification
handle: the routine's correctness can be property-tested exhaustively
across the f32 boundary classes without invoking any function library.

Modes
-----

- ``"NE"``  IEEE 754 roundTiesToEven (the default for most code)
- ``"RNA"`` IEEE 754 roundTiesToAway
- ``"RZ"``  IEEE 754 roundTowardZero (truncate)
- ``"RP"``  IEEE 754 roundTowardPositive (ceiling)
- ``"RM"``  IEEE 754 roundTowardNegative (floor)

Special values
--------------

A bracket can include the rationals representing ``+0`` (the rational
``0``; ``f32`` distinguishes signed zero but the rational does not),
or extend below the negative-finite range and above the positive-finite
range. Values that round to ``+inf`` / ``-inf`` are returned with the
corresponding bit patterns. The caller is responsible for representing
``NaN`` results out-of-band (a bracket cannot include NaN; if the true
function value is NaN, the worker reports NaN directly, not via this
routine).
"""

from fractions import Fraction
from typing import Optional

# IEEE 754 binary32 layout constants.
F32_MANT_BITS = 23
F32_BIAS = 127
F32_MIN_NORMAL_EXP = -126  # smallest normal exponent (unbiased)
F32_MAX_NORMAL_EXP = 127  # largest normal exponent (unbiased)
F32_SUBNORMAL_EXP = -149  # exponent of the smallest positive subnormal

POS_ZERO_BITS = 0x0000_0000
NEG_ZERO_BITS = 0x8000_0000
POS_INF_BITS = 0x7F80_0000
NEG_INF_BITS = 0xFF80_0000

# Largest finite f32 magnitude = (2 - 2^-23) * 2^127.
# Stored as exact Fraction for boundary comparisons.
F32_MAX_FINITE = Fraction((2**24 - 1) * (1 << 104), 1)
# Smallest positive subnormal f32 = 2^-149.
F32_SMALLEST_SUBNORMAL = Fraction(1, 1 << 149)

# Direction tokens for the underlying directed-rounding primitive.
# These describe rounding direction in the SIGNED rational sense:
# ``DOWN`` rounds toward -inf, ``UP`` rounds toward +inf, etc.
DOWN = "down"
UP = "up"
ZERO = "zero"
AWAY = "away"
NEAREST_EVEN = "nearest_even"
NEAREST_AWAY = "nearest_away"


def _pack_f32(sign: int, exp_field: int, mant_field: int) -> int:
    """Pack sign | exponent field | mantissa field into a 32-bit pattern.
    All three arguments are assumed to be in range."""
    return (sign << 31) | (exp_field << 23) | mant_field


def _f32_bits_of_rational(x: Fraction, direction: str) -> int:
    """Round exact rational ``x`` to f32, returning the bit pattern.

    ``direction`` is one of ``DOWN``, ``UP``, ``ZERO``, ``AWAY``,
    ``NEAREST_EVEN``, ``NEAREST_AWAY``.

    The routine works in exact arithmetic throughout. The algorithm:

    1. Handle zero (exact).
    2. Extract sign; work with magnitude.
    3. Find the binary exponent of the magnitude (the unique ``e``
       such that ``2^e <= |x| < 2^(e+1)``).
    4. Decide normal vs subnormal regime based on ``e`` vs
       ``F32_MIN_NORMAL_EXP``.
    5. Scale the magnitude into an exact rational mantissa
       (``[2^23, 2^24)`` for normals, ``[0, 2^23)`` for subnormals).
    6. Round the rational mantissa to an integer under the requested
       direction, accounting for sign on the asymmetric modes.
    7. Handle mantissa overflow (carry into the next f32 grid step).
    8. Handle exponent overflow (saturate to ``+inf`` / ``-inf``).
    9. Pack sign, exponent field, mantissa field into the bit pattern.
    """
    if x == 0:
        # The rational ``0`` has no sign; we return ``+0``. The caller
        # that needs signed zero distinguishes via bracket structure
        # (a bracket with lower < 0 < upper resolves to ``None`` for
        # most modes, except RZ which collapses both signs to ``+0``
        # only when lower == upper == 0; signed zero discrimination
        # lives at the caller, not here).
        return POS_ZERO_BITS

    sign = 1 if x < 0 else 0
    ax = abs(x)

    # Find the binary exponent: 2^e <= ax < 2^(e+1).
    # Work with the integer log2 of the ratio numerator/denominator.
    e = _floor_log2_fraction(ax)

    # Handle magnitudes above f32_max: saturate to the appropriate inf
    # for any directed mode that would round there.
    if ax > F32_MAX_FINITE:
        # Magnitude exceeds the largest finite f32. The behavior under
        # each mode (relative to f32_max and the original sign):
        # - DOWN: round toward -inf; for positive ax, this is
        #   f32_max (not inf); for negative -ax, this is -inf.
        # - UP: symmetric.
        # - ZERO: round toward 0; for positive ax, this is f32_max;
        #   for negative -ax, this is -f32_max.
        # - AWAY: round away from 0; positive -> +inf, negative -> -inf.
        # - NEAREST_*: depends on which side of the
        #   "halfway-between f32_max and +inf" boundary the value lies.
        #   For both NE and NA, IEEE 754 specifies overflow to inf
        #   (the halfway is exactly between f32_max and 2^128, and any
        #   finite ax > f32_max + 2^103 rounds to inf).
        return _saturate_above_max(sign, ax, direction)

    # Determine regime. F32 normal values have e in [F32_MIN_NORMAL_EXP,
    # F32_MAX_NORMAL_EXP]; subnormals have e < F32_MIN_NORMAL_EXP with
    # the mantissa fixed at scale 2^-149.
    if e >= F32_MIN_NORMAL_EXP:
        # Normal regime: mantissa carries 24 bits (1 implicit + 23
        # explicit). Scale ax into [2^23, 2^24) by dividing by 2^(e-23).
        mant_scale = e - F32_MANT_BITS
    else:
        # Subnormal regime: mantissa is at fixed scale 2^-149. Below
        # 2^-149 the value rounds to either 0 or 2^-149 depending on
        # the direction; that case is handled by the scaling +
        # rounding logic uniformly (mant_int can be 0 here, which we
        # encode as either +0 or -0 based on direction-and-sign).
        mant_scale = F32_SUBNORMAL_EXP

    # mant_exact = ax * 2^(-mant_scale)
    if mant_scale >= 0:
        mant_exact = ax / Fraction(1 << mant_scale, 1)
    else:
        mant_exact = ax * Fraction(1 << (-mant_scale), 1)

    # Round mant_exact (non-negative Fraction) to integer under the
    # requested direction. Direction interpretation depends on sign:
    # DOWN of the SIGNED value = if sign=0: floor(mant_exact);
    #                            if sign=1: ceil(mant_exact)
    # (because DOWN of a negative is more negative i.e. larger
    # magnitude).
    mant_int = _round_magnitude_to_int(mant_exact, direction, sign)

    # Handle subnormal rounded to zero: the bit pattern is +0 or -0
    # depending on the original sign. This is the only place a
    # nonzero input can produce a zero bit pattern.
    if mant_int == 0 and e < F32_MIN_NORMAL_EXP:
        return NEG_ZERO_BITS if sign else POS_ZERO_BITS

    # Handle mantissa overflow due to rounding up at the top of the
    # mantissa range.
    if e >= F32_MIN_NORMAL_EXP:
        # Normal regime: mant_int should be in [2^23, 2^24].
        if mant_int >= (1 << 24):
            # Round-up overflow into the next binade. Reset to the
            # next exponent's implicit leading bit and bump e.
            mant_int = 1 << 23
            e += 1
    else:
        # Subnormal regime: mant_int in [0, 2^23]. If it overflows,
        # the value has crossed into the smallest normal (2^-126).
        if mant_int >= (1 << 23):
            return _pack_f32(sign, 1, 0)

    # Check for exponent overflow to +/- inf.
    if e > F32_MAX_NORMAL_EXP:
        return NEG_INF_BITS if sign else POS_INF_BITS

    # Pack the final bit pattern.
    if e >= F32_MIN_NORMAL_EXP:
        exp_field = e + F32_BIAS
        mant_field = mant_int - (1 << 23)
    else:
        exp_field = 0
        mant_field = mant_int

    return _pack_f32(sign, exp_field, mant_field)


def _floor_log2_fraction(x: Fraction) -> int:
    """Return the unique integer ``e`` such that ``2^e <= x < 2^(e+1)``
    for positive Fraction ``x``."""
    # Use bit_length on numerator and denominator; floor(log2(n/d)) =
    # n.bit_length() - 1 - (d.bit_length() - 1) - (1 if 2^that > n/d
    # else 0).
    n = x.numerator
    d = x.denominator
    # Initial estimate of e.
    e = n.bit_length() - d.bit_length()
    # 2^e <= n/d iff n >= 2^e * d iff n * 1 >= d << e (for e >= 0)
    # or n << -e >= d (for e < 0).
    # Verify the inequality 2^e <= x < 2^(e+1) and adjust e by at
    # most one in either direction.
    while _scaled_lt(x, e):
        e -= 1
    while _scaled_lt(x, e + 1) is False and _scaled_lt(x, e + 1) is not True:
        # Defensive: shouldn't loop; bit_length already gives us
        # within ±1 of the correct e. Break if no adjustment needed.
        break
    while not _scaled_lt(x, e + 1):
        e += 1
    return e


def _scaled_lt(x: Fraction, e: int) -> bool:
    """Return ``x < 2^e`` for positive Fraction ``x`` and integer ``e``."""
    # x < 2^e iff x.numerator < x.denominator * 2^e (for e >= 0)
    #         iff x.numerator * 2^-e < x.denominator (for e < 0)
    if e >= 0:
        return x.numerator < x.denominator << e
    else:
        return (x.numerator << (-e)) < x.denominator


def _round_magnitude_to_int(m: Fraction, direction: str, sign: int) -> int:
    """Round non-negative Fraction ``m`` to integer under direction.

    Direction interpretation depends on the original signed value's
    sign (``sign=0`` positive, ``sign=1`` negative):

    - ``DOWN`` (toward -inf): positive -> floor, negative -> ceil.
    - ``UP`` (toward +inf): positive -> ceil, negative -> floor.
    - ``ZERO`` (toward 0): always floor of magnitude (sign-symmetric).
    - ``AWAY`` (away from 0): always ceil of magnitude
      (sign-symmetric).
    - ``NEAREST_EVEN``: nearest integer; ties to even LSB
      (sign-symmetric: the LSB is the same regardless of sign).
    - ``NEAREST_AWAY``: nearest integer; ties to ceil of magnitude
      (i.e. away from zero in the signed sense; for the magnitude
      this means ceil regardless of sign).
    """
    floor = m.numerator // m.denominator
    has_frac = (m.numerator % m.denominator) != 0
    ceil = floor + (1 if has_frac else 0)

    if direction == DOWN:
        return floor if sign == 0 else ceil
    if direction == UP:
        return ceil if sign == 0 else floor
    if direction == ZERO:
        return floor
    if direction == AWAY:
        return ceil

    # Nearest modes: compare m to floor + 1/2.
    twice_remainder_num = 2 * (m.numerator - floor * m.denominator)
    # twice_remainder_num / m.denominator = 2 * (m - floor) which is
    # in [0, 2). Compare to 1 (the half-bit boundary).
    if twice_remainder_num < m.denominator:
        # Strictly below half: round down to floor.
        return floor
    if twice_remainder_num > m.denominator:
        # Strictly above half: round up to ceil.
        return ceil
    # Exactly at half: tie.
    if direction == NEAREST_EVEN:
        return floor if (floor % 2 == 0) else ceil
    if direction == NEAREST_AWAY:
        return ceil
    raise ValueError(f"unknown direction: {direction!r}")


def _saturate_above_max(sign: int, ax: Fraction, direction: str) -> int:
    """Determine the f32 bit pattern when ``|x| > f32_max`` for each
    rounding direction. The boundary is at ``(f32_max + 2^128) / 2 =
    f32_max + 2^103`` for the nearest modes; above that, NE/NA
    saturate to inf, otherwise to f32_max."""
    # f32_max bit pattern (positive): 0x7F7FFFFF.
    f32_max_bits = 0x7F7F_FFFF
    f32_max_signed = (sign << 31) | f32_max_bits
    inf_signed = NEG_INF_BITS if sign else POS_INF_BITS

    if direction == DOWN:
        # Toward -inf: positive saturates to f32_max; negative goes to -inf.
        return f32_max_signed if sign == 0 else NEG_INF_BITS
    if direction == UP:
        # Toward +inf: positive goes to +inf; negative saturates to -f32_max.
        return POS_INF_BITS if sign == 0 else f32_max_signed
    if direction == ZERO:
        return f32_max_signed
    if direction == AWAY:
        return inf_signed

    # Nearest modes: the half-bit boundary above f32_max is at
    # f32_max + 2^103 = (2^24 - 1) * 2^104 + 2^103. Above that, NE/NA
    # round to inf; below or equal (NE ties go to even = inf because
    # f32_max's LSB is odd) for NE, NA ties go away to inf.
    half_boundary = F32_MAX_FINITE + Fraction(1 << 103, 1)
    if ax > half_boundary:
        return inf_signed
    if ax < half_boundary:
        return f32_max_signed
    # Exactly at half_boundary.
    if direction == NEAREST_EVEN:
        # f32_max has mantissa 0x7FFFFF, LSB = 1 (odd).
        # The next neighbor up is +inf which is "even" in the sense
        # that the IEEE 754 spec says infinity is selected on tie.
        # Per the spec: when the magnitude rounds at the overflow
        # boundary, NE rounds to inf.
        return inf_signed
    if direction == NEAREST_AWAY:
        return inf_signed
    raise ValueError(f"unknown direction: {direction!r}")


# ---------------------------------------------------------------------------
# The certified-rounding routine.
# ---------------------------------------------------------------------------


def certified_round_f32(
    lower: Fraction, upper: Fraction, mode: str
) -> Optional[int]:
    """The contract function. Returns the f32 bit pattern every value
    in ``[lower, upper]`` rounds to under ``mode``, or ``None`` if the
    bracket straddles a rounding boundary.

    Each supported mode is monotone (non-decreasing as the input
    increases), so the check reduces to: round the lower endpoint
    under the mode, round the upper endpoint under the mode, and
    compare. If they agree, every intermediate value also rounds to
    the same f32 (by monotonicity).

    A subtle point: rounding modes like NE technically have
    non-monotonic behavior AT exact ties (NE picks even, RNA picks
    away), but as a function the result is constant on a "basin"
    around each f32 grid point and the basin boundaries align with
    the half-bit midpoints. The lower-and-upper-endpoint check works
    for NE / RNA because both endpoints landing in the same basin
    implies the entire bracket lies in that basin.

    Raises ``ValueError`` on ``lower > upper`` (invalid bracket).
    """
    if lower > upper:
        raise ValueError(f"invalid bracket: lower={lower} > upper={upper}")

    mode_to_direction = {
        "NE": NEAREST_EVEN,
        "RNA": NEAREST_AWAY,
        "RZ": ZERO,
        "RP": UP,
        "RM": DOWN,
    }
    if mode not in mode_to_direction:
        raise ValueError(f"unknown mode: {mode!r}")
    direction = mode_to_direction[mode]

    lo_bits = _f32_bits_of_rational(lower, direction)
    hi_bits = _f32_bits_of_rational(upper, direction)
    if lo_bits == hi_bits:
        return lo_bits

    # Endpoints round to different f32s. The bracket straddles a
    # rounding boundary for this mode; report None so the caller can
    # widen working precision (Ziv-at-oracle) or declare inconclusive.
    return None
