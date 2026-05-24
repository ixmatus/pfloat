#!/usr/bin/env python3
"""Tests for ``certified_rounding``.

Three test layers:

1. **Boundary-class tests**: for every f32 boundary class (exact
   normals, exact subnormals, exact midpoints, exact ties, the
   normal/subnormal transition, signed zero, the f32_max boundary,
   the overflow boundary), construct synthetic rational inputs whose
   rounding answer is hand-derived and verify the routine matches.

2. **Differential cross-check against Python f32 (struct.pack)**:
   for rational inputs that are exactly f32-representable, or for
   inputs constructed by perturbing an exact f32 by a known
   sub-ULP offset, verify the routine agrees with what Python's
   ``struct.pack('>f', float_value)`` returns under NE (the only
   mode Python's f64->f32 cast supports). This catches gross errors
   in the routine but is weak for directed modes.

3. **Property tests via random rational generation**: generate
   thousands of random rationals at varied magnitudes; for each,
   compute the "true f32 rounding" by an independent algorithm
   (find the two adjacent f32 grid points by exact integer
   arithmetic, then decide which side based on the mode's rule)
   and verify the routine matches.

4. **Certified-rounding bracket tests**: synthesize brackets at
   known positions relative to the f32 grid (entirely inside one
   bucket; straddling a midpoint; spanning multiple grid points)
   and verify ``certified_round_f32`` returns ``Some(f)`` for the
   former and ``None`` for the latter.

Standalone runner: ``python3 test_certified_rounding.py``. Exits 0
on all passing, 1 on first failure with a diagnostic line.

No third-party deps. Uses only ``fractions`` and ``struct``.
"""

import random
import struct
import sys
from fractions import Fraction

# Import the module under test from the same directory (the runner
# can be invoked as either ``python3 test_certified_rounding.py`` or
# ``python3 -m scripts.oracle_workers.test_certified_rounding``).
try:
    from certified_rounding import (
        certified_round_f32,
        _f32_bits_of_rational,
        NEAREST_EVEN,
        NEAREST_AWAY,
        DOWN,
        UP,
        ZERO,
        AWAY,
        POS_ZERO_BITS,
        NEG_ZERO_BITS,
        POS_INF_BITS,
        NEG_INF_BITS,
        F32_MAX_FINITE,
        F32_SMALLEST_SUBNORMAL,
    )
except ImportError:
    from .certified_rounding import (  # type: ignore[no-redef]
        certified_round_f32,
        _f32_bits_of_rational,
        NEAREST_EVEN,
        NEAREST_AWAY,
        DOWN,
        UP,
        ZERO,
        AWAY,
        POS_ZERO_BITS,
        NEG_ZERO_BITS,
        POS_INF_BITS,
        NEG_INF_BITS,
        F32_MAX_FINITE,
        F32_SMALLEST_SUBNORMAL,
    )


# ---------------------------------------------------------------------------
# Test infrastructure.
# ---------------------------------------------------------------------------

FAILURES = []


def check(condition, label, *details):
    """Record a failure if condition is false. Always returns the
    condition value so callers can chain checks."""
    if not condition:
        FAILURES.append((label, details))
        print(f"FAIL: {label}")
        for d in details:
            print(f"      {d}")
    return condition


def f32_to_rational(bits: int) -> Fraction:
    """Convert an f32 bit pattern to its exact Fraction value.
    NaN and inf return None (they have no rational representation;
    NaN is filtered before calling, inf is handled separately)."""
    sign = (bits >> 31) & 1
    exp_field = (bits >> 23) & 0xFF
    mant = bits & 0x7FFFFF
    if exp_field == 0xFF:
        return None  # inf or NaN; caller handles separately
    if exp_field == 0:
        # Subnormal: value = mant * 2^-149
        if mant == 0:
            return Fraction(0)
        v = Fraction(mant, 1 << 149)
    else:
        # Normal: value = (1.mant) * 2^(exp_field - 127)
        #               = (2^23 | mant) * 2^(exp_field - 150)
        mantissa_int = (1 << 23) | mant
        scale = exp_field - 127 - 23
        if scale >= 0:
            v = Fraction(mantissa_int * (1 << scale), 1)
        else:
            v = Fraction(mantissa_int, 1 << (-scale))
    return -v if sign else v


def reference_f32_round_ne(x: Fraction) -> int:
    """Independent reference implementation of f32 NE rounding.

    Algorithm: convert through Python's f64 (which IS IEEE NE for f64),
    then cast f64 to f32 via struct.pack('>f', ...) (which IS IEEE NE
    for f64->f32). This gives the correct NE f32 for ALMOST every
    rational: the only failure mode is "double rounding" where the
    f64 intermediate happens to land at an f32 midpoint when the
    rational was just off-midpoint and the second rounding ties
    incorrectly. We avoid that by checking whether the f64
    representation is exact at f32-grid precision (in which case
    no double rounding occurs) or sub-ULP (in which case we
    re-derive from the rational).

    Returns the f32 bit pattern.
    """
    # Convert via float for fast path; check if the conversion is
    # exact (Fraction-to-float is lossless when the fraction is
    # exactly representable in f64; otherwise we must check for
    # double-rounding hazard).
    f64 = float(x)
    f64_as_rational = Fraction(f64)
    if f64_as_rational == x:
        # f64 represents x exactly; struct.pack to f32 NE is correct.
        return struct.unpack(">I", struct.pack(">f", f64))[0]
    # f64 was rounded. Determine if the rounding could double-round
    # to a different f32 than direct rational->f32 NE would. The
    # safe path: re-derive from the rational using the routine's
    # SAME algorithm, but we cannot use the routine itself as the
    # reference. Instead, use a completely independent algorithm:
    # find the two f32 neighbors of x by integer mantissa arithmetic,
    # determine which is closer.
    return _independent_round_to_nearest_f32(x, ties_even=True)


def _independent_round_to_nearest_f32(x: Fraction, ties_even: bool) -> int:
    """Independent algorithm: round Fraction x to f32 NE or NA.

    1. Find the largest f32 value f_below such that f_below <= x.
    2. Find the smallest f32 value f_above such that f_above >= x.
    3. If f_below == f_above, return that. (Exact f32.)
    4. Otherwise compute |x - f_below| vs |f_above - x| as Fractions.
       Closer side wins. On tie, ties_even -> even LSB, else away from zero.
    """
    if x == 0:
        return POS_ZERO_BITS

    # Find f_below by binary search over the f32 bit patterns
    # ordered as signed magnitude. Easier: convert through float and
    # adjust by one ULP if needed.
    f64 = float(x)
    f64_rat = Fraction(f64)

    if f64_rat == x:
        # x is exactly representable in f64; check if it's exactly
        # representable in f32 too.
        f32 = struct.unpack(">f", struct.pack(">f", f64))[0]
        if Fraction(f32) == x:
            return struct.unpack(">I", struct.pack(">f", f32))[0]

    # x is not exactly an f32. Find f_below and f_above by starting
    # near struct.pack rounding and stepping by one ULP.
    nearest_f32_bits = struct.unpack(">I", struct.pack(">f", f64))[0]
    nearest_f32 = struct.unpack(">f", struct.pack(">I", nearest_f32_bits))[0]
    nearest_rat = f32_to_rational(nearest_f32_bits)
    if nearest_rat is None:
        # Nearest is inf or NaN; the routine should saturate, but
        # for the reference we treat this as "x is outside finite range".
        if x > 0:
            return POS_INF_BITS if x > F32_MAX_FINITE else nearest_f32_bits
        return NEG_INF_BITS if -x > F32_MAX_FINITE else nearest_f32_bits

    if nearest_rat == x:
        return nearest_f32_bits
    if nearest_rat < x:
        below_bits = nearest_f32_bits
        above_bits = _f32_ulp_step(nearest_f32_bits, +1)
    else:
        above_bits = nearest_f32_bits
        below_bits = _f32_ulp_step(nearest_f32_bits, -1)

    below_rat = f32_to_rational(below_bits)
    above_rat = f32_to_rational(above_bits)
    if above_rat is None or below_rat is None:
        # Hit inf boundary. Fall back to whichever is finite.
        return below_bits if above_rat is None else above_bits

    d_below = x - below_rat  # >= 0
    d_above = above_rat - x  # >= 0

    if d_below < d_above:
        return below_bits
    if d_below > d_above:
        return above_bits
    # Exact tie.
    if ties_even:
        # Pick the one whose mantissa LSB is even.
        if (below_bits & 1) == 0:
            return below_bits
        return above_bits
    else:
        # Ties away from zero: pick the one with larger magnitude.
        if abs(below_rat) > abs(above_rat):
            return below_bits
        return above_bits


def _f32_ulp_step(bits: int, direction: int) -> int:
    """Step bits by +1 or -1 ULP in the f32 sense. Handles the sign
    transition (incrementing/decrementing across zero). Direction is
    +1 (toward +inf) or -1 (toward -inf)."""
    if direction == +1:
        if bits == 0x8000_0000:
            return 0x0000_0000  # -0 -> +0
        if (bits & 0x8000_0000) == 0:
            # positive: increment toward +inf
            return bits + 1
        else:
            # negative: decrement magnitude toward -0 (which is
            # decrementing the absolute value, which is bits - 1
            # in the negative half)
            return bits - 1
    elif direction == -1:
        if bits == 0x0000_0000:
            return 0x8000_0000  # +0 -> -0
        if (bits & 0x8000_0000) == 0:
            return bits - 1
        else:
            return bits + 1
    raise ValueError(f"direction must be +1 or -1, got {direction}")


# ---------------------------------------------------------------------------
# Layer 1: boundary-class tests.
# ---------------------------------------------------------------------------


def test_exact_zero():
    """Zero (any sign) rounds to +0 under every direction."""
    for direction in (NEAREST_EVEN, NEAREST_AWAY, DOWN, UP, ZERO, AWAY):
        result = _f32_bits_of_rational(Fraction(0), direction)
        check(
            result == POS_ZERO_BITS,
            f"zero under {direction}",
            f"got {result:#010x}, expected {POS_ZERO_BITS:#010x}",
        )


def test_exact_f32_values_roundtrip():
    """Every f32 value should round to itself under every direction."""
    # Sample across normal range, subnormal range, and special points.
    sample_bits = (
        [0x0000_0001, 0x0000_0042, 0x007F_FFFF]  # subnormals
        + [0x0080_0000, 0x3F80_0000, 0x4000_0000, 0x4248_F5C3]  # normals
        + [0x7F7F_FFFF]  # f32_max
        + [0x8000_0001, 0xBF80_0000, 0xFF7F_FFFF]  # negative versions
    )
    for bits in sample_bits:
        x = f32_to_rational(bits)
        for direction in (NEAREST_EVEN, NEAREST_AWAY, DOWN, UP, ZERO, AWAY):
            result = _f32_bits_of_rational(x, direction)
            check(
                result == bits,
                f"exact f32 {bits:#010x} under {direction}",
                f"got {result:#010x}, expected {bits:#010x}",
            )


def test_subnormal_midpoint_smallest():
    """The hand-derived I1 case: 2^-150 is exactly the midpoint between
    f32 +0 and the smallest positive subnormal (mantissa 1 = 2^-149)."""
    midpoint = Fraction(1, 1 << 150)  # 2^-150
    just_above = midpoint + Fraction(1, 1 << 500)
    just_below = midpoint - Fraction(1, 1 << 500)

    # NE: midpoint exactly ties to even. Mantissa 0 has LSB 0 (even),
    # mantissa 1 has LSB 1 (odd). NE picks mantissa 0 = +0.
    check(
        _f32_bits_of_rational(midpoint, NEAREST_EVEN) == POS_ZERO_BITS,
        "NE of 2^-150 (smallest-subnormal midpoint) -> +0 (tie-to-even)",
        f"got {_f32_bits_of_rational(midpoint, NEAREST_EVEN):#010x}",
    )

    # NE: just above midpoint -> mantissa 1.
    check(
        _f32_bits_of_rational(just_above, NEAREST_EVEN) == 0x0000_0001,
        "NE of 2^-150 + epsilon -> 0x00000001",
        f"got {_f32_bits_of_rational(just_above, NEAREST_EVEN):#010x}",
    )

    # NE: just below midpoint -> mantissa 0 (= +0).
    check(
        _f32_bits_of_rational(just_below, NEAREST_EVEN) == POS_ZERO_BITS,
        "NE of 2^-150 - epsilon -> +0",
        f"got {_f32_bits_of_rational(just_below, NEAREST_EVEN):#010x}",
    )

    # RNA: midpoint -> away from zero -> mantissa 1.
    check(
        _f32_bits_of_rational(midpoint, NEAREST_AWAY) == 0x0000_0001,
        "RNA of 2^-150 (tie) -> 0x00000001 (away from zero)",
        f"got {_f32_bits_of_rational(midpoint, NEAREST_AWAY):#010x}",
    )

    # DOWN (toward -inf) of positive value: floor -> mantissa 0.
    check(
        _f32_bits_of_rational(midpoint, DOWN) == POS_ZERO_BITS,
        "DOWN of positive 2^-150 -> +0",
        f"got {_f32_bits_of_rational(midpoint, DOWN):#010x}",
    )

    # UP (toward +inf) of positive value: ceil -> mantissa 1.
    check(
        _f32_bits_of_rational(midpoint, UP) == 0x0000_0001,
        "UP of positive 2^-150 -> 0x00000001",
        f"got {_f32_bits_of_rational(midpoint, UP):#010x}",
    )


def test_i1_truth_case():
    """The literal pf-6a4e diagnostic case: I1(2^-149) = 2^-150 +
    2^-451. NE-rounds to mantissa 1 because it's slightly above the
    midpoint."""
    truth = Fraction(1, 1 << 150) + Fraction(1, 1 << 451)
    check(
        _f32_bits_of_rational(truth, NEAREST_EVEN) == 0x0000_0001,
        "NE of I1(2^-149) truth = 2^-150 + 2^-451 -> 0x00000001",
        f"got {_f32_bits_of_rational(truth, NEAREST_EVEN):#010x}",
    )


def test_normal_subnormal_transition():
    """The transition between subnormals and normals at 2^-126."""
    smallest_normal = Fraction(1, 1 << 126)  # 2^-126
    largest_subnormal = Fraction((1 << 23) - 1, 1 << 149)
    midpoint = (smallest_normal + largest_subnormal) / 2

    # 2^-126 is exactly representable as smallest normal: exp_field=1,
    # mantissa=0, sign=0 -> bits 0x00800000.
    check(
        _f32_bits_of_rational(smallest_normal, NEAREST_EVEN) == 0x0080_0000,
        "NE of 2^-126 (smallest normal) -> 0x00800000",
        f"got {_f32_bits_of_rational(smallest_normal, NEAREST_EVEN):#010x}",
    )

    # Largest subnormal: exp_field=0, mantissa=2^23-1.
    check(
        _f32_bits_of_rational(largest_subnormal, NEAREST_EVEN) == 0x007F_FFFF,
        "NE of largest subnormal -> 0x007FFFFF",
        f"got {_f32_bits_of_rational(largest_subnormal, NEAREST_EVEN):#010x}",
    )

    # Midpoint between largest subnormal and smallest normal.
    # NE tie: 0x007FFFFF has LSB 1 (odd); 0x00800000 has LSB 0 (even).
    # Tie to even -> 0x00800000.
    check(
        _f32_bits_of_rational(midpoint, NEAREST_EVEN) == 0x0080_0000,
        "NE of midpoint(largest_subnormal, smallest_normal) -> 0x00800000",
        f"got {_f32_bits_of_rational(midpoint, NEAREST_EVEN):#010x}",
    )


def test_f32_max_overflow():
    """Values above f32_max saturate to inf under appropriate modes."""
    # Value above f32_max but below the overflow boundary
    # (f32_max + 2^103): NE rounds down to f32_max, AWAY rounds to inf.
    just_above_max = F32_MAX_FINITE + Fraction(1, 1 << 100)
    f32_max_bits = 0x7F7F_FFFF
    check(
        _f32_bits_of_rational(just_above_max, NEAREST_EVEN) == f32_max_bits,
        "NE of f32_max + tiny -> f32_max",
        f"got {_f32_bits_of_rational(just_above_max, NEAREST_EVEN):#010x}",
    )
    check(
        _f32_bits_of_rational(just_above_max, AWAY) == POS_INF_BITS,
        "AWAY of f32_max + tiny -> +inf",
        f"got {_f32_bits_of_rational(just_above_max, AWAY):#010x}",
    )
    # Value above the overflow boundary: NE rounds to inf.
    well_above_max = F32_MAX_FINITE + Fraction(1 << 110, 1)
    check(
        _f32_bits_of_rational(well_above_max, NEAREST_EVEN) == POS_INF_BITS,
        "NE of f32_max + huge -> +inf",
        f"got {_f32_bits_of_rational(well_above_max, NEAREST_EVEN):#010x}",
    )


def test_negative_values():
    """Negative values round correctly under each direction; signs flip
    the DOWN/UP relationship."""
    # -2^-150 is the midpoint between -0 and the smallest negative
    # subnormal (-1 * 2^-149 = 0x80000001).
    neg_midpoint = -Fraction(1, 1 << 150)
    # NE: tie to even -> mantissa 0 (= -0 = 0x80000000).
    check(
        _f32_bits_of_rational(neg_midpoint, NEAREST_EVEN) == NEG_ZERO_BITS,
        "NE of -2^-150 -> -0 (tie-to-even)",
        f"got {_f32_bits_of_rational(neg_midpoint, NEAREST_EVEN):#010x}",
    )
    # DOWN (toward -inf) of negative: ceil of magnitude (= mantissa 1)
    # with sign 1.
    check(
        _f32_bits_of_rational(neg_midpoint, DOWN) == 0x8000_0001,
        "DOWN of -2^-150 -> -smallest-subnormal (0x80000001)",
        f"got {_f32_bits_of_rational(neg_midpoint, DOWN):#010x}",
    )
    # UP (toward +inf) of negative: floor of magnitude (= 0) with
    # sign 1 -> -0.
    check(
        _f32_bits_of_rational(neg_midpoint, UP) == NEG_ZERO_BITS,
        "UP of -2^-150 -> -0",
        f"got {_f32_bits_of_rational(neg_midpoint, UP):#010x}",
    )


# ---------------------------------------------------------------------------
# Layer 2: differential cross-check against Python f32 (NE only).
# ---------------------------------------------------------------------------


def test_ne_agrees_with_struct_pack_on_exact_f32():
    """For exact f32 inputs, NE should agree with struct.pack."""
    # Many exact f32 bit patterns; convert each to rational, round
    # under NE, expect the same bits back.
    random.seed(42)
    for _ in range(200):
        # Random f32 normal value
        bits = random.randint(0x0080_0000, 0x7F7F_FFFF)
        if random.random() < 0.5:
            bits |= 0x8000_0000  # make negative
        x = f32_to_rational(bits)
        result = _f32_bits_of_rational(x, NEAREST_EVEN)
        check(
            result == bits,
            f"exact normal {bits:#010x} NE round-trip",
            f"got {result:#010x}",
        )


def test_ne_agrees_with_struct_pack_on_perturbed_exact():
    """For x = exact_f32 + small_offset, NE should round to one of the
    two adjacent f32 grid points; compare with reference."""
    random.seed(43)
    for _ in range(200):
        bits = random.randint(0x0080_0000, 0x7F7F_FFFF)
        x_base = f32_to_rational(bits)
        # Perturb by a sub-ULP fraction.
        offset = Fraction(random.randint(-1000, 1000), 1 << 200)
        x = x_base + offset
        if x <= 0:
            continue
        result = _f32_bits_of_rational(x, NEAREST_EVEN)
        reference = _independent_round_to_nearest_f32(x, ties_even=True)
        check(
            result == reference,
            f"NE perturbed exact {bits:#010x} + {offset}",
            f"got {result:#010x}, reference {reference:#010x}, x={float(x)}",
        )


# ---------------------------------------------------------------------------
# Layer 3: random rational property tests.
# ---------------------------------------------------------------------------


def test_random_rationals_ne():
    """Random rationals at varied magnitudes; NE result must match
    the independent reference."""
    random.seed(44)
    for _ in range(500):
        # Magnitude spans normal range, subnormal range, near f32_max.
        # Generate by picking a random binary exponent and a random
        # 30-bit mantissa, then assembling.
        # Range chosen to stay strictly inside the finite f32 range
        # so the reference (which routes through Python f64 -> f32
        # struct.pack) does not overflow. Overflow handling is tested
        # separately in test_f32_max_overflow.
        exp = random.randint(-180, 90)
        mant = random.randint(1, (1 << 30) - 1)
        sign = random.choice([-1, 1])
        x = Fraction(mant * sign, 1) * (Fraction(1 << exp, 1) if exp >= 0 else Fraction(1, 1 << -exp))
        if x == 0:
            continue
        result = _f32_bits_of_rational(x, NEAREST_EVEN)
        reference = _independent_round_to_nearest_f32(x, ties_even=True)
        check(
            result == reference,
            f"random NE: sign={sign} mant={mant} exp={exp}",
            f"x={float(x)}, got {result:#010x}, reference {reference:#010x}",
        )


def test_random_rationals_rna():
    """Random rationals; RNA result must match independent reference."""
    random.seed(45)
    for _ in range(500):
        # Range chosen to stay strictly inside the finite f32 range
        # so the reference (which routes through Python f64 -> f32
        # struct.pack) does not overflow. Overflow handling is tested
        # separately in test_f32_max_overflow.
        exp = random.randint(-180, 90)
        mant = random.randint(1, (1 << 30) - 1)
        sign = random.choice([-1, 1])
        x = Fraction(mant * sign, 1) * (Fraction(1 << exp, 1) if exp >= 0 else Fraction(1, 1 << -exp))
        if x == 0:
            continue
        result = _f32_bits_of_rational(x, NEAREST_AWAY)
        reference = _independent_round_to_nearest_f32(x, ties_even=False)
        check(
            result == reference,
            f"random RNA: sign={sign} mant={mant} exp={exp}",
            f"x={float(x)}, got {result:#010x}, reference {reference:#010x}",
        )


# ---------------------------------------------------------------------------
# Layer 4: certified-rounding bracket tests.
# ---------------------------------------------------------------------------


def test_certify_bracket_inside_single_bucket():
    """A bracket entirely inside one f32 NE bucket should certify
    that bucket's f32 value."""
    # f32 normal 1.0 = 0x3F800000. The NE bucket around 1.0 is
    # (1 - ulp/2, 1 + ulp/2] where ulp = 2^-23 at this magnitude.
    one = Fraction(1)
    ulp_half = Fraction(1, 1 << 24)
    lower = one - ulp_half / 2  # well inside the bucket below 1.0
    upper = one + ulp_half / 2  # well inside the bucket above 1.0
    result = certified_round_f32(lower, upper, "NE")
    check(
        result == 0x3F80_0000,
        "NE certify bracket entirely inside 1.0's bucket",
        f"got {result if result is not None else 'None'}",
    )


def test_certify_bracket_straddles_midpoint():
    """A bracket straddling the NE midpoint between 1.0 and 1+ulp
    should return None."""
    one = Fraction(1)
    ulp = Fraction(1, 1 << 23)  # ULP at 1.0
    midpoint = one + ulp / 2
    lower = midpoint - Fraction(1, 1 << 50)
    upper = midpoint + Fraction(1, 1 << 50)
    result = certified_round_f32(lower, upper, "NE")
    check(
        result is None,
        "NE bracket straddling midpoint -> None",
        f"got {result:#010x}" if result is not None else "got None (ok)",
    )


def test_certify_bracket_exact_zero_to_subnormal():
    """The I1(2^-149) case: a tight bracket around 2^-150 + 2^-451
    should certify mantissa 1 because both endpoints round to 1
    (the truth is above the midpoint, and any sufficiently tight
    bracket around the truth lies entirely above the midpoint)."""
    truth = Fraction(1, 1 << 150) + Fraction(1, 1 << 451)
    rad = Fraction(1, 1 << 500)
    lower = truth - rad
    upper = truth + rad
    result = certified_round_f32(lower, upper, "NE")
    check(
        result == 0x0000_0001,
        "NE certify tight bracket around I1(2^-149) truth -> 0x00000001",
        f"got {result if result is not None else 'None'}",
    )


def test_certify_bracket_wide_around_midpoint():
    """The pf-6a4e oracle-bug shape: a WIDE bracket centered on the
    f32 midpoint should return None (the bracket spans across the
    midpoint and resolves to different f32s on either side)."""
    midpoint = Fraction(1, 1 << 150)
    rad = Fraction(1, 1 << 60)  # WIDER than 2^-451, straddles midpoint
    lower = midpoint - rad
    upper = midpoint + rad
    result = certified_round_f32(lower, upper, "NE")
    check(
        result is None,
        "NE wide bracket around midpoint -> None",
        f"got {result:#010x}" if result is not None else "got None (ok)",
    )


def test_certify_bracket_directed_mode():
    """RZ certifies a positive bracket as long as both endpoints
    floor to the same f32."""
    one = Fraction(1)
    ulp = Fraction(1, 1 << 23)
    # Bracket entirely in (1.0, 1.0 + ulp): both endpoints floor to 1.0
    # under RZ (positive).
    lower = one + ulp / 4
    upper = one + ulp / 2 - Fraction(1, 1 << 50)
    result = certified_round_f32(lower, upper, "RZ")
    check(
        result == 0x3F80_0000,
        "RZ certify positive bracket entirely above 1.0 below next grid",
        f"got {result if result is not None else 'None'}",
    )


def test_certify_invalid_bracket():
    """lower > upper raises ValueError."""
    try:
        certified_round_f32(Fraction(2), Fraction(1), "NE")
    except ValueError:
        return  # expected
    check(False, "lower > upper should raise ValueError", "no exception raised")


# ---------------------------------------------------------------------------
# Runner.
# ---------------------------------------------------------------------------


ALL_TESTS = [
    # Layer 1: boundary classes
    test_exact_zero,
    test_exact_f32_values_roundtrip,
    test_subnormal_midpoint_smallest,
    test_i1_truth_case,
    test_normal_subnormal_transition,
    test_f32_max_overflow,
    test_negative_values,
    # Layer 2: differential against struct.pack
    test_ne_agrees_with_struct_pack_on_exact_f32,
    test_ne_agrees_with_struct_pack_on_perturbed_exact,
    # Layer 3: random property tests
    test_random_rationals_ne,
    test_random_rationals_rna,
    # Layer 4: certified-rounding bracket tests
    test_certify_bracket_inside_single_bucket,
    test_certify_bracket_straddles_midpoint,
    test_certify_bracket_exact_zero_to_subnormal,
    test_certify_bracket_wide_around_midpoint,
    test_certify_bracket_directed_mode,
    test_certify_invalid_bracket,
]


def main() -> int:
    for t in ALL_TESTS:
        t()
    n_failures = len(FAILURES)
    print()
    if n_failures == 0:
        print(f"OK: {len(ALL_TESTS)} test groups passed.")
        return 0
    print(f"FAIL: {n_failures} check(s) failed across {len(ALL_TESTS)} test groups.")
    return 1


if __name__ == "__main__":
    sys.exit(main())
