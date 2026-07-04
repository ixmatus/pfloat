//! [`Ball<T>`]: a rigorous midpoint-radius real enclosure.
//!
//! See [`crate::spec`] for the enclosure laws this type and its
//! operations uphold. ADR-0076.

use core::cmp::Ordering;

use pfloat::RoundingMode;

use crate::mag::Mag;
use crate::scalar::RealScalar;

/// A rigorous real enclosure `[mid − rad, mid + rad]`: a full-precision
/// pfloat midpoint `mid` and an upward-rounded radius `rad`.
///
/// The midpoint is always finite (the constructors reject NaN and `±∞`);
/// the radius is a [`Mag`], non-negative and never inward-rounded by
/// construction. `rad = 0` is an *exact* ball denoting the single point
/// `{mid}`; `rad = +∞` denotes the whole real line.
///
/// Equality is structural (same `mid`, same `rad`), matching pfloat's
/// scalar equality: `+0` and `−0` midpoints are distinct, as are two
/// balls denoting the same interval through different `(mid, rad)` pairs.
///
/// `Deserialize` is hand-written (see below) rather than derived: the
/// derive would admit a non-finite midpoint (NaN or `±∞`), violating the
/// ball's core invariant that a midpoint is a finite real. The custom
/// impl reconstructs through [`Ball::new`], which enforces it, and the
/// `rad` field's own validating [`Mag`] deserialize guards the radius.
#[derive(Clone, Debug, PartialEq, Eq, Hash)]
#[cfg_attr(feature = "serde", derive(serde::Serialize))]
pub struct Ball<T> {
    mid: T,
    rad: Mag,
}

/// Why constructing a [`Ball`] failed.
#[derive(Clone, Copy, Debug, PartialEq, Eq, Hash)]
#[cfg_attr(feature = "serde", derive(serde::Serialize, serde::Deserialize))]
pub enum BallError {
    /// The midpoint was NaN or `±∞`. A ball denotes a real interval, so
    /// its centre must be a finite real.
    NonFiniteMidpoint,
    /// An endpoint passed to [`Ball::from_interval`] was NaN or `±∞`.
    /// (Unbounded and half-bounded intervals are the IEEE 1788 interval
    /// face, a separate later crate.)
    NonFiniteEndpoint,
    /// `from_interval(lo, hi)` was called with `lo > hi` (or an
    /// unordered NaN comparison), which denotes no real interval.
    ReversedInterval,
}

impl<T: RealScalar> Ball<T> {
    /// Builds a ball from a finite midpoint and a radius.
    ///
    /// Returns [`BallError::NonFiniteMidpoint`] if `mid` is NaN or `±∞`.
    /// Radius non-negativity needs no check — it is a [`Mag`] type fact.
    pub fn new(mid: T, rad: Mag) -> Result<Self, BallError> {
        if mid.is_finite() {
            Ok(Self { mid, rad })
        } else {
            Err(BallError::NonFiniteMidpoint)
        }
    }

    /// Builds an *exact* ball (`rad = 0`) denoting the single point
    /// `{value}`.
    ///
    /// Returns [`BallError::NonFiniteMidpoint`] if `value` is not finite.
    pub fn point(value: T) -> Result<Self, BallError> {
        Self::new(value, Mag::ZERO)
    }

    /// Internal: assemble from parts when the midpoint is already known
    /// finite (the arithmetic kernels guarantee it).
    #[inline]
    pub(crate) fn from_parts(mid: T, rad: Mag) -> Self {
        debug_assert!(mid.is_finite(), "from_parts midpoint must be finite");
        Self { mid, rad }
    }

    /// Internal: the entire real line `[0 ± ∞]`, the sound result of an
    /// operation that leaves the reals (division by a zero-containing
    /// ball, sqrt of a wholly-negative ball).
    #[inline]
    pub(crate) fn entire(precision: u32) -> Self {
        Self {
            mid: T::zero(precision),
            rad: Mag::INFINITY,
        }
    }

    /// The midpoint.
    #[inline]
    pub fn midpoint(&self) -> &T {
        &self.mid
    }

    /// The radius.
    #[inline]
    pub fn radius(&self) -> Mag {
        self.rad
    }

    /// Precision of the midpoint, in bits.
    #[inline]
    pub fn precision(&self) -> u32 {
        self.mid.precision()
    }

    /// `true` when the ball is exact (`rad = 0`): it denotes the single
    /// point `{mid}`.
    #[inline]
    pub fn is_exact(&self) -> bool {
        self.rad.is_zero()
    }

    /// `true` when the radius is unbounded (`rad = +∞`): the ball denotes
    /// the whole real line.
    #[inline]
    pub fn is_entire(&self) -> bool {
        self.rad.is_infinite()
    }

    /// The lower endpoint `mid − rad`, rounded toward `−∞` — the tightest
    /// representable lower bound of the enclosure (Law 4: ball-to-endpoints
    /// is exact). An entire ball yields `−∞`; an exact ball yields `mid`.
    #[must_use]
    pub fn lower(&self) -> T {
        let rad = T::radius_to_scalar(self.rad);
        self.mid.sub(&rad, RoundingMode::TowardNegative).0
    }

    /// The upper endpoint `mid + rad`, rounded toward `+∞` — the tightest
    /// representable upper bound (Law 4). An entire ball yields `+∞`; an
    /// exact ball yields `mid`.
    #[must_use]
    pub fn upper(&self) -> T {
        let rad = T::radius_to_scalar(self.rad);
        self.mid.add(&rad, RoundingMode::TowardPositive).0
    }

    /// Builds the smallest sound ball enclosing the interval `[lo, hi]`
    /// (Law 4: endpoints-to-ball is sound but inflating).
    ///
    /// The midpoint `(lo + hi)/2` is computed at working precision and
    /// may round, so the radius is set by the **never-assume-centered**
    /// formula `rad ≥ round_up(max(mid − lo, hi − mid))`, which contains
    /// both endpoints unconditionally.
    ///
    /// Returns [`BallError::NonFiniteEndpoint`] if either endpoint is not
    /// finite, or [`BallError::ReversedInterval`] if `lo > hi`.
    pub fn from_interval(lo: &T, hi: &T) -> Result<Self, BallError> {
        if !lo.is_finite() || !hi.is_finite() {
            return Err(BallError::NonFiniteEndpoint);
        }
        match lo.compare(hi).0 {
            Some(Ordering::Less) | Some(Ordering::Equal) => {}
            _ => return Err(BallError::ReversedInterval),
        }

        // mid = (lo + hi) / 2 at the natural (max) precision. The add may
        // round; halving by scale_by_pow2(-1) is exact, so all the
        // rounding lives in the add and is covered by the radius below.
        let (sum, _) = lo.add(hi, RoundingMode::NearestEven);
        let (mid, _) = sum.scale_by_pow2(-1);

        // mid ∈ [lo, hi] up to the add's rounding, so mid − lo and
        // hi − mid can each be slightly negative under rounding; bounding
        // each ABOVE toward +∞ and taking the magnitude up to Mag keeps
        // the result a sound upper bound regardless.
        let (left, _) = mid.sub(lo, RoundingMode::TowardPositive);
        let (right, _) = hi.sub(&mid, RoundingMode::TowardPositive);
        let rad = left.magnitude_to_mag().max(right.magnitude_to_mag());

        // mid is finite (sum of two finite values, halved); new() still
        // guards defensively.
        Self::new(mid, rad)
    }
}

/// Validating [`Deserialize`](serde::Deserialize) for [`Ball<T>`].
///
/// Deserialize is a trust boundary. The derived impl would admit a NaN or
/// `±∞` midpoint, breaking the finite-midpoint invariant; this impl routes
/// the deserialized parts through [`Ball::new`], which rejects a
/// non-finite midpoint, while the `rad` field is a [`Mag`] whose own
/// validating deserialize (see `mag.rs`) revalidates the radius canonical
/// form. Malformed input is rejected with a serde error, never silently
/// coerced. Mirrors pfloat's `BigFloat` deserialize discipline
/// (ADR-0068); ADR-0116.
#[cfg(feature = "serde")]
impl<'de, T> serde::Deserialize<'de> for Ball<T>
where
    T: RealScalar + serde::Deserialize<'de>,
{
    fn deserialize<D>(deserializer: D) -> Result<Self, D::Error>
    where
        D: serde::Deserializer<'de>,
    {
        // A shadow struct with the identical wire form (same field names),
        // deserialized without invariants and then validated through the
        // normal constructor.
        #[derive(serde::Deserialize)]
        struct RawBall<T> {
            mid: T,
            rad: Mag,
        }
        let raw = RawBall::<T>::deserialize(deserializer)?;
        // The only construction failure is a non-finite midpoint; the
        // radius is already a validated `Mag`.
        Ball::new(raw.mid, raw.rad).map_err(|_| {
            serde::de::Error::custom(
                "pfloat-ball: Ball midpoint must be finite (NaN or ±∞ rejected)",
            )
        })
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use pfloat::BigFloat;

    fn bf(n: i64, p: u32) -> BigFloat {
        BigFloat::try_from_i64_exact(n, p).unwrap()
    }

    fn eq(a: &BigFloat, b: &BigFloat) -> bool {
        a.partial_cmp(b).0 == Some(Ordering::Equal)
    }

    #[test]
    fn point_is_exact_and_endpoints_collapse() {
        let b = Ball::point(bf(5, 53)).unwrap();
        assert!(b.is_exact());
        assert!(eq(&b.lower(), &bf(5, 53)));
        assert!(eq(&b.upper(), &bf(5, 53)));
    }

    #[test]
    fn rejects_non_finite_midpoint() {
        let inf = BigFloat::try_new_infinity(pfloat::Sign::Positive, 53).unwrap();
        assert_eq!(Ball::point(inf).unwrap_err(), BallError::NonFiniteMidpoint);
        let nan = BigFloat::try_new_quiet_nan(pfloat::Sign::Positive, 53, &[]).unwrap();
        assert_eq!(
            Ball::new(nan, Mag::ZERO).unwrap_err(),
            BallError::NonFiniteMidpoint
        );
    }

    #[test]
    fn endpoints_bracket_the_midpoint() {
        // A ball [4 ± 1]: lower 3, upper 5.
        let one = Mag::from_pow2(0); // 2^0 = 1
        let b = Ball::new(bf(4, 53), one).unwrap();
        assert!(eq(&b.lower(), &bf(3, 53)));
        assert!(eq(&b.upper(), &bf(5, 53)));
    }

    #[test]
    fn entire_ball_has_infinite_endpoints() {
        let b = Ball::new(bf(0, 53), Mag::INFINITY).unwrap();
        assert!(b.is_entire());
        assert!(b.lower().is_infinite() && b.lower().is_sign_negative());
        assert!(b.upper().is_infinite() && !b.upper().is_sign_negative());
    }

    #[test]
    fn from_interval_contains_both_endpoints() {
        // [lo, hi] = [3, 7]: mid 5, rad >= 2.
        let lo = bf(3, 53);
        let hi = bf(7, 53);
        let b = Ball::from_interval(&lo, &hi).unwrap();
        // lower <= lo and upper >= hi (containment).
        assert!(b.lower().partial_cmp(&lo).0 != Some(Ordering::Greater));
        assert!(b.upper().partial_cmp(&hi).0 != Some(Ordering::Less));
    }

    #[test]
    fn from_interval_off_center_is_sound() {
        // An interval whose midpoint is not exactly representable forces
        // the never-assume-centered radius to inflate. [0, 1] at p=2:
        // mid = 0.5 exact here, but the containment must hold regardless.
        let lo = bf(0, 8);
        let (hi, _) = bf(1, 8).scale_by_pow2(0);
        let b = Ball::from_interval(&lo, &hi).unwrap();
        assert!(b.lower().partial_cmp(&lo).0 != Some(Ordering::Greater));
        assert!(b.upper().partial_cmp(&hi).0 != Some(Ordering::Less));
    }

    #[test]
    fn from_interval_rejects_reversed_and_non_finite() {
        let lo = bf(7, 53);
        let hi = bf(3, 53);
        assert_eq!(
            Ball::from_interval(&lo, &hi).unwrap_err(),
            BallError::ReversedInterval
        );
        let inf = BigFloat::try_new_infinity(pfloat::Sign::Positive, 53).unwrap();
        assert_eq!(
            Ball::from_interval(&bf(0, 53), &inf).unwrap_err(),
            BallError::NonFiniteEndpoint
        );
    }

    #[test]
    fn from_interval_degenerate_point() {
        // lo == hi: an exact (or near-exact) ball at that point.
        let p = bf(42, 53);
        let b = Ball::from_interval(&p, &p).unwrap();
        assert!(b.lower().partial_cmp(&p).0 != Some(Ordering::Greater));
        assert!(b.upper().partial_cmp(&p).0 != Some(Ordering::Less));
        // For lo == hi exactly representable, the radius is zero.
        assert!(b.is_exact());
    }
}
