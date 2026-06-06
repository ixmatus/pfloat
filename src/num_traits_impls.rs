//! `num-traits` impls for `FixedFloat<PREC>`.
//!
//! Gated behind the `num-traits` feature (which implies `fixed` and
//! `ops`). The impls are deliberately scoped to `FixedFloat<PREC>` and
//! not [`BigFloat`](crate::big::BigFloat): `num-traits`' constructors
//! carry no precision argument (`zero`, `one`, `from_str_radix`,
//! `from_i64`), so they are well defined only when the precision is
//! fixed by the type. `BigFloat`'s dynamic precision would force a
//! hidden default, which this crate declines to introduce; dynamic
//! precision code uses `BigFloat`'s explicit methods instead.
//!
//! Every operation runs at the type's precision `PREC`. The arithmetic
//! traits (`Num`'s `Add`/`Sub`/`Mul`/`Div`/`Rem`, `Inv`'s `Div`) reuse
//! the `core::ops` overloads from the `ops` feature, which round to
//! `PREC` under `NearestEven` and discard the `Status`. The numeric
//! conversions (`ToPrimitive`, `NumCast`) route through `f64`, so a
//! value with more than 53 significant bits loses precision converting
//! to an integer or float, the standard lossy contract of those
//! traits. `Float` and `PartialOrd` are intentionally absent: pfloat's
//! comparison returns `(Option<Ordering>, Status)` per IEEE 754-2019
//! §5.11, and `BigFloat`/`FixedFloat` carry no fixed associated
//! constants for `Float`. ADR-0070.

use core::cmp::Ordering;

use num_traits::{FromPrimitive, Inv, Num, NumCast, One, Signed, ToPrimitive, Zero};

use crate::big::BigFloat;
use crate::fixed::FixedFloat;
use crate::mantissa::limbs_for;
use crate::parse::ParseError;
use crate::rounding::RoundingMode;

const NE: RoundingMode = RoundingMode::NearestEven;

/// Error returned by `<FixedFloat<PREC> as Num>::from_str_radix`. Only
/// radix 10 is supported (pfloat's parser is decimal); other radixes
/// are rejected rather than silently mishandled.
#[derive(Clone, Debug, Eq, PartialEq)]
pub enum RadixParseError {
    /// A radix other than 10 was requested.
    UnsupportedRadix(u32),
    /// The decimal parse failed.
    Parse(ParseError),
}

impl core::fmt::Display for RadixParseError {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        match self {
            RadixParseError::UnsupportedRadix(r) => {
                write!(f, "pfloat: from_str_radix supports radix 10 only, got {r}")
            }
            RadixParseError::Parse(e) => write!(f, "pfloat: {e:?}"),
        }
    }
}

#[cfg(feature = "std")]
impl std::error::Error for RadixParseError {}

impl<const PREC: u32> Zero for FixedFloat<PREC>
where
    [(); limbs_for(PREC)]:,
{
    #[inline]
    fn zero() -> Self {
        FixedFloat::try_from_i64_exact(0).expect("0 is exact at any precision")
    }

    #[inline]
    fn is_zero(&self) -> bool {
        // Route through the predicate, not `== zero()`: the derived
        // structural equality would miss negative zero.
        self.to_big().is_zero()
    }
}

impl<const PREC: u32> One for FixedFloat<PREC>
where
    [(); limbs_for(PREC)]:,
{
    #[inline]
    fn one() -> Self {
        FixedFloat::try_from_i64_exact(1).expect("1 is exact at precision >= 1")
    }
    // `is_one`'s default (`*self == Self::one()`) is correct here: the
    // precision is fixed, so 1.0 has a single canonical representation
    // and the derived equality is value equality for it.
}

impl<const PREC: u32> Num for FixedFloat<PREC>
where
    [(); limbs_for(PREC)]:,
{
    type FromStrRadixErr = RadixParseError;

    fn from_str_radix(s: &str, radix: u32) -> Result<Self, Self::FromStrRadixErr> {
        if radix != 10 {
            return Err(RadixParseError::UnsupportedRadix(radix));
        }
        FixedFloat::parse_str(s, NE)
            .map(|(v, _status)| v)
            .map_err(RadixParseError::Parse)
    }
}

impl<const PREC: u32> Signed for FixedFloat<PREC>
where
    [(); limbs_for(PREC)]:,
{
    #[inline]
    fn abs(&self) -> Self {
        (*self).abs()
    }

    #[inline]
    fn abs_sub(&self, other: &Self) -> Self {
        // `max(self - other, 0)`, the (deprecated) positive difference.
        if matches!(self.partial_cmp(other).0, Some(Ordering::Greater)) {
            *self - *other
        } else {
            Self::zero()
        }
    }

    #[inline]
    fn signum(&self) -> Self {
        (*self).signum()
    }

    #[inline]
    fn is_positive(&self) -> bool {
        matches!(self.partial_cmp(&Self::zero()).0, Some(Ordering::Greater))
    }

    #[inline]
    fn is_negative(&self) -> bool {
        matches!(self.partial_cmp(&Self::zero()).0, Some(Ordering::Less))
    }
}

impl<const PREC: u32> FromPrimitive for FixedFloat<PREC>
where
    [(); limbs_for(PREC)]:,
{
    #[inline]
    fn from_i64(n: i64) -> Option<Self> {
        Some(FixedFloat::try_from_i64_round(n, NE).0)
    }

    fn from_u64(n: u64) -> Option<Self> {
        if let Ok(i) = i64::try_from(n) {
            Some(FixedFloat::try_from_i64_round(i, NE).0)
        } else {
            // n >= 2^63: split as (n - 2^63) + 2^63 so both pieces are
            // i64-constructible. Exact when PREC >= 64, rounded to PREC
            // otherwise — the same rounding `from_i64` would apply.
            let low = (n - (1u64 << 63)) as i64;
            let a = FixedFloat::try_from_i64_round(low, NE).0;
            let half = FixedFloat::try_from_i64_round(1i64 << 62, NE).0;
            let two = FixedFloat::try_from_i64_round(2, NE).0;
            Some(a + half * two)
        }
    }

    #[inline]
    fn from_f64(x: f64) -> Option<Self> {
        Some(FixedFloat::try_from_big_round(&BigFloat::from_f64(x), NE).0)
    }

    #[inline]
    fn from_f32(x: f32) -> Option<Self> {
        Some(FixedFloat::try_from_big_round(&BigFloat::from_f32(x), NE).0)
    }
}

impl<const PREC: u32> ToPrimitive for FixedFloat<PREC>
where
    [(); limbs_for(PREC)]:,
{
    fn to_i64(&self) -> Option<i64> {
        let f = self.to_big().to_f64();
        // 2^63 is the smallest f64 above i64::MAX, so the half-open
        // bound rejects exactly the values that would overflow.
        if f.is_finite() && (-9_223_372_036_854_775_808.0..9_223_372_036_854_775_808.0).contains(&f)
        {
            Some(f.trunc() as i64)
        } else {
            None
        }
    }

    fn to_u64(&self) -> Option<u64> {
        let f = self.to_big().to_f64();
        if f.is_finite() && (0.0..18_446_744_073_709_551_616.0).contains(&f) {
            Some(f.trunc() as u64)
        } else {
            None
        }
    }

    #[inline]
    fn to_f64(&self) -> Option<f64> {
        Some(self.to_big().to_f64())
    }

    #[inline]
    fn to_f32(&self) -> Option<f32> {
        Some(self.to_big().to_f32())
    }
}

impl<const PREC: u32> NumCast for FixedFloat<PREC>
where
    [(); limbs_for(PREC)]:,
{
    #[inline]
    fn from<T: ToPrimitive>(n: T) -> Option<Self> {
        n.to_f64().and_then(<Self as FromPrimitive>::from_f64)
    }
}

impl<const PREC: u32> Inv for FixedFloat<PREC>
where
    [(); limbs_for(PREC)]:,
{
    type Output = Self;

    #[inline]
    fn inv(self) -> Self {
        Self::one() / self
    }
}
