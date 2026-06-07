//! Hand-written serde impls for [`BigFloat`] and `FixedFloat<PREC>`.
//!
//! Gated behind the `serde` feature. The simple public types (`Sign`,
//! `Status`, `RoundingMode`, `IeeeClass`, `BuildError`, `ParseError`)
//! get a plain `#[derive]` at their definition; the two precision
//! profiles are hand-written here because their canonical form carries
//! invariants that a derive cannot enforce on deserialize.
//!
//! Wire form: a `{ precision, class }` struct where `class` is an
//! externally-tagged enum mirroring the internal [`Class`] variants
//! (`Zero`, `Infinity`, `Nan`, `Normal`). The form is bit-exact by
//! construction (raw sign, exponent, and mantissa or payload limbs),
//! so a round-trip recovers the value exactly under any serde format,
//! human-readable or compact. Serialization borrows the limbs (no
//! clone); deserialization owns them.
//!
//! Deserialization is a trust boundary (it may run on attacker-supplied
//! bytes), so [`big_from_repr`] revalidates the canonical form rather
//! than trusting the input: precision at least one bit, the mantissa
//! limb count matching `limbs_for(precision)`, the top bit of the
//! most-significant limb set, and the storage-padding bits below the
//! precision clear. A `FixedFloat<PREC>` additionally routes through
//! [`FixedFloat::try_from_big_exact`], which rejects any value whose
//! precision is not `PREC`. Malformed input is rejected with a serde
//! error, never silently coerced. ADR-0068.

#[cfg(feature = "big")]
extern crate alloc;

#[cfg(feature = "big")]
use alloc::vec::Vec;

#[cfg(feature = "big")]
use serde::{Deserialize, Deserializer, Serialize, Serializer};

#[cfg(feature = "big")]
use crate::big::BigFloat;
#[cfg(feature = "big")]
use crate::class::Class;
#[cfg(feature = "big")]
use crate::mantissa::limbs_for;
#[cfg(feature = "big")]
use crate::sign::Sign;

// ----- wire representations -----
//
// A borrowed pair for serialization (no limb clone) and an owned pair
// for deserialization. Both derive the same field and variant names, so
// they share one wire form.

#[cfg(feature = "big")]
#[derive(Serialize)]
struct BigFloatRef<'a> {
    precision: u32,
    class: ClassRef<'a>,
}

#[cfg(feature = "big")]
#[derive(Serialize)]
enum ClassRef<'a> {
    Zero {
        sign: Sign,
    },
    Infinity {
        sign: Sign,
    },
    Nan {
        quiet: bool,
        sign: Sign,
        payload: &'a [u64],
    },
    Normal {
        sign: Sign,
        exponent: i64,
        mantissa: &'a [u64],
    },
}

#[cfg(feature = "big")]
#[derive(Deserialize)]
struct BigFloatOwned {
    precision: u32,
    class: ClassOwned,
}

#[cfg(feature = "big")]
#[derive(Deserialize)]
enum ClassOwned {
    Zero {
        sign: Sign,
    },
    Infinity {
        sign: Sign,
    },
    Nan {
        quiet: bool,
        sign: Sign,
        payload: Vec<u64>,
    },
    Normal {
        sign: Sign,
        exponent: i64,
        mantissa: Vec<u64>,
    },
}

// ----- BigFloat -----

#[cfg(feature = "big")]
impl Serialize for BigFloat {
    fn serialize<S: Serializer>(&self, serializer: S) -> Result<S::Ok, S::Error> {
        let class = match &self.class {
            Class::Zero { sign } => ClassRef::Zero { sign: *sign },
            Class::Infinity { sign } => ClassRef::Infinity { sign: *sign },
            Class::Nan {
                quiet,
                sign,
                payload,
            } => ClassRef::Nan {
                quiet: *quiet,
                sign: *sign,
                payload,
            },
            Class::Normal {
                sign,
                exponent,
                mantissa,
            } => ClassRef::Normal {
                sign: *sign,
                exponent: *exponent,
                mantissa,
            },
        };
        BigFloatRef {
            precision: self.precision,
            class,
        }
        .serialize(serializer)
    }
}

#[cfg(feature = "big")]
impl<'de> Deserialize<'de> for BigFloat {
    fn deserialize<D: Deserializer<'de>>(deserializer: D) -> Result<Self, D::Error> {
        let repr = BigFloatOwned::deserialize(deserializer)?;
        big_from_repr(repr).map_err(serde::de::Error::custom)
    }
}

/// Reconstruct a [`BigFloat`] from a deserialized representation,
/// revalidating every canonical-form invariant. Untrusted input is
/// rejected rather than coerced.
#[cfg(feature = "big")]
fn big_from_repr(repr: BigFloatOwned) -> Result<BigFloat, &'static str> {
    let precision = repr.precision;
    if precision == 0 {
        return Err("pfloat: precision must be at least 1 bit");
    }
    let class = match repr.class {
        ClassOwned::Zero { sign } => Class::Zero { sign },
        ClassOwned::Infinity { sign } => Class::Infinity { sign },
        ClassOwned::Nan {
            quiet,
            sign,
            payload,
        } => {
            if payload.len() != limbs_for(precision) {
                return Err("pfloat: NaN payload limb count does not match precision");
            }
            Class::Nan {
                quiet,
                sign,
                payload,
            }
        }
        ClassOwned::Normal {
            sign,
            exponent,
            mantissa,
        } => {
            validate_normal(&mantissa, precision)?;
            Class::Normal {
                sign,
                exponent,
                mantissa,
            }
        }
    };
    Ok(BigFloat { class, precision })
}

/// Validate a `Normal` mantissa against the canonical-form rules
/// (ADR-0001, ADR-0002): the limb count matches `limbs_for(precision)`,
/// the most-significant bit of the most-significant limb is set, and
/// the storage-padding bits below the precision are clear.
#[cfg(feature = "big")]
fn validate_normal(mantissa: &[u64], precision: u32) -> Result<(), &'static str> {
    let limbs = limbs_for(precision);
    if mantissa.len() != limbs {
        return Err("pfloat: Normal mantissa limb count does not match precision");
    }
    // `limbs >= 1` because `precision >= 1` (checked by the caller).
    if mantissa[limbs - 1] >> 63 != 1 {
        return Err("pfloat: Normal mantissa is not normalized (top bit clear)");
    }
    // Bits below the precision live in the least-significant limb; the
    // padding width is `limbs * 64 - precision`, always in `0..64`. The
    // product is computed in u64 so a near-`u32::MAX` precision (where
    // `limbs * 64` reaches 2^32) cannot overflow before the subtraction
    // brings it back into range; this site takes attacker bytes.
    let pad = ((limbs as u64) * 64 - u64::from(precision)) as u32;
    if pad > 0 && mantissa[0] & ((1u64 << pad) - 1) != 0 {
        return Err("pfloat: Normal mantissa has nonzero bits below the precision");
    }
    Ok(())
}

// ----- FixedFloat<PREC> -----

#[cfg(feature = "fixed")]
use crate::fixed::FixedFloat;

#[cfg(feature = "fixed")]
impl<const PREC: u32> Serialize for FixedFloat<PREC>
where
    [(); limbs_for(PREC)]:,
{
    fn serialize<S: Serializer>(&self, serializer: S) -> Result<S::Ok, S::Error> {
        // Serialize through the BigFloat view: a FixedFloat<PREC> is a
        // BigFloat pinned at precision PREC, so the wire form coincides.
        self.to_big().serialize(serializer)
    }
}

#[cfg(feature = "fixed")]
impl<'de, const PREC: u32> Deserialize<'de> for FixedFloat<PREC>
where
    [(); limbs_for(PREC)]:,
{
    fn deserialize<D: Deserializer<'de>>(deserializer: D) -> Result<Self, D::Error> {
        let big = BigFloat::deserialize(deserializer)?;
        // try_from_big_exact rejects any precision other than PREC, so a
        // mismatched serialized precision is an error rather than a
        // silent re-round.
        FixedFloat::try_from_big_exact(big).map_err(|_| {
            serde::de::Error::custom(
                "pfloat: deserialized precision does not match FixedFloat<PREC>",
            )
        })
    }
}
