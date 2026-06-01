//! Arithmetic kernels (slice 1c onward).
//!
//! Each kernel routes through the rounding pipeline in
//! [`crate::rounding::round_finite_to_precision`] for IEEE 754-2019
//! §6.5 rounding behavior. Special-case dispatch (NaN, infinity,
//! signed zero) lives at the top of each kernel and matches the
//! IEEE 754 §6.2 / §6.3 rules.
//!
//! Slice 1c ships [`addsub`] (add, sub). Slices 1d–1f add mul, div,
//! sqrt, fma in their own files.

#[cfg(feature = "big")]
pub(crate) mod addsub;
#[cfg(feature = "big")]
pub(crate) mod cbrt;
#[cfg(feature = "big")]
pub(crate) mod div;
#[cfg(feature = "big")]
pub(crate) mod fma;
#[cfg(feature = "big")]
pub(crate) mod limbs;
#[cfg(feature = "big")]
pub(crate) mod mul;
// `rootn` uses the Ziv driver (gated `exp-log`), so it cannot live under
// bare `big`; `cbrt` (exact-integer, no Ziv) does. ADR-0056.
#[cfg(feature = "exp-log")]
pub(crate) mod rootn;
#[cfg(feature = "big")]
pub(crate) mod sqrt;

#[cfg(feature = "big")]
use crate::big::BigFloat;
#[cfg(feature = "big")]
use crate::class::Class;
#[cfg(feature = "big")]
use crate::sign::Sign;
#[cfg(feature = "big")]
use crate::status::Status;

/// Propagate NaN per IEEE 754-2019 §6.2.3 for a two-operand op.
///
/// Returns:
/// - `Some((value, Status::INVALID))` if either operand is a
///   signaling NaN (the result is a quiet NaN; INVALID is raised).
/// - `Some((value, Status::OK))` if exactly one operand is a quiet
///   NaN (that NaN propagates).
/// - `Some((value, Status::OK))` if both are quiet NaNs (`a` wins
///   by convention).
/// - `None` if neither operand is a NaN; the caller continues with
///   the operation.
///
/// The propagated NaN is re-emitted at `target_precision` so the
/// caller's chosen result precision applies even for NaN passthrough.
/// IEEE 754-2019 §6.2.3 leaves the payload-propagation policy
/// implementation-defined within the constraints of §6.2.1's
/// quiet/signaling rules; pfloat's policy matches ferrodec's:
/// preserve the first NaN's payload (zero-padded or truncated to
/// fit the new precision), or use empty payload when an sNaN was
/// converted to qNaN.
#[cfg(feature = "big")]
#[allow(dead_code)] // consumed by slice 1c+ kernels
pub(crate) fn propagate_nan2(
    a: &BigFloat,
    b: &BigFloat,
    target_precision: u32,
) -> Option<(BigFloat, Status)> {
    if a.is_signaling_nan() || b.is_signaling_nan() {
        let nan = BigFloat::try_new_quiet_nan(Sign::Positive, target_precision, &[])
            .expect("BigFloat invariant: precision >= 1");
        crate::status::auto_raise(Status::INVALID);
        return Some((nan, Status::INVALID));
    }
    if let Class::Nan {
        quiet,
        sign,
        payload,
    } = &a.class
    {
        let nan = BigFloat::try_new_quiet_nan(*sign, target_precision, payload)
            .expect("BigFloat invariant: precision >= 1");
        let _ = quiet; // already gated above
        return Some((nan, Status::OK));
    }
    if let Class::Nan {
        quiet: _,
        sign,
        payload,
    } = &b.class
    {
        let nan = BigFloat::try_new_quiet_nan(*sign, target_precision, payload)
            .expect("BigFloat invariant: precision >= 1");
        return Some((nan, Status::OK));
    }
    None
}

/// Propagate NaN per IEEE 754-2019 §6.2.3 for a three-operand op
/// (FMA). Same shape as [`propagate_nan2`] extended to three
/// operands; sNaN priority is left-to-right.
#[cfg(feature = "big")]
#[allow(dead_code)] // consumed by slice 1f FMA kernel
pub(crate) fn propagate_nan3(
    a: &BigFloat,
    b: &BigFloat,
    c: &BigFloat,
    target_precision: u32,
) -> Option<(BigFloat, Status)> {
    if a.is_signaling_nan() || b.is_signaling_nan() || c.is_signaling_nan() {
        let nan = BigFloat::try_new_quiet_nan(Sign::Positive, target_precision, &[])
            .expect("BigFloat invariant: precision >= 1");
        crate::status::auto_raise(Status::INVALID);
        return Some((nan, Status::INVALID));
    }
    for op in [a, b, c] {
        if let Class::Nan { sign, payload, .. } = &op.class {
            let nan = BigFloat::try_new_quiet_nan(*sign, target_precision, payload)
                .expect("BigFloat invariant: precision >= 1");
            return Some((nan, Status::OK));
        }
    }
    None
}
