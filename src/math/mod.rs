//! Elementary transcendentals (Phase 3).
//!
//! Each function uses range reduction to bring the argument into a
//! Taylor-friendly window, evaluates the reduced argument's series,
//! then composes the result. Slice 3a ships [`exp`](crate::big::BigFloat::exp);
//! subsequent slices fill in `ln`, trig, hyperbolic, and `pow`.
//!
//! All functions are gated behind the `exp-log` or `trig` cluster
//! features so embedded users can compile in only what they need.

#[cfg(feature = "exp-log")]
pub(crate) mod exp;
