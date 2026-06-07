//! The in-tree ball enclosure specification.
//!
//! This module is the written-down contract every `pfloat-ball`
//! operation must uphold. It is prose and a few `const` law statements,
//! not runtime code: the laws are the obligations the implementation and
//! its verification lanes discharge, kept in tree (per the "promote
//! durable rationale" discipline) rather than deferred to a reviewer's
//! memory or to Arb's runtime behaviour. ADR-0076.
//!
//! # What a ball denotes
//!
//! A ball `Ball { mid, rad }` with finite midpoint `mid` and radius
//! `rad` (a [`Mag`](crate::Mag)) denotes the closed real interval
//!
//! ```text
//!     [mid − rad, mid + rad]
//! ```
//!
//! when `rad` is finite, and the whole real line when `rad = +∞`. A ball
//! with `rad = 0` denotes the single point `{mid}` (an *exact* ball).
//!
//! # Law 1: enclosure soundness (the Fundamental Theorem of Interval
//! Arithmetic, FTIA)
//!
//! For a ball operation implementing the real function `f` (unary or
//! binary), the output ball must contain `f` applied to every point of
//! the input ball(s):
//!
//! ```text
//!     for all x in [a],  f(x) ∈ op([a])          (unary)
//!     for all x in [a], y in [b],  f(x, y) ∈ op([a], [b])   (binary)
//! ```
//!
//! This is the one property that makes a ball *rigorous*: the result is
//! never a claim narrower than the truth. Soundness is one-directional —
//! the radius may over-estimate (a wider-than-necessary enclosure is
//! still correct), but it must **never** under-estimate. An
//! under-estimating radius is the single defect that turns every
//! downstream enclosure into a falsehood, which is why
//! [`Mag`](crate::Mag) makes an inward-rounded radius unrepresentable and
//! the verification lanes weight radius soundness above all else.
//!
//! # Law 2: the directed-pair radius-soundness law (v1.0 route)
//!
//! The v1.0 radius for an arithmetic operation is built from the
//! *directed pair*: compute the result toward `−∞` (`lo`) and toward
//! `+∞` (`hi`) with the correctly-rounded scalar kernel, take the
//! midpoint as the round-to-nearest result, and set
//!
//! ```text
//!     rad ≥ round_up( max(mid − lo, hi − mid) )  +  propagated input error
//! ```
//!
//! under [`Mag`](crate::Mag)'s unconditional round-up. Soundness follows
//! directly from the directed kernels' correctness: the true result lies
//! in `[lo, hi]`, so `rad` bounding half that spread (rounded up) bounds
//! the kernel's own residual. **There is no separate Ziv-residual term**,
//! because this route never invokes the Ziv driver — a point that must
//! not be conflated with the surfaced-half-width route (deferred), whose
//! soundness *does* require adding the Ziv residual even on convergence.
//!
//! # Law 3: exact-in-exact-out
//!
//! When a scalar kernel returns its result exactly (the directed pair
//! coincides: `lo == hi`), the ball op emits a **zero-radius** ball.
//! Exactness in must produce exactness out; a spurious positive radius on
//! an exact result is a (sound but) avoidable loss of information, and
//! this law forbids it. Reuses pfloat's exact-dispatch discipline.
//!
//! # Law 4: the conversion boundary is asymmetric
//!
//! Ball-to-endpoints is **exact**: `lower = mid ⊖ rad` toward `−∞` and
//! `upper = mid ⊕ rad` toward `+∞` are computed by the directed kernels
//! and are the tightest representable endpoints containing the ball, with
//! no slack beyond one outward rounding.
//!
//! Endpoints-to-ball is **sound but inflating**: `from_interval(lo, hi)`
//! must never assume the midpoint is centred. It sets
//!
//! ```text
//!     rad ≥ round_up( max(mid − lo, hi − mid) )
//! ```
//!
//! which contains both endpoints unconditionally even though `mid`
//! rounded. Reserve the word "lossless" for the ball-to-endpoints
//! direction only; the reverse direction is "sound" or "containing."
//!
//! # Law 5: status is the secondary channel
//!
//! `Status` composes through the OR-monoid: `INEXACT` if any component
//! rounded, `INVALID` for a NaN component, `DIV_BY_ZERO` for a zero
//! divisor. On a ball, `INEXACT` is the *normal* correct outcome (a small
//! positive radius), so the radius — surfaced as `rel_accuracy_bits` — is
//! the primary accuracy channel and `Status` is the secondary IEEE-flag
//! channel. `OVERFLOW`/`UNDERFLOW` are essentially unreachable at the
//! `i64` exponent scale, so radius blow-up is read off the radius
//! directly, as Arb does.
//!
//! # The one-directional rule, in one line
//!
//! Every radius rounding in this crate is *outward*: radii round up, and
//! endpoints round away from the midpoint. A future change that rounds
//! any radius inward, on any path, is a soundness regression regardless
//! of how small the effect looks.
