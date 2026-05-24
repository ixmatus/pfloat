//! `MetaOracle`: per-`FnId` dispatcher across the MPFR and Arb
//! backends. The Phase 1 surface splits into 33 MPFR-primary
//! `FnId`s and 12 Arb-primary ones; the runner uses one
//! [`OracleBackend`] handle so the verification core does not have
//! to switch backends per call. ADR-0034 calls this the "single
//! `FnId` enum dispatch" shape.
//!
//! Compilation: this module is gated under `differential-mpfr` (the
//! base feature; MPFR backend is always present in the oracle
//! harness). The Arb backend is opt-in via `differential-arb`; when
//! absent, the Arb-primary `FnId`s route through a NaN-enclosure
//! fallback that propagates as `Verdict::OracleInconclusive`
//! downstream (the verifier's `certified_round_f32` returns
//! `Some(f32::NAN)` only when *both* endpoints are NaN, which the
//! pfloat-side `kernel(f, input, mode)` will not match unless
//! pfloat also returns NaN — i.e. the inconclusive verdict is
//! honest about the absence of an oracle for the `FnId`).

#![cfg(all(unix, feature = "differential-mpfr"))]

#[cfg(feature = "differential-arb")]
use super::arb::{ArbError, ArbOracle};
use super::mpfr::MpfrOracle;
use super::types::{Enclosure, FnId, OracleBackend};

#[cfg(not(feature = "differential-arb"))]
use rug::float::Special;
#[cfg(not(feature = "differential-arb"))]
use rug::Float;

/// Errors `MetaOracle::new` can surface. Today the only failure
/// mode is the Arb backend's setup error; when the
/// `differential-arb` feature is not enabled the constructor is
/// infallible and this enum has no constructible variants.
#[derive(Debug)]
pub enum MetaError {
    #[cfg(feature = "differential-arb")]
    Arb(ArbError),
}

impl std::fmt::Display for MetaError {
    #[cfg(feature = "differential-arb")]
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        match self {
            Self::Arb(e) => write!(f, "MetaOracle Arb backend setup: {e}"),
        }
    }

    // Without `differential-arb` the enum has no constructible
    // variants, so this branch is unreachable at runtime; the type
    // system enforces the absence of any value to format and the
    // `match *self {}` body discharges the function with no code.
    #[cfg(not(feature = "differential-arb"))]
    fn fmt(&self, _f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        match *self {}
    }
}

impl std::error::Error for MetaError {}

#[cfg(feature = "differential-arb")]
impl From<ArbError> for MetaError {
    fn from(e: ArbError) -> Self {
        Self::Arb(e)
    }
}

/// The per-`FnId` dispatching oracle. Holds one [`MpfrOracle`] (a
/// unit struct, free) and, under `differential-arb`, one
/// [`ArbOracle`] (owns a Python subprocess).
pub struct MetaOracle {
    mpfr: MpfrOracle,
    #[cfg(feature = "differential-arb")]
    arb: ArbOracle,
}

impl MetaOracle {
    /// Construct the `MetaOracle`. Under `differential-arb` this
    /// spawns the python-flint worker subprocess; without that
    /// feature the constructor is essentially free. The `Result`
    /// shape is preserved across feature flavors so callers do not
    /// have to branch; without `differential-arb` the constructor
    /// is infallible but the wrap stays for signature stability.
    #[allow(clippy::unnecessary_wraps)]
    pub fn new() -> Result<Self, MetaError> {
        let mpfr = MpfrOracle;
        #[cfg(feature = "differential-arb")]
        {
            let arb = ArbOracle::new()?;
            Ok(Self { mpfr, arb })
        }
        #[cfg(not(feature = "differential-arb"))]
        {
            Ok(Self { mpfr })
        }
    }

    /// Dispatch to the Arb backend when available, otherwise return
    /// a NaN enclosure (the verifier will then report
    /// `Verdict::OracleInconclusive` since the pfloat kernel's
    /// output will not match the certified `Some(f32::NAN)` unless
    /// the pfloat kernel also returns NaN). The `&self` argument
    /// stays for signature parity with the `differential-arb` arm
    /// where it dispatches to the owned `ArbOracle`.
    #[cfg(feature = "differential-arb")]
    fn enclose_arb(
        &self,
        f: FnId,
        input: u32,
        mode: pfloat::RoundingMode,
        working_prec: u32,
    ) -> Enclosure {
        self.arb.enclose(f, input, mode, working_prec)
    }

    #[cfg(not(feature = "differential-arb"))]
    #[allow(clippy::unused_self)]
    fn enclose_arb(
        &self,
        _f: FnId,
        _input: u32,
        _mode: pfloat::RoundingMode,
        working_prec: u32,
    ) -> Enclosure {
        let nan = Float::with_val(working_prec, Special::Nan);
        Enclosure {
            lo: nan.clone(),
            hi: nan,
        }
    }
}

impl OracleBackend for MetaOracle {
    fn enclose(
        &self,
        f: FnId,
        input: u32,
        mode: pfloat::RoundingMode,
        working_prec: u32,
    ) -> Enclosure {
        if is_arb_primary(f) {
            self.enclose_arb(f, input, mode, working_prec)
        } else {
            self.mpfr.enclose(f, input, mode, working_prec)
        }
    }

    /// Authoritative iff the dispatched backend would be
    /// authoritative for the `FnId` being queried. Without
    /// information about WHICH `FnId` is being asked, default to
    /// `false` so the verifier runs its Ziv loop, which is correct
    /// for MPFR-primary inputs; the Arb-primary path then short
    /// circuits at the first iteration because the Arb worker
    /// returns a single-point enclosure that certifies immediately.
    ///
    /// Note: this is a per-instance default; once
    /// [`OracleBackend`] gains a per-call authority signal (a
    /// follow-up slice; ADR-0035 calls this out as an
    /// optimization) the dispatcher can per-FnId-route the
    /// authority flag.
    fn is_authoritative(&self) -> bool {
        false
    }

    fn name(&self) -> &'static str {
        // The MetaOracle is a dispatcher; per-row status emission
        // should record the actual backend that produced the
        // enclosure, not "Meta". Callers that need the per-row
        // oracle name (the runner does, for the TOML row) use
        // `oracle_name_for(f)` below.
        "Meta"
    }
}

/// Returns `true` for the twelve Arb-primary `FnId`s (the ones MPFR
/// cannot cover). The list mirrors the dispatch in
/// [`crate::oracle::arb::fnid_to_worker_args`].
pub fn is_arb_primary(f: FnId) -> bool {
    matches!(
        f,
        FnId::Si
            | FnId::Ci
            | FnId::Li
            | FnId::Bi
            | FnId::AiPrime
            | FnId::BiPrime
            | FnId::BesselI0
            | FnId::BesselI1
            | FnId::BesselIn(_)
            | FnId::BesselK0
            | FnId::BesselK1
            | FnId::BesselKn(_)
    )
}

/// The backend name that should appear in the status row's
/// `oracle` field for a given `FnId`. The runner uses this to set
/// `StatusRow::oracle` per row rather than recording the
/// dispatcher's name.
pub fn oracle_name_for(f: FnId) -> &'static str {
    if is_arb_primary(f) {
        "Arb"
    } else {
        "MPFR"
    }
}
