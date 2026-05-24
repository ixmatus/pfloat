//! ADR-0035 three-way agreement: Arb + mpmath + Maxima all certify
//! the same `f32` for every entry in the pinned corpus. This is
//! the load-bearing cross-check: three independent libraries with
//! no shared code lineage agreeing on every sampled input is the
//! strongest empirical evidence short of formal proof that the
//! certified value is correct.
//!
//! Use mode: sampling layer per ADR-0035 Tier 6. Runs only on the
//! pinned corpus (~25 entries currently) rather than the full f32
//! sweep, because Maxima's per-request cost (~500ms-1s including
//! the nix-shell + Maxima startup) makes full-sweep coverage
//! impractical. The full sweep relies on Arb + mpmath two-oracle
//! agreement (`tests/oracle_arb_mpmath_agreement.rs`); Maxima
//! corroborates the high-value cases the pinned corpus enumerates.
//!
//! Maxima coverage gaps the worker handles by INC: very small
//! subnormals on `bessel_i` trigger "Exceeded maximum allowed
//! fpprec" in Maxima's internal hypergeometric. The test treats
//! Maxima INC as "abstain" rather than "disagree" so Arb + mpmath
//! agreement still passes when Maxima can't reach a verdict; only
//! a genuine `OK <bits>` from Maxima participates in the
//! agreement check.
//!
//! Cost: ~25 Maxima requests at ~1s each = ~25 seconds. Skipped
//! from per-push (gated on `differential-arb`); runs at
//! slice-close cadence.

#![cfg(all(unix, feature = "differential-arb"))]

#[path = "oracle/mod.rs"]
mod oracle;

use std::path::{Path, PathBuf};

use oracle::{ArbOracle, Enclosure, FnId, MaximaOracle, MpmathOracle, OracleBackend};
use pfloat::RoundingMode;
use rug::float::Round;

fn fnid_from_filename(name: &str) -> Option<FnId> {
    let stem = name.strip_suffix(".toml")?;
    match stem {
        "Si" => Some(FnId::Si),
        "Ci" => Some(FnId::Ci),
        "li" => Some(FnId::Li),
        "Bi" => Some(FnId::Bi),
        "Ai_prime" => Some(FnId::AiPrime),
        "Bi_prime" => Some(FnId::BiPrime),
        "BesselI0" => Some(FnId::BesselI0),
        "BesselI1" => Some(FnId::BesselI1),
        "BesselK0" => Some(FnId::BesselK0),
        "BesselK1" => Some(FnId::BesselK1),
        _ => None,
    }
}

fn parse_mode(s: &str) -> Option<RoundingMode> {
    match s {
        "NE" => Some(RoundingMode::NearestEven),
        "RNA" => Some(RoundingMode::NearestAway),
        "RZ" => Some(RoundingMode::TowardZero),
        "RP" => Some(RoundingMode::TowardPositive),
        "RM" => Some(RoundingMode::TowardNegative),
        _ => None,
    }
}

fn parse_hex_u32(s: &str) -> Option<u32> {
    let stripped = s.strip_prefix("0x").or_else(|| s.strip_prefix("0X"))?;
    u32::from_str_radix(stripped, 16).ok()
}

#[derive(Debug)]
struct PinEntry {
    input_bits: u32,
    mode: RoundingMode,
}

fn parse_pin_file(path: &Path) -> Vec<PinEntry> {
    let text = std::fs::read_to_string(path).expect("read pin file");
    let mut entries = Vec::new();
    let mut current: Option<(Option<u32>, Option<RoundingMode>)> = None;
    for line in text.lines() {
        let line = line.trim();
        if line.is_empty() || line.starts_with('#') {
            continue;
        }
        if line.starts_with("[[entry]]") {
            if let Some((Some(i), Some(m))) = current {
                entries.push(PinEntry {
                    input_bits: i,
                    mode: m,
                });
            }
            current = Some((None, None));
            continue;
        }
        let parts: Vec<&str> = line.splitn(2, '=').collect();
        if parts.len() != 2 {
            continue;
        }
        let key = parts[0].trim();
        let val = parts[1]
            .trim()
            .strip_prefix('"')
            .and_then(|s| s.strip_suffix('"'))
            .unwrap_or(parts[1].trim());
        match key {
            "input_bits" => {
                if let Some(ref mut c) = current {
                    c.0 = parse_hex_u32(val);
                }
            }
            "mode" => {
                if let Some(ref mut c) = current {
                    c.1 = parse_mode(val);
                }
            }
            _ => {}
        }
    }
    if let Some((Some(i), Some(m))) = current {
        entries.push(PinEntry {
            input_bits: i,
            mode: m,
        });
    }
    entries
}

/// Result of an authoritative oracle call: either the certified
/// `f32` bit pattern, or `None` if the worker returned `INC`
/// (NaN enclosure).
fn certified_bits_or_none(enc: &Enclosure) -> Option<u32> {
    if enc.lo.is_nan() && enc.hi.is_nan() {
        return None;
    }
    let lo_f32 = enc.lo.to_f32_round(Round::Nearest);
    let hi_f32 = enc.hi.to_f32_round(Round::Nearest);
    assert_eq!(lo_f32.to_bits(), hi_f32.to_bits());
    Some(lo_f32.to_bits())
}

#[test]
#[ignore = "slow: Maxima nix-shell startup adds ~1s per request; run at slice-close cadence"]
fn arb_mpmath_maxima_three_way_agreement_on_pinned_corpus() {
    let pinned_dir = PathBuf::from(env!("CARGO_MANIFEST_DIR")).join("tests/oracle/pinned");
    let arb = ArbOracle::new().expect("ArbOracle::new");
    let mpm = MpmathOracle::new().expect("MpmathOracle::new");
    let maxima = MaximaOracle::new().expect("MaximaOracle::new (needs nix-shell + maxima)");

    let mut total_entries = 0u32;
    let mut maxima_inc_count = 0u32;
    let mut three_way_ok = 0u32;
    let mut divergences: Vec<(FnId, u32, RoundingMode, u32, u32, Option<u32>)> = Vec::new();

    let read_dir = std::fs::read_dir(&pinned_dir).expect("read pinned dir");
    for entry in read_dir {
        let path = entry.expect("dir entry").path();
        let name = match path.file_name().and_then(|n| n.to_str()) {
            Some(n) => n,
            None => continue,
        };
        if !name.to_ascii_lowercase().ends_with(".toml") {
            continue;
        }
        let Some(f) = fnid_from_filename(name) else {
            continue;
        };
        for pin in parse_pin_file(&path) {
            total_entries += 1;
            let arb_bits = certified_bits_or_none(&arb.enclose(f, pin.input_bits, pin.mode, 64))
                .expect("Arb worker must certify (pinned corpus is curated)");
            let mpm_bits = certified_bits_or_none(&mpm.enclose(f, pin.input_bits, pin.mode, 64))
                .expect("mpmath worker must certify (pinned corpus is curated)");
            let maxima_bits =
                certified_bits_or_none(&maxima.enclose(f, pin.input_bits, pin.mode, 64));

            if maxima_bits.is_none() {
                maxima_inc_count += 1;
                // Maxima abstained; Arb + mpmath agreement still
                // counts, but only as two-way at this entry.
                if arb_bits != mpm_bits {
                    divergences.push((
                        f,
                        pin.input_bits,
                        pin.mode,
                        arb_bits,
                        mpm_bits,
                        maxima_bits,
                    ));
                }
                continue;
            }
            let m_bits = maxima_bits.unwrap();
            if arb_bits == mpm_bits && mpm_bits == m_bits {
                three_way_ok += 1;
            } else {
                divergences.push((
                    f,
                    pin.input_bits,
                    pin.mode,
                    arb_bits,
                    mpm_bits,
                    Some(m_bits),
                ));
            }
        }
    }

    eprintln!(
        "[three-way] checked {total_entries} pinned entries; \
         three-way agreement on {three_way_ok}; \
         Maxima INC abstentions: {maxima_inc_count}; \
         divergences: {}",
        divergences.len()
    );
    for (f, input, mode, arb_bits, mpm_bits, maxima_bits) in &divergences {
        eprintln!(
            "[three-way] DIVERGENCE: {f:?} input={input:#010x} mode={mode:?} \
             arb={arb_bits:#010x} mpmath={mpm_bits:#010x} maxima={maxima_bits:?}"
        );
    }

    assert!(
        divergences.is_empty(),
        "{} pinned entries showed oracle divergence; see eprintln above",
        divergences.len()
    );
    assert!(
        total_entries > 0,
        "no pinned entries found; check tests/oracle/pinned/ layout"
    );
}
