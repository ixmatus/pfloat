//! ADR-0035 Tier 5: pinned worker-output corpus diff check.
//!
//! Reads every TOML file in `tests/oracle/pinned/` and verifies the
//! live Arb worker emits the certified `f32` bit pattern each entry
//! pins. Any divergence halts with a diagnostic listing the
//! `(FnId, input_bits, mode, pinned_bits, live_bits)` of the
//! divergence.
//!
//! The pin file is the durable contract. A push that intentionally
//! changes worker output (an Arb upgrade, a routine fix, a
//! function-coverage extension) must update the pin file with a
//! provenance note explaining the change; without the pin update
//! the gate fails and the push is blocked. This converts "trust
//! the worker" into "review the diff."
//!
//! Cost: ~25 worker calls (one per pinned entry); ~1 second debug.
//! Cheap enough for per-slice; could plausibly fit in per-push CI
//! if the venv were a CI standard, but currently gated on
//! `differential-arb` like the other Arb-using tests.
//!
//! Adding entries: edit the appropriate `<FnId>.toml`. The
//! provenance line must explain how the certified bit pattern was
//! derived (hand-derived from first principles, two-oracle
//! agreement, Maxima triple-check, etc.). See
//! `tests/oracle/pinned/README.md` for the file format and
//! regeneration ritual.

#![cfg(all(unix, feature = "differential-arb"))]

#[path = "oracle/mod.rs"]
mod oracle;

use std::path::{Path, PathBuf};

use oracle::{ArbOracle, Enclosure, FnId, OracleBackend};
use pfloat::RoundingMode;
use rug::float::Round;

/// Parse a `<FnId>.toml` filename to the corresponding `FnId`. The
/// file basename matches the function's display name (the
/// `FnId::name()` value), e.g. `BesselI1.toml` for `FnId::BesselI1`.
/// Parametric Bessel orders are not pinned in this slice (they
/// would need per-order files; deferred to a follow-up).
fn fnid_from_filename(name: &str) -> Option<FnId> {
    // Strip ".toml" suffix.
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
    certified_bits: u32,
}

/// Parse a TOML file into a vector of pin entries. Minimal
/// hand-written parser scoped to the fixed schema we use; avoids
/// pulling in a TOML dep just for this test.
fn parse_pin_file(path: &Path) -> Result<Vec<PinEntry>, String> {
    let text =
        std::fs::read_to_string(path).map_err(|e| format!("read {}: {e}", path.display()))?;
    let mut entries = Vec::new();
    let mut current: Option<(Option<u32>, Option<RoundingMode>, Option<u32>)> = None;

    for (line_no, raw_line) in text.lines().enumerate() {
        let line = raw_line.trim();
        // Skip comments, blanks, and the """ multiline boundary lines.
        if line.is_empty() || line.starts_with('#') {
            continue;
        }
        if line.starts_with("[[entry]]") {
            // Close any open entry, then open a new one.
            if let Some((Some(i), Some(m), Some(c))) = current {
                entries.push(PinEntry {
                    input_bits: i,
                    mode: m,
                    certified_bits: c,
                });
            }
            current = Some((None, None, None));
            continue;
        }
        // Skip multi-line strings: we only handle the simple
        // `key = "value"` form. The `provenance = """..."""`
        // block is multi-line; we recognize lines that don't
        // start with one of our known keys and skip them.
        let parts: Vec<&str> = line.splitn(2, '=').collect();
        if parts.len() != 2 {
            continue;
        }
        let key = parts[0].trim();
        let val_raw = parts[1].trim();
        // Strip surrounding quotes if present (single-line "...")
        let val = val_raw
            .strip_prefix('"')
            .and_then(|s| s.strip_suffix('"'))
            .unwrap_or(val_raw);
        match key {
            "input_bits" => {
                if let Some(b) = parse_hex_u32(val) {
                    if let Some(ref mut c) = current {
                        c.0 = Some(b);
                    }
                } else {
                    return Err(format!(
                        "{}:{}: invalid input_bits `{val}`",
                        path.display(),
                        line_no + 1
                    ));
                }
            }
            "mode" => {
                if let Some(m) = parse_mode(val) {
                    if let Some(ref mut c) = current {
                        c.1 = Some(m);
                    }
                } else {
                    return Err(format!(
                        "{}:{}: invalid mode `{val}`",
                        path.display(),
                        line_no + 1
                    ));
                }
            }
            "certified_bits" => {
                if let Some(b) = parse_hex_u32(val) {
                    if let Some(ref mut c) = current {
                        c.2 = Some(b);
                    }
                } else {
                    return Err(format!(
                        "{}:{}: invalid certified_bits `{val}`",
                        path.display(),
                        line_no + 1
                    ));
                }
            }
            _ => {
                // Unknown key (likely `provenance` or a multi-line
                // continuation); ignore.
            }
        }
    }
    // Close the final entry.
    if let Some((Some(i), Some(m), Some(c))) = current {
        entries.push(PinEntry {
            input_bits: i,
            mode: m,
            certified_bits: c,
        });
    }
    Ok(entries)
}

/// Extract the f32 bits an authoritative single-point enclosure
/// certifies. Mirrors `oracle_arb_mpmath_agreement.rs`.
fn extract_certified_f32(enc: &Enclosure) -> u32 {
    if enc.lo.is_nan() && enc.hi.is_nan() {
        return f32::NAN.to_bits();
    }
    let lo_f32 = enc.lo.to_f32_round(Round::Nearest);
    let hi_f32 = enc.hi.to_f32_round(Round::Nearest);
    assert_eq!(
        lo_f32.to_bits(),
        hi_f32.to_bits(),
        "authoritative enclosure endpoints disagree on f32: lo={lo_f32}, hi={hi_f32}"
    );
    lo_f32.to_bits()
}

#[test]
fn arb_worker_output_matches_pinned_corpus() {
    let pinned_dir = PathBuf::from(env!("CARGO_MANIFEST_DIR")).join("tests/oracle/pinned");
    let arb = ArbOracle::new()
        .expect("ArbOracle::new (requires the python-flint venv; run scripts/setup_arb_oracle.sh)");

    let mut total_entries = 0u32;
    let mut divergences: Vec<(FnId, u32, RoundingMode, u32, u32, PathBuf)> = Vec::new();

    let read_dir = std::fs::read_dir(&pinned_dir).expect("read tests/oracle/pinned directory");
    for entry in read_dir {
        let entry = entry.expect("read pinned dir entry");
        let path = entry.path();
        let name = path
            .file_name()
            .and_then(|n| n.to_str())
            .expect("filename utf-8");
        // Skip non-TOML files (e.g. README.md).
        if !name.to_ascii_lowercase().ends_with(".toml") {
            continue;
        }
        let Some(f) = fnid_from_filename(name) else {
            continue;
        };
        let pins = match parse_pin_file(&path) {
            Ok(p) => p,
            Err(e) => panic!("parse {}: {e}", path.display()),
        };
        for pin in pins {
            total_entries += 1;
            let enc = arb.enclose(f, pin.input_bits, pin.mode, 64);
            let live_bits = extract_certified_f32(&enc);
            if live_bits != pin.certified_bits {
                divergences.push((
                    f,
                    pin.input_bits,
                    pin.mode,
                    pin.certified_bits,
                    live_bits,
                    path.clone(),
                ));
            }
        }
    }

    if !divergences.is_empty() {
        for (f, input, mode, pinned, live, path) in &divergences {
            eprintln!(
                "[pinned] {f:?} input={input:#010x} mode={mode:?} \
                 pinned={pinned:#010x} live={live:#010x}  ({})",
                path.display()
            );
        }
        panic!(
            "Arb worker output diverged from pin on {} of {} entries; \
             see eprintln above. If the change is intentional, update \
             the corresponding pin file with a provenance note.",
            divergences.len(),
            total_entries
        );
    }

    eprintln!("[pinned] {total_entries} entries checked; all match");
    assert!(
        total_entries > 0,
        "no pin entries found; check tests/oracle/pinned/ layout"
    );
}
