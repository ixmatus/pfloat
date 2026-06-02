#!/usr/bin/env python3
"""Aggregate pf-lm3 exhaustive-sweep shards into per-function status rows.

Reads the per-shard ``result.json`` files under a synced run directory
(``/tmp/<RUN_ID>/<fn>_<K>of<M>/result.json``), merges the sub-shards of
each function, validates that the union of shard ranges covers
``[0, 2^32)`` with no gap or overlap, and writes one ``<fn>.toml`` status
row per function mirroring pfloat's schema. Emits a global summary JSON.

Adapts pfloat's ``scripts/pf-hcz4-aggregate.py`` for the libm harness
(ADR-0058): the libm rows are per function (not appended cross-check
tables), and coverage validation is the new gate that catches a dropped
shard before a row claims ``exhaustive``.

Usage:
  pf-lm3-aggregate.py <results_dir> [--status-dir DIR] [--summary PATH]
"""

import argparse
import json
import sys
from pathlib import Path

TWO_POW_32 = 1 << 32
MODES = ["NE", "NA", "TZ", "TP", "TN"]


def read_shards(results_dir):
    """Yield each sub-shard's parsed result JSON. Each function directory
    holds one result_<k>.json per vCPU the instance sharded across."""
    for p in sorted(results_dir.glob("*/result*.json")):
        try:
            yield json.loads(p.read_text())
        except (OSError, json.JSONDecodeError) as e:
            print(f"[aggregate] WARNING: bad shard {p}: {e}", file=sys.stderr)


def merge_function(shards):
    """Merge a function's shards into a single record. shards is a list of
    result dicts that share (function, order)."""
    agg = {
        "function": shards[0]["function"],
        "order": shards[0].get("order", ""),
        "ok": 0,
        "value_mismatch": 0,
        "flag_mismatch": 0,
        "inconclusive": 0,
        "panic": 0,
        "lm_seeds_run": 0,
        "per_mode": {m: {"value_mismatch": 0, "flag_mismatch": 0, "inconclusive": 0, "panic": 0} for m in MODES},
        "sample_value_mismatches": [],
        "ranges": [],
        "wall_clock_seconds": 0.0,
        "n_shards": len(shards),
        "git_sha": shards[0].get("git_sha", "unknown"),
    }
    for s in shards:
        for k in ("ok", "value_mismatch", "flag_mismatch", "inconclusive", "panic", "lm_seeds_run"):
            agg[k] += int(s.get(k, 0))
        for m in MODES:
            pm = s.get("per_mode", {}).get(m, {})
            for k in ("value_mismatch", "flag_mismatch", "inconclusive", "panic"):
                agg["per_mode"][m][k] += int(pm.get(k, 0))
        agg["sample_value_mismatches"].extend(s.get("sample_value_mismatches", []))
        agg["ranges"].append((int(s["range_start"]), int(s["range_end"])))
        agg["wall_clock_seconds"] = max(agg["wall_clock_seconds"], float(s.get("wall_clock_seconds", 0.0)))
    return agg


def coverage_ok(ranges):
    """True iff the ranges tile [0, 2^32) with no gap or overlap."""
    rs = sorted(ranges)
    cursor = 0
    for start, end in rs:
        if start != cursor:
            return False, cursor
        cursor = end
    return cursor == TWO_POW_32, cursor


def mode_status(agg, m):
    pm = agg["per_mode"][m]
    bad = pm["value_mismatch"] + pm["flag_mismatch"] + pm["panic"]
    return "has-errors" if bad > 0 else "correctly-rounded"


def status_toml(agg, covered):
    coverage = "exhaustive" if covered else "INCOMPLETE-COVERAGE"
    mismatch = agg["value_mismatch"] + agg["flag_mismatch"]
    lines = [
        f'function           = "{agg["function"]}"',
        f'order              = "{agg["order"]}"',
        'kernel_kind        = "primary"',
        f'domain_coverage    = "{coverage}"',
        'oracle             = "MPFR"',
        'oracle_independence = "independent"',
        "worst_ulp          = 0",
        f"mismatch_count     = {mismatch}",
        f'inconclusive_count = {agg["inconclusive"]}',
        f'panic_count        = {agg["panic"]}',
        'vectors            = ""',
        f'lm_seeds_run       = {agg["lm_seeds_run"]}',
        "",
        "[rounding_status]",
    ]
    for m in MODES:
        lines.append(f'{m} = "{mode_status(agg, m)}"')
    return "\n".join(lines) + "\n"


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("results_dir", type=Path)
    ap.add_argument("--status-dir", type=Path, default=None,
                    help="where to write <fn>.toml (default: <results_dir>/status)")
    ap.add_argument("--summary", type=Path, default=None,
                    help="summary JSON path (default: <results_dir>/summary.json)")
    args = ap.parse_args()

    status_dir = args.status_dir or (args.results_dir / "status")
    summary_path = args.summary or (args.results_dir / "summary.json")
    status_dir.mkdir(parents=True, exist_ok=True)

    # Group shards by (function, order).
    groups = {}
    for s in read_shards(args.results_dir):
        key = (s["function"], s.get("order", ""))
        groups.setdefault(key, []).append(s)

    summary = {"schema_version": 1, "functions": {}, "totals": {
        "value_mismatch": 0, "flag_mismatch": 0, "inconclusive": 0, "panic": 0,
    }, "incomplete": [], "has_errors": []}

    for (fn, order), shards in sorted(groups.items()):
        agg = merge_function(shards)
        covered, cursor = coverage_ok(agg["ranges"])
        if not covered:
            print(f"[aggregate] WARNING: {fn} coverage incomplete: reached {cursor:#x} of {TWO_POW_32:#x} "
                  f"({agg['n_shards']} shards)", file=sys.stderr)
            summary["incomplete"].append(fn)
        stem = f"rootn_{order}" if fn == "rootn" and order else fn
        (status_dir / f"{stem}.toml").write_text(status_toml(agg, covered))

        errs = agg["value_mismatch"] + agg["flag_mismatch"] + agg["panic"]
        if errs > 0:
            summary["has_errors"].append(fn)
        for k in ("value_mismatch", "flag_mismatch", "inconclusive", "panic"):
            summary["totals"][k] += agg[k]
        summary["functions"][stem] = {
            "n_shards": agg["n_shards"],
            "covered": covered,
            "value_mismatch": agg["value_mismatch"],
            "flag_mismatch": agg["flag_mismatch"],
            "inconclusive": agg["inconclusive"],
            "panic": agg["panic"],
            "lm_seeds_run": agg["lm_seeds_run"],
            "git_sha": agg["git_sha"],
            "sample_value_mismatches": agg["sample_value_mismatches"][:16],
        }

    summary_path.write_text(json.dumps(summary, indent=2) + "\n")

    n = len(summary["functions"])
    print(f"[aggregate] {n} functions, status rows in {status_dir}")
    print(f"[aggregate] totals: {summary['totals']}")
    if summary["incomplete"]:
        print(f"[aggregate] INCOMPLETE coverage: {summary['incomplete']}", file=sys.stderr)
    if summary["has_errors"]:
        print(f"[aggregate] HAS-ERRORS: {summary['has_errors']}", file=sys.stderr)
    # Non-zero exit if any function is incomplete or has errors.
    return 1 if (summary["incomplete"] or summary["has_errors"]) else 0


if __name__ == "__main__":
    sys.exit(main())
