#!/usr/bin/env python3
"""pf-hcz4 result aggregator. ADR-0049.

Reads per-shard `result.json` from the synced S3 directory tree, emits:

  1. A global summary JSON (`/tmp/<RUN_ID>-summary.json` by default)
     with totals + per-FnId violation counts + top-N highest
     ratio_log2 cells across all shards.

  2. Per-FnId `[cross_check]` table appended to existing
     `tests/oracle/status/<fn>.toml` files. This is the v1.0 baseline
     durable artifact per bead pf-hcz4 acceptance #4. Append-only;
     does not modify any existing row in the TOML.

Usage:
    pf-hcz4-aggregate.py <RESULTS_DIR>
        [--emit-summary <PATH>]
        [--emit-status-baseline <STATUS_DIR>]
        [--top-n 10]
        [--no-status-write]

Layout of <RESULTS_DIR>:
    <RESULTS_DIR>/<shard>/result.json
    <RESULTS_DIR>/<shard>/sweep.log
    <RESULTS_DIR>/<shard>/_DONE  (or _FAILED_PREFLIGHT)
"""

from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path


SCHEMA_VERSION = 1


def shard_to_toml_basename(fn_id: str, order):
    """fn_id "Yn" + order 5  →  "Yn_5"; fn_id "exp" + order None → "exp"."""
    if order is None:
        return fn_id
    return f"{fn_id}_{order}"


def read_shard_results(results_dir: Path):
    """Yield (shard_name, payload_dict) for each result.json found."""
    for shard_dir in sorted(results_dir.iterdir()):
        if not shard_dir.is_dir():
            continue
        result_path = shard_dir / "result.json"
        if not result_path.exists():
            print(f"[aggregate] WARNING: {shard_dir.name} has no result.json", file=sys.stderr)
            continue
        with result_path.open() as f:
            payload = json.load(f)
        if payload.get("schema_version") != SCHEMA_VERSION:
            print(
                f"[aggregate] WARNING: {shard_dir.name} schema_version "
                f"{payload.get('schema_version')} != {SCHEMA_VERSION}",
                file=sys.stderr,
            )
        yield shard_dir.name, payload


def aggregate(results_dir: Path, top_n: int):
    """Walk the shard results and build the summary structure."""
    shards = list(read_shard_results(results_dir))
    summary = {
        "schema_version": SCHEMA_VERSION,
        "n_shards": len(shards),
        "totals": {
            "passes": 0,
            "skipped_no_ziv_path": 0,
            "skipped_no_midpoint": 0,
            "skipped_non_finite": 0,
            "skipped_trace_not_final": 0,
            "violations": 0,
        },
        "per_shard": {},
        "top_violations": [],
        "arb_midpoint_calls": 0,
        "mpfr_midpoint_calls": 0,
        "wall_clock_seconds_max": 0.0,
        "wall_clock_seconds_sum": 0.0,
        "git_sha": None,
    }
    all_violations = []
    for name, payload in shards:
        t = payload["totals"]
        summary["totals"]["passes"] += t["passes"]
        summary["totals"]["skipped_no_ziv_path"] += t["skipped_no_ziv_path"]
        summary["totals"]["skipped_no_midpoint"] += t["skipped_no_midpoint"]
        summary["totals"]["skipped_non_finite"] += t["skipped_non_finite"]
        summary["totals"]["skipped_trace_not_final"] += t.get("skipped_trace_not_final", 0)
        summary["totals"]["violations"] += t["violations"]
        summary["arb_midpoint_calls"] += payload.get("arb_midpoint_calls", 0)
        summary["mpfr_midpoint_calls"] += payload.get("mpfr_midpoint_calls", 0)
        w = float(payload.get("wall_clock_seconds", 0.0))
        summary["wall_clock_seconds_sum"] += w
        summary["wall_clock_seconds_max"] = max(summary["wall_clock_seconds_max"], w)
        if summary["git_sha"] is None:
            summary["git_sha"] = payload.get("git_sha")
        summary["per_shard"][name] = {
            "fn_id": payload["fn_id"],
            "order": payload["order"],
            "passes": t["passes"],
            "violations": t["violations"],
            "wall_clock_seconds": w,
        }
        all_violations.extend((name, v) for v in payload.get("violations", []))

    # Top-N violations by ratio_log2 descending.
    all_violations.sort(key=lambda nv: nv[1].get("ratio_log2", 0.0), reverse=True)
    summary["top_violations"] = [
        {"shard": name, **v} for name, v in all_violations[:top_n]
    ]
    return summary


def emit_status_baseline(results_dir: Path, status_dir: Path):
    """For each shard, append a [cross_check] table to the
    matching tests/oracle/status/<basename>.toml. Append-only:
    existing tables are untouched. Idempotent in the sense that
    re-running adds a duplicate [cross_check] block, so a guard
    rejects re-application if the marker is already present."""
    count_added = 0
    count_skipped = 0
    for shard_name, payload in read_shard_results(results_dir):
        basename = shard_to_toml_basename(payload["fn_id"], payload["order"])
        toml_path = status_dir / f"{basename}.toml"
        if not toml_path.exists():
            print(f"[aggregate] WARNING: no {toml_path}; skipping", file=sys.stderr)
            continue
        existing = toml_path.read_text()
        if "[cross_check]" in existing:
            print(f"[aggregate] {basename}: [cross_check] already present, skipping", file=sys.stderr)
            count_skipped += 1
            continue
        table = render_cross_check_table(payload)
        with toml_path.open("a") as f:
            if not existing.endswith("\n"):
                f.write("\n")
            f.write("\n")
            f.write(table)
        count_added += 1

    print(f"[aggregate] {count_added} status TOMLs appended; {count_skipped} skipped (already had [cross_check])")


def render_cross_check_table(payload):
    """Render the [cross_check] TOML block from a shard's result.json
    payload. Captures the verifiable v1.0 baseline data: pass / skip
    / violation counts, the error_guard constant, the run identifying
    metadata."""
    t = payload["totals"]
    n_violations = t["violations"]
    lines = [
        "[cross_check]",
        f'# pf-hcz4 cross-check baseline. ADR-0049.',
        f'git_sha = "{payload.get("git_sha", "unknown")}"',
        f'error_guard_const = {payload.get("error_guard_const", 24)}',
        f'sample = {payload.get("sample", 65536)}',
        f'lm_seeds_used = {payload.get("lm_seeds_used", 0)}',
        f'passes = {t["passes"]}',
        f'skipped_no_ziv_path = {t["skipped_no_ziv_path"]}',
        f'skipped_no_midpoint = {t["skipped_no_midpoint"]}',
        f'skipped_non_finite = {t["skipped_non_finite"]}',
        f'skipped_trace_not_final = {t.get("skipped_trace_not_final", 0)}',
        f'violations = {n_violations}',
        f'wall_clock_seconds = {payload.get("wall_clock_seconds", 0.0)}',
        f'instance_arch = "{payload.get("instance_arch", "unknown")}"',
        f'instance_type = "{payload.get("instance_type", "unknown")}"',
    ]
    # If any violations, embed the top 3 ratio_log2 entries inline as
    # commented diagnostics; the full per-mode violation list lives
    # in the per-shard result.json (preserved in tests/vectors/ or
    # the S3 archive).
    if n_violations > 0:
        violations = sorted(
            payload.get("violations", []),
            key=lambda v: v.get("ratio_log2", 0.0),
            reverse=True,
        )[:3]
        lines.append("# Top 3 violations (full list in shard result.json):")
        for v in violations:
            lines.append(
                f'# input=0x{v["input_u32"]:08x} mode={v["mode"]} '
                f'working={v["working_prec"]} ratio_log2={v["ratio_log2"]:.3f}'
            )
    return "\n".join(lines) + "\n"


def main():
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument("results_dir", type=Path)
    ap.add_argument("--emit-summary", type=Path, default=None)
    ap.add_argument("--emit-status-baseline", type=Path, default=None)
    ap.add_argument("--top-n", type=int, default=10)
    ap.add_argument("--no-status-write", action="store_true")
    args = ap.parse_args()

    if not args.results_dir.is_dir():
        print(f"error: {args.results_dir} not a directory", file=sys.stderr)
        sys.exit(1)

    summary = aggregate(args.results_dir, args.top_n)

    # Default summary path.
    summary_path = args.emit_summary
    if summary_path is None:
        summary_path = args.results_dir.parent / f"{args.results_dir.name}-summary.json"
    with summary_path.open("w") as f:
        json.dump(summary, f, indent=2)
    print(f"[aggregate] summary → {summary_path}")
    print(
        f"[aggregate] totals: passes={summary['totals']['passes']:,} "
        f"violations={summary['totals']['violations']:,} "
        f"n_shards={summary['n_shards']}/63 "
        f"arb_calls={summary['arb_midpoint_calls']:,} "
        f"mpfr_calls={summary['mpfr_midpoint_calls']:,}"
    )
    if summary["top_violations"]:
        print(f"[aggregate] top violations (ratio_log2 desc):")
        for v in summary["top_violations"]:
            print(
                f"  {v['shard']:<14} mode={v['mode']:<14} "
                f"input=0x{v['input_u32']:08x} "
                f"ratio_log2={v['ratio_log2']:+.3f}"
            )

    if args.emit_status_baseline and not args.no_status_write:
        emit_status_baseline(args.results_dir, args.emit_status_baseline)


if __name__ == "__main__":
    main()
