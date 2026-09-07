#!/usr/bin/env python3
"""Run the Fast-SCC module-overhead experiment and write averaged CSV data."""

from __future__ import annotations

import argparse
import csv
import os
import re
import statistics
import subprocess
import sys
from pathlib import Path


SCRIPT_DIR = Path(__file__).resolve().parent
ROOT = SCRIPT_DIR.parent
FAST_SRC = ROOT / "Fast-SCC" / "src"
PROFILE_SOURCE = FAST_SRC / "scc_paper_swap_overhead.cpp"
PROFILE_BINARY = FAST_SRC / "scc_paper_swap_overhead"
BASELINE_BINARY = ROOT / "GBBS-VGC" / "src" / "scc"
DATA_DIR = ROOT / "data"
RESULTS_DIR = ROOT / "results"
DEFAULT_CSV = RESULTS_DIR / "overhead.csv"
LOG_DIR = RESULTS_DIR / "overhead-logs"

GRAPH_FILES = {
    "CH5": "CHEM_5.bin",
    "GL2": "GeoLifeNoScale_2.bin",
    "GL5": "GeoLifeNoScale_5.bin",
    "GL10": "GeoLifeNoScale_10.bin",
    "GL15": "GeoLifeNoScale_15.bin",
    "GL20": "GeoLifeNoScale_20.bin",
}
GRAPH_ORDER = tuple(GRAPH_FILES)

CSV_FIELDS = (
    "graph",
    "gbbs_vgc_seconds",
    "fast_scc_seconds",
    "fast_core_seconds",
    "fake_link_seconds",
    "pivot_selection_seconds",
    "fake_link_scan_seconds",
    "fake_link_scope_seconds",
    "fake_link_partition_seconds",
    "fake_link_identification_calls",
    "fake_link_scope_updates",
    "fake_link_edges_processed",
    "fake_link_profile_samples",
    "pivot_candidates_examined",
    "pivots_selected",
    "scc_count",
    "correctness",
)


def build_profiler() -> None:
    command = [
        "g++",
        "-std=c++17",
        "-Wall",
        "-Wextra",
        "-Werror",
        "-O3",
        "-mcx16",
        "-march=native",
        "-pthread",
        "-I../parlaylib/include/",
        PROFILE_SOURCE.name,
        "-o",
        PROFILE_BINARY.name,
    ]
    print("[build]", " ".join(command), flush=True)
    subprocess.run(command, cwd=FAST_SRC, check=True)


def run_and_log(
    command: list[str], log_path: Path, environment: dict[str, str]
) -> str:
    print("[run]", " ".join(command), flush=True)
    completed = subprocess.run(
        command,
        cwd=ROOT,
        env=environment,
        text=True,
        stdout=subprocess.PIPE,
        stderr=subprocess.STDOUT,
        check=False,
    )
    output = completed.stdout
    log_path.write_text(
        "$ " + " ".join(command) + "\n\n" + output,
        encoding="utf-8",
    )
    if completed.returncode != 0:
        tail = "\n".join(output.splitlines()[-20:])
        raise RuntimeError(
            f"command failed with exit code {completed.returncode}; "
            f"see {log_path}\n{tail}"
        )
    return output


def parse_baseline(output: str, expected_runs: int) -> tuple[list[float], int]:
    costs = [
        float(value)
        for value in re.findall(
            r"^scc cost:\s*([0-9.eE+-]+)\s*$", output, flags=re.MULTILINE
        )
    ]
    if len(costs) != expected_runs:
        raise ValueError(
            f"expected {expected_runs} GBBS-VGC timings, found {len(costs)}"
        )
    count_matches = re.findall(r"^n_scc\s*=\s*(\d+)\s*$", output, re.MULTILINE)
    if not count_matches:
        raise ValueError("GBBS-VGC output does not contain an SCC count")
    return costs, int(count_matches[-1])


def parse_key_values(line: str) -> dict[str, str]:
    values: dict[str, str] = {}
    for item in line.split()[1:]:
        if "=" in item:
            key, value = item.split("=", 1)
            values[key] = value
    return values


def parse_profiler(output: str, expected_runs: int) -> list[dict[str, str]]:
    records = [
        parse_key_values(line)
        for line in output.splitlines()
        if line.startswith("OVERHEAD_RESULT ")
    ]
    if len(records) != expected_runs:
        raise ValueError(
            f"expected {expected_runs} Fast-SCC profiles, found {len(records)}"
        )
    required = {
        "fast_total_seconds",
        "fast_core_seconds",
        "fake_link_seconds",
        "fake_link_scan_seconds",
        "fake_link_scope_seconds",
        "fake_link_partition_seconds",
        "pivot_selection_seconds",
        "fake_link_identification_calls",
        "fake_link_scope_updates",
        "fake_link_edges_processed",
        "fake_link_profile_samples",
        "pivot_candidates_examined",
        "pivots_selected",
        "scc_count",
        "correctness",
    }
    for record in records:
        missing = required.difference(record)
        if missing:
            raise ValueError(f"profile record is missing: {sorted(missing)}")
    return records


def mean_field(records: list[dict[str, str]], key: str) -> float:
    return statistics.fmean(float(record[key]) for record in records)


def make_row(
    graph: str,
    baseline_costs: list[float],
    baseline_scc_count: int,
    profiles: list[dict[str, str]],
) -> dict[str, object]:
    profile_counts = {int(record["scc_count"]) for record in profiles}
    internal_pass = all(record["correctness"] == "PASS" for record in profiles)
    baseline_pass = profile_counts == {baseline_scc_count}
    correctness = "PASS" if internal_pass and baseline_pass else "FAIL"

    fast_total = mean_field(profiles, "fast_total_seconds")
    core = mean_field(profiles, "fast_core_seconds")
    fake_link = mean_field(profiles, "fake_link_seconds")
    pivot = mean_field(profiles, "pivot_selection_seconds")
    if fake_link + pivot > fast_total * 1.02:
        raise ValueError(
            f"{graph}: measured module overhead exceeds Fast-SCC total time; "
            "the profile is not reliable"
        )

    return {
        "graph": graph,
        "gbbs_vgc_seconds": statistics.fmean(baseline_costs),
        "fast_scc_seconds": fast_total,
        "fast_core_seconds": core,
        "fake_link_seconds": fake_link,
        "pivot_selection_seconds": pivot,
        "fake_link_scan_seconds": mean_field(
            profiles, "fake_link_scan_seconds"
        ),
        "fake_link_scope_seconds": mean_field(
            profiles, "fake_link_scope_seconds"
        ),
        "fake_link_partition_seconds": mean_field(
            profiles, "fake_link_partition_seconds"
        ),
        "fake_link_identification_calls": mean_field(
            profiles, "fake_link_identification_calls"
        ),
        "fake_link_scope_updates": mean_field(
            profiles, "fake_link_scope_updates"
        ),
        "fake_link_edges_processed": mean_field(
            profiles, "fake_link_edges_processed"
        ),
        "fake_link_profile_samples": mean_field(
            profiles, "fake_link_profile_samples"
        ),
        "pivot_candidates_examined": mean_field(
            profiles, "pivot_candidates_examined"
        ),
        "pivots_selected": mean_field(profiles, "pivots_selected"),
        "scc_count": baseline_scc_count,
        "correctness": correctness,
    }


def write_csv(path: Path, rows: list[dict[str, object]]) -> None:
    temporary = path.with_suffix(path.suffix + ".tmp")
    with temporary.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=CSV_FIELDS)
        writer.writeheader()
        for row in rows:
            formatted = dict(row)
            for field in (
                "gbbs_vgc_seconds",
                "fast_scc_seconds",
                "fast_core_seconds",
                "fake_link_seconds",
                "pivot_selection_seconds",
                "fake_link_scan_seconds",
                "fake_link_scope_seconds",
                "fake_link_partition_seconds",
            ):
                formatted[field] = f"{float(row[field]):.9f}"
            for field in (
                "fake_link_identification_calls",
                "fake_link_scope_updates",
                "fake_link_edges_processed",
                "fake_link_profile_samples",
                "pivot_candidates_examined",
                "pivots_selected",
            ):
                formatted[field] = f"{float(row[field]):.3f}"
            writer.writerow(formatted)
    temporary.replace(path)


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Measure Fake-link and Pivot-selection overhead."
    )
    parser.add_argument("--runs", type=int, default=10)
    parser.add_argument("--threads", type=int, default=24)
    parser.add_argument(
        "--graphs", nargs="+", choices=GRAPH_ORDER, default=list(GRAPH_ORDER)
    )
    parser.add_argument("--output", type=Path, default=DEFAULT_CSV)
    parser.add_argument("--skip-build", action="store_true")
    return parser.parse_args()


def main() -> int:
    args = parse_args()
    if args.runs < 1 or args.threads < 1:
        raise ValueError("--runs and --threads must be positive")

    args.output.parent.mkdir(parents=True, exist_ok=True)
    LOG_DIR.mkdir(parents=True, exist_ok=True)
    if not args.skip_build:
        build_profiler()

    required = [PROFILE_BINARY, BASELINE_BINARY]
    required.extend(DATA_DIR / GRAPH_FILES[graph] for graph in args.graphs)
    missing = [str(path) for path in required if not path.exists()]
    if missing:
        raise FileNotFoundError("missing experiment inputs:\n" + "\n".join(missing))

    environment = os.environ.copy()
    environment["PARLAY_NUM_THREADS"] = str(args.threads)
    rows: list[dict[str, object]] = []

    for graph in args.graphs:
        graph_path = DATA_DIR / GRAPH_FILES[graph]
        common = [str(graph_path), "-t", str(args.runs), "-local_reach", "-local_scc"]
        baseline_output = run_and_log(
            [str(BASELINE_BINARY), *common, "-status"],
            LOG_DIR / f"{graph}-gbbs-vgc.log",
            environment,
        )
        profile_output = run_and_log(
            [str(PROFILE_BINARY), *common],
            LOG_DIR / f"{graph}-fast-scc-profile.log",
            environment,
        )

        baseline_costs, baseline_count = parse_baseline(
            baseline_output, args.runs
        )
        profiles = parse_profiler(profile_output, args.runs)
        row = make_row(graph, baseline_costs, baseline_count, profiles)
        rows.append(row)
        write_csv(args.output, rows)
        print(
            f"[{graph}] GBBS-VGC={row['gbbs_vgc_seconds']:.6f}s, "
            f"Fast-SCC={row['fast_scc_seconds']:.6f}s, "
            f"FL={row['fake_link_seconds']:.6f}s, "
            f"PS={row['pivot_selection_seconds']:.6f}s, "
            f"FL edges={row['fake_link_edges_processed']:.0f}, "
            f"pivot candidates={row['pivot_candidates_examined']:.0f}, "
            f"correctness={row['correctness']}",
            flush=True,
        )
        if row["correctness"] != "PASS":
            raise RuntimeError(f"{graph}: SCC result does not match GBBS-VGC")

    print(f"[done] wrote {args.output}")
    return 0


if __name__ == "__main__":
    try:
        raise SystemExit(main())
    except (FileNotFoundError, RuntimeError, ValueError) as error:
        print(f"error: {error}", file=sys.stderr)
        raise SystemExit(1)
