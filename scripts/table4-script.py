#!/usr/bin/env python3
"""Run Table 4 experiments and fill average cost values into table4.csv.

Default usage on the server:
  cd /home/thinker/fhz/Fast-SCC-revision
  python3 scripts/table4-script.py

The script runs datasets row by row. For each dataset it runs the four
methods in Table 4 order, extracts "average cost" from the program output
or averages GBBS "time per iter" / "Running Time" lines, keeps three decimals,
calculates Fast-SCC's rank among the four methods, and writes partial results
to results/table4.csv after every successful run.
"""

from __future__ import annotations

import argparse
import csv
import os
import re
import shlex
import subprocess
import sys
from decimal import Decimal, InvalidOperation, ROUND_HALF_UP
from pathlib import Path


METHODS = ("Fast-SCC", "GBBS-VGC", "GBBS", "Multistep")

DATASETS = (
    ("LJ", "soc-LiveJournal1.bin"),
    ("HH5", "Household.lines_5.bin"),
    ("CH5", "CHEM_5.bin"),
    ("GL2", "GeoLifeNoScale_2.bin"),
    ("GL5", "GeoLifeNoScale_5.bin"),
    ("GL10", "GeoLifeNoScale_10.bin"),
    ("GL15", "GeoLifeNoScale_15.bin"),
    ("GL20", "GeoLifeNoScale_20.bin"),
    ("SQR", "grid_4000_4000.bin"),
    ("REC", "grid_1000_10000.bin"),
)

AVERAGE_COST_RE = re.compile(
    r"\baverage\s+cost\b\s*[:=]?\s*([-+]?(?:\d+(?:\.\d*)?|\.\d+)(?:[eE][-+]?\d+)?)",
    re.IGNORECASE,
)
RUNNING_TIME_RE = re.compile(
    r"\brunning\s+time\b\s*[:=]?\s*([-+]?(?:\d+(?:\.\d*)?|\.\d+)(?:[eE][-+]?\d+)?)",
    re.IGNORECASE,
)
TIME_PER_ITER_RE = re.compile(
    r"\btime\s+per\s+iter\b\s*[:=]?\s*([-+]?(?:\d+(?:\.\d*)?|\.\d+)(?:[eE][-+]?\d+)?)",
    re.IGNORECASE,
)


def command_for(method: str, graph_path: str) -> list[str]:
    if method == "Fast-SCC":
        return [
            # "Fast-SCC/src/scc",
            "Fast-SCC/src/scc_paper_swap",
            graph_path,
            "-local_reach",
            "-local_scc",
            "-status",
            "-t",
            "100",
        ]
    if method == "GBBS-VGC":
        return [
            "GBBS-VGC/src/scc",
            graph_path,
            "-local_reach",
            "-local_scc",
            "-status",
            "-t",
            "100",
        ]
    if method == "GBBS":
        return [
            "GBBS-VGC/baselines/gbbs/benchmarks/StronglyConnectedComponents/"
            "RandomGreedyBGSS16/StronglyConnectedComponents",
            "-b",
            "-beta",
            "1.5",
            "-rounds",
            "100",
            graph_path,
        ]
    if method == "Multistep":
        return ["MultiStep/multistep/scc", graph_path, "-t", "100"]
    raise ValueError(f"unknown method: {method}")


def load_existing_csv(path: Path) -> dict[str, dict[str, str]]:
    rows = blank_rows()
    if not path.exists():
        return rows

    with path.open("r", newline="", encoding="utf-8") as f:
        reader = csv.DictReader(f)
        for row in reader:
            dataset = row.get("Graph", "") or row.get("Dataset", "")
            if dataset in rows:
                for method in METHODS:
                    rows[dataset][method] = format_average_cost(row.get(method, ""))
    return rows


def blank_rows() -> dict[str, dict[str, str]]:
    return {name: {method: "" for method in METHODS} for name, _ in DATASETS}


def write_csv(path: Path, rows: dict[str, dict[str, str]]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="", encoding="utf-8") as f:
        writer = csv.writer(f)
        writer.writerow(("Graph", *METHODS, "Our rank"))
        for dataset, _ in DATASETS:
            writer.writerow(
                (
                    dataset,
                    *(rows[dataset][method] for method in METHODS),
                    calculate_our_rank(rows[dataset]),
                )
            )


def calculate_our_rank(row: dict[str, str]) -> str:
    if any(not row[method] for method in METHODS):
        return ""

    costs = {method: Decimal(row[method]) for method in METHODS}
    fast_scc_cost = costs["Fast-SCC"]
    return str(1 + sum(cost < fast_scc_cost for cost in costs.values()))


def parse_result_value(output: str, method: str) -> str:
    matches = AVERAGE_COST_RE.findall(output)
    if matches:
        return format_average_cost(matches[-1])

    if method == "GBBS":
        gbbs_times = [
            Decimal(value)
            for value in (*TIME_PER_ITER_RE.findall(output), *RUNNING_TIME_RE.findall(output))
        ]
        if gbbs_times:
            average = sum(gbbs_times) / Decimal(len(gbbs_times))
            return format_average_cost(str(average))

    raise ValueError(
        'could not find "average cost" or GBBS "time per iter"/"Running Time" in output'
    )


def format_average_cost(value: str) -> str:
    value = value.strip()
    if not value:
        return ""

    try:
        return str(Decimal(value).quantize(Decimal("0.001"), rounding=ROUND_HALF_UP))
    except InvalidOperation as exc:
        raise ValueError(f"invalid average cost value: {value}") from exc


def run_one(
    project_root: Path,
    method: str,
    dataset: str,
    graph_file: str,
    threads: int,
    use_numactl: bool,
    log_dir: Path,
    dry_run: bool,
) -> str:
    graph_path = f"./data/{graph_file}"
    cmd = command_for(method, graph_path)
    if use_numactl:
        cmd = ["numactl", "-i", "all", *cmd]

    printable = " ".join(shlex.quote(part) for part in cmd)
    print(f"[RUN] {dataset} / {method}: {printable}", flush=True)
    if dry_run:
        return ""

    env = os.environ.copy()
    env["PARLAY_NUM_THREADS"] = str(threads)
    env["OMP_NUM_THREADS"] = str(threads)

    completed = subprocess.run(
        cmd,
        cwd=project_root,
        env=env,
        text=True,
        stdout=subprocess.PIPE,
        stderr=subprocess.STDOUT,
    )

    log_dir.mkdir(parents=True, exist_ok=True)
    log_path = log_dir / f"{dataset}-{method}.log"
    log_path.write_text(completed.stdout, encoding="utf-8", errors="replace")

    if completed.returncode != 0:
        raise RuntimeError(
            f"command failed with exit code {completed.returncode}; see {log_path}"
        )

    value = parse_result_value(completed.stdout, method)
    print(f"[OK] {dataset} / {method}: average cost = {value}", flush=True)
    return value


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Run Table 4 experiments and write average cost to CSV."
    )
    parser.add_argument(
        "--project-root",
        type=Path,
        default=Path("/home/thinker/fhz/Fast-SCC-revision"),
        help="Fast-SCC-revision project root.",
    )
    parser.add_argument(
        "--output",
        type=Path,
        default=None,
        help="CSV output path. Default: <project-root>/results/table4.csv",
    )
    parser.add_argument(
        "--threads",
        type=int,
        default=24,
        help="Thread count for PARLAY_NUM_THREADS and OMP_NUM_THREADS.",
    )
    parser.add_argument(
        "--methods",
        nargs="+",
        choices=METHODS,
        default=list(METHODS),
        help="Methods to run. A subset preserves the other existing CSV columns.",
    )
    parser.add_argument(
        "--resume",
        action="store_true",
        help="Skip CSV cells that already have values.",
    )
    parser.add_argument(
        "--keep-going",
        action="store_true",
        help="Continue with later cells when one command fails.",
    )
    parser.add_argument(
        "--no-numactl",
        action="store_true",
        help="Run commands without the 'numactl -i all' prefix.",
    )
    parser.add_argument(
        "--dry-run",
        action="store_true",
        help="Print commands and create the CSV layout without running experiments.",
    )
    return parser.parse_args()


def main() -> int:
    args = parse_args()
    project_root = args.project_root.resolve()
    output_path = args.output or (project_root / "results" / "table4.csv")
    output_path = output_path.resolve()
    log_dir = output_path.parent / "table4-logs"

    selected_methods = tuple(args.methods)
    partial_run = set(selected_methods) != set(METHODS)
    rows = load_existing_csv(output_path) if args.resume or partial_run else blank_rows()
    write_csv(output_path, rows)

    print(f"[INFO] project root: {project_root}", flush=True)
    print(f"[INFO] output csv: {output_path}", flush=True)
    print(f"[INFO] threads: {args.threads}", flush=True)
    print(
        f"[INFO] PARLAY_NUM_THREADS={args.threads}, OMP_NUM_THREADS={args.threads}",
        flush=True,
    )

    failed = False
    for dataset, graph_file in DATASETS:
        for method in selected_methods:
            if args.resume and rows[dataset][method]:
                print(f"[SKIP] {dataset} / {method}: {rows[dataset][method]}", flush=True)
                continue

            try:
                rows[dataset][method] = run_one(
                    project_root=project_root,
                    method=method,
                    dataset=dataset,
                    graph_file=graph_file,
                    threads=args.threads,
                    use_numactl=not args.no_numactl,
                    log_dir=log_dir,
                    dry_run=args.dry_run,
                )
                write_csv(output_path, rows)
            except Exception as exc:
                failed = True
                print(f"[ERROR] {dataset} / {method}: {exc}", file=sys.stderr, flush=True)
                write_csv(output_path, rows)
                if not args.keep_going:
                    return 1

    return 1 if failed else 0


if __name__ == "__main__":
    raise SystemExit(main())
