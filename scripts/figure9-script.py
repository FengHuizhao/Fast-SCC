#!/usr/bin/env python3
"""Run Figure 9 self-speedup experiments and write figure9.csv.

Default usage on the server:
  cd /home/thinker/fhz/Fast-SCC-revision
  python3 scripts/figure9-script.py

The script tests Fast-SCC on threads 1/2/4/8/12/24 for the datasets shown in
Figure 9. Self-speedup is computed as:
  self-speedup(x threads) = time(1 thread) / time(x threads)
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


PROJECT_ROOT = Path("/home/thinker/fhz/Fast-SCC-revision")
THREADS = (1, 2, 4, 8, 12, 24)
DATASETS = (
    ("LJ", "soc-LiveJournal1.bin"),
    ("HH5", "Household.lines_5.bin"),
    ("CH5", "CHEM_5.bin"),
    ("GL5", "GeoLifeNoScale_5.bin"),
    ("GL10", "GeoLifeNoScale_10.bin"),
    ("REC", "grid_1000_10000.bin"),
)

AVERAGE_COST_RE = re.compile(
    r"\baverage\s+cost\b\s*[:=]?\s*([-+]?(?:\d+(?:\.\d*)?|\.\d+)(?:[eE][-+]?\d+)?)",
    re.IGNORECASE,
)


def command_for(graph_path: str) -> list[str]:
    return [
        "Fast-SCC/src/scc_paper_swap",
        graph_path,
        "-local_reach",
        "-local_scc",
        "-status",
        "-t",
        "10",
    ]


def blank_rows() -> dict[tuple[str, int], dict[str, str]]:
    return {
        (graph, threads): {"Fast-SCC": "", "Self-speedup": ""}
        for graph, _ in DATASETS
        for threads in THREADS
    }


def format_decimal(value: str | Decimal, places: str = "0.001") -> str:
    text = str(value).strip()
    if not text:
        return ""
    try:
        return str(Decimal(text).quantize(Decimal(places), rounding=ROUND_HALF_UP))
    except InvalidOperation as exc:
        raise ValueError(f"invalid numeric value: {value}") from exc


def parse_average_cost(output: str) -> str:
    matches = AVERAGE_COST_RE.findall(output)
    if not matches:
        raise ValueError('could not find "average cost" in command output')
    return format_decimal(matches[-1])


def compute_self_speedups(rows: dict[tuple[str, int], dict[str, str]]) -> None:
    for graph, _ in DATASETS:
        baseline = rows[(graph, 1)]["Fast-SCC"]
        if not baseline:
            for threads in THREADS:
                rows[(graph, threads)]["Self-speedup"] = ""
            continue

        baseline_value = Decimal(baseline)
        for threads in THREADS:
            value = rows[(graph, threads)]["Fast-SCC"]
            rows[(graph, threads)]["Self-speedup"] = (
                format_decimal(baseline_value / Decimal(value)) if value else ""
            )


def load_existing_csv(path: Path) -> dict[tuple[str, int], dict[str, str]]:
    rows = blank_rows()
    if not path.exists():
        return rows

    with path.open("r", newline="", encoding="utf-8-sig") as f:
        reader = csv.DictReader(f)
        for raw in reader:
            graph = (raw.get("Graph") or raw.get("Dataset") or "").strip()
            threads_text = (raw.get("Threads") or raw.get("Thread") or "").strip()
            if not graph or not threads_text:
                continue
            try:
                threads = int(threads_text)
            except ValueError:
                continue
            key = (graph, threads)
            if key in rows:
                rows[key]["Fast-SCC"] = format_decimal(raw.get("Fast-SCC", ""))

    compute_self_speedups(rows)
    return rows


def write_csv(path: Path, rows: dict[tuple[str, int], dict[str, str]]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    compute_self_speedups(rows)
    with path.open("w", newline="", encoding="utf-8") as f:
        writer = csv.writer(f)
        writer.writerow(("Graph", "Threads", "Fast-SCC", "Self-speedup"))
        for graph, _ in DATASETS:
            for threads in THREADS:
                row = rows[(graph, threads)]
                writer.writerow((graph, threads, row["Fast-SCC"], row["Self-speedup"]))


def run_one(
    project_root: Path,
    graph: str,
    graph_file: str,
    threads: int,
    use_numactl: bool,
    log_dir: Path,
    dry_run: bool,
) -> str:
    graph_path = f"./data/{graph_file}"
    cmd = command_for(graph_path)
    if use_numactl:
        cmd = ["numactl", "-i", "all", *cmd]

    printable = " ".join(shlex.quote(part) for part in cmd)
    print(f"[RUN] {graph} / {threads} threads: {printable}", flush=True)
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
    log_path = log_dir / f"{graph}-t{threads}-Fast-SCC.log"
    log_path.write_text(completed.stdout, encoding="utf-8", errors="replace")

    if completed.returncode != 0:
        raise RuntimeError(
            f"command failed with exit code {completed.returncode}; see {log_path}"
        )

    value = parse_average_cost(completed.stdout)
    print(f"[OK] {graph} / {threads} threads: time = {value}", flush=True)
    return value


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Run Figure 9 Fast-SCC self-speedup experiments."
    )
    parser.add_argument(
        "--project-root",
        type=Path,
        default=PROJECT_ROOT,
        help="Fast-SCC-revision project root.",
    )
    parser.add_argument(
        "--output",
        type=Path,
        default=None,
        help="CSV output path. Default: <project-root>/results/figure9.csv",
    )
    parser.add_argument(
        "--resume",
        action="store_true",
        help="Skip cells that already have values.",
    )
    parser.add_argument(
        "--keep-going",
        action="store_true",
        help="Continue with later runs when one command fails.",
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
    output_path = args.output or (project_root / "results" / "figure9.csv")
    output_path = output_path.resolve()
    log_dir = output_path.parent / "figure9-logs"

    rows = load_existing_csv(output_path) if args.resume else blank_rows()
    write_csv(output_path, rows)

    print(f"[INFO] project root: {project_root}", flush=True)
    print(f"[INFO] output csv: {output_path}", flush=True)
    print(f"[INFO] threads: {', '.join(map(str, THREADS))}", flush=True)

    failed = False
    for graph, graph_file in DATASETS:
        for threads in THREADS:
            if args.resume and rows[(graph, threads)]["Fast-SCC"]:
                print(
                    f"[SKIP] {graph} / {threads} threads: {rows[(graph, threads)]['Fast-SCC']}",
                    flush=True,
                )
                continue

            print(
                f"[INFO] PARLAY_NUM_THREADS={threads}, OMP_NUM_THREADS={threads}",
                flush=True,
            )
            try:
                rows[(graph, threads)]["Fast-SCC"] = run_one(
                    project_root=project_root,
                    graph=graph,
                    graph_file=graph_file,
                    threads=threads,
                    use_numactl=not args.no_numactl,
                    log_dir=log_dir,
                    dry_run=args.dry_run,
                )
                write_csv(output_path, rows)
            except Exception as exc:
                failed = True
                print(f"[ERROR] {graph} / {threads} threads: {exc}", file=sys.stderr, flush=True)
                write_csv(output_path, rows)
                if not args.keep_going:
                    return 1

    return 1 if failed else 0


if __name__ == "__main__":
    raise SystemExit(main())
