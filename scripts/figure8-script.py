#!/usr/bin/env python3
"""Run Figure 8 thread-scaling experiments and write figure8.csv.

Default usage on the server:
  cd /home/thinker/fhz/Fast-SCC-revision
  python3 scripts/figure8-script.py

The script tests threads 1/2/4/8/12/24 on LJ, HH5, GL10, and REC.
For each graph/thread pair it runs the four methods, records running time,
and computes speedup(method) = GBBS time / method time.
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
METHODS = ("Fast-SCC", "GBBS-VGC", "GBBS", "Multistep")
SPEEDUP_COLUMNS = ("Ours-speedup", "GBBS-VGC-speedup", "GBBS-speedup", "Multistep-speedup")

DATASETS = (
    ("LJ", "soc-LiveJournal1.bin"),
    ("HH5", "Household.lines_5.bin"),
    ("GL10", "GeoLifeNoScale_10.bin"),
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
            "Fast-SCC/src/scc_paper_swap",
            graph_path,
            "-local_reach",
            "-local_scc",
            "-status",
            "-t",
            "10",
        ]
    if method == "GBBS-VGC":
        return [
            "GBBS-VGC/src/scc",
            graph_path,
            "-local_reach",
            "-local_scc",
            "-status",
            "-t",
            "10",
        ]
    if method == "GBBS":
        return [
            "GBBS-VGC/baselines/gbbs/benchmarks/StronglyConnectedComponents/"
            "RandomGreedyBGSS16/StronglyConnectedComponents",
            "-b",
            "-beta",
            "1.5",
            "-rounds",
            "10",
            graph_path,
        ]
    if method == "Multistep":
        return ["MultiStep/multistep/scc", graph_path, "-t", "10"]
    raise ValueError(f"unknown method: {method}")


def blank_rows() -> dict[tuple[str, int], dict[str, str]]:
    rows: dict[tuple[str, int], dict[str, str]] = {}
    for graph, _ in DATASETS:
        for threads in THREADS:
            rows[(graph, threads)] = {
                **{method: "" for method in METHODS},
                **{column: "" for column in SPEEDUP_COLUMNS},
            }
    return rows


def format_decimal(value: str | Decimal, places: str = "0.001") -> str:
    text = str(value).strip()
    if not text:
        return ""
    try:
        return str(Decimal(text).quantize(Decimal(places), rounding=ROUND_HALF_UP))
    except InvalidOperation as exc:
        raise ValueError(f"invalid numeric value: {value}") from exc


def parse_result_value(output: str, method: str) -> str:
    matches = AVERAGE_COST_RE.findall(output)
    if matches:
        return format_decimal(matches[-1])

    if method == "GBBS":
        gbbs_times = [
            Decimal(value)
            for value in (*TIME_PER_ITER_RE.findall(output), *RUNNING_TIME_RE.findall(output))
        ]
        if gbbs_times:
            average = sum(gbbs_times) / Decimal(len(gbbs_times))
            return format_decimal(average)

    raise ValueError(
        'could not find "average cost" or GBBS "time per iter"/"Running Time" in output'
    )


def compute_speedups(row: dict[str, str]) -> None:
    if any(not row[method] for method in METHODS):
        for column in SPEEDUP_COLUMNS:
            row[column] = ""
        return

    gbbs = Decimal(row["GBBS"])
    row["Ours-speedup"] = format_decimal(gbbs / Decimal(row["Fast-SCC"]))
    row["GBBS-VGC-speedup"] = format_decimal(gbbs / Decimal(row["GBBS-VGC"]))
    row["GBBS-speedup"] = format_decimal(Decimal("1"))
    row["Multistep-speedup"] = format_decimal(gbbs / Decimal(row["Multistep"]))


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
            if key not in rows:
                continue
            for method in METHODS:
                rows[key][method] = format_decimal(raw.get(method, ""))
            compute_speedups(rows[key])
    return rows


def write_csv(path: Path, rows: dict[tuple[str, int], dict[str, str]]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="", encoding="utf-8") as f:
        writer = csv.writer(f)
        writer.writerow(("Graph", "Threads", *METHODS, *SPEEDUP_COLUMNS))
        for graph, _ in DATASETS:
            for threads in THREADS:
                row = rows[(graph, threads)]
                compute_speedups(row)
                writer.writerow(
                    (
                        graph,
                        threads,
                        *(row[method] for method in METHODS),
                        *(row[column] for column in SPEEDUP_COLUMNS),
                    )
                )


def run_one(
    project_root: Path,
    method: str,
    graph: str,
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
    print(f"[RUN] {graph} / {threads} threads / {method}: {printable}", flush=True)
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
    log_path = log_dir / f"{graph}-t{threads}-{method}.log"
    log_path.write_text(completed.stdout, encoding="utf-8", errors="replace")

    if completed.returncode != 0:
        raise RuntimeError(
            f"command failed with exit code {completed.returncode}; see {log_path}"
        )

    value = parse_result_value(completed.stdout, method)
    print(f"[OK] {graph} / {threads} threads / {method}: time = {value}", flush=True)
    return value


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Run Figure 8 experiments and write timing/speedup CSV."
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
        help="CSV output path. Default: <project-root>/results/figure8.csv",
    )
    parser.add_argument(
        "--resume",
        action="store_true",
        help="Skip method cells that already have values.",
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
    output_path = args.output or (project_root / "results" / "figure8.csv")
    output_path = output_path.resolve()
    log_dir = output_path.parent / "figure8-logs"

    rows = load_existing_csv(output_path) if args.resume else blank_rows()
    write_csv(output_path, rows)

    print(f"[INFO] project root: {project_root}", flush=True)
    print(f"[INFO] output csv: {output_path}", flush=True)
    print(f"[INFO] threads: {', '.join(map(str, THREADS))}", flush=True)

    failed = False
    for graph, graph_file in DATASETS:
        for threads in THREADS:
            key = (graph, threads)
            print(
                f"[INFO] PARLAY_NUM_THREADS={threads}, OMP_NUM_THREADS={threads}",
                flush=True,
            )
            for method in METHODS:
                if args.resume and rows[key][method]:
                    print(
                        f"[SKIP] {graph} / {threads} threads / {method}: {rows[key][method]}",
                        flush=True,
                    )
                    continue
                try:
                    rows[key][method] = run_one(
                        project_root=project_root,
                        method=method,
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
                    print(
                        f"[ERROR] {graph} / {threads} threads / {method}: {exc}",
                        file=sys.stderr,
                        flush=True,
                    )
                    write_csv(output_path, rows)
                    if not args.keep_going:
                        return 1

    return 1 if failed else 0


if __name__ == "__main__":
    raise SystemExit(main())
