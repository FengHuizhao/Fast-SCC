#!/usr/bin/env python3
"""Run the 2x2 Fast-SCC ablation study and write averaged results to CSV.

Default usage on the server:
  cd /home/thinker/fhz/Fast-SCC-revision
  python3 scripts/ablation.py

Datasets are processed in paper order. For each dataset, the four variants
are run in ablation order. Each executable invocation performs ten measured
runs; this script averages the ten ABLATION_RESULT records and writes the
partial CSV after every completed configuration.
"""

from __future__ import annotations

import argparse
import csv
import os
import shlex
import subprocess
import sys
from decimal import Decimal, InvalidOperation, ROUND_HALF_UP
from pathlib import Path


VERSIONS = ("GBBS-VGC", "Fast-SCC-FL", "Fast-SCC-PS", "Fast-SCC")

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

CSV_FIELDS = (
    "graph",
    "version",
    "runtime_seconds",
    "edge_checks_total",
    "fake_link_checks",
    "global_search_rounds",
    "executed_batches",
    "pivots_selected",
    "scc_count",
    "correctness",
)

AVERAGED_FIELDS = (
    "runtime_seconds",
    "edge_checks_total",
    "fake_link_checks",
    "global_search_rounds",
    "executed_batches",
    "pivots_selected",
)


def blank_rows() -> dict[tuple[str, str], dict[str, str]]:
    return {
        (graph, version): {
            field: (graph if field == "graph" else version if field == "version" else "")
            for field in CSV_FIELDS
        }
        for graph, _ in DATASETS
        for version in VERSIONS
    }


def load_existing_csv(path: Path) -> dict[tuple[str, str], dict[str, str]]:
    rows = blank_rows()
    if not path.exists():
        return rows

    with path.open("r", newline="", encoding="utf-8") as stream:
        for row in csv.DictReader(stream):
            key = (row.get("graph", ""), row.get("version", ""))
            if key in rows:
                rows[key].update({field: row.get(field, "") for field in CSV_FIELDS})
    return rows


def write_csv(path: Path, rows: dict[tuple[str, str], dict[str, str]]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    temporary = path.with_suffix(path.suffix + ".tmp")
    with temporary.open("w", newline="", encoding="utf-8") as stream:
        writer = csv.DictWriter(stream, fieldnames=CSV_FIELDS)
        writer.writeheader()
        for graph, _ in DATASETS:
            for version in VERSIONS:
                writer.writerow(rows[(graph, version)])
    temporary.replace(path)


def has_current_log(log_dir: Path, graph: str, version: str) -> bool:
    log_path = log_dir / f"{graph}-{version}.log"
    if not log_path.is_file():
        return False
    with log_path.open("r", encoding="utf-8", errors="replace") as stream:
        config_line = next(
            (line.strip() for line in stream if line.startswith("ABLATION_CONFIG ")),
            "",
        )
    expected_compatibility = "ON" if version == "GBBS-VGC" else "OFF"
    return (
        "revision=2" in config_line
        and f"method={version}" in config_line
        and f"gbbs_vgc_compatibility={expected_compatibility}" in config_line
    )


def parse_result_line(line: str) -> dict[str, str] | None:
    stripped = line.strip()
    if not stripped.startswith("ABLATION_RESULT "):
        return None

    record: dict[str, str] = {}
    for token in stripped.split()[1:]:
        if "=" in token:
            key, value = token.split("=", 1)
            record[key] = value
    return record


def decimal_value(record: dict[str, str], field: str) -> Decimal:
    try:
        return Decimal(record[field])
    except KeyError as exc:
        raise ValueError(f"ABLATION_RESULT is missing {field}") from exc
    except InvalidOperation as exc:
        raise ValueError(f"invalid {field} value: {record.get(field)!r}") from exc


def format_average(value: Decimal, places: int) -> str:
    quantum = Decimal(1).scaleb(-places)
    rounded = value.quantize(quantum, rounding=ROUND_HALF_UP)
    if rounded == rounded.to_integral_value():
        return str(rounded.to_integral_value())
    return format(rounded, "f")


def summarize(
    graph: str,
    version: str,
    records: list[dict[str, str]],
    runs: int,
) -> dict[str, str]:
    if len(records) != runs:
        raise ValueError(
            f"expected {runs} ABLATION_RESULT records, found {len(records)}"
        )

    for index, record in enumerate(records, start=1):
        if record.get("method") != version:
            raise ValueError(
                f"run {index} reports method={record.get('method')!r}, expected {version}"
            )
        if record.get("correctness") != "PASS":
            raise ValueError(
                f"run {index} failed SCC correctness: {record.get('correctness')!r}"
            )

    scc_counts = {record.get("scc_count") for record in records}
    if None in scc_counts or len(scc_counts) != 1:
        values = sorted(str(value) for value in scc_counts)
        raise ValueError(f"inconsistent scc_count values: {values}")

    row = {field: "" for field in CSV_FIELDS}
    row["graph"] = graph
    row["version"] = version
    divisor = Decimal(runs)
    for field in AVERAGED_FIELDS:
        average = sum(decimal_value(record, field) for record in records) / divisor
        row[field] = format_average(average, 9 if field == "runtime_seconds" else 3)
    row["scc_count"] = next(iter(scc_counts)) or ""
    row["correctness"] = "PASS"
    return row


def run_configuration(
    project_root: Path,
    binary: Path,
    graph: str,
    graph_file: str,
    version: str,
    runs: int,
    threads: int | None,
    log_dir: Path,
) -> dict[str, str]:
    graph_path = project_root / "data" / graph_file
    command = [
        str(binary),
        version,
        str(graph_path),
        "-t",
        str(runs),
        "-local_reach",
        "-local_scc",
    ]
    print(
        f"[RUN] {graph} / {version}: "
        + " ".join(shlex.quote(part) for part in command),
        flush=True,
    )

    environment = os.environ.copy()
    if threads is not None:
        environment["PARLAY_NUM_THREADS"] = str(threads)
        environment["OMP_NUM_THREADS"] = str(threads)

    log_dir.mkdir(parents=True, exist_ok=True)
    log_path = log_dir / f"{graph}-{version}.log"
    records: list[dict[str, str]] = []
    with log_path.open("w", encoding="utf-8") as log_stream:
        process = subprocess.Popen(
            command,
            cwd=project_root,
            env=environment,
            text=True,
            stdout=subprocess.PIPE,
            stderr=subprocess.STDOUT,
            bufsize=1,
        )
        assert process.stdout is not None
        for line in process.stdout:
            sys.stdout.write(line)
            sys.stdout.flush()
            log_stream.write(line)
            parsed = parse_result_line(line)
            if parsed is not None:
                records.append(parsed)
        return_code = process.wait()

    if return_code != 0:
        raise RuntimeError(
            f"command exited with status {return_code}; see {log_path}"
        )

    row = summarize(graph, version, records, runs)
    print(
        f"[OK] {graph} / {version}: average runtime "
        f"{row['runtime_seconds']} s, correctness={row['correctness']}",
        flush=True,
    )
    return row


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Run the Fast-SCC 2x2 ablation study.")
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
        help="CSV path. Default: <project-root>/results/ablation.csv",
    )
    parser.add_argument(
        "--runs",
        type=int,
        default=10,
        help="Measured runs per graph/version invocation (default: 10).",
    )
    parser.add_argument(
        "--threads",
        type=int,
        default=24,
        help="PARLAY_NUM_THREADS and OMP_NUM_THREADS (default: 24).",
    )
    parser.add_argument(
        "--resume",
        action="store_true",
        help="Keep completed PASS rows in an existing CSV and skip them.",
    )
    parser.add_argument(
        "--keep-going",
        action="store_true",
        help="Continue with later configurations after a failure.",
    )
    return parser.parse_args()


def main() -> int:
    args = parse_args()
    if args.runs <= 0:
        raise ValueError("--runs must be greater than zero")
    if args.threads <= 0:
        raise ValueError("--threads must be greater than zero")

    project_root = args.project_root.resolve()
    binary = project_root / "Fast-SCC" / "src" / "scc_paper_swap_ablation"
    output = args.output or project_root / "results" / "ablation.csv"
    output = output.resolve()
    log_dir = project_root / "results" / "ablation-logs"

    if not binary.is_file():
        raise FileNotFoundError(f"ablation executable not found: {binary}")
    sources = (
        project_root / "Fast-SCC" / "src" / "scc_paper_swap_ablation.cpp",
        project_root / "Fast-SCC" / "src" / "scc_paper_swap_ablation.hpp",
    )
    if any(source.stat().st_mtime_ns > binary.stat().st_mtime_ns for source in sources):
        raise ValueError(
            "ablation executable is older than its corrected sources; rebuild it with:\n"
            "  cd /home/thinker/fhz/Fast-SCC-revision/Fast-SCC/src\n"
            "  g++ -std=c++17 -Wall -Wextra -Werror -pthread -O3 -mcx16 "
            "-march=native -I../parlaylib/include/ "
            "scc_paper_swap_ablation.cpp -o scc_paper_swap_ablation"
        )
    missing_graphs = [
        str(project_root / "data" / graph_file)
        for _, graph_file in DATASETS
        if not (project_root / "data" / graph_file).is_file()
    ]
    if missing_graphs:
        raise FileNotFoundError("missing graph files:\n  " + "\n  ".join(missing_graphs))

    rows = load_existing_csv(output) if args.resume else blank_rows()
    write_csv(output, rows)
    failures: list[str] = []

    for graph, graph_file in DATASETS:
        for version in VERSIONS:
            key = (graph, version)
            if (
                args.resume
                and rows[key].get("correctness") == "PASS"
                and has_current_log(log_dir, graph, version)
            ):
                print(f"[SKIP] {graph} / {version}: completed in {output}", flush=True)
                continue
            try:
                rows[key] = run_configuration(
                    project_root,
                    binary,
                    graph,
                    graph_file,
                    version,
                    args.runs,
                    args.threads,
                    log_dir,
                )
                write_csv(output, rows)
            except (OSError, RuntimeError, ValueError) as exc:
                message = f"{graph} / {version}: {exc}"
                failures.append(message)
                print(f"[ERROR] {message}", file=sys.stderr, flush=True)
                write_csv(output, rows)
                if not args.keep_going:
                    return 1

    if failures:
        print("Completed with failures:", file=sys.stderr)
        for failure in failures:
            print(f"  {failure}", file=sys.stderr)
        return 1

    print(f"[DONE] Results written to {output}", flush=True)
    return 0


if __name__ == "__main__":
    try:
        raise SystemExit(main())
    except (FileNotFoundError, ValueError) as exc:
        print(f"error: {exc}", file=sys.stderr)
        raise SystemExit(2)
