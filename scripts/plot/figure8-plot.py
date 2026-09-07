#!/usr/bin/env python3
"""Plot Figure 8 speedup curves from figure8.csv.

Default usage on the server:
  cd /home/thinker/fhz/Fast-SCC-revision
  python3 scripts/plot/figure8-plot.py

Input:
  results/figure8.csv

Outputs:
  figures/figure8-LJ.pdf
  figures/figure8-HH5.pdf
  figures/figure8-GL10.pdf
  figures/figure8-REC.pdf
"""

from __future__ import annotations

import argparse
import csv
import math
from pathlib import Path


try:
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
except ModuleNotFoundError as exc:
    raise ModuleNotFoundError(
        "figure8-plot.py needs matplotlib. Install it on the server, for example: "
        "python3 -m pip install matplotlib"
    ) from exc


PROJECT_ROOT = Path("/home/thinker/fhz/Fast-SCC-revision")
THREADS = (1, 2, 4, 8, 12, 24)
GRAPHS = ("LJ", "HH5", "GL10", "REC")
METHODS = (
    ("OURS", "Ours-speedup", "#f26d6d", "o"),
    ("GBBS-VGC", "GBBS-VGC-speedup", "#9fb88d", "v"),
    ("GBBS", "GBBS-speedup", "#8e65d3", "s"),
    ("Multistep", "Multistep-speedup", "#4f91d9", "*"),
)


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Plot Figure 8 speedup PDFs.")
    parser.add_argument(
        "--project-root",
        type=Path,
        default=PROJECT_ROOT,
        help="Fast-SCC-revision project root.",
    )
    parser.add_argument(
        "--input",
        type=Path,
        default=None,
        help="Input CSV. Default: <project-root>/results/figure8.csv.",
    )
    parser.add_argument(
        "--output-dir",
        type=Path,
        default=None,
        help="Output directory. Default: <project-root>/figures.",
    )
    return parser.parse_args()


def parse_float(value: str, graph: str, threads: int, column: str) -> float:
    text = str(value).strip()
    if not text:
        raise ValueError(f"missing value for {graph} / {threads} / {column}")
    try:
        return float(text)
    except ValueError as exc:
        raise ValueError(f"invalid value for {graph} / {threads} / {column}: {text}") from exc


def load_speedups(csv_path: Path) -> dict[str, dict[str, list[float]]]:
    if not csv_path.exists():
        raise FileNotFoundError(f"input CSV not found: {csv_path}")

    data = {graph: {label: [] for label, _, _, _ in METHODS} for graph in GRAPHS}
    by_key: dict[tuple[str, int], dict[str, str]] = {}

    with csv_path.open("r", newline="", encoding="utf-8-sig") as f:
        reader = csv.DictReader(f)
        for row in reader:
            graph = (row.get("Graph") or row.get("Dataset") or "").strip()
            if graph not in GRAPHS:
                continue
            threads_text = (row.get("Threads") or row.get("Thread") or "").strip()
            if not threads_text:
                continue
            by_key[(graph, int(threads_text))] = row

    for graph in GRAPHS:
        for threads in THREADS:
            row = by_key.get((graph, threads))
            if row is None:
                raise ValueError(f"missing row for {graph} / {threads} threads")
            for label, speedup_column, _, _ in METHODS:
                data[graph][label].append(parse_float(row.get(speedup_column, ""), graph, threads, speedup_column))

    return data


def plot_one(graph: str, series: dict[str, list[float]], output_path: Path) -> None:
    fig, ax = plt.subplots(figsize=(4.25, 3.05))
    fig.subplots_adjust(top=0.80)

    for label, _, color, marker in METHODS:
        ax.plot(
            THREADS,
            series[label],
            label=label,
            color=color,
            marker=marker,
            linewidth=2.0,
            markersize=7.0 if marker != "*" else 9.0,
            alpha=0.9,
        )

    ax.set_xlabel("Number of Threads", fontsize=13)
    ax.set_ylabel("Speedup", fontsize=13)
    ax.set_xticks(THREADS)
    ax.set_xticklabels([str(thread) for thread in THREADS], rotation=0, fontsize=11)
    ax.tick_params(axis="y", labelsize=11)
    ax.grid(False)

    ymax = max(max(values) for values in series.values())
    y_top = max(1, math.ceil(ymax * 1.12))
    ax.set_ylim(bottom=0, top=y_top)
    ax.set_yticks(range(0, y_top + 1, 1))

    for spine in ax.spines.values():
        spine.set_linewidth(1.1)

    ax.legend(
        loc="upper center",
        bbox_to_anchor=(0.5, 1.20),
        ncol=4,
        frameon=True,
        fontsize=10,
        columnspacing=1.1,
        handlelength=1.8,
    )

    output_path.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(output_path, format="pdf", bbox_inches="tight", pad_inches=0.03)
    plt.close(fig)


def main() -> int:
    args = parse_args()
    project_root = args.project_root.resolve()
    input_path = (args.input or (project_root / "results" / "figure8.csv")).resolve()
    output_dir = (args.output_dir or (project_root / "figures")).resolve()

    print(f"[INFO] input csv: {input_path}")
    print(f"[INFO] output dir: {output_dir}")

    data = load_speedups(input_path)
    for graph in GRAPHS:
        output_path = output_dir / f"figure8-{graph}.pdf"
        plot_one(graph, data[graph], output_path)
        print(f"[OK] saved {output_path}")

    return 0


if __name__ == "__main__":
    raise SystemExit(main())
