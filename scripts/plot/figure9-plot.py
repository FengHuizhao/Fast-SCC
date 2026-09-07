#!/usr/bin/env python3
"""Plot Figure 9 self-speedup curves from figure9.csv.

Default usage on the server:
  cd /home/thinker/fhz/Fast-SCC-revision
  python3 scripts/plot/figure9-plot.py

Input:
  results/figure9.csv

Output:
  figures/figure9.pdf
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
        "figure9-plot.py needs matplotlib. Install it on the server, for example: "
        "python3 -m pip install matplotlib"
    ) from exc


PROJECT_ROOT = Path("/home/thinker/fhz/Fast-SCC-revision")
THREADS = (1, 2, 4, 8, 12, 24)
GRAPHS = (
    ("LJ", "#f26d6d", "o"),
    ("HH5", "#9fb88d", "v"),
    ("CH5", "#8e65d3", "s"),
    ("GL5", "#f3a0c3", "*"),
    ("GL10", "#9cc6ef", "d"),
    ("REC", "#f2bd93", "X"),
)


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Plot Figure 9 self-speedup PDF.")
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
        help="Input CSV. Default: <project-root>/results/figure9.csv.",
    )
    parser.add_argument(
        "--output",
        type=Path,
        default=None,
        help="Output PDF. Default: <project-root>/figures/figure9.pdf.",
    )
    return parser.parse_args()


def parse_float(value: str, graph: str, threads: int) -> float:
    text = str(value).strip()
    if not text:
        raise ValueError(f"missing self-speedup for {graph} / {threads}")
    try:
        return float(text)
    except ValueError as exc:
        raise ValueError(f"invalid self-speedup for {graph} / {threads}: {text}") from exc


def load_self_speedups(csv_path: Path) -> dict[str, list[float]]:
    if not csv_path.exists():
        raise FileNotFoundError(f"input CSV not found: {csv_path}")

    graph_names = {graph for graph, _, _ in GRAPHS}
    by_key: dict[tuple[str, int], dict[str, str]] = {}

    with csv_path.open("r", newline="", encoding="utf-8-sig") as f:
        reader = csv.DictReader(f)
        for row in reader:
            graph = (row.get("Graph") or row.get("Dataset") or "").strip()
            if graph not in graph_names:
                continue
            threads_text = (row.get("Threads") or row.get("Thread") or "").strip()
            if not threads_text:
                continue
            by_key[(graph, int(threads_text))] = row

    data: dict[str, list[float]] = {}
    for graph, _, _ in GRAPHS:
        values = []
        for threads in THREADS:
            row = by_key.get((graph, threads))
            if row is None:
                raise ValueError(f"missing row for {graph} / {threads} threads")
            values.append(parse_float(row.get("Self-speedup", ""), graph, threads))
        data[graph] = values

    return data


def plot(data: dict[str, list[float]], output_path: Path) -> None:
    fig, ax = plt.subplots(figsize=(6.25, 3.55))

    for graph, color, marker in GRAPHS:
        ax.plot(
            THREADS,
            data[graph],
            label=graph,
            color=color,
            marker=marker,
            linewidth=2.0,
            markersize=7.0 if marker not in ("*", "X") else 8.5,
            alpha=0.85,
        )

    ax.set_xlabel("Number of Threads", fontsize=14)
    ax.set_ylabel("Self-Speedup", fontsize=16)
    ax.set_xticks(THREADS)
    ax.set_xticklabels([str(thread) for thread in THREADS], rotation=0, fontsize=12)

    ymax = max(max(values) for values in data.values())
    y_top = max(1, math.ceil(ymax * 1.08))
    ax.set_ylim(bottom=0.5, top=y_top)
    ax.set_yticks([1, *range(2, y_top + 1, 2)])
    ax.tick_params(axis="y", labelsize=12)

    for spine in ax.spines.values():
        spine.set_linewidth(1.0)

    ax.legend(
        loc="upper center",
        bbox_to_anchor=(0.5, 1.12),
        ncol=6,
        frameon=False,
        fontsize=11,
        columnspacing=1.2,
        handlelength=1.8,
    )

    output_path.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(output_path, format="pdf", bbox_inches="tight", pad_inches=0.04)
    plt.close(fig)


def main() -> int:
    args = parse_args()
    project_root = args.project_root.resolve()
    input_path = (args.input or (project_root / "results" / "figure9.csv")).resolve()
    output_path = (args.output or (project_root / "figures" / "figure9.pdf")).resolve()

    print(f"[INFO] input csv: {input_path}")
    print(f"[INFO] output pdf: {output_path}")
    data = load_self_speedups(input_path)
    plot(data, output_path)
    print("[OK] saved Figure 9 PDF.")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
