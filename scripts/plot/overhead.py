#!/usr/bin/env python3
"""Plot GBBS-VGC runtime and the profiled Fast-SCC module overhead."""

from __future__ import annotations

import argparse
import csv
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt


PLOT_DIR = Path(__file__).resolve().parent
ROOT = PLOT_DIR.parents[1]
DEFAULT_INPUT = ROOT / "results" / "overhead.csv"
DEFAULT_OUTPUT = ROOT / "figures" / "overhead"
GRAPH_ORDER = ("CH5", "GL2", "GL5", "GL10", "GL15", "GL20")


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Draw the Fast-SCC overhead stacked-bar figure."
    )
    parser.add_argument("--input", type=Path, default=DEFAULT_INPUT)
    parser.add_argument("--output", type=Path, default=DEFAULT_OUTPUT)
    return parser.parse_args()


def read_rows(path: Path) -> list[dict[str, str]]:
    with path.open(newline="", encoding="utf-8") as handle:
        by_graph = {row["graph"]: row for row in csv.DictReader(handle)}
    missing = [graph for graph in GRAPH_ORDER if graph not in by_graph]
    if missing:
        raise ValueError(f"CSV is missing graphs: {', '.join(missing)}")
    rows = [by_graph[graph] for graph in GRAPH_ORDER]
    failed = [row["graph"] for row in rows if row["correctness"] != "PASS"]
    if failed:
        raise ValueError(f"SCC correctness failed for: {', '.join(failed)}")
    return rows


def values(rows: list[dict[str, str]], field: str) -> list[float]:
    return [float(row[field]) for row in rows]


def main() -> None:
    args = parse_args()
    rows = read_rows(args.input)
    labels = [row["graph"] for row in rows]
    baseline = values(rows, "gbbs_vgc_seconds")
    core = values(rows, "fast_core_seconds")
    fake_link = values(rows, "fake_link_seconds")
    pivot = values(rows, "pivot_selection_seconds")

    x = list(range(len(labels)))
    width = 0.34
    left = [value - width / 2 for value in x]
    right = [value + width / 2 for value in x]

    plt.rcParams.update(
        {
            "font.family": "serif",
            "font.size": 16,
            "axes.labelsize": 18,
            "axes.linewidth": 0.8,
            "hatch.linewidth": 0.7,
            "legend.fontsize": 12,
            "pdf.fonttype": 42,
            "ps.fonttype": 42,
        }
    )
    fig, ax = plt.subplots(figsize=(7.4, 4.15), constrained_layout=True)

    edge = "#303030"
    ax.bar(
        left,
        baseline,
        width,
        label="GBBS-VGC",
        color="#9fb88d",
        edgecolor=edge,
        linewidth=0.8,
        hatch="//",
        zorder=3,
    )
    ax.bar(
        right,
        core,
        width,
        label="Fast-SCC",
        color="#D99058",
        edgecolor=edge,
        linewidth=0.8,
        zorder=3,
    )
    ax.bar(
        right,
        fake_link,
        width,
        bottom=core,
        label="Fake-link Identification",
        color="#ECC8AC",
        edgecolor=edge,
        linewidth=0.8,
        hatch="xx",
        zorder=3,
    )
    second_bottom = [a + b for a, b in zip(core, fake_link)]
    ax.bar(
        right,
        pivot,
        width,
        bottom=second_bottom,
        label="Bounded Pivot Selection",
        color="#F8E9DE",
        edgecolor=edge,
        linewidth=0.8,
        hatch="...",
        zorder=3,
    )

    ax.set_ylabel("Running time (s)")
    ax.set_xlabel("Graph")
    ax.set_xticks(x, labels)
    ax.set_xlim(-0.65, len(labels) - 0.35)
    ax.set_ylim(bottom=0)
    ax.grid(axis="y", linestyle="--", linewidth=0.65, color="#B7B7B7", alpha=0.8)
    ax.set_axisbelow(True)
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    ax.legend(
        loc="upper left",
        frameon=True,
        ncol=2,
        columnspacing=1.2,
        handlelength=1.8,
    )

    args.output.parent.mkdir(parents=True, exist_ok=True)
    pdf_path = args.output.with_suffix(".pdf")
    # png_path = args.output.with_suffix(".png")
    fig.savefig(pdf_path, bbox_inches="tight")
    # fig.savefig(png_path, dpi=300, bbox_inches="tight")
    plt.close(fig)
    print(f"wrote {pdf_path}")
    # print(f"wrote {png_path}")


if __name__ == "__main__":
    main()
