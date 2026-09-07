#!/usr/bin/env python3
"""Plot Figure 7 heatmap from Figure 4/Table 4 timing CSV.

Default usage on the server:
  cd /home/thinker/fhz/Fast-SCC-revision
  python3 scripts/plot/figure7-plot.py

The script reads results/figure4.csv if it exists, otherwise results/table4.csv.
It computes speedup as:
  speedup(method) = GBBS time / method time
and saves figures/figure7.pdf.
"""

from __future__ import annotations

import argparse
import csv
import math
from pathlib import Path

try:
    import matplotlib.pyplot as plt
    from matplotlib.colors import LinearSegmentedColormap, Normalize

    HAS_MATPLOTLIB = True
except ModuleNotFoundError:
    HAS_MATPLOTLIB = False


PROJECT_ROOT = Path("/home/thinker/fhz/Fast-SCC-revision")
DATASETS = ("LJ", "HH5", "CH5", "GL2", "GL5", "GL10", "GL15", "GL20", "SQR", "REC")
PLOT_ROWS = (
    ("Ours", ("Fast-SCC", "Ours")),
    ("GBBS-VGC", ("GBBS-VGC",)),
    ("Multistep", ("Multistep", "MultiStep")),
    ("GBBS", ("GBBS",)),
)
GROUPS = (
    ("Social", 0, 1),
    ("KNN", 1, 8),
    ("Lattice", 8, 10),
    ("Average", 10, 11),
)


def default_input_path(project_root: Path) -> Path:
    figure4 = project_root / "results" / "figure4.csv"
    table4 = project_root / "results" / "table4.csv"
    return figure4 if figure4.exists() else table4


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Plot Figure 7 speedup heatmap.")
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
        help="Input timing CSV. Default: results/figure4.csv, falling back to results/table4.csv.",
    )
    parser.add_argument(
        "--output",
        type=Path,
        default=None,
        help="Output PDF path. Default: <project-root>/figures/figure7.pdf.",
    )
    return parser.parse_args()


def first_present(row: dict[str, str], names: tuple[str, ...]) -> str:
    normalized = {key.strip().lower(): key for key in row}
    for name in names:
        key = normalized.get(name.strip().lower())
        if key is not None:
            return row[key]
    raise KeyError(f"missing column, expected one of: {', '.join(names)}")


def parse_float(value: str, dataset: str, method: str) -> float:
    value = str(value).strip()
    if not value:
        raise ValueError(f"empty value for {dataset} / {method}")
    try:
        parsed = float(value)
    except ValueError as exc:
        raise ValueError(f"invalid value for {dataset} / {method}: {value}") from exc
    if parsed <= 0:
        raise ValueError(f"value must be positive for {dataset} / {method}: {value}")
    return parsed


def load_times(csv_path: Path) -> dict[str, dict[str, float]]:
    if not csv_path.exists():
        raise FileNotFoundError(f"input CSV not found: {csv_path}")

    rows: dict[str, dict[str, float]] = {}
    with csv_path.open("r", newline="", encoding="utf-8-sig") as f:
        reader = csv.DictReader(f)
        for raw_row in reader:
            graph = (
                raw_row.get("Graph")
                or raw_row.get("Dataset")
                or raw_row.get("graph")
                or raw_row.get("dataset")
            )
            if not graph:
                continue

            graph = graph.strip()
            if graph not in DATASETS:
                continue

            rows[graph] = {}
            for label, aliases in PLOT_ROWS:
                value = first_present(raw_row, aliases)
                rows[graph][label] = parse_float(value, graph, label)

    missing = [dataset for dataset in DATASETS if dataset not in rows]
    if missing:
        raise ValueError(f"input CSV is missing dataset rows: {', '.join(missing)}")

    return rows


def geometric_mean(values: list[float]) -> float:
    if any(value <= 0 for value in values):
        raise ValueError("geometric mean requires positive values")
    return math.exp(sum(math.log(value) for value in values) / len(values))


def compute_speedups(times: dict[str, dict[str, float]]) -> tuple[list[str], list[str], list[list[float]]]:
    columns = [*DATASETS, ""]
    row_labels = [label for label, _ in PLOT_ROWS]
    matrix: list[list[float]] = []

    for label, _ in PLOT_ROWS:
        speedups = []
        for dataset in DATASETS:
            gbbs_time = times[dataset]["GBBS"]
            method_time = times[dataset][label]
            speedups.append(gbbs_time / method_time)
        speedups.append(geometric_mean(speedups))
        matrix.append(speedups)

    return row_labels, columns, matrix


def format_cell(value: float, row_label: str) -> str:
    if row_label == "GBBS" and abs(value - 1.0) < 0.0005:
        return "1"
    return f"{value:.2f}"


def plot_heatmap(row_labels: list[str], columns: list[str], matrix: list[list[float]], output_path: Path) -> None:
    if not HAS_MATPLOTLIB:
        plot_heatmap_reportlab(row_labels, columns, matrix, output_path)
        return

    nrows = len(row_labels)
    ncols = len(columns)

    cmap = LinearSegmentedColormap.from_list(
        "figure7_orange",
        ["#fff6f0", "#f5b27c", "#d96a13"],
    )
    norm = Normalize(vmin=0, vmax=6, clip=True)

    fig = plt.figure(figsize=(9.6, 3.15))
    ax = fig.add_axes([0.08, 0.14, 0.82, 0.78])
    cax = fig.add_axes([0.92, 0.14, 0.016, 0.78])

    ax.imshow(matrix, cmap=cmap, norm=norm, aspect="auto", extent=[0, ncols, nrows, 0])

    # Cell values.
    for i, row_label in enumerate(row_labels):
        for j in range(ncols):
            ax.text(
                j + 0.5,
                i + 0.5,
                format_cell(matrix[i][j], row_label),
                ha="center",
                va="center",
                fontsize=11,
                fontweight="bold",
                color="black",
            )

    # Column labels and group labels above the heatmap.
    for j, column in enumerate(columns):
        if column:
            ax.text(j + 0.5, -0.36, column, ha="center", va="center", fontsize=11, fontweight="bold")

    for label, start, end in GROUPS:
        group_fontsize = 10 if label == "Average" else 12
        ax.text(
            (start + end) / 2,
            -1.28,
            label,
            ha="center",
            va="center",
            fontsize=group_fontsize,
            fontweight="bold",
        )

    ax.set_yticks([i + 0.5 for i in range(nrows)])
    ax.set_yticklabels(row_labels, fontsize=12, fontweight="bold")
    ax.set_xticks([])
    ax.tick_params(axis="y", length=0)
    ax.set_xlim(-0.02, ncols + 0.02)
    ax.set_ylim(nrows, -1.62)

    # Table-like rules to match the reference figure.
    for x in (0, 1, 8, 10, 11):
        ax.plot([x, x], [-1.55, nrows], color="black", linewidth=0.9)
    for y in (0, 1, 2, 3, 4):
        ax.plot([0, ncols], [y, y], color="black", linewidth=0.45)
    ax.plot([0, ncols], [0, 0], color="black", linewidth=0.9)

    for spine in ax.spines.values():
        spine.set_visible(False)

    colorbar = fig.colorbar(
        plt.cm.ScalarMappable(norm=norm, cmap=cmap),
        cax=cax,
        ticks=[0, 3, 6],
    )
    colorbar.ax.set_yticklabels(["0", "3", ">6"], fontsize=16, fontweight="bold")
    colorbar.outline.set_visible(False)
    colorbar.ax.tick_params(length=8, width=1.4, direction="in", pad=2)

    output_path.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(output_path, format="pdf", bbox_inches="tight", pad_inches=0.04)
    plt.close(fig)


def interpolate_color(value: float) -> tuple[float, float, float]:
    value = max(0.0, min(6.0, value))
    stops = (
        (0.0, (255, 246, 240)),
        (3.0, (245, 178, 124)),
        (6.0, (217, 106, 19)),
    )
    if value <= 3.0:
        left, right = stops[0], stops[1]
    else:
        left, right = stops[1], stops[2]

    span = right[0] - left[0]
    t = 0.0 if span == 0 else (value - left[0]) / span
    rgb = tuple((left[1][idx] + (right[1][idx] - left[1][idx]) * t) / 255 for idx in range(3))
    return rgb


def plot_heatmap_reportlab(
    row_labels: list[str], columns: list[str], matrix: list[list[float]], output_path: Path
) -> None:
    try:
        from reportlab.lib.colors import Color, black
        from reportlab.lib.pagesizes import landscape
        from reportlab.pdfgen import canvas
    except ModuleNotFoundError as exc:
        raise ModuleNotFoundError(
            "figure7-plot.py needs matplotlib or reportlab to create the PDF. "
            "Install one of them, for example: python3 -m pip install matplotlib"
        ) from exc

    nrows = len(row_labels)
    ncols = len(columns)
    page_w, page_h = landscape((250, 720))
    grid_x = 92
    grid_y = 30
    cell_w = 48
    cell_h = 23
    grid_w = ncols * cell_w
    grid_h = nrows * cell_h
    grid_top = grid_y + grid_h

    output_path.parent.mkdir(parents=True, exist_ok=True)
    c = canvas.Canvas(str(output_path), pagesize=(page_w, page_h))

    def set_font(size: int, bold: bool = True) -> None:
        c.setFont("Helvetica-Bold" if bold else "Helvetica", size)

    # Heatmap cells and values.
    for i, row_label in enumerate(row_labels):
        y = grid_top - (i + 1) * cell_h
        for j in range(ncols):
            x = grid_x + j * cell_w
            r, g, b = interpolate_color(matrix[i][j])
            c.setFillColor(Color(r, g, b))
            c.rect(x, y, cell_w, cell_h, stroke=0, fill=1)
            c.setFillColor(black)
            set_font(10)
            c.drawCentredString(x + cell_w / 2, y + 7, format_cell(matrix[i][j], row_label))

    # Row labels.
    set_font(11)
    for i, row_label in enumerate(row_labels):
        y = grid_top - (i + 0.5) * cell_h - 4
        c.drawRightString(grid_x - 12, y, row_label)

    # Column labels.
    set_font(10)
    for j, column in enumerate(columns):
        x = grid_x + j * cell_w + cell_w / 2
        if not column:
            continue
        label = column.replace("\n", " ")
        if "\n" in column:
            parts = column.split("\n")
            c.drawCentredString(x, grid_top + 14, parts[0])
            c.drawCentredString(x, grid_top + 2, parts[1])
        else:
            c.drawCentredString(x, grid_top + 7, label)

    # Group labels.
    for label, start, end in GROUPS:
        set_font(9.5 if label == "Average" else 11)
        x = grid_x + (start + end) * cell_w / 2
        if "\n" in label:
            top, bottom = label.split("\n")
            c.drawCentredString(x, grid_top + 37, top)
            c.drawCentredString(x, grid_top + 24, bottom)
        else:
            c.drawCentredString(x, grid_top + 31, label)

    # Table rules.
    c.setStrokeColor(black)
    c.setLineWidth(0.8)
    for boundary in (0, 1, 8, 10, 11):
        x = grid_x + boundary * cell_w
        c.line(x, grid_y, x, grid_top + 42)
    c.line(grid_x, grid_top, grid_x + grid_w, grid_top)
    c.line(grid_x, grid_y, grid_x + grid_w, grid_y)
    c.setLineWidth(0.4)
    for i in range(1, nrows):
        y = grid_top - i * cell_h
        c.line(grid_x, y, grid_x + grid_w, y)

    # Colorbar.
    bar_x = grid_x + grid_w + 31
    bar_y = grid_y
    bar_w = 13
    bar_h = grid_h + 42
    steps = 80
    for step in range(steps):
        value = 6 * step / (steps - 1)
        r, g, b = interpolate_color(value)
        c.setFillColor(Color(r, g, b))
        y = bar_y + step * bar_h / steps
        c.rect(bar_x, y, bar_w, bar_h / steps + 0.5, stroke=0, fill=1)

    c.setStrokeColor(black)
    c.setLineWidth(1.1)
    for tick, label in ((0, "0"), (3, "3"), (6, ">6")):
        y = bar_y + tick / 6 * bar_h
        c.line(bar_x - 7, y, bar_x, y)
        set_font(14)
        c.drawRightString(bar_x - 10, y - 5, label)

    c.save()


def main() -> int:
    args = parse_args()
    project_root = args.project_root.resolve()
    input_path = (args.input or default_input_path(project_root)).resolve()
    output_path = (args.output or (project_root / "figures" / "figure7.pdf")).resolve()

    print(f"[INFO] input csv: {input_path}")
    print(f"[INFO] output pdf: {output_path}")

    times = load_times(input_path)
    row_labels, columns, matrix = compute_speedups(times)
    plot_heatmap(row_labels, columns, matrix, output_path)

    print("[OK] Figure 7 heatmap saved.")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
