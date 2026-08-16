"""
Figure 4. Distribution of between-subgroup heterogeneity (I-squared) across the
meta-analyzed SNPs.

Panel A gives the count of SNPs in each I-squared category; panel B gives the
corresponding proportions. Categories follow the convention used in the paper:
exactly zero, 0-25%, 25-50%, 50-75%, and above 75%.

Input : Meta_analysis_results_all.csv (column I2), selected on the command line
        or supplied with --input
Output: Figure4_I2_Distribution.png (600 dpi) and .pdf, written to the current
        working directory

Palette: colorblind-safe ColorBrewer RdBu, matching Figures 1-3 and
Supplementary Figure S3. Sequential blues encode increasing heterogeneity; red
is reserved for the highest category.

Usage:  python Generate_Figure4_I2.py --input path/to/Meta_analysis_results_all.csv
"""
import argparse
import csv
import os
import sys

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

BLUE_D = "#2166ac"
RED = "#b2182b"
FILLS = ["#f0f6fb", "#c6dbef", "#92c5de", "#4393c3", RED]
LABELS = [
    "I\u00b2 = 0%",
    "0% < I\u00b2 \u2264 25%",
    "25% < I\u00b2 \u2264 50%",
    "50% < I\u00b2 \u2264 75%",
    "I\u00b2 > 75%",
]


def read_i2(path):
    values = []
    with open(path, newline="") as handle:
        reader = csv.DictReader(handle)
        if "I2" not in reader.fieldnames:
            sys.exit("Column 'I2' not found in %s" % path)
        for row in reader:
            raw = row["I2"]
            if raw in (None, "", "NA"):
                continue
            values.append(float(raw))
    if not values:
        sys.exit("No usable I2 values in %s" % path)
    return values


def bin_counts(values):
    return [
        sum(1 for v in values if v == 0),
        sum(1 for v in values if 0 < v <= 25),
        sum(1 for v in values if 25 < v <= 50),
        sum(1 for v in values if 50 < v <= 75),
        sum(1 for v in values if v > 75),
    ]


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--input", help="Meta_analysis_results_all.csv")
    parser.add_argument("--outdir", default=".")
    args = parser.parse_args()

    path = args.input
    if not path:
        path = input("Full path to Meta_analysis_results_all.csv: ").strip().strip('"')
    if not os.path.exists(path):
        sys.exit("File not found: %s" % path)

    values = read_i2(path)
    counts = bin_counts(values)
    total = len(values)
    print("SNPs read: %d" % total)
    for label, count in zip(LABELS, counts):
        print("  %-16s %3d (%.1f%%)" % (label, count, 100 * count / total))

    plt.rcParams.update({"font.family": "DejaVu Sans", "font.size": 10})
    fig, (ax_bar, ax_pie) = plt.subplots(
        1, 2, figsize=(12.5, 5.2), gridspec_kw={"width_ratios": [1.35, 1]}
    )

    bars = ax_bar.bar(range(5), counts, color=FILLS, edgecolor=BLUE_D,
                      linewidth=1.2, width=0.68)
    for count, bar in zip(counts, bars):
        ax_bar.text(bar.get_x() + bar.get_width() / 2, count + 0.8,
                    "%d\n(%.1f%%)" % (count, 100 * count / total),
                    ha="center", va="bottom", fontsize=10)
    ax_bar.set_xticks(range(5))
    ax_bar.set_xticklabels(LABELS, rotation=28, ha="right", fontsize=10)
    ax_bar.set_ylabel("Number of SNPs", fontsize=11.5)
    ax_bar.set_ylim(0, max(counts) * 1.18)
    ax_bar.spines[["top", "right"]].set_visible(False)
    ax_bar.text(-0.12, 1.04, "A", transform=ax_bar.transAxes,
                fontsize=16, fontweight="bold")

    shown = [(l, c, f) for l, c, f in zip(LABELS, counts, FILLS) if c > 0]
    wedges, _, autotexts = ax_pie.pie(
        [c for _, c, _ in shown],
        colors=[f for _, _, f in shown],
        startangle=90, counterclock=False,
        autopct=lambda pct: "%.1f%%" % pct, pctdistance=0.72,
        wedgeprops=dict(edgecolor=BLUE_D, linewidth=1.2),
    )
    for text in autotexts:
        text.set_fontsize(9.5)
    ax_pie.legend(wedges, ["%s (n = %d)" % (l, c) for l, c, _ in shown],
                  loc="center left", bbox_to_anchor=(0.97, 0.5),
                  fontsize=9.5, frameon=False)
    ax_pie.text(-0.05, 1.04, "B", transform=ax_pie.transAxes,
                fontsize=16, fontweight="bold")

    plt.tight_layout()
    png = os.path.join(args.outdir, "Figure4_I2_Distribution.png")
    pdf = os.path.join(args.outdir, "Figure4_I2_Distribution.pdf")
    fig.savefig(png, dpi=600, bbox_inches="tight", facecolor="white")
    fig.savefig(pdf, bbox_inches="tight", facecolor="white")
    plt.close(fig)
    print("Written:\n  %s\n  %s" % (png, pdf))


if __name__ == "__main__":
    main()
