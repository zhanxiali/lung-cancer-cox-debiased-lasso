"""
Figure 1. Analytic workflow  (CEBP revision, v3)

Change from v2: the free-floating "77 / 76 / 65" tally numbers that sat in the
left margin have been removed. The counts remain inside the boxes themselves
(77 candidate SNPs -> 76 SNPs present in >=1 dataset -> 65 SNPs in >=2 subgroups),
which is what Reviewer 1 asked for.

Palette: colorblind-safe ColorBrewer RdBu.
  blues  #2166ac / #92c5de / #d1e5f0 / #f0f6fb  -> data, process, and method boxes
  red    #b2182b                                -> reserved for the primary finding
This matches Figures 2-4 and Supplementary Figure S3.

Output (written to the current working directory):
  Figure1_Workflow.png   600 dpi raster
  Figure1_Workflow.pdf   vector (preferred for submission)

Requires: matplotlib.  Run with:  python Generate_Figure1_Workflow_v3.py
"""

import os

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import FancyBboxPatch, FancyArrowPatch

# ---------------------------------------------------------------- palette ----
BLUE_D = "#2166ac"   # borders, arrows, emphasis
BLUE_L = "#d1e5f0"   # light blue fill: key nodes and method boxes
BLUE_XL = "#f0f6fb"  # extra light fill: data and process boxes
RED = "#b2182b"      # accent: headline result only
GREY = "#4d4d4d"

plt.rcParams.update({"font.family": "DejaVu Sans", "font.size": 10})

fig, ax = plt.subplots(figsize=(10.5, 13.2))
ax.set_xlim(0, 10)
ax.set_ylim(0, 15.6)
ax.axis("off")


def box(x, y, w, h, text, fc=BLUE_XL, ec=BLUE_D, fs=10,
        bold=False, tc="#1a1a1a", lw=1.6):
    """Rounded box with centered, optionally bold label."""
    ax.add_patch(FancyBboxPatch(
        (x, y), w, h,
        boxstyle="round,pad=0.06,rounding_size=0.14",
        fc=fc, ec=ec, lw=lw))
    ax.text(x + w / 2, y + h / 2, text, ha="center", va="center", fontsize=fs,
            fontweight="bold" if bold else "normal", color=tc, linespacing=1.5)


def arrow(x1, y1, x2, y2, color=BLUE_D):
    ax.add_patch(FancyArrowPatch(
        (x1, y1), (x2, y2), arrowstyle="-|>",
        mutation_scale=16, lw=1.6, color=color, shrinkA=2, shrinkB=2))


# ------------------------------------------------- row 1: 77 candidate SNPs --
box(2.6, 14.55, 4.8, 0.85,
    "77 candidate SNPs\n(Kachuri et al. 2020 + pleiotropic loci)",
    fc=BLUE_L, bold=True, fs=10.5)

# --------------------------------------- row 2: 76 available (R1-1 request) --
arrow(5.0, 14.49, 5.0, 13.85)
ax.text(5.18, 14.16, "1 variant not genotyped\non any platform",
        fontsize=8.2, color=GREY, ha="left", va="center", style="italic")
box(2.6, 13.0, 4.8, 0.85, "76 SNPs present in \u22651 dataset",
    fc=BLUE_L, bold=True, fs=10.5)

# ------------------------------------------------------ row 3: the cohorts --
cohx = [0.35, 3.55, 6.75]
cw = 2.9
cohorts = [
    "HSPH\nn = 979\nIllumina arrays\n68 SNPs available",
    "OncoArray\nn = 2,322\nIllumina OncoArray\n45 SNPs available",
    "MGH\nn = 1,088\nLow-pass WGS\n62 SNPs available",
]
for cx, lab in zip(cohx, cohorts):
    box(cx, 11.15, cw, 1.45, lab, fc=BLUE_XL, fs=9.5)
for cx in cohx:
    arrow(5.0, 12.94, cx + cw / 2, 12.66)

# ------------------------------------------------------ QC and merge bands --
for cx in cohx:
    arrow(cx + cw / 2, 11.09, cx + cw / 2, 10.72)
box(0.35, 9.95, 9.3, 0.72,
    "Genotype QC (call rate, MAF, HWE)  \u2192  PCA (PC1\u2013PC3)  \u2192  "
    "allele alignment (HSPH as reference)",
    fc=BLUE_XL, fs=9.5)

arrow(5.0, 9.89, 5.0, 9.52)
box(0.35, 8.75, 9.3, 0.72,
    "Merge SNPs + clinical covariates "
    "(age, sex, smoking, stage, treatment, PC1\u2013PC3)",
    fc=BLUE_XL, fs=9.5)

# ------------------------------------------- Cox debiased lasso per subgroup --
for cx in cohx:
    arrow(cx + cw / 2, 8.69, cx + cw / 2, 8.32)
methods = [
    "Cox debiased lasso\n(HSPH)\nL1 penalty \u2192 CV\nscore-based debiasing",
    "Cox debiased lasso\n(OncoArray)\nL1 penalty \u2192 CV\nscore-based debiasing",
    "Cox debiased lasso\n(MGH)\nL1 penalty \u2192 CV\nscore-based debiasing",
]
for cx, lab in zip(cohx, methods):
    box(cx, 6.85, cw, 1.42, lab, fc=BLUE_L, fs=9.3)

for cx in cohx:
    arrow(cx + cw / 2, 6.79, cx + cw / 2, 6.44)
outputs = [
    "\u03b2, SE, P\n(79 covariates)",
    "\u03b2, SE, P\n(56 covariates)",
    "\u03b2, SE, P\n(73 covariates)",
]
for cx, lab in zip(cohx, outputs):
    box(cx, 5.6, cw, 0.8, lab, fc=BLUE_XL, fs=9.3)

# ------------------------------------------------------ allele harmonization --
for cx in cohx:
    arrow(cx + cw / 2, 5.54, 5.0, 5.2)
box(0.35, 4.45, 9.3, 0.72,
    "Allele harmonization: EAF concordance check + sign correction (28 SNPs) "
    "\u2192 unified effect direction",
    fc=BLUE_XL, fs=9.5)

# --------------------------------------------------- meta-analysis (65 SNPs) --
arrow(5.0, 4.39, 5.0, 4.02)
box(1.3, 3.05, 7.4, 0.92,
    "Fixed-effect meta-analysis: 65 SNPs available in \u22652 subgroups\n"
    "IVW pooling + Cochran Q + I\u00b2 + HKSJ random-effects sensitivity",
    fc=BLUE_L, bold=True, fs=9.8)

arrow(5.0, 2.99, 5.0, 2.62)
box(1.3, 1.55, 7.4, 1.02,
    "FDR correction (BH, across 65 SNPs) + evidence hierarchy\n"
    "Tier 1: FDR q < 0.05 and available in all 3 subgroups\n"
    "Tier 2: P < 0.05 not meeting Tier 1",
    fc=BLUE_XL, fs=9.3)

# ---------------------------------------------------------- primary finding --
arrow(5.0, 1.49, 5.0, 1.12)
box(1.3, 0.18, 7.4, 0.9,
    "rs11022690 (Tier 1): HR = 0.927, P = 7.79\u00d710\u207b\u2074, FDR = 0.050\n"
    "10 additional Tier 2 signals",
    fc=RED, ec=RED, bold=True, tc="white", fs=10)

# NOTE (v3): the left-margin "77 / 76 / 65" tally labels that appeared here in
# v2 have been removed; the counts are already stated inside the boxes.

plt.tight_layout()
out_png = os.path.join(os.getcwd(), "Figure1_Workflow.png")
out_pdf = os.path.join(os.getcwd(), "Figure1_Workflow.pdf")
fig.savefig(out_png, dpi=600, bbox_inches="tight", facecolor="white")
fig.savefig(out_pdf, bbox_inches="tight", facecolor="white")
plt.close(fig)

print("Written:")
print("  " + out_png)
print("  " + out_pdf)
