"""Regenerate README figures with: python docs/generate_input_processing_figures.py."""

from pathlib import Path
import sys
import warnings

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle
import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))
from gritic import dataloader, sampletools

OUT = Path(__file__).resolve().parent / "figures"
OUT.mkdir(parents=True, exist_ok=True)
plt.rcParams.update({
    "font.family": "DejaVu Sans", "font.size": 11,
    "axes.spines.top": False, "axes.spines.right": False,
    "axes.titleweight": "bold",
    "savefig.facecolor": "white",
})
BLUE = "#2476a8"
RED = "#ba4949"
INK = "#203147"
GREEN = "#25856c"

cn = pd.DataFrame({
    "Chromosome": ["1"] * 3,
    "Segment_Start": [0, 100, 200], "Segment_End": [100, 200, 300],
    "Major_CN": [2, 0, 2], "Minor_CN": [1, 0, 1],
})
mutations = pd.DataFrame({
    "Chromosome": ["1"] * 3, "Position": [50, 150, 250],
    "Tumor_Ref_Count": [23] * 3, "Tumor_Alt_Count": [7] * 3,
})
with warnings.catch_warnings():
    warnings.simplefilter("ignore", UserWarning)
    explicit_zero = dataloader.prepare_sample_inputs(cn, mutations)
uncovered_gap = dataloader.prepare_sample_inputs(cn.iloc[[0, 2]], mutations)
assert explicit_zero.mutation_table.Position.tolist() == [50, 250]
assert uncovered_gap.mutation_table.Position.tolist() == [50, 150, 250]
assert explicit_zero.copy_number_table.Segment_ID.tolist() == ["1-0-300"]

fig, axes = plt.subplots(2, 1, figsize=(11, 7.2))
fig.subplots_adjust(left=.10, right=.97, top=.81, bottom=.12, hspace=.70)
fig.suptitle("A deleted interval and an uncovered gap behave differently", x=.10,
             y=.97, ha="left", fontsize=18, fontweight="bold", color=INK)
fig.text(.10, .91, "Coordinate assignment • default merging • all three SNVs pass read-count filters",
         fontsize=11, color="#526479", va="top")
for ax, zero, prepared, title in zip(
    axes, [True, False], [explicit_zero, uncovered_gap],
    ["A  Explicit 0+0 interval: its SNV is removed before merging",
     "B  No copy-number row for the gap: its SNV joins the merged interval"],
):
    ax.set_title(title, loc="left", fontsize=12, pad=12, color=INK)
    for start, end in [(0, 100), (200, 300)]:
        ax.add_patch(Rectangle((start, .69), end-start, .24, color=BLUE))
        ax.text((start+end)/2, .81, "2+1", ha="center", va="center", color="white", weight="bold")
    ax.add_patch(Rectangle((100, .69), 100, .24, facecolor="#f6dddd" if zero else "#f4f6f9",
                           edgecolor=RED if zero else "#a4afbd", linestyle="-" if zero else "--"))
    ax.text(150, .81, "0+0" if zero else "uncovered", ha="center", va="center", color=INK)
    ax.add_patch(Rectangle((0, .16), 300, .24, color=BLUE, alpha=.16))
    ax.text(150, .28, "2+1  ·  Segment_ID = 1-0-300", ha="center", va="center", color=INK)
    ax.text(-12, .81, "Input CN", ha="right", va="center", fontsize=10)
    ax.text(-12, .28, "Final CN", ha="right", va="center", fontsize=10)
    kept = set(prepared.mutation_table.Position)
    for pos in mutations.Position:
        ax.plot(pos, 1.08, "o", color=INK, markersize=5)
        ax.annotate("", xy=(pos, .46), xytext=(pos, .65),
                    arrowprops=dict(arrowstyle="->", color=GREEN if pos in kept else RED, lw=1.4))
        ax.plot(pos, .05, "o" if pos in kept else "x", color=GREEN if pos in kept else RED,
                markersize=7, markeredgewidth=2)
        ax.text(pos, -.12, "retained" if pos in kept else "removed", ha="center", va="top",
                fontsize=10, color=GREEN if pos in kept else RED)
    ax.set_xlim(-4, 304)
    ax.set_ylim(-.32, 1.15)
    ax.set_yticks([])
    ax.set_xticks([0, 100, 200, 300])
    ax.tick_params(axis="x", labelsize=10)
    ax.spines[["left", "bottom"]].set_visible(False)
fig.text(.10, .035, "Intervals are [start, end). Final width is 300 bases in both cases. "
         "--max-merge-gap 0 would keep the flanks separate.", fontsize=10, color="#526479")
fig.savefig(OUT / "segment_assignment.png", dpi=180)
plt.close(fig)

purity, total_cn, normal_cn = .6, 3, 2
state_vafs = purity * np.array([1, 2]) / (purity * total_cn + (1-purity) * normal_cn)
alt = np.arange(31)
likelihoods = sampletools._normalised_count_likelihoods(30-alt, alt, state_vafs)
fig, ax = plt.subplots(figsize=(10, 5.4))
fig.subplots_adjust(left=.10, right=.97, top=.73, bottom=.23)
fig.suptitle("Read counts give a likelihood over multiplicities", x=.10, y=.96,
             ha="left", fontsize=17, fontweight="bold", color=INK)
fig.text(.10, .865, "Example: autosomal 2+1 segment • purity 60% • depth 30 • no subclones\n"
         "Expected VAF: 23.1% at multiplicity 1; 46.2% at multiplicity 2", color="#526479", fontsize=11, va="top")
ax.axvspan(-.4, 2.5, color="#eee9e6")
ax.plot(alt, likelihoods[0], color=BLUE, lw=2.5, label="Multiplicity 1")
ax.plot(alt, likelihoods[1], color=GREEN, lw=2.5, label="Multiplicity 2")
for k in [7, 14]:
    ax.axvline(k, color="#8996a7", lw=1, linestyle=":")
    for row, color in zip(likelihoods, [BLUE, GREEN]):
        ax.plot(k, row[k], "o", color=color, markersize=6)
ax.set(xlim=(-.4, 30.4), ylim=(-.03, 1.05), xlabel="Alternate reads out of 30 total reads",
       ylabel="Normalized binomial likelihood")
ax.set_xticks([0, 3, 7, 10, 14, 20, 25, 30])
ax.legend(loc="center right", frameon=False)
ax.grid(axis="y", alpha=.15)
fig.text(.10, .035, "Shaded counts fail the default alternate-read filter (<3).\n"
         "These curves are read evidence, before route weights and detection correction.",
         fontsize=10, color="#526479")
fig.savefig(OUT / "multiplicity_likelihood.png", dpi=180)
plt.close(fig)

print("Verified merged-interval examples and regenerated both README figures.")
