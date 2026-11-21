"""Redraw Supplementary Figure S7 from the frozen shared-budget comparison table."""
from __future__ import annotations

from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import pandas as pd


ROOT = Path(__file__).resolve().parents[1]
SOURCE = ROOT / "analysis" / "rq4_outputs" / "run_20261114_231925" / "rq4c_comparisons.csv"
OUT = ROOT / "publication" / "manuscript" / "figures"
OUT.mkdir(parents=True, exist_ok=True)

NAVY = "#183b56"
BLUE = "#2c7fb8"
ORANGE = "#e69f00"
LABELS = {
    "primary_six": "Primary, six targets",
    "primary_five": "Primary, five targets (excluding AAH-I)",
    "strict_six": "Strict, six targets",
}


def main() -> None:
    comparisons = pd.read_csv(SOURCE)
    pilot = comparisons[comparisons["comparison"].eq("pilot_I_minus_raw_ipTM")].copy()

    plt.rcParams.update(
        {
            "font.family": "DejaVu Sans",
            "font.size": 9,
            "axes.spines.top": False,
            "axes.spines.right": False,
            "axes.edgecolor": "#555555",
            "axes.linewidth": 0.8,
            "legend.frameon": False,
            "figure.facecolor": "white",
            "savefig.facecolor": "white",
        }
    )
    fig, axes = plt.subplots(1, 3, figsize=(11.2, 4.5), sharey=True)
    fig.suptitle(
        "Supplementary Figure S7 | Retrospective pilot-allocation effects depend on budget and pilot size",
        x=0.02,
        ha="left",
        fontweight="bold",
        color=NAVY,
    )
    for axis, analysis_set in zip(axes, LABELS):
        frame = pilot[pilot["analysis_set"].eq(analysis_set)]
        for budget, color in ((120, BLUE), (300, ORANGE)):
            subset = frame[frame["B"].eq(budget)].sort_values("m_per_target")
            x = subset["m_per_target"].to_numpy(float)
            mean = subset["mean_hit_difference"].to_numpy(float)
            low = subset["difference_lo"].to_numpy(float)
            high = subset["difference_hi"].to_numpy(float)
            axis.errorbar(
                x,
                mean,
                yerr=[mean - low, high - mean],
                fmt="o-",
                capsize=3,
                markersize=4.5,
                linewidth=1.2,
                color=color,
                label=f"B={budget}",
            )
            for x_value, y_value in zip(x, mean):
                axis.annotate(
                    f"{y_value:+.2f}",
                    (x_value, y_value),
                    xytext=(0, 6),
                    textcoords="offset points",
                    ha="center",
                    fontsize=7,
                    color=color,
                )
        axis.axhline(0, color=NAVY, linestyle="--", linewidth=0.8)
        axis.set_xticks([10, 20, 40])
        axis.set_xlabel("Pilot size m per target")
        axis.set_title(LABELS[analysis_set], fontsize=9, fontweight="bold")
        axis.grid(axis="y", alpha=0.2)

    axes[0].set_ylabel("Mean hit difference: pilot I - raw ipTM")
    axes[-1].legend(fontsize=8, loc="lower left")
    fig.text(
        0.02,
        0.035,
        "Post hoc exploratory observed-pool simulation. Ranges are not confidence intervals for future targets; no universal pilot size is implied.",
        ha="left",
        va="bottom",
        fontsize=7.5,
        color="#555555",
    )
    fig.subplots_adjust(left=0.08, right=0.98, top=0.82, bottom=0.22, wspace=0.18)

    for suffix in (".png", ".svg", ".pdf"):
        fig.savefig(OUT / f"Supplementary_Figure_S7{suffix}", bbox_inches="tight", dpi=320)
    plt.close(fig)


if __name__ == "__main__":
    main()
