"""Redraw post-review manuscript figures from frozen numerical tables."""
from __future__ import annotations

import hashlib
import json
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import FancyArrowPatch, FancyBboxPatch
import numpy as np
import pandas as pd


ROOT = Path(__file__).resolve().parents[3]
OUT = ROOT / "analysis" / "postreview_revision_2026-11-21" / "figures"
OUT.mkdir(parents=True, exist_ok=True)
MAIN_FIG = ROOT / "publication" / "figures" / "main"
SUPP_FIG = ROOT / "publication" / "figures" / "supplementary"
MAIN_FIG.mkdir(parents=True, exist_ok=True)
SUPP_FIG.mkdir(parents=True, exist_ok=True)

TARGET_ORDER = [
    "aah-i",
    "human-vegf110",
    "nipah-glycoprotein-g",
    "pd-l1",
    "sars-cov-2-rbd",
    "trem2",
]
TARGET_LABEL = {
    "aah-i": "AAH-I",
    "human-vegf110": "VEGF110",
    "nipah-glycoprotein-g": "Nipah G",
    "pd-l1": "PD-L1",
    "sars-cov-2-rbd": "SARS-CoV-2\nRBD",
    "trem2": "TREM2",
}
COLORS = {
    "navy": "#244a73",
    "blue": "#4c78a8",
    "teal": "#1b9e77",
    "purple": "#6f4c9b",
    "orange": "#e69f00",
    "red": "#c44e52",
    "gray": "#8c8c8c",
    "light_gray": "#d9d9d9",
    "ink": "#202020",
}
LAYER_COLOR = {"U": COLORS["gray"], "I": COLORS["teal"], "F": COLORS["purple"]}
LAYER_LABEL = {
    "U": "U (no target labels)",
    "I": "I (local intercept)",
    "F": "F (local intercept + slope)",
}


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def style() -> None:
    plt.rcParams.update(
        {
            "font.family": "DejaVu Sans",
            "font.size": 9,
            "axes.titlesize": 11,
            "axes.labelsize": 9,
            "axes.spines.top": False,
            "axes.spines.right": False,
            "axes.edgecolor": "#555555",
            "axes.linewidth": 0.8,
            "xtick.major.width": 0.7,
            "ytick.major.width": 0.7,
            "legend.frameon": False,
            "figure.facecolor": "white",
            "savefig.facecolor": "white",
        }
    )


def save(fig: plt.Figure, stem: str, destination: Path) -> dict:
    destination.mkdir(parents=True, exist_ok=True)
    files = {}
    for suffix, kwargs in (
        (".png", {"dpi": 320}),
        (".pdf", {}),
        (".svg", {}),
    ):
        path = destination / f"{stem}{suffix}"
        fig.savefig(path, bbox_inches="tight", **kwargs)
        files[suffix] = {
            "path": str(path.relative_to(ROOT)),
            "bytes": path.stat().st_size,
            "sha256": sha256(path),
        }
    plt.close(fig)
    return files


def panel_label(ax: plt.Axes, text: str) -> None:
    ax.text(
        -0.11,
        1.05,
        text,
        transform=ax.transAxes,
        fontsize=13,
        fontweight="bold",
        va="top",
        ha="left",
        color=COLORS["ink"],
    )


def add_box(
    ax: plt.Axes,
    x: float,
    y: float,
    width: float,
    height: float,
    text: str,
    *,
    face: str,
    edge: str,
) -> None:
    patch = FancyBboxPatch(
        (x - width / 2, y - height / 2),
        width,
        height,
        boxstyle="round,pad=0.02,rounding_size=0.02",
        facecolor=face,
        edgecolor=edge,
        linewidth=1.0,
    )
    ax.add_patch(patch)
    ax.text(
        x,
        y,
        text,
        ha="center",
        va="center",
        fontsize=8.2,
        color=COLORS["ink"],
        wrap=True,
    )


def arrow(ax: plt.Axes, x1: float, y1: float, x2: float, y2: float) -> None:
    ax.add_patch(
        FancyArrowPatch(
            (x1, y1),
            (x2, y2),
            arrowstyle="-|>",
            mutation_scale=11,
            linewidth=0.9,
            color="#666666",
        )
    )


def figure1() -> dict:
    fig, ax = plt.subplots(figsize=(12.2, 6.2))
    ax.set_xlim(0, 1)
    ax.set_ylim(0, 1)
    ax.axis("off")
    ax.text(0.5, 0.98, "Study cohorts and evidence branches", ha="center", va="top", fontsize=14, fontweight="bold")

    ax.text(0.24, 0.90, "Proteinbase main source", ha="center", weight="bold", color=COLORS["navy"])
    add_box(ax, 0.24, 0.79, 0.34, 0.105, "4,280 designs\n4,943 design-target pairs", face="#eaf1f8", edge=COLORS["navy"])
    add_box(ax, 0.24, 0.63, 0.34, 0.10, "Six prespecified targets\n2,185 design-target pairs", face="#eaf1f8", edge=COLORS["navy"])
    add_box(ax, 0.24, 0.47, 0.34, 0.10, "2,084 with released ipTM\n101 without ipTM (excluded)", face="#dbe8f4", edge=COLORS["navy"])
    add_box(ax, 0.24, 0.29, 0.40, 0.11, "Assayed design pool\nOutcome is experimental binding after upstream selection", face="#f5f5f5", edge="#777777")
    add_box(ax, 0.24, 0.14, 0.34, 0.07, "2,074 complete cases for\nPinc / ipSAE sensitivity", face="#f7f1e3", edge=COLORS["orange"])
    arrow(ax, 0.24, 0.737, 0.24, 0.683)
    arrow(ax, 0.24, 0.578, 0.24, 0.523)
    arrow(ax, 0.24, 0.418, 0.24, 0.348)
    arrow(ax, 0.24, 0.233, 0.24, 0.177)

    ax.text(0.75, 0.90, "Claude binder-design source", ha="center", weight="bold", color=COLORS["navy"])
    add_box(ax, 0.75, 0.79, 0.34, 0.105, "1,440 designs / 16 targets\nindependent campaign", face="#eaf1f8", edge=COLORS["navy"])
    add_box(ax, 0.75, 0.63, 0.34, 0.10, "Exclude 120 Mature GDF-8\n1,320 interpretable designs", face="#eaf1f8", edge=COLORS["navy"])
    add_box(ax, 0.75, 0.47, 0.34, 0.10, "≥80 usable pairs and\n≥15 positives per target", face="#eaf1f8", edge=COLORS["navy"])
    add_box(ax, 0.75, 0.30, 0.34, 0.10, "720 pairs / 8 targets\nfrozen external evaluation", face="#dbe8f4", edge=COLORS["navy"])
    add_box(ax, 0.75, 0.14, 0.34, 0.08, "521 pairs / 6 low-event targets\ndescriptive stress test", face="#f3f3f3", edge="#777777")
    arrow(ax, 0.75, 0.737, 0.75, 0.683)
    arrow(ax, 0.75, 0.578, 0.75, 0.523)
    arrow(ax, 0.75, 0.418, 0.75, 0.352)
    arrow(ax, 0.75, 0.248, 0.75, 0.182)

    ax.plot([0.54, 0.54], [0.10, 0.84], color="#bdbdbd", linewidth=0.8, linestyle=(0, (4, 4)))
    ax.text(
        0.24,
        0.06,
        "Other Proteinbase targets: 2,758 pairs\n809 with ipTM (outside six-target scope); 1,949 without ipTM",
        ha="center",
        va="center",
        fontsize=7.5,
        color="#555555",
    )
    ax.text(
        0.75,
        0.06,
        "79 not analyzed: 60 Latent GDF-8 (single-vendor protocol);\n19 without Adaptyv-valid outcome (BBF-14 11, Cas9 8)",
        ha="center",
        va="center",
        fontsize=7.5,
        color="#555555",
    )
    ax.text(0.525, 0.005, "Blue: primary/confirmatory cohorts   Orange: post hoc complete-case sensitivity", ha="center", fontsize=7.5, color="#555555")
    return save(fig, "Figure_1_study_design_postreview", MAIN_FIG)


def figure2() -> dict:
    pair = pd.read_csv(ROOT / "analysis" / "v3_outputs" / "pair_master_v3.csv", low_memory=False)
    metrics = pd.read_csv(ROOT / "analysis" / "rq12_outputs" / "rq12_target_metrics.csv", low_memory=False)
    metrics = metrics[
        metrics["analysis_set"].eq("primary_six")
        & metrics["metric"].isin(["roc_auc", "average_precision"])
    ]
    fig = plt.figure(figsize=(12.4, 10.0))
    gs = fig.add_gridspec(2, 1, height_ratios=(1.05, 1.0), hspace=0.38)
    ax1 = fig.add_subplot(gs[0])
    ax2 = fig.add_subplot(gs[1])
    positions = np.arange(len(TARGET_ORDER))
    for index, target in enumerate(TARGET_ORDER):
        for offset, outcome, color, label in (
            (-0.18, 0, COLORS["gray"], "Negative"),
            (0.18, 1, COLORS["teal"], "Positive"),
        ):
            values = pair[
                pair["target"].eq(target)
                & pair["primary_outcome"].astype(int).eq(outcome)
            ]["boltz2_iptm"].to_numpy(dtype=float)
            bp = ax1.boxplot(
                [values],
                positions=[index + offset],
                widths=0.28,
                patch_artist=True,
                showfliers=False,
                medianprops={"color": "white", "linewidth": 1.1},
                whiskerprops={"color": color, "linewidth": 0.8},
                capprops={"color": color, "linewidth": 0.8},
            )
            bp["boxes"][0].set(facecolor=color, edgecolor=color, alpha=0.75, label=label if index == 0 else None)
    ax1.set_xticks(positions, [TARGET_LABEL[t] for t in TARGET_ORDER])
    ax1.set_ylabel("Boltz-2 ipTM")
    ax1.set_ylim(0.25, 1.0)
    ax1.set_title("Raw ipTM distributions by target and observed primary outcome")
    ax1.legend(
        ncol=2,
        loc="upper right",
        bbox_to_anchor=(1.0, -0.14),
        borderaxespad=0,
    )
    panel_label(ax1, "A")

    width = 0.34
    auc = metrics[metrics["metric"].eq("roc_auc")].set_index("target").loc[TARGET_ORDER]
    ap = metrics[metrics["metric"].eq("average_precision")].set_index("target").loc[TARGET_ORDER]
    y = np.arange(len(TARGET_ORDER))
    for values, offset, color, label in (
        (auc, -width / 2, COLORS["navy"], "AUC"),
        (ap, width / 2, COLORS["orange"], "AP"),
    ):
        errors = np.vstack(
            [
                values["estimate"] - values["ci_low"],
                values["ci_high"] - values["estimate"],
            ]
        )
        ax2.barh(y + offset, values["estimate"], height=width, color=color, alpha=0.86, label=label)
        ax2.errorbar(
            values["estimate"],
            y + offset,
            xerr=errors,
            fmt="none",
            ecolor=COLORS["ink"],
            elinewidth=0.8,
            capsize=2,
        )
    prevalence = pair.groupby("target")["primary_outcome"].mean().loc[TARGET_ORDER]
    for index, target in enumerate(TARGET_ORDER):
        ax2.scatter(prevalence[target], index + width / 2, marker="D", s=22, color=COLORS["red"], zorder=4)
    ax2.scatter([], [], marker="D", s=22, color=COLORS["red"], label="Target prevalence")
    ax2.set_yticks(y, [TARGET_LABEL[t] for t in TARGET_ORDER])
    ax2.set_xlim(0, 1)
    ax2.set_xlabel("Metric value")
    ax2.set_title("Within-target discrimination with 95% sequence-cluster intervals")
    ax2.legend(
        loc="upper left",
        bbox_to_anchor=(0.0, -0.18),
        ncol=3,
        borderaxespad=0,
    )
    panel_label(ax2, "B")
    fig.suptitle("ipTM ranking signal varies strongly by target", fontsize=14, fontweight="bold", y=0.98)
    return save(fig, "Figure_2_ipTM_discrimination_postreview", MAIN_FIG)


def figure3() -> dict:
    reliability = pd.read_csv(ROOT / "analysis" / "rq12_outputs" / "rq2_reliability_bins.csv", low_memory=False)
    reliability = reliability[reliability["analysis_set"].eq("primary_six")]
    hierarchy = pd.read_csv(ROOT / "analysis" / "rq12_outputs" / "rq2_hierarchical_targets.csv", low_memory=False)
    hierarchy = hierarchy[hierarchy["analysis_set"].eq("primary_six")].set_index("target").loc[TARGET_ORDER]
    fig = plt.figure(figsize=(12.4, 9.0))
    gs = fig.add_gridspec(
        2,
        2,
        height_ratios=(1.15, 0.85),
        hspace=0.38,
        wspace=0.28,
    )
    axes = [
        fig.add_subplot(gs[0, :]),
        fig.add_subplot(gs[1, 0]),
        fig.add_subplot(gs[1, 1]),
    ]
    for target in TARGET_ORDER:
        frame = reliability[reliability["target"].eq(target) & reliability["n"].gt(0)]
        axes[0].plot(frame["mean_iptm"], frame["observed_rate"], marker="o", markersize=3.5, linewidth=1.2, label=TARGET_LABEL[target].replace("\n", " "))
    axes[0].plot([0, 1], [0, 1], linestyle="--", color="#999999", linewidth=0.8)
    axes[0].set(xlim=(0, 1), ylim=(0, 1), xlabel="Mean ipTM in bin", ylabel="Observed binding rate")
    axes[0].set_title("Reliability by target")
    axes[0].legend(fontsize=7, loc="upper left")
    panel_label(axes[0], "A")

    y = np.arange(len(TARGET_ORDER))
    for idx, target in enumerate(TARGET_ORDER):
        row = hierarchy.loc[target]
        axes[1].plot([row.alpha_unpooled, row.alpha_partial], [idx, idx], color="#bbbbbb", linewidth=1.3)
        axes[1].scatter(row.alpha_unpooled, idx, color=COLORS["orange"], s=24, label="Unpooled" if idx == 0 else None)
        axes[1].scatter(row.alpha_partial, idx, color=COLORS["navy"], s=24, label="Partially pooled" if idx == 0 else None)
    axes[1].set_yticks(y, [TARGET_LABEL[t].replace("\n", " ") for t in TARGET_ORDER])
    axes[1].set_xlabel("Intercept on clipped logit-ipTM")
    axes[1].set_title("Intercept heterogeneity")
    axes[1].legend(loc="lower right", fontsize=7)
    panel_label(axes[1], "B")

    for idx, target in enumerate(TARGET_ORDER):
        row = hierarchy.loc[target]
        axes[2].plot([row.beta_unpooled, row.beta_partial], [idx, idx], color="#bbbbbb", linewidth=1.3)
        axes[2].scatter(row.beta_unpooled, idx, color=COLORS["orange"], s=24, label="Unpooled" if idx == 0 else None)
        axes[2].scatter(row.beta_partial, idx, color=COLORS["navy"], s=24, label="Partially pooled" if idx == 0 else None)
    axes[2].set_yticks(y, [TARGET_LABEL[t].replace("\n", " ") for t in TARGET_ORDER])
    axes[2].set_xlabel("Slope on clipped logit-ipTM")
    axes[2].set_title("Slope heterogeneity")
    axes[2].legend(loc="lower right", fontsize=7)
    panel_label(axes[2], "C")
    fig.suptitle("Raw ipTM is not calibrated and calibration parameters vary by target", fontsize=14, fontweight="bold", y=1.02)
    return save(fig, "Figure_3_raw_calibration_postreview", MAIN_FIG)


def figure4() -> dict:
    main = pd.read_csv(
        ROOT / "analysis" / "codex_outputs" / "RQ3_F_correction_2026-11-13" / "run_20261113_182311" / "rq3_target_adaptation_metrics.csv",
        low_memory=False,
    )
    main = main[
        main["analysis_set"].eq("primary_six")
        & main["cluster_rule"].eq("sequence_cluster90")
        & main["source_weighting"].eq("target_equal")
    ].copy()
    main["layer"] = main["model_layer_label"].map(
        {
            "U_transport_zero_label": "U",
            "predictive_I": "I",
            "predictive_F": "F",
        }
    )
    external = pd.read_csv(ROOT / "analysis" / "rq3e_if_outputs" / "run_20261114_174403" / "target_metrics.csv")
    intervals = pd.read_csv(
        ROOT / "analysis" / "postreview_revision_2026-11-21" / "bootstrap" / "bootstrap_macro_intervals.csv"
    )
    rel = intervals[intervals["metric"].eq("relative_macro_U_to_I_brier_reduction")].copy()
    rel["label"] = rel["dataset"].map(
        {
            "main_primary_six": "Main six",
            "main_primary_five": "Main five",
            "external_primary_weak_plus": "External eight",
        }
    )
    rel["plot_order"] = rel["dataset"].map(
        {
            "main_primary_six": 0,
            "main_primary_five": 1,
            "external_primary_weak_plus": 2,
        }
    )
    rel = rel.sort_values("plot_order").reset_index(drop=True)
    fig = plt.figure(figsize=(12.4, 9.2))
    gs = fig.add_gridspec(
        2,
        2,
        height_ratios=(1.0, 0.78),
        hspace=0.38,
        wspace=0.24,
    )
    axes = [
        fig.add_subplot(gs[0, 0]),
        fig.add_subplot(gs[0, 1]),
        fig.add_subplot(gs[1, :]),
    ]
    width = 0.24
    x = np.arange(len(TARGET_ORDER))
    for offset, layer in zip((-width, 0, width), ("U", "I", "F")):
        frame = main[main["layer"].eq(layer)].set_index("held_out_target").loc[TARGET_ORDER]
        axes[0].bar(
            x + offset,
            frame["bss_transport"],
            width=width,
            color=LAYER_COLOR[layer],
            label=LAYER_LABEL[layer],
        )
    axes[0].axhline(0, color="#444444", linewidth=0.8)
    axes[0].set_xticks(x, [TARGET_LABEL[t] for t in TARGET_ORDER], rotation=25, ha="right")
    axes[0].set_ylabel("BSS-transport")
    axes[0].set_title("Main cohort: held-out targets")
    axes[0].legend(fontsize=7, ncol=1, loc="lower left")
    panel_label(axes[0], "A")

    ext_targets = ["BHRF1", "IL-7Ra", "Nipah-G", "PD-L1", "RBX1", "TREM2", "TrkA", "VEGF-A"]
    x2 = np.arange(len(ext_targets))
    for offset, layer in zip((-width, 0, width), ("U", "I", "F")):
        frame = external[external["layer"].eq(layer)].set_index("target").loc[ext_targets]
        axes[1].bar(
            x2 + offset,
            frame["bss_transport"],
            width=width,
            color=LAYER_COLOR[layer],
        )
    axes[1].axhline(0, color="#444444", linewidth=0.8)
    axes[1].set_xticks(x2, ext_targets, rotation=30, ha="right")
    axes[1].set_ylabel("BSS-transport")
    axes[1].set_title("External cohort: target/campaign context")
    panel_label(axes[1], "B")

    y = np.arange(len(rel))
    axes[2].errorbar(
        rel["point_estimate"],
        y,
        xerr=np.vstack(
            [
                rel["point_estimate"] - rel["bootstrap_95_low"],
                rel["bootstrap_95_high"] - rel["point_estimate"],
            ]
        ),
        fmt="o",
        color=COLORS["navy"],
        ecolor=COLORS["navy"],
        capsize=3,
        markersize=5,
    )
    axes[2].set_yticks(y, rel["label"])
    axes[2].invert_yaxis()
    axes[2].set_xlim(0.1, 0.55)
    axes[2].set_xlabel("Relative macro-Brier reduction, U to I")
    axes[2].set_title("Intercept gain with 95% interval")
    panel_label(axes[2], "C")
    fig.suptitle("Local intercept updating improves held-out probability error, but transport remains context-dependent", fontsize=13.5, fontweight="bold", y=1.02)
    return save(fig, "Figure_4_transport_adaptation_postreview", MAIN_FIG)


def figure5() -> dict:
    metrics = pd.read_csv(
        ROOT / "analysis" / "rq4_outputs" / "run_20261114_231925" / "rq4ab_target_metrics.csv",
        low_memory=False,
    )
    fraction = metrics[
        metrics["analysis_set"].eq("primary_six")
        & metrics["analysis_family"].eq("RQ4-A")
        & metrics["budget_kind"].eq("fraction")
    ].copy()
    threshold = metrics[
        metrics["analysis_set"].eq("primary_six")
        & metrics["analysis_family"].eq("RQ4-B")
        & metrics["budget_kind"].eq("threshold")
    ].copy()
    comparisons = pd.read_csv(
        ROOT / "analysis" / "rq4_outputs" / "run_20261114_231925" / "rq4c_comparisons.csv"
    )
    comparisons = comparisons[
        comparisons["analysis_set"].eq("primary_six")
        & comparisons["comparison"].isin(["pilot_I_minus_raw_ipTM", "pilot_F_minus_raw_ipTM"])
    ].copy()
    fig = plt.figure(figsize=(12.4, 9.2))
    gs = fig.add_gridspec(
        2,
        2,
        height_ratios=(1.0, 0.82),
        hspace=0.40,
        wspace=0.30,
    )
    axes = [
        fig.add_subplot(gs[0, 0]),
        fig.add_subplot(gs[0, 1]),
        fig.add_subplot(gs[1, :]),
    ]
    budgets = [0.05, 0.10, 0.20]
    budget_colors = ["#9ecae1", "#4292c6", "#08519c"]
    x = np.arange(len(TARGET_ORDER))
    width = 0.23
    for offset, budget, color in zip((-width, 0, width), budgets, budget_colors):
        frame = fraction[fraction["budget_value"].eq(budget)].set_index("target").loc[TARGET_ORDER]
        axes[0].bar(x + offset, frame["enrichment"], width=width, color=color, label=f"Top {int(budget*100)}%")
    axes[0].axhline(1, linestyle="--", color="#666666", linewidth=0.8)
    axes[0].set_xticks(x, [TARGET_LABEL[t] for t in TARGET_ORDER], rotation=25, ha="right")
    axes[0].set_ylabel("Enrichment over same-target random")
    axes[0].set_title("Proportional selection")
    axes[0].legend(fontsize=7)
    panel_label(axes[0], "A")

    thresholds = [0.80, 0.90, 0.95]
    matrix = np.full((len(TARGET_ORDER), len(thresholds)), np.nan)
    for row, target in enumerate(TARGET_ORDER):
        for column, value in enumerate(thresholds):
            subset = threshold[
                threshold["target"].eq(target) & threshold["budget_value"].eq(value)
            ]
            if len(subset):
                matrix[row, column] = float(subset.iloc[0]["selected_fraction"])
    image = axes[1].imshow(matrix, cmap="Blues", vmin=0, vmax=1, aspect="auto")
    for row in range(matrix.shape[0]):
        for column in range(matrix.shape[1]):
            value = matrix[row, column]
            axes[1].text(
                column,
                row,
                f"{100*value:.0f}%",
                ha="center",
                va="center",
                fontsize=7,
                color="white" if value > 0.55 else COLORS["ink"],
            )
    axes[1].set_xticks(range(len(thresholds)), ["0.80", "0.90", "0.95"])
    axes[1].set_yticks(range(len(TARGET_ORDER)), [TARGET_LABEL[t].replace("\n", " ") for t in TARGET_ORDER])
    axes[1].set_xlabel("Raw ipTM threshold")
    axes[1].set_title("Selected fraction")
    fig.colorbar(image, ax=axes[1], fraction=0.05, pad=0.04, label="Selected fraction")
    panel_label(axes[1], "B")

    scenario_rows = (
        comparisons[["B", "m_per_target"]]
        .drop_duplicates()
        .sort_values(["B", "m_per_target"])
        .reset_index(drop=True)
    )
    scenario_rows["x"] = np.arange(len(scenario_rows))
    comparisons = comparisons.merge(
        scenario_rows,
        on=["B", "m_per_target"],
        how="left",
        validate="many_to_one",
    )
    for offset, comparison, color, label in (
        (-0.10, "pilot_I_minus_raw_ipTM", COLORS["teal"], "Pilot I - raw"),
        (0.10, "pilot_F_minus_raw_ipTM", COLORS["purple"], "Pilot F - raw"),
    ):
        frame = comparisons[comparisons["comparison"].eq(comparison)].sort_values("x")
        yerr = np.vstack(
            [
                frame["mean_hit_difference"] - frame["difference_lo"],
                frame["difference_hi"] - frame["mean_hit_difference"],
            ]
        )
        axes[2].errorbar(
            frame["x"].to_numpy() + offset,
            frame["mean_hit_difference"],
            yerr=yerr,
            fmt="o",
            markersize=4,
            capsize=2,
            color=color,
            label=label,
        )
    axes[2].axhline(0, color="#444444", linewidth=0.8)
    axes[2].set_xticks(
        scenario_rows["x"],
        [
            f"B={int(b)}, m={int(m)}"
            for b, m in zip(scenario_rows["B"], scenario_rows["m_per_target"])
        ],
        rotation=30,
        ha="right",
    )
    axes[2].set_ylabel("Total-hit difference")
    axes[2].set_title("Shared-budget pilot scenarios")
    axes[2].legend(fontsize=7, bbox_to_anchor=(1.0, 1.0), loc="upper left")
    panel_label(axes[2], "C")
    fig.suptitle("Retrospective pilot-allocation effects depend on budget and pilot size", fontsize=13.5, fontweight="bold", y=1.02)
    return save(fig, "Figure_5_experimental_selection_postreview", MAIN_FIG)


def figure6() -> dict:
    raw = pd.read_csv(
        ROOT / "analysis" / "codex_outputs" / "RQ5_execution_2026-11-15" / "run_20261115_111720" / "rq5_score_macro_summary.csv"
    )
    raw = raw[
        raw["analysis_set"].eq("primary_six")
        & raw["cluster_rule"].eq("sequence_cluster90")
        & raw["summary_type"].eq("target_equal_macro")
        & raw["metric"].isin(["roc_auc", "ap"])
    ]
    pairwise = pd.read_csv(
        ROOT / "analysis" / "codex_outputs" / "RQ5_execution_2026-11-15" / "run_20261115_111720" / "rq5_pairwise_comparisons.csv"
    )
    pairwise = pairwise[
        pairwise["analysis_set"].eq("primary_six")
        & pairwise["cluster_rule"].eq("sequence_cluster90")
        & pairwise["source_weighting"].eq("target_equal")
        & pairwise["comparison"].eq("ipTM_plus_Pinc_minus_ipTM_only")
        & pairwise["summary_type"].eq("target_equal_macro")
        & pairwise["layer"].eq("U_transport_zero_label")
        & pairwise["metric"].isin(["brier", "log_loss", "bss_transport"])
    ]
    fig, axes = plt.subplots(1, 2, figsize=(10.8, 4.9), gridspec_kw={"width_ratios": (1.0, 1.25)})
    metric_order = ["roc_auc", "ap"]
    x = np.arange(len(metric_order))
    for offset, score, color, label in (
        (-0.16, "ipTM", COLORS["navy"], "ipTM"),
        (0.16, "Pinc", COLORS["orange"], "Pinc"),
    ):
        frame = raw[raw["score"].eq(score)].set_index("metric").loc[metric_order]
        axes[0].bar(x + offset, frame["mean"], width=0.30, color=color, label=label)
    axes[0].set_xticks(x, ["AUC", "AP"])
    axes[0].set_ylim(0, 0.8)
    axes[0].set_ylabel("Target-equal mean")
    axes[0].set_title("Raw discrimination")
    axes[0].legend()
    panel_label(axes[0], "A")

    labels = {"brier": "Brier", "log_loss": "Log loss", "bss_transport": "BSS-transport"}
    values = pairwise.set_index("metric").loc[["brier", "log_loss", "bss_transport"]]
    colors = [
        COLORS["red"] if value > 0 else COLORS["teal"]
        for value in values["mean_difference"]
    ]
    axes[1].barh(np.arange(3), values["mean_difference"], color=colors, alpha=0.85)
    axes[1].axvline(0, color="#444444", linewidth=0.8)
    axes[1].set_yticks(np.arange(3), [labels[k] for k in ["brier", "log_loss", "bss_transport"]])
    axes[1].set_xlabel("ipTM+Pinc minus ipTM-only")
    axes[1].set_title("U-layer incremental change")
    axes[1].text(
        0.98,
        0.03,
        "Lower Brier/log loss is better;\nhigher BSS-transport is better",
        transform=axes[1].transAxes,
        ha="right",
        va="bottom",
        fontsize=7,
        color="#555555",
    )
    panel_label(axes[1], "B")
    fig.suptitle("Pinc changes some probability summaries without a stable ranking gain", fontsize=13.5, fontweight="bold", y=1.02)
    return save(fig, "Figure_6_Pinc_increment_postreview", MAIN_FIG)


def figure7() -> dict:
    summary = pd.read_csv(
        ROOT / "analysis" / "codex_outputs" / "RQ6_execution_2026-11-15" / "run_20261115_151147" / "rq6_target_equal_summary.csv"
    )
    summary = summary[
        summary["metric"].eq("auc")
        & summary["outcome"].isin(["adaptyv_weak_plus", "adaptyv_medium_plus"])
    ].copy()
    predictor_order = ["rf3", "afm3", "ptxv2", "boltz2", "chai1", "of3", "ef2fast", "ef2full", "odde", "af3of3"]
    predictor_label = {
        "boltz2": "Boltz-2",
        "afm3": "AFM3",
        "chai1": "Chai-1",
        "ef2fast": "EF2-fast",
        "ef2full": "EF2-full",
        "ptxv2": "PTXv2",
        "of3": "OF3",
        "odde": "ODDE",
        "rf3": "RF3",
        "af3of3": "AF3/OF3",
    }
    fig, ax = plt.subplots(figsize=(10.8, 4.8))
    y = np.arange(len(predictor_order))
    for offset, outcome, color, label in (
        (-0.15, "adaptyv_weak_plus", COLORS["teal"], "Adaptyv Weak+"),
        (0.15, "adaptyv_medium_plus", COLORS["purple"], "Adaptyv Medium+"),
    ):
        frame = summary[summary["outcome"].eq(outcome)].set_index("predictor").loc[predictor_order]
        xerr = np.vstack(
            [
                frame["target_equal_mean"] - frame["target_equal_min"],
                frame["target_equal_max"] - frame["target_equal_mean"],
            ]
        )
        ax.errorbar(
            frame["target_equal_mean"],
            y + offset,
            xerr=xerr,
            fmt="o",
            markersize=4,
            capsize=2,
            linewidth=1,
            color=color,
            label=label,
        )
    ax.set_yticks(y, [predictor_label[p] for p in predictor_order])
    ax.set_xlabel("Target-equal AUC (mean and observed target range)")
    ax.set_xlim(0.55, 0.80)
    ax.legend(ncol=2, loc="upper left", bbox_to_anchor=(0.0, -0.16))
    ax.set_title("No predictor has a uniformly best ranking across targets and outcomes", fontsize=13, fontweight="bold")
    return save(fig, "Figure_7_predictor_pipelines_postreview", MAIN_FIG)


def supplementary_s11() -> dict:
    ipsae = pd.read_csv(
        ROOT / "analysis" / "postreview_revision_2026-11-21" / "ipSAE" / "ipSAE_macro_metrics.csv"
    )
    combo = pd.read_csv(
        ROOT / "analysis" / "postreview_revision_2026-11-21" / "ipSAE_combo" / "combo_paired_comparisons.csv"
    )
    intervals = pd.read_csv(
        ROOT / "analysis" / "postreview_revision_2026-11-21" / "bootstrap" / "bootstrap_macro_intervals.csv"
    )
    fig, axes = plt.subplots(2, 2, figsize=(12.0, 8.0), gridspec_kw={"wspace": 0.34, "hspace": 0.48})
    axes = axes.ravel()

    raw = ipsae[
        ipsae["source_weighting"].eq("target_equal") & ipsae["layer"].eq("raw")
    ].set_index("model")
    metric_order = ["roc_auc", "ap"]
    x = np.arange(2)
    for offset, model, color, label in (
        (-0.16, "ipTM_only", COLORS["navy"], "ipTM"),
        (0.16, "ipSAE_only", COLORS["orange"], "ipSAE"),
    ):
        axes[0].bar(x + offset, raw.loc[model, metric_order], width=0.30, color=color, label=label)
    axes[0].set_xticks(x, ["AUC", "AP"])
    axes[0].set_ylim(0, 0.8)
    axes[0].set_ylabel("Target-equal mean")
    axes[0].set_title("Raw-score discrimination")
    axes[0].legend()
    panel_label(axes[0], "A")

    deltas = combo[
        combo["summary_type"].eq("target_equal_macro")
        & combo["source_weighting"].eq("target_equal")
        & combo["comparison"].eq("ipTM_plus_ipSAE_vs_ipTM_only")
    ].set_index("layer")
    layers = ["U", "I", "F"]
    x2 = np.arange(len(layers))
    width = 0.32
    for offset, metric, label, color in zip(
        (-width / 2, width / 2),
        [
            "auc_combination_minus_comparator",
            "ap_combination_minus_comparator",
        ],
        ["Delta AUC", "Delta AP"],
        [COLORS["navy"], COLORS["orange"]],
    ):
        axes[1].bar(x2 + offset, deltas.loc[layers, metric], width=width, color=color, label=label)
    axes[1].axhline(0, color="#444444", linewidth=0.8)
    axes[1].set_xticks(x2, ["U", "I", "F"])
    axes[1].set_title("Discrimination difference")
    axes[1].set_ylabel("ipTM+ipSAE minus ipTM-only")
    axes[1].legend(fontsize=7)
    panel_label(axes[1], "B")

    for offset, metric, label, color in zip(
        (-width / 2, width / 2),
        [
            "brier_combination_minus_comparator",
            "log_loss_combination_minus_comparator",
        ],
        ["Delta Brier", "Delta log loss"],
        [COLORS["teal"], COLORS["red"]],
    ):
        axes[2].bar(x2 + offset, deltas.loc[layers, metric], width=width, color=color, label=label)
    axes[2].axhline(0, color="#444444", linewidth=0.8)
    axes[2].set_xticks(x2, ["U", "I", "F"])
    axes[2].set_title("Probability-error difference")
    axes[2].set_ylabel("ipTM+ipSAE minus ipTM-only (lower is better)")
    axes[2].legend(fontsize=7)
    panel_label(axes[2], "C")

    rel = intervals[intervals["metric"].eq("relative_macro_U_to_I_brier_reduction")].copy()
    rel["label"] = rel["dataset"].map(
        {
            "main_primary_six": "Main six",
            "main_primary_five": "Main five",
            "external_primary_weak_plus": "External eight",
        }
    )
    rel["plot_order"] = rel["dataset"].map(
        {
            "main_primary_six": 0,
            "main_primary_five": 1,
            "external_primary_weak_plus": 2,
        }
    )
    rel = rel.sort_values("plot_order").reset_index(drop=True)
    y = np.arange(len(rel))
    axes[3].errorbar(
        rel["point_estimate"],
        y,
        xerr=np.vstack(
            [
                rel["point_estimate"] - rel["bootstrap_95_low"],
                rel["bootstrap_95_high"] - rel["point_estimate"],
            ]
        ),
        fmt="o",
        capsize=3,
        color=COLORS["teal"],
        ecolor=COLORS["teal"],
    )
    axes[3].set_yticks(y, rel["label"])
    axes[3].invert_yaxis()
    axes[3].set_xlim(0.1, 0.55)
    axes[3].set_xlabel("Relative macro-Brier reduction, U to I")
    axes[3].set_title("Paired cluster bootstrap")
    panel_label(axes[3], "D")
    fig.suptitle("ipSAE and ipTM+ipSAE did not show a stable increment over ipTM", fontsize=13.5, fontweight="bold", y=1.02)
    return save(fig, "Supplementary_Figure_S11_ipSAE_postreview", SUPP_FIG)


def main() -> None:
    style()
    source_paths = {
        "pair_master_v3": ROOT / "analysis" / "v3_outputs" / "pair_master_v3.csv",
        "rq12_target_metrics": ROOT / "analysis" / "rq12_outputs" / "rq12_target_metrics.csv",
        "rq12_reliability_bins": ROOT / "analysis" / "rq12_outputs" / "rq2_reliability_bins.csv",
        "rq12_hierarchical_targets": ROOT / "analysis" / "rq12_outputs" / "rq2_hierarchical_targets.csv",
        "rq3_target_adaptation_metrics": ROOT / "analysis" / "codex_outputs" / "RQ3_F_correction_2026-11-13" / "run_20261113_182311" / "rq3_target_adaptation_metrics.csv",
        "rq3e_target_metrics": ROOT / "analysis" / "rq3e_if_outputs" / "run_20261114_174403" / "target_metrics.csv",
        "rq4_target_metrics": ROOT / "analysis" / "rq4_outputs" / "run_20261114_231925" / "rq4ab_target_metrics.csv",
        "rq4_comparisons": ROOT / "analysis" / "rq4_outputs" / "run_20261114_231925" / "rq4c_comparisons.csv",
        "rq5_score_macro": ROOT / "analysis" / "codex_outputs" / "RQ5_execution_2026-11-15" / "run_20261115_111720" / "rq5_score_macro_summary.csv",
        "rq5_pairwise": ROOT / "analysis" / "codex_outputs" / "RQ5_execution_2026-11-15" / "run_20261115_111720" / "rq5_pairwise_comparisons.csv",
        "rq6_target_equal": ROOT / "analysis" / "codex_outputs" / "RQ6_execution_2026-11-15" / "run_20261115_151147" / "rq6_target_equal_summary.csv",
        "ipSAE_macro": ROOT / "analysis" / "postreview_revision_2026-11-21" / "ipSAE" / "ipSAE_macro_metrics.csv",
        "ipSAE_combo": ROOT / "analysis" / "postreview_revision_2026-11-21" / "ipSAE_combo" / "combo_paired_comparisons.csv",
        "bootstrap_macro": ROOT / "analysis" / "postreview_revision_2026-11-21" / "bootstrap" / "bootstrap_macro_intervals.csv",
    }
    outputs = {
        "Figure 1": figure1(),
        "Figure 2": figure2(),
        "Figure 3": figure3(),
        "Figure 4": figure4(),
        "Figure 5": figure5(),
        "Figure 6": figure6(),
        "Figure 7": figure7(),
        "Supplementary Figure S11": supplementary_s11(),
    }
    manifest = {
        "status": "complete",
        "figures": outputs,
        "source_inputs": {
            name: {
                "path": str(path.relative_to(ROOT)),
                "sha256": sha256(path),
            }
            for name, path in source_paths.items()
        },
    }
    (OUT / "figure_manifest_postreview_2026-11-21.json").write_text(
        json.dumps(manifest, indent=2, ensure_ascii=False) + "\n",
        encoding="utf-8",
    )


if __name__ == "__main__":
    main()
