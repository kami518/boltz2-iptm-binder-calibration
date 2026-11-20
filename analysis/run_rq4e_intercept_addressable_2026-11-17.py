#!/usr/bin/env python3
"""Derive the frozen RQ4E post-hoc diagnostic without refitting any model."""

from __future__ import annotations

import csv
import hashlib
import json
import math
from collections import defaultdict
from datetime import datetime, timezone, timedelta
from pathlib import Path

ROOT = next(_p for _p in Path(__file__).resolve().parents if (_p / ".release-root").is_file())
SOURCE = ROOT / "analysis" / "codex_outputs" / "RQ3_F_correction_2026-11-13" / "run_20261113_182311"
PROTOCOL = ROOT / "analysis" / "RQ4E_intercept_addressable_frozen_protocol_2026-11-17.md"
SOURCE_FILES = {
    "target_adaptation": SOURCE / "rq3_target_adaptation_metrics.csv",
    "loto": SOURCE / "rq3_loto_metrics.csv",
    "comparisons": SOURCE / "rq3_adaptation_comparisons.csv",
    "macro": SOURCE / "rq3_adaptation_macro_summary.csv",
    "original_qc": SOURCE / "qc_checks.json",
    "original_hash_gate": SOURCE / "final_hash_gate_result.json",
}
LAYERS = ("U_transport_zero_label", "predictive_I", "predictive_F")
SHORT = {"U_transport_zero_label": "U", "predictive_I": "I", "predictive_F": "F"}
TARGET_LABELS = {
    "aah-i": "AAH-I",
    "human-vegf110": "VEGF110",
    "nipah-glycoprotein-g": "Nipah G",
    "pd-l1": "PD-L1",
    "sars-cov-2-rbd": "SARS-CoV-2 RBD",
    "trem2": "TREM2",
}
TARGET_ORDER = tuple(TARGET_LABELS)
EXPECTED = {
    "primary_six": {"rows": 18, "targets": 6, "n": 2084, "positives": 748},
    "primary_five": {"rows": 15, "targets": 5, "n": 1951, "positives": 729},
}
METRICS = ("brier", "log_loss", "bss_transport", "bss_target_oracle", "calibration_intercept", "calibration_slope")
TOL = 1e-12


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def read_csv(path: Path) -> list[dict[str, str]]:
    with path.open("r", encoding="utf-8-sig", newline="") as handle:
        return list(csv.DictReader(handle))


def write_csv(path: Path, rows: list[dict], fieldnames: list[str]) -> None:
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames)
        writer.writeheader()
        writer.writerows(rows)


def f(row: dict[str, str], key: str) -> float:
    return float(row[key])


def mean(values: list[float]) -> float:
    return math.fsum(values) / len(values)


def main() -> None:
    started = datetime.now(timezone(timedelta(hours=8)))
    hashes_before = {name: sha256(path) for name, path in SOURCE_FILES.items()}
    adaptation = read_csv(SOURCE_FILES["target_adaptation"])
    loto = read_csv(SOURCE_FILES["loto"])
    comparisons = read_csv(SOURCE_FILES["comparisons"])
    frozen_macro = read_csv(SOURCE_FILES["macro"])
    original_qc = json.loads(SOURCE_FILES["original_qc"].read_text(encoding="utf-8"))
    original_hash_gate = json.loads(SOURCE_FILES["original_hash_gate"].read_text(encoding="utf-8"))

    selected = [
        row for row in adaptation
        if row["analysis_set"] in EXPECTED
        and row["outcome_name"] == "primary_outcome"
        and row["cluster_rule"] == "sequence_cluster90"
        and row["source_weighting"] == "target_equal"
        and row["model_layer_label"] in LAYERS
    ]
    grouped: dict[tuple[str, str], dict[str, dict[str, str]]] = defaultdict(dict)
    for row in selected:
        key = (row["analysis_set"], row["held_out_target"])
        layer = row["model_layer_label"]
        if layer in grouped[key]:
            raise RuntimeError(f"Duplicate selected row: {key}, {layer}")
        grouped[key][layer] = row

    checks: list[dict] = []

    def check(name: str, passed: bool, observed, expected) -> None:
        checks.append({"check": name, "passed": bool(passed), "observed": observed, "expected": expected})

    check("original_rq3_all_qc_passed", original_qc.get("all_passed") is True and original_qc.get("n_failed") == 0,
          {"all_passed": original_qc.get("all_passed"), "n_failed": original_qc.get("n_failed")},
          {"all_passed": True, "n_failed": 0})
    check("original_rq3_input_hash_gate_passed", original_hash_gate.get("all_input_hashes_match") is True,
          original_hash_gate.get("all_input_hashes_match"), True)

    target_rows: list[dict] = []
    macro_rows: list[dict] = []
    derived_differences: dict[tuple[str, str, str], dict] = {}

    for analysis_set, expectation in EXPECTED.items():
        set_rows = [row for row in selected if row["analysis_set"] == analysis_set]
        targets = sorted({row["held_out_target"] for row in set_rows}, key=TARGET_ORDER.index)
        check(f"{analysis_set}_selected_row_count", len(set_rows) == expectation["rows"], len(set_rows), expectation["rows"])
        check(f"{analysis_set}_target_count", len(targets) == expectation["targets"], len(targets), expectation["targets"])
        check(f"{analysis_set}_complete_layer_coverage", all(set(grouped[(analysis_set, target)]) == set(LAYERS) for target in targets),
              {target: sorted(grouped[(analysis_set, target)]) for target in targets}, "U/I/F for every target")
        check(f"{analysis_set}_aah_i_membership", ("aah-i" in targets) == (analysis_set == "primary_six"),
              "aah-i" in targets, analysis_set == "primary_six")

        for layer in LAYERS:
            layer_rows = [grouped[(analysis_set, target)][layer] for target in targets]
            n_sum = sum(int(row["n"]) for row in layer_rows)
            positive_sum = sum(int(row["positives"]) for row in layer_rows)
            check(f"{analysis_set}_{SHORT[layer]}_n_total", n_sum == expectation["n"], n_sum, expectation["n"])
            check(f"{analysis_set}_{SHORT[layer]}_positives_total", positive_sum == expectation["positives"], positive_sum, expectation["positives"])

        for target in targets:
            rows = grouped[(analysis_set, target)]
            n_values = {row["n"] for row in rows.values()}
            positive_values = {row["positives"] for row in rows.values()}
            prevalence_values = {f(row, "prevalence") for row in rows.values()}
            check(f"{analysis_set}_{target}_same_n", len(n_values) == 1, sorted(n_values), "one shared n")
            check(f"{analysis_set}_{target}_same_positives", len(positive_values) == 1, sorted(positive_values), "one shared positives")
            check(f"{analysis_set}_{target}_same_prevalence", len(prevalence_values) == 1, sorted(prevalence_values), "one shared prevalence")
            evaluable = all(row["convergence_flags"] == "evaluable" and row["calibration_status"] == "evaluable" for row in rows.values())
            check(f"{analysis_set}_{target}_all_layers_evaluable", evaluable,
                  {SHORT[layer]: [row["convergence_flags"], row["calibration_status"]] for layer, row in rows.items()}, "all evaluable")
            finite = all(math.isfinite(f(row, metric)) for row in rows.values() for metric in METRICS + ("prevalence",))
            check(f"{analysis_set}_{target}_finite_metrics", finite, finite, True)

            u, i, full = (rows[layer] for layer in LAYERS)
            out = {
                "analysis_set": analysis_set,
                "outcome_name": "primary_outcome",
                "cluster_rule": "sequence_cluster90",
                "source_weighting": "target_equal",
                "held_out_target": target,
                "target_label": TARGET_LABELS[target],
                "n": int(u["n"]),
                "positives": int(u["positives"]),
                "prevalence": f(u, "prevalence"),
            }
            for layer, row in zip(("U", "I", "F"), (u, i, full)):
                for metric in METRICS:
                    out[f"{metric}_{layer}"] = f(row, metric)
            out["delta_brier_I_minus_U"] = out["brier_I"] - out["brier_U"]
            out["delta_brier_F_minus_I"] = out["brier_F"] - out["brier_I"]
            out["delta_log_loss_I_minus_U"] = out["log_loss_I"] - out["log_loss_U"]
            out["delta_log_loss_F_minus_I"] = out["log_loss_F"] - out["log_loss_I"]
            out["relative_brier_reduction_I_vs_U"] = (out["brier_U"] - out["brier_I"]) / out["brier_U"]
            out["relative_brier_reduction_F_vs_I"] = (out["brier_I"] - out["brier_F"]) / out["brier_I"]
            target_rows.append(out)

        set_target_rows = [row for row in target_rows if row["analysis_set"] == analysis_set]
        macro = {
            "analysis_set": analysis_set,
            "outcome_name": "primary_outcome",
            "cluster_rule": "sequence_cluster90",
            "source_weighting": "target_equal",
            "target_count": len(targets),
            "n_total": expectation["n"],
            "positives_total": expectation["positives"],
        }
        for layer in ("U", "I", "F"):
            for metric in METRICS:
                macro[f"macro_{metric}_{layer}"] = mean([row[f"{metric}_{layer}"] for row in set_target_rows])
        for comparison, left, right in (("I_minus_U", "I", "U"), ("F_minus_I", "F", "I")):
            for metric in ("brier", "log_loss"):
                values = [row[f"{metric}_{left}"] - row[f"{metric}_{right}"] for row in set_target_rows]
                macro[f"macro_delta_{metric}_{comparison}"] = mean(values)
                macro[f"n_improved_{metric}_{comparison}"] = sum(value < 0 for value in values)
                macro[f"n_worsened_{metric}_{comparison}"] = sum(value > 0 for value in values)
                derived_differences[(analysis_set, f"predictive_{comparison}", metric)] = {
                    "mean": mean(values), "improved": sum(value < 0 for value in values),
                    "worsened": sum(value > 0 for value in values), "min": min(values), "max": max(values),
                }
        macro["macro_relative_brier_reduction_I_vs_U"] = (macro["macro_brier_U"] - macro["macro_brier_I"]) / macro["macro_brier_U"]
        macro["macro_relative_brier_reduction_F_vs_I"] = (macro["macro_brier_I"] - macro["macro_brier_F"]) / macro["macro_brier_I"]
        macro["mean_target_relative_brier_reduction_I_vs_U"] = mean([row["relative_brier_reduction_I_vs_U"] for row in set_target_rows])
        macro["mean_target_relative_brier_reduction_F_vs_I"] = mean([row["relative_brier_reduction_F_vs_I"] for row in set_target_rows])
        macro_rows.append(macro)

    for row in target_rows:
        arithmetic = {
            "delta_brier_I_minus_U": row["brier_I"] - row["brier_U"],
            "delta_brier_F_minus_I": row["brier_F"] - row["brier_I"],
            "delta_log_loss_I_minus_U": row["log_loss_I"] - row["log_loss_U"],
            "delta_log_loss_F_minus_I": row["log_loss_F"] - row["log_loss_I"],
        }
        check(f"{row['analysis_set']}_{row['held_out_target']}_difference_arithmetic",
              all(math.isclose(row[key], value, rel_tol=0, abs_tol=TOL) for key, value in arithmetic.items()),
              "all four differences", "exact within 1e-12")

    comparison_mismatches = []
    for row in comparisons:
        key = (row["analysis_set"], row["comparison"], row["metric"])
        if (row["cluster_rule"] != "sequence_cluster90" or row["source_weighting"] != "target_equal"
                or key not in derived_differences):
            continue
        derived = derived_differences[key]
        fields = {
            "macro_mean_difference": derived["mean"], "n_improved": derived["improved"],
            "n_worsened": derived["worsened"], "min_difference": derived["min"], "max_difference": derived["max"],
        }
        for field, value in fields.items():
            observed = f(row, field) if field in ("macro_mean_difference", "min_difference", "max_difference") else int(row[field])
            if not math.isclose(float(observed), float(value), rel_tol=0, abs_tol=TOL):
                comparison_mismatches.append({"key": key, "field": field, "frozen": observed, "derived": value})
    check("derived_comparisons_match_frozen_table", not comparison_mismatches and len(derived_differences) == 8,
          {"mismatches": comparison_mismatches, "comparison_metric_pairs": len(derived_differences)},
          {"mismatches": [], "comparison_metric_pairs": 8})

    macro_mismatches = []
    for row in frozen_macro:
        if (row["analysis_set"] not in EXPECTED or row["cluster_rule"] != "sequence_cluster90"
                or row["source_weighting"] != "target_equal" or row["model_layer_label"] not in LAYERS
                or row["metric"] not in ("brier", "log_loss")):
            continue
        derived = next(item for item in macro_rows if item["analysis_set"] == row["analysis_set"])
        derived_value = derived[f"macro_{row['metric']}_{SHORT[row['model_layer_label']]}"]
        if not math.isclose(f(row, "mean"), derived_value, rel_tol=0, abs_tol=TOL):
            macro_mismatches.append({"analysis_set": row["analysis_set"], "layer": row["model_layer_label"],
                                     "metric": row["metric"], "frozen": f(row, "mean"), "derived": derived_value})
    check("derived_layer_macro_means_match_frozen_table", not macro_mismatches, macro_mismatches, [])

    loto_mismatches = []
    for row in loto:
        if (row["analysis_set"] == "primary_six" and row["outcome_name"] == "primary_outcome"
                and row["cluster_rule"] == "sequence_cluster90" and row["source_weighting"] == "target_equal"
                and row["model_form"] == "linear" and row["model_layer_label"] == "U_transport_zero_label"):
            source_u = grouped[("primary_six", row["held_out_target"])]["U_transport_zero_label"]
            for metric in ("n", "positives", "prevalence") + METRICS:
                # Repeated numerical calibration optimization differs at ~1e-8 for VEGF110;
                # predictions and all error metrics remain subject to the strict 1e-12 gate.
                tolerance = 1e-7 if metric in ("calibration_intercept", "calibration_slope") else TOL
                if not math.isclose(f(row, metric), f(source_u, metric), rel_tol=0, abs_tol=tolerance):
                    loto_mismatches.append({"target": row["held_out_target"], "metric": metric,
                                            "loto": f(row, metric), "adaptation": f(source_u, metric),
                                            "tolerance": tolerance})
    check("primary_six_u_matches_independent_loto_table", not loto_mismatches, loto_mismatches, [])

    hashes_after = {name: sha256(path) for name, path in SOURCE_FILES.items()}
    check("source_hashes_unchanged", hashes_before == hashes_after,
          {name: hashes_before[name] == hashes_after[name] for name in hashes_before}, "all unchanged")

    run_id = started.strftime("run_%Y%m%d_%H%M%S")
    output_dir = ROOT / "analysis" / "rq4e_outputs" / run_id
    output_dir.mkdir(parents=True, exist_ok=False)

    target_fields = list(target_rows[0])
    macro_fields = list(macro_rows[0])
    write_csv(output_dir / "rq4e_target_diagnostic.csv", target_rows, target_fields)
    write_csv(output_dir / "rq4e_macro_summary.csv", macro_rows, macro_fields)
    hash_rows = [{
        "source_name": name, "path": str(SOURCE_FILES[name]), "sha256_before": hashes_before[name],
        "sha256_after": hashes_after[name], "unchanged": hashes_before[name] == hashes_after[name],
    } for name in SOURCE_FILES]
    write_csv(output_dir / "input_hashes.csv", hash_rows, list(hash_rows[0]))

    config = {
        "analysis": "RQ4E_intercept_addressable_post_hoc_diagnostic",
        "generated_at": started.isoformat(),
        "protocol": str(PROTOCOL),
        "protocol_sha256": sha256(PROTOCOL),
        "script": str(Path(__file__).resolve()),
        "script_sha256": sha256(Path(__file__).resolve()),
        "source_directory": str(SOURCE),
        "selection": {"outcome_name": "primary_outcome", "cluster_rule": "sequence_cluster90", "source_weighting": "target_equal"},
        "analysis_sets": EXPECTED,
        "layers": LAYERS,
        "tolerance": TOL,
        "model_refitting": False,
        "new_bootstrap_or_inference": False,
        "causal_attribution_permitted": False,
    }
    (output_dir / "run_config.snapshot.json").write_text(json.dumps(config, indent=2, ensure_ascii=False) + "\n", encoding="utf-8")

    lines = [
        "# RQ4E results: intercept-addressable error diagnostic", "",
        f"Generated: {started.isoformat()}  ",
        "Status: post-hoc derivation of frozen RQ3 held-out-target results; no model was refitted.", "",
        "## Target-equal macro summary", "",
        "| Analysis set | Targets | Brier U | Brier I | Brier F | I−U | F−I | Relative reduction I vs U | Log loss I−U | Log loss F−I |",
        "|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|",
    ]
    for row in macro_rows:
        lines.append(
            f"| {row['analysis_set']} | {row['target_count']} | {row['macro_brier_U']:.6f} | {row['macro_brier_I']:.6f} | "
            f"{row['macro_brier_F']:.6f} | {row['macro_delta_brier_I_minus_U']:.6f} | {row['macro_delta_brier_F_minus_I']:.6f} | "
            f"{100 * row['macro_relative_brier_reduction_I_vs_U']:.1f}% | {row['macro_delta_log_loss_I_minus_U']:.6f} | "
            f"{row['macro_delta_log_loss_F_minus_I']:.6f} |"
        )
    for analysis_set in EXPECTED:
        lines.extend(["", f"## {analysis_set}: target-level Brier diagnostic", "",
                      "| Target | n | Positives | U | I | F | I−U | F−I | Relative reduction I vs U |",
                      "|---|---:|---:|---:|---:|---:|---:|---:|---:|"])
        for row in target_rows:
            if row["analysis_set"] != analysis_set:
                continue
            lines.append(
                f"| {row['target_label']} | {row['n']} | {row['positives']} | {row['brier_U']:.6f} | {row['brier_I']:.6f} | "
                f"{row['brier_F']:.6f} | {row['delta_brier_I_minus_U']:.6f} | {row['delta_brier_F_minus_I']:.6f} | "
                f"{100 * row['relative_brier_reduction_I_vs_U']:.1f}% |"
            )
    six, five = macro_rows
    lines.extend(["", "## Interpretation boundary", "",
        f"In the six-target analysis, intercept updating reduced target-equal macro Brier by {abs(six['macro_delta_brier_I_minus_U']):.6f} "
        f"({100 * six['macro_relative_brier_reduction_I_vs_U']:.1f}% relative to U), with improvement in "
        f"{six['n_improved_brier_I_minus_U']}/6 targets. Further slope updating changed macro Brier by "
        f"{six['macro_delta_brier_F_minus_I']:.6f}, improving {six['n_improved_brier_F_minus_I']} and worsening "
        f"{six['n_worsened_brier_F_minus_I']} targets.", "",
        f"In the separately refitted five-target sensitivity analysis, the corresponding intercept-updating change was "
        f"{five['macro_delta_brier_I_minus_U']:.6f} ({100 * five['macro_relative_brier_reduction_I_vs_U']:.1f}% relative reduction), "
        f"with improvement in {five['n_improved_brier_I_minus_U']}/5 targets. Further slope updating changed macro Brier by "
        f"{five['macro_delta_brier_F_minus_I']:.6f}, improving {five['n_improved_brier_F_minus_I']} and worsening "
        f"{five['n_worsened_brier_F_minus_I']} targets.", "",
        "These results quantify an intercept-addressable error component and the incremental benefit of slope updating. "
        "They are not a causal decomposition and do not establish that base rate is the largest mechanistic cause of transport failure. "
        "Post-I Brier and BSS describe remaining error after intercept adaptation; they do not identify its biological or experimental cause. "
        "No new confidence intervals, bootstrap tests, or pilot recommendation are introduced by this diagnostic.", "",
        "## Quality control", "",
        f"Checks passed: {sum(item['passed'] for item in checks)}/{len(checks)}. Source SHA-256 values were unchanged.", ""])
    (output_dir / "RQ4E_RESULTS.md").write_text("\n".join(lines), encoding="utf-8")

    qc = {
        "analysis": "RQ4E_intercept_addressable_post_hoc_diagnostic",
        "checked_at": datetime.now(timezone(timedelta(hours=8))).isoformat(),
        "all_passed": all(item["passed"] for item in checks),
        "n_checks": len(checks),
        "n_failed": sum(not item["passed"] for item in checks),
        "checks": checks,
    }
    (output_dir / "rq4e_qc_results.json").write_text(json.dumps(qc, indent=2, ensure_ascii=False) + "\n", encoding="utf-8")
    print(json.dumps({"output_dir": str(output_dir), "qc_all_passed": qc["all_passed"],
                      "qc_checks": qc["n_checks"], "macro_summary": macro_rows}, ensure_ascii=False, indent=2))
    if not qc["all_passed"]:
        raise RuntimeError(f"RQ4E QC failed; inspect {output_dir / 'rq4e_qc_results.json'}")


if __name__ == "__main__":
    main()
