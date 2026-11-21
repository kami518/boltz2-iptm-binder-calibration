"""Preflight or execute the frozen derived reproduction chain.

The restored legacy RQ3/RQ4 scripts are kept as accepted historical evidence,
but are not called by this entry point. Use their COMMANDS.md records when
reconstructing those earlier analyses from the documented upstream inputs.
"""
from __future__ import annotations

import argparse
import subprocess
import sys
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
REQUIRED = [
    "analysis/v3_outputs/pair_master_v3.csv",
    "analysis/audit_outputs/pair_audit.csv",
    "analysis/pinc_batch/run_20261110_181500/pinc_pair_extension.csv",
    "analysis/rq12_outputs/sequence_clusters_rq12.csv",
    "analysis/codex_outputs/RQ3_prep_2026-11-13/run_20261113_105301/loto_pair_assignments.csv",
    "analysis/codex_outputs/RQ3_prep_2026-11-13/run_20261113_105301/target_cv5_pair_assignments.csv",
    "analysis/codex_outputs/RQ3_F_correction_2026-11-13/run_20261113_182311/rq3_adaptation_pair_predictions.csv",
    "analysis/rq3e_if_outputs/run_20261114_174403/oof_predictions.csv",
    "analysis/postreview_revision_2026-11-21/ipSAE/ipSAE_pair_extension.csv",
    "analysis/postreview_revision_2026-11-21/ipSAE/ipSAE_target_metrics.csv",
    "publication/manuscript/main_article.md",
    "publication/manuscript/supplementary_materials.md",
]
COMMANDS = [
    ["analysis/postreview_revision_2026-11-21/scripts/run_c01b_postreview.py"],
    ["analysis/postreview_revision_2026-11-21/scripts/run_c01_c02_postreview.py", "--bootstrap-only"],
    ["analysis/postreview_revision_2026-11-21/scripts/make_postreview_figures.py"],
    ["scripts/redraw_supplementary_figure_s7.py"],
    ["publication/manuscript/build_review_pdf.py"],
    ["publication/manuscript/qc_review_pdf.py"],
    ["analysis/postreview_revision_2026-11-21/scripts/qc_postreview_delivery.py"],
]


def check_inputs() -> list[str]:
    return [relative for relative in REQUIRED if not (ROOT / relative).is_file()]


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--check", action="store_true", help="Check inputs without running analyses.")
    parser.add_argument("--execute", action="store_true", help="Run the frozen derived reproduction chain.")
    args = parser.parse_args()
    if not args.check and not args.execute:
        args.check = True
    missing = check_inputs()
    if missing:
        raise SystemExit("Missing frozen inputs: " + ", ".join(missing))
    print(f"frozen derived workflow inputs verified: {len(REQUIRED)}")
    if not args.execute:
        print("preflight only; rerun with --execute to run C01B, C02, figures, PDF, and QC")
        return
    for command in COMMANDS:
        path = ROOT / command[0]
        if not path.is_file():
            raise SystemExit(f"Missing reproduction script: {command[0]}")
        print("RUN", " ".join(command), flush=True)
        subprocess.run([sys.executable, str(path), *command[1:]], cwd=ROOT, check=True)
    print("frozen derived workflow complete")


if __name__ == "__main__":
    main()
