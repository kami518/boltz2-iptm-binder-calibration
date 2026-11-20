"""RQ3 stage 0: frozen-input, fold, and environment QC.

This stage writes the run status and an exact configuration snapshot before
performing a complete structural audit of the accepted 44 LOTO and 220 CV5
fold assignments.
"""
from __future__ import annotations

import shutil
import sys

sys.dont_write_bytecode = True

import pandas as pd

from rq3_common import (
    CONFIG_PATH,
    FAILURE_PATH,
    RUN_DIR,
    atomic_csv,
    atomic_json,
    environment_snapshot,
    load_config,
    load_joined_data,
    now_iso,
    require_hash_gate,
    set_stage_status,
    sha256_file,
    update_status,
    validate_fold_structure,
)

STAGE = "stage_0_input_and_fold_qc"


def main() -> None:
    config = load_config()
    update_status(
        overall_status="running",
        current_stage=STAGE,
        started_at=now_iso(),
        config_frozen_before_any_rq3_result=bool(
            config["config_frozen_before_any_rq3_result"]
        ),
    )
    set_stage_status(STAGE, "running")
    snapshot_path = RUN_DIR / "run_config.snapshot.json"
    shutil.copyfile(CONFIG_PATH, snapshot_path)
    original_hash = sha256_file(CONFIG_PATH)
    snapshot_hash = sha256_file(snapshot_path)
    if snapshot_hash != original_hash:
        raise RuntimeError("Configuration snapshot hash differs from source")
    update_status(
        config_source_path=str(CONFIG_PATH),
        config_source_sha256=original_hash,
        config_snapshot_path=str(snapshot_path),
        config_snapshot_sha256=snapshot_hash,
    )

    initial_hashes = require_hash_gate(config)
    atomic_csv(RUN_DIR / "input_hash_initial.csv", initial_hashes)
    environment = environment_snapshot()
    atomic_json(RUN_DIR / "environment.json", environment)

    joined = load_joined_data(config)
    loto, loto_manifest, cv5, cv5_manifest, exclusions, qc = (
        validate_fold_structure(config, joined)
    )
    atomic_json(RUN_DIR / "stage0_qc.json", qc)

    fold_summary = {
        "loto_fold_count": int(loto_manifest["fold_id"].nunique()),
        "loto_assignment_rows": int(len(loto)),
        "loto_train_rows": int(loto["split_status"].eq("train").sum()),
        "loto_test_rows": int(loto["split_status"].eq("test").sum()),
        "loto_removed_rows": int(
            loto["split_status"].eq("removed_cross_target_cluster").sum()
        ),
        "cv5_fold_count": int(cv5_manifest["fold_id"].nunique()),
        "cv5_assignment_rows": int(len(cv5)),
        "cv5_train_rows": int(cv5["split"].eq("train").sum()),
        "cv5_test_rows": int(cv5["split"].eq("test").sum()),
        "cross_target_exclusions": int(len(exclusions)),
        "unique_pair_ids": int(joined["pair_id"].nunique()),
        "targets": {
            str(target): int(count)
            for target, count in joined.groupby("target").size().items()
        },
    }
    atomic_json(RUN_DIR / "stage0_fold_summary.json", fold_summary)
    if not FAILURE_PATH.exists():
        atomic_csv(
            FAILURE_PATH,
            pd.DataFrame(
                columns=[
                    "timestamp",
                    "stage",
                    "item_id",
                    "failure_type",
                    "detail",
                    "context_json",
                ]
            ),
        )
    set_stage_status(
        STAGE,
        "complete",
        checks_passed=qc["n_checks"],
        failures=0,
        initial_hash_gate_passed=True,
        config_snapshot_sha256=snapshot_hash,
    )
    print(
        f"Stage 0 complete: {qc['n_checks']}/{qc['n_checks']} structural checks "
        f"passed; {fold_summary['loto_fold_count']} LOTO folds and "
        f"{fold_summary['cv5_fold_count']} CV5 folds reused."
    )


if __name__ == "__main__":
    main()
