"""Final immutable-style hash manifest after all execution logs are closed."""
from __future__ import annotations

import hashlib
import json
import sys

sys.dont_write_bytecode = True

from pathlib import Path

import pandas as pd

from rq3_common import (
    RUN_DIR,
    atomic_csv,
    atomic_json,
    load_config,
    now_iso,
    require_hash_gate,
    sha256_file,
)


def main() -> None:
    config = load_config()
    final_inputs = require_hash_gate(config)
    atomic_csv(RUN_DIR / "final_input_hashes.csv", final_inputs)
    gate = {
        "checked_at": now_iso(),
        "all_input_hashes_match": bool(final_inputs["sha256_matches"].all()),
        "input_count": int(len(final_inputs)),
        "mismatches": final_inputs.loc[
            ~final_inputs["sha256_matches"], "input_key"
        ].astype(str).tolist(),
    }
    atomic_json(RUN_DIR / "final_hash_gate_result.json", gate)
    if not gate["all_input_hashes_match"]:
        raise RuntimeError("Final frozen-input hash gate failed")

    handoff_path = RUN_DIR / "handoff_result.json"
    handoff = json.loads(handoff_path.read_text(encoding="utf-8"))
    handoff["final_hash_manifest_pending_separate_finalizer"] = False
    handoff["final_hash_manifest"] = {
        "path": "final_hash_manifest.csv",
        "manifest_sha256_sidecar": "manifest_sha256.txt",
        "excludes_itself_and_sidecar": True,
    }
    atomic_json(handoff_path, handoff)

    excluded = {
        "final_hash_manifest.csv",
        "manifest_sha256.txt",
    }
    rows = []
    for path in sorted(RUN_DIR.rglob("*")):
        if not path.is_file() or path.name in excluded:
            continue
        rows.append(
            {
                "relative_path": str(path.relative_to(RUN_DIR)),
                "bytes": int(path.stat().st_size),
                "sha256": sha256_file(path),
            }
        )
    manifest = pd.DataFrame(rows)
    atomic_csv(RUN_DIR / "final_hash_manifest.csv", manifest)
    manifest_hash = sha256_file(RUN_DIR / "final_hash_manifest.csv")
    (RUN_DIR / "manifest_sha256.txt").write_text(
        f"{manifest_hash}  final_hash_manifest.csv\n",
        encoding="utf-8",
        newline="\n",
    )
    print(
        f"Final hash manifest complete: {len(manifest)} files; "
        f"manifest_sha256={manifest_hash}"
    )


if __name__ == "__main__":
    main()
