"""Finalize RQ4 hashes after report-only clarification and accepted summaries."""
from pathlib import Path
import hashlib, json, shutil
import pandas as pd
ROOT = next(_p for _p in Path(__file__).resolve().parents if (_p / ".release-root").is_file())
RUN = ROOT / "analysis/rq4_outputs/run_20261114_231925"
CLAR = ROOT / "analysis/RQ4_protocol_clarification_2026-11-14.md"

def sha(p):
    h = hashlib.sha256()
    with p.open("rb") as f:
        for b in iter(lambda: f.read(1024*1024), b""): h.update(b)
    return h.hexdigest()

old_manifest = sha(RUN / "final_hash_manifest.csv")
shutil.copy2(CLAR, RUN / CLAR.name)
note = {
    "reason": "Report-only finalization: clarify empty-selection recall=0 and append accepted sensitivity/verification summaries. No source data, selected sets, hit counts, model parameters, intervals, figures, or executed-code snapshot changed.",
    "previous_manifest_sha256": old_manifest,
    "run_status": json.loads((RUN / "RUN_STATUS.json").read_text(encoding="utf-8"))["status"],
}
(RUN / "finalization_note.json").write_text(json.dumps(note, ensure_ascii=False, indent=2)+"\n", encoding="utf-8")
files = sorted(p for p in RUN.iterdir() if p.is_file() and p.name not in {"final_hash_manifest.csv", "manifest_sha256.txt"})
manifest = pd.DataFrame({"file": [p.name for p in files], "bytes": [p.stat().st_size for p in files], "sha256": [sha(p) for p in files]})
manifest.to_csv(RUN / "final_hash_manifest.csv", index=False)
new_hash = sha(RUN / "final_hash_manifest.csv")
(RUN / "manifest_sha256.txt").write_text(new_hash+"\n", encoding="utf-8")
print(json.dumps({"files": len(files), "manifest_sha256": new_hash}, indent=2))
