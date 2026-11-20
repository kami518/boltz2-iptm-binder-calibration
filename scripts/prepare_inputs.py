from __future__ import annotations

import argparse
import json
import sys
import urllib.request
from pathlib import Path
from typing import Any

ROOT = next(path for path in Path(__file__).resolve().parents if (path / ".release-root").is_file())
ANTHROPIC_REVISION = "9e1b81696da46835e9e9cde9a3da976e0abc92ab"
PROTEINBASE_REVISION = "2d2825070da5ba77e0c7e1a0942146e8307e9282"
HF = "https://huggingface.co/datasets"

ANTHROPIC_FILES = {
    "data/tables/design_summary.parquet": "data/external_target_audit_2026-11-13/anthropic_design_summary.parquet",
    "data/tables/wetlab/summary.parquet": "data/external_target_audit_2026-11-13/anthropic_wetlab_summary.parquet",
    "data/tables/insilico/cofold_predictions.parquet": "data/external_target_audit_2026-11-13/anthropic_cofold_predictions.parquet",
}


def download(url: str, target: Path, *, force: bool) -> None:
    target.parent.mkdir(parents=True, exist_ok=True)
    if target.exists() and not force:
        print(f"present: {target.relative_to(ROOT)}")
        return
    temporary = target.with_suffix(target.suffix + ".part")
    print(f"download: {url}")
    request = urllib.request.Request(url, headers={"User-Agent": "boltz2-iptm-reproducibility/0.1"})
    with urllib.request.urlopen(request) as response, temporary.open("wb") as stream:
        while chunk := response.read(1024 * 1024):
            stream.write(chunk)
    temporary.replace(target)
    print(f"saved: {target.relative_to(ROOT)}")


def download_anthropic(force: bool) -> None:
    for upstream, relative in ANTHROPIC_FILES.items():
        url = f"{HF}/Anthropic/claude-protein-binder-design/resolve/{ANTHROPIC_REVISION}/{upstream}"
        download(url, ROOT / relative, force=force)


def _jsonable(value: Any) -> Any:
    if value is None or isinstance(value, (str, bool, int, float)):
        return value
    if isinstance(value, dict):
        return {str(k): _jsonable(v) for k, v in value.items()}
    if isinstance(value, (list, tuple)):
        return [_jsonable(v) for v in value]
    if hasattr(value, "item"):
        return value.item()
    return str(value)


def prepare_proteinbase(force: bool) -> None:
    parquet = ROOT / "data/proteinbase_latest.parquet"
    url = f"{HF}/adaptyvbio/proteinbase/resolve/{PROTEINBASE_REVISION}/data/proteins-00000-of-00001.parquet"
    download(url, parquet, force=force)
    output = ROOT / "data/proteinbase_latest.csv"
    if output.exists() and not force:
        print(f"present: {output.relative_to(ROOT)}")
        return
    try:
        import duckdb
    except ImportError as exc:
        raise SystemExit("DuckDB is required. Install dependencies from requirements.txt.") from exc
    frame = duckdb.connect().execute("SELECT * FROM read_parquet(?)", [str(parquet)]).df()
    frame["evaluations"] = frame["evaluations"].map(
        lambda values: json.dumps([_jsonable(dict(value)) for value in values], ensure_ascii=False, separators=(",", ":"))
    )
    frame.to_csv(output, index=False)
    print(f"saved: {output.relative_to(ROOT)}")
    print("For exact byte-level comparison, record the source CSV SHA-256 documented in docs/DATA_SOURCES.md.")


def main() -> None:
    parser = argparse.ArgumentParser(description="Download pinned public upstream tables.")
    parser.add_argument("--anthropic", action="store_true", help="download the three upstream Anthropic tables")
    parser.add_argument("--proteinbase-parquet", action="store_true", help="download and convert Proteinbase v2026-11-01")
    parser.add_argument("--force", action="store_true", help="replace existing local inputs")
    args = parser.parse_args()
    if not args.anthropic and not args.proteinbase_parquet:
        parser.error("select at least one of --anthropic or --proteinbase-parquet")
    if args.anthropic:
        download_anthropic(args.force)
    if args.proteinbase_parquet:
        prepare_proteinbase(args.force)


if __name__ == "__main__":
    main()
