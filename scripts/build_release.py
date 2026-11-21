"""Build a clean publication candidate from an explicit allowlist."""
from __future__ import annotations

import argparse
import csv
import hashlib
import re
import shutil
from pathlib import Path


TEXT_SUFFIXES = {".py", ".md", ".txt", ".json", ".csv", ".yaml", ".yml", ".toml", ".cff"}
ROOT_FILES = [
    ".gitattributes",
    ".gitignore",
    ".release-root",
    "CITATION.cff",
    "DATA_LICENSE.md",
    "LICENSE",
    "README.md",
    "THIRD_PARTY_NOTICES.md",
    "ZENODO_METADATA.json",
    "Proteinbase校准研究_实验流程_v3_2026-11-06.md",
    "pyproject.toml",
    "requirements.txt",
]
PUBLICATION_FILES = [
    "publication/RELEASE_ASSETS.md",
    "publication/manuscript/main_article.md",
    "publication/manuscript/supplementary_materials.md",
    "publication/manuscript/build_review_pdf.py",
    "publication/manuscript/qc_review_pdf.py",
    "publication/manuscript/review_pdf_build.json",
    "publication/manuscript/review_pdf_qc_results.json",
    "publication/manuscript/supplement_source_manifest.md",
    "publication/manuscript/supplement_build_manifest.json",
]
DATA_FILES = [
    "data/postreview_dependency_manifest.json",
    "data/release_assets_manifest.json",
]
ANALYSIS_FILES = [
    "analysis/audit_outputs/pair_audit.csv",
    "analysis/pinc_batch/run_20261110_181500/pinc_pair_extension.csv",
    "analysis/rq12_outputs/sequence_clusters_rq12.csv",
    "analysis/rq12_outputs/rq12_target_metrics.csv",
    "analysis/rq12_outputs/rq2_reliability_bins.csv",
    "analysis/rq12_outputs/rq2_hierarchical_targets.csv",
    "analysis/codex_outputs/RQ3_prep_2026-11-13/run_20261113_105301/loto_pair_assignments.csv",
    "analysis/codex_outputs/RQ3_prep_2026-11-13/run_20261113_105301/loto_fold_manifest.csv",
    "analysis/codex_outputs/RQ3_prep_2026-11-13/run_20261113_105301/target_cv5_pair_assignments.csv",
    "analysis/codex_outputs/RQ3_prep_2026-11-13/run_20261113_105301/target_cv5_fold_manifest.csv",
    "analysis/codex_outputs/RQ3_prep_2026-11-13/run_20261113_105301/cross_target_cluster_exclusions.csv",
    "analysis/codex_outputs/RQ3_prep_2026-11-13/run_20261113_105301/fold_qc.json",
    "analysis/codex_outputs/RQ3_prep_2026-11-13/run_20261113_105301/phase_ab.py",
    "analysis/codex_outputs/RQ3_F_correction_2026-11-13/run_20261113_182311/rq3_target_adaptation_metrics.csv",
    "analysis/codex_outputs/RQ3_F_correction_2026-11-13/run_20261113_182311/rq3_adaptation_pair_predictions.csv",
    "analysis/codex_outputs/RQ3_F_correction_2026-11-13/run_20261113_182311/rq3_common.py",
    "analysis/rq3e_if_outputs/run_20261114_174403/target_metrics.csv",
    "analysis/rq3e_if_outputs/run_20261114_174403/oof_predictions.csv",
    "analysis/rq4_outputs/run_20261114_231925/rq4ab_target_metrics.csv",
    "analysis/rq4_outputs/run_20261114_231925/rq4c_comparisons.csv",
    "analysis/codex_outputs/RQ5_execution_2026-11-15/run_20261115_111720/rq5_score_macro_summary.csv",
    "analysis/codex_outputs/RQ5_execution_2026-11-15/run_20261115_111720/rq5_pairwise_comparisons.csv",
    "analysis/codex_outputs/RQ6_execution_2026-11-15/run_20261115_151147/rq6_target_equal_summary.csv",
    "analysis/v3_outputs/pair_master_v3.csv",
]
POSTREVIEW_ROOT_FILES = [
    "analysis/postreview_revision_2026-11-21/C01_C02_RESULTS.md",
    "analysis/postreview_revision_2026-11-21/C03_EVIDENCE_LOCK_2026-11-21.md",
    "analysis/postreview_revision_2026-11-21/RUN_STATUS.json",
    "analysis/postreview_revision_2026-11-21/POSTREVIEW_DELIVERY_QC_2026-11-21.md",
    "analysis/postreview_revision_2026-11-21/POSTREVIEW_DELIVERY_QC_2026-11-21.json",
]
POSTREVIEW_DIRS = {
    "ipSAE": {".csv", ".json", ".md"},
    "ipSAE_combo": {".csv", ".json", ".md"},
    "bootstrap": {".csv", ".json", ".md"},
    "figures": {".json"},
    "scripts": {".py"},
}
MAIN_FIGURE_SUFFIXES = {".png", ".svg", ".pdf"}
ANALYSIS_SUFFIXES = {".py", ".md", ".json", ".csv", ".txt", ".png"}
ANALYSIS_EXCLUDED_PARTS = {"__pycache__", "sample_qc", "vendor", "batch_work", "sources"}


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        for chunk in iter(lambda: stream.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def sanitize(text: str) -> str:
    text = re.sub(r"(?i)\b[A-Z]:[\\/]+[^\s\"']+", "<ABSOLUTE_PATH_REDACTED>", text)
    text = re.sub(r"(?i)\b[A-Z]:[\\/]+diyipian", "<REPRODUCTION_ROOT>", text)
    text = re.sub(r"(?i)\b[A-Z]:[\\/]+Users[\\/]+[^\\/\\\s\"']+", "<LOCAL_USER>", text)
    text = text.replace("/home/", "<HOME>/").replace("/Users/", "<USERS>/")
    return text


def source_files(source_root: Path) -> list[Path]:
    selected: set[Path] = set()
    for rel in ROOT_FILES + PUBLICATION_FILES + DATA_FILES + ANALYSIS_FILES + POSTREVIEW_ROOT_FILES:
        path = source_root / rel
        if path.is_file():
            selected.add(path)
    analysis_root = source_root / "analysis"
    if analysis_root.is_dir():
        for path in analysis_root.rglob("*"):
            if not path.is_file() or path.suffix.lower() not in ANALYSIS_SUFFIXES:
                continue
            relative_parts = path.relative_to(source_root).parts
            if any(part in ANALYSIS_EXCLUDED_PARTS for part in relative_parts):
                continue
            if path.name == "rq4c_pilot_parameters.csv":
                continue
            selected.add(path)
    for rel_root, suffixes in POSTREVIEW_DIRS.items():
        directory = source_root / "analysis" / "postreview_revision_2026-11-21" / rel_root
        if directory.is_dir():
            selected.update(path for path in directory.iterdir() if path.is_file() and path.suffix.lower() in suffixes)
    docs = source_root / "docs"
    if docs.is_dir():
        selected.update(docs.glob("*.md"))
    scripts = source_root / "scripts"
    if scripts.is_dir():
        selected.update(path for path in scripts.glob("*.py") if path.is_file())
    main_figures = source_root / "publication" / "figures" / "main"
    if main_figures.is_dir():
        for path in main_figures.iterdir():
            if path.is_file() and "_postreview" in path.stem and path.suffix.lower() in MAIN_FIGURE_SUFFIXES:
                selected.add(path)
    supplementary_figures = source_root / "publication" / "figures" / "supplementary"
    if supplementary_figures.is_dir():
        for path in supplementary_figures.iterdir():
            name = path.name
            if re.fullmatch(r"Supplementary_Figure_S(?:[1-9]|10)\.png", name):
                selected.add(path)
            elif re.fullmatch(r"Supplementary_Figure_S11_ipSAE_postreview\.(?:png|svg|pdf)", name):
                selected.add(path)
    return sorted(selected)


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--source-root", type=Path, required=True)
    parser.add_argument("--candidate-root", type=Path, required=True)
    args = parser.parse_args()
    source_root = args.source_root.resolve()
    candidate_root = args.candidate_root.resolve()
    if candidate_root.exists() and any(candidate_root.iterdir()):
        raise SystemExit(f"Candidate root must be absent or empty: {candidate_root}")
    candidate_root.mkdir(parents=True, exist_ok=True)

    rows: list[dict[str, object]] = []
    for source in source_files(source_root):
        rel = source.relative_to(source_root)
        target = candidate_root / rel
        target.parent.mkdir(parents=True, exist_ok=True)
        shutil.copy2(source, target)
        transformed = False
        non_sanitized_tools = {
            (source_root / "scripts" / "build_release.py").resolve(),
            (source_root / "scripts" / "verify_publication.py").resolve(),
        }
        if source.suffix.lower() in TEXT_SUFFIXES and source.resolve() not in non_sanitized_tools:
            original = target.read_text(encoding="utf-8", errors="replace")
            cleaned = sanitize(original)
            if cleaned != original:
                target.write_text(cleaned, encoding="utf-8")
                transformed = True
        rows.append(
            {
                "source_path": rel.as_posix(),
                "published_path": rel.as_posix(),
                "source_sha256": sha256(source),
                "published_sha256": sha256(target),
                "published_bytes": target.stat().st_size,
                "category": "figure" if source.suffix.lower() in MAIN_FIGURE_SUFFIXES else "code_or_documentation",
                "transform": "text-path-redacted" if transformed else "verbatim",
            }
        )

    manifest = candidate_root / "data" / "PUBLICATION_MANIFEST.csv"
    manifest.parent.mkdir(parents=True, exist_ok=True)
    with manifest.open("w", encoding="utf-8", newline="") as stream:
        writer = csv.DictWriter(
            stream,
            fieldnames=[
                "source_path",
                "published_path",
                "source_sha256",
                "published_sha256",
                "published_bytes",
                "category",
                "transform",
            ],
        )
        writer.writeheader()
        writer.writerows(sorted(rows, key=lambda row: str(row["published_path"])))

    published = {
        path.relative_to(candidate_root).as_posix()
        for path in candidate_root.rglob("*")
        if path.is_file()
    }
    expected = {str(row["published_path"]) for row in rows} | {"data/PUBLICATION_MANIFEST.csv"}
    if published != expected:
        raise SystemExit(
            "Reverse coverage mismatch: "
            f"missing={sorted(expected - published)}, extra={sorted(published - expected)}"
        )
    print(f"candidate: {candidate_root}")
    print(f"files: {len(rows) + 1}")


if __name__ == "__main__":
    main()
