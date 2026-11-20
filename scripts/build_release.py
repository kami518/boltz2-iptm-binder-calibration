from __future__ import annotations

import csv
import hashlib
import re
import shutil
from pathlib import Path

ROOT = Path(r"E:/diyipian")
RELEASE = ROOT / "release/boltz2-iptm-binder-calibration"
MANIFEST = RELEASE / "data/PUBLICATION_MANIFEST.csv"


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        for chunk in iter(lambda: stream.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def relative(path: Path) -> str:
    return path.relative_to(ROOT).as_posix()


def public_relative(path: Path) -> str:
    return path.relative_to(RELEASE).as_posix()


def manifest_entries() -> list[str]:
    text = (ROOT / "manuscript/supplement/supplement_source_manifest_2026-11-16.md").read_text(encoding="utf-8")
    return re.findall(r"^\| `([^`]+)` \| \d+ \| `[0-9a-f]{64}` \|$", text, flags=re.MULTILINE)


def source_manifest_relative_paths() -> list[str]:
    paths = {path for path in manifest_entries() if path.startswith("analysis/")}
    accepted_roots = [
        "analysis/codex_outputs/RQ3_F_correction_2026-11-13/run_20261113_182311",
        "analysis/rq3e_outputs/run_20261114_164950",
        "analysis/rq3e_if_outputs/run_20261114_174403",
        "analysis/rq4_outputs/run_20261114_231925",
        "analysis/rq4e_outputs/run_20261117_162929",
        "analysis/codex_outputs/RQ5_execution_2026-11-15/run_20261115_111720",
        "analysis/codex_outputs/RQ6_execution_2026-11-15/run_20261115_151147",
        "analysis/codex_outputs/RQ5RF_model_complexity_2026-11-17/run_20261117_220855",
    ]
    for rel_root in accepted_roots:
        base = ROOT / rel_root
        for path in base.rglob("*"):
            if not path.is_file():
                continue
            if path.suffix.lower() in {".py", ".md", ".json"} or path.name in {
                "COMMANDS.md", "RUN_STATUS.json", "qc_checks.json", "final_hash_manifest.csv",
                "final_input_hashes.csv", "manifest_sha256.txt",
            }:
                paths.add(relative(path))
    return sorted(paths)


CODE_FILES = [
    "analysis/build_pair_master_v3.py",
    "analysis/audit_proteinbase.py",
    "analysis/audit_external_proteinbase_2026-11-13.py",
    "analysis/audit_anthropic_external_2026-11-13.py",
    "analysis/audit_anthropic_sequence_similarity_2026-11-13.py",
    "analysis/build_sequence_clusters_rq12_2026-11-09.py",
    "analysis/check_rq12_adaptive_modes.py",
    "analysis/check_rq12_quadrature_convergence.py",
    "analysis/debug_rq12_hierarchy.py",
    "analysis/rq12_adaptive.py",
    "analysis/rq12_core.py",
    "analysis/run_rq12_analysis_2026-11-09.py",
    "analysis/setup_mmseqs_rq12_2026-11-09.py",
    "analysis/summarize_rq12_results_2026-11-09.py",
    "analysis/test_rq12_core.py",
    "analysis/trace_rq12_adaptive_mode.py",
    "analysis/validate_rq3e_low_event_stress_2026-11-14.py",
    "analysis/verify_rq4_decision_value_2026-11-14.py",
    "analysis/finalize_rq4_run_2026-11-14.py",
    "analysis/run_rq3e_confirmatory_720_2026-11-14.py",
    "analysis/run_rq3e_if_external_adaptation_2026-11-14.py",
    "analysis/run_rq3e_low_event_stress_2026-11-14.py",
    "analysis/run_rq4_decision_value_2026-11-14.py",
    "analysis/run_rq4e_intercept_addressable_2026-11-17.py",
    "analysis/RQ3E_external_transport_frozen_protocol_2026-11-14.md",
    "analysis/RQ4E_intercept_addressable_frozen_protocol_2026-11-17.md",
    "analysis/RQ5RF_model_complexity_frozen_protocol_2026-11-17.md",
    "manuscript/assemble_figures_2026-11-15.py",
    "manuscript/build_review_pdf_v2_2026-11-18.py",
    "manuscript/build_supplement_2026-11-16.py",
    "manuscript/integrate_rq5rf_manuscript_2026-11-18.py",
    "manuscript/qc_manuscript_v2_2026-11-17.py",
    "manuscript/qc_rq5rf_integration_2026-11-18.py",
    "manuscript/supplement/build_supplement_s16_2026-11-18.py",
    "manuscript/review/qc_review_pdf_v2_2026-11-18.py",
    "Proteinbase校准研究_实验流程_v3_2026-11-06.md",
]

EXACT_FILES = [
    ("manuscript/build_public_release_2026-11-19.py", "scripts/build_release.py"),
    ("manuscript/PEDS_Research_Article_English_draft_v2_2026-11-17.md", "publication/manuscript/main_article.md"),
    ("manuscript/supplement/Supplementary_Materials_PEDS_2026-11-16.md", "publication/manuscript/supplementary_materials.md"),
    ("manuscript/review/review_pdf_qc_results_v2_2026-11-18.json", "publication/manuscript/review_pdf_qc_results.json"),
    ("manuscript/DATA_SOURCE_VERIFICATION_2026-11-18.md", "docs/DATA_SOURCE_VERIFICATION.md"),
    ("manuscript/supplement/supplement_source_manifest_2026-11-16.md", "data/provenance/supplement_source_manifest.md"),
    ("manuscript/supplement/supplement_build_manifest.json", "data/provenance/supplement_build_manifest.json"),
    ("manuscript/review/RQ5RF_results_and_supplement_note_for_review_2026-11-18.md", "docs/RQ5RF_REVIEW_NOTE.md"),
]

LARGE_RELEASE_ASSET_NAMES = {
    "PEDS_Main_and_Supplement_Review_Draft_v2_2026-11-18.pdf",
    "Supplementary_Tables_S1_S16.xlsx",
    "rq4c_pilot_parameters.csv",
}

FIGURE_GLOBS = [
    ("manuscript/figures/Figure_*.png", "publication/figures/main"),
    ("manuscript/supplement/figures/Supplementary_Figure_S*.png", "publication/figures/supplementary"),
]


def should_exclude(rel: str) -> bool:
    lower = rel.lower()
    return any(
        token in lower
        for token in (
            "overath",
            "pinc_pair_extension",
            "pair_processing_manifest",
            "pilot_pair_processing_manifest",
            "pilot_pinc_pair_extension",
        )
    )


def sanitize_text(text: str) -> str:
    root_patterns = [
        r'ROOT\s*=\s*Path\(r?["\']E:[\\/]diyipian["\']\)',
        r'ROOT\s*=\s*Path\(__file__\)\.resolve\(\)\.parents\[1\]',
    ]
    root_line = 'ROOT = next(_p for _p in Path(__file__).resolve().parents if (_p / ".release-root").is_file())'
    for pattern in root_patterns:
        text = re.sub(pattern, root_line, text)
    replacements = {
        "E:\\diyipian": "<REPRODUCTION_ROOT>",
        "E:/diyipian": "<REPRODUCTION_ROOT>",
        "D:\\KimiData\\daimon-share\\daimon\\runtime\\python\\.venv\\Scripts\\python.exe": "<EXECUTION_PYTHON>",
        "C:\\Users\\leo": "<LOCAL_USER>",
    }
    for old, new in replacements.items():
        text = text.replace(old, new)
    return text


def copy_file(source_rel: str, target_rel: str, *, text_mode: bool, category: str) -> dict[str, object]:
    source = ROOT / source_rel
    target = RELEASE / target_rel
    if not source.is_file():
        raise FileNotFoundError(source)
    target.parent.mkdir(parents=True, exist_ok=True)
    shutil.copy2(source, target)
    transformed = False
    if text_mode:
        original = target.read_text(encoding="utf-8", errors="replace")
        sanitized = sanitize_text(original)
        if sanitized != original:
            target.write_text(sanitized, encoding="utf-8")
            transformed = True
    return {
        "source_path": source_rel,
        "published_path": target_rel,
        "source_sha256": sha256(source),
        "published_sha256": sha256(target),
        "published_bytes": target.stat().st_size,
        "category": category,
        "transform": "text-path-redacted" if transformed else "verbatim",
    }


def main() -> None:
    RELEASE.mkdir(parents=True, exist_ok=True)
    entries: list[dict[str, object]] = []
    seen_targets: set[str] = set()

    candidate_paths = sorted(set(source_manifest_relative_paths()) | set(CODE_FILES))
    for source_rel in candidate_paths:
        if should_exclude(source_rel) or Path(source_rel).name in LARGE_RELEASE_ASSET_NAMES:
            continue
        source = ROOT / source_rel
        if not source.is_file():
            raise FileNotFoundError(f"Selected source is missing: {source_rel}")
        target_rel = source_rel
        text_mode = source.suffix.lower() in {".py", ".md", ".json", ".txt", ".csv", ".yaml", ".yml"}
        category = "code" if source.suffix == ".py" else "derived_data_or_provenance"
        if target_rel not in seen_targets:
            entries.append(copy_file(source_rel, target_rel, text_mode=text_mode, category=category))
            seen_targets.add(target_rel)

    for source_rel, target_rel in EXACT_FILES:
        if target_rel in seen_targets:
            continue
        text_mode = Path(source_rel).suffix.lower() in {".md", ".json", ".txt", ".csv"}
        entries.append(copy_file(source_rel, target_rel, text_mode=text_mode, category="publication_output"))
        seen_targets.add(target_rel)

    for pattern, target_dir in FIGURE_GLOBS:
        for source in sorted(ROOT.glob(pattern)):
            target_rel = f"{target_dir}/{source.name}"
            if target_rel in seen_targets:
                continue
            entries.append(copy_file(relative(source), target_rel, text_mode=False, category="figure"))
            seen_targets.add(target_rel)

    # Authored repository files are generated in place; include them in the integrity manifest too.
    for path in sorted(RELEASE.rglob("*")):
        if not path.is_file() or path == MANIFEST:
            continue
        relative_parts = path.relative_to(RELEASE).parts
        if ".git" in relative_parts or "__pycache__" in relative_parts or path.suffix in {".pyc", ".pyo"}:
            continue
        target_rel = public_relative(path)
        if target_rel in seen_targets:
            continue
        entries.append(
            {
                "source_path": "<repository>",
                "published_path": target_rel,
                "source_sha256": "",
                "published_sha256": sha256(path),
                "published_bytes": path.stat().st_size,
                "category": "repository_document",
                "transform": "authored",
            }
        )
        seen_targets.add(target_rel)

    MANIFEST.parent.mkdir(parents=True, exist_ok=True)
    with MANIFEST.open("w", encoding="utf-8", newline="") as stream:
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
        writer.writerows(sorted(entries, key=lambda row: str(row["published_path"])))

    total_bytes = sum(int(row["published_bytes"]) for row in entries)
    print(f"release: {RELEASE}")
    print(f"files: {len(entries)}")
    print(f"bytes: {total_bytes}")
    print(f"manifest: {MANIFEST}")


if __name__ == "__main__":
    main()
