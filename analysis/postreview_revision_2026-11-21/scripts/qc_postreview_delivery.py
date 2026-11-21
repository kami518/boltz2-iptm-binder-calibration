"""Final consistency checks for the post-review manuscript package."""
from __future__ import annotations

import hashlib
import json
import re
from pathlib import Path

import pandas as pd
from pypdf import PdfReader


ROOT = Path(__file__).resolve().parents[3]
OUT = ROOT / "analysis" / "postreview_revision_2026-11-21"
MAIN = ROOT / "publication" / "manuscript" / "main_article.md"
SUPP = ROOT / "publication" / "manuscript" / "supplementary_materials.md"
PDF = ROOT / "publication" / "manuscript" / "PEDS_Main_and_Supplement_Review_Draft_v1.1.0.pdf"
PDF_QC = ROOT / "publication" / "manuscript" / "review_pdf_qc_results.json"
FIGURE_MANIFEST = OUT / "figures" / "figure_manifest_postreview_2026-11-21.json"
RESULT = OUT / "POSTREVIEW_DELIVERY_QC_2026-11-21.json"


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def abtest(name: str, condition: bool, detail) -> tuple[str, dict]:
    return name, {"passed": bool(condition), "detail": detail}


main = MAIN.read_text(encoding="utf-8")
supp = SUPP.read_text(encoding="utf-8")
c01 = json.loads((OUT / "ipSAE" / "run_manifest.json").read_text(encoding="utf-8"))
c01b = json.loads((OUT / "ipSAE_combo" / "run_manifest.json").read_text(encoding="utf-8"))
c02 = json.loads((OUT / "bootstrap" / "run_manifest.json").read_text(encoding="utf-8"))
fig_manifest = json.loads(FIGURE_MANIFEST.read_text(encoding="utf-8"))
pdf_qc = json.loads(PDF_QC.read_text(encoding="utf-8"))
reader = PdfReader(str(PDF))

refs_section = main.split("## References", 1)[1].split("## Main tables", 1)[0]
refs = [part.strip() for part in refs_section.split("\n\n") if part.strip()]
abstract = re.search(
    r"^## Abstract\s*$(.*?)(?=^\*\*Keywords)",
    main,
    re.M | re.S,
).group(1)
abstract_sentences = [
    item for item in re.split(r"(?<=[.!?])\s+", abstract.strip()) if item
]
abstract_words = len(
    re.findall(r"[A-Za-z0-9]+(?:[-'][A-Za-z0-9]+)*", abstract)
)

figure_files = []
for figure in fig_manifest["figures"].values():
    for item in figure.values():
        figure_files.append(ROOT / item["path"])

official = pd.read_csv(OUT / "ipSAE" / "ipsae_official_pair_metrics.csv", low_memory=False)
combo = pd.read_csv(OUT / "ipSAE_combo" / "combo_target_metrics.csv", low_memory=False)
bootstrap_macro = pd.read_csv(OUT / "bootstrap" / "bootstrap_macro_intervals.csv", low_memory=False)

checks = dict(
    [
        abtest(
            "reference_count_and_uniqueness",
            len(refs) == 41
            and len(re.findall(r"https://doi.org/", refs_section)) == 35
            and len(set(re.findall(r"https://doi.org/([^\s)]+)", refs_section))) == 35,
            {
                "references": len(refs),
                "doi_mentions": len(re.findall(r"https://doi.org/", refs_section)),
                "unique_dois": len(set(re.findall(r"https://doi.org/([^\s)]+)", refs_section))),
            },
        ),
        abtest(
            "references_alphabetical",
            all(refs[index].casefold() <= refs[index + 1].casefold() for index in range(len(refs) - 1)),
            "Reference entries are alphabetically ordered.",
        ),
        abtest(
            "abstract_length",
            abstract_words <= 250 and 7 <= len(abstract_sentences) <= 9,
            {"words": abstract_words, "sentences": len(abstract_sentences)},
        ),
        abtest(
            "title_consistent",
            main.splitlines()[0].lstrip("# ").strip()
            == supp.splitlines()[2].lstrip("# ").strip()
            == "Context-dependent calibration and experimental decision value of Boltz-2 ipTM in binder design",
            {
                "main": main.splitlines()[0],
                "supplement": supp.splitlines()[2],
            },
        ),
        abtest(
            "rejected_wording_absent",
            "prespecified two-feature" not in main.casefold()
            and "universal selection threshold" not in main.casefold(),
            "Rejected wording is absent from the main manuscript.",
        ),
        abtest(
            "canonical_c01b_sentence",
            "neither ipSAE alone nor a post hoc two-feature logistic model combining ipSAE with ipTM" in main,
            "Canonical post hoc wording is present.",
        ),
        abtest(
            "bootstrap_intervals_locked",
            all(
                value in main and value in supp
                for value in (
                    "45.1%",
                    "42.9%",
                    "25.9%",
                    "39.5% to 49.8%",
                    "38.4% to 46.8%",
                    "19.3% to 31.9%",
                )
            ),
            "Primary, five-target, and external intervals are present in both documents.",
        ),
        abtest(
            "anthropic_report_cited",
            "Claude Science and Shanehsazzadeh, 2026" in main
            and "Autonomous de novo protein binder design with Claude" in refs_section,
            "Technical report and data release are cited distinctly.",
        ),
        abtest(
            "author_metadata_and_declarations_complete",
            all(
                value in main
                for value in (
                    "KangMing Chen",
                    "0009-0006-9634-5397",
                    "LeoChen070809@outlook.com",
                    "This research received no specific grant from any funding agency.",
                    "The author declares no conflicts of interest.",
                    "The author thanks Kai Li for supervision, guidance, and support during this work.",
                )
            )
            and "To be supplied" not in main,
            "Confirmed author metadata and declarations are complete.",
        ),
        abtest(
            "c01_official_batch_complete",
            c01["status"] == "complete"
            and len(official) == 2074
            and official["status"].eq("passed").all(),
            {
                "manifest_status": c01["status"],
                "rows": len(official),
                "passed": int(official["status"].eq("passed").sum()),
            },
        ),
        abtest(
            "c01b_complete_no_rf",
            c01b["status"] == "complete"
            and not c01b["random_forest_run"]
            and not c01b["other_features_added"]
            and len(combo) == 36,
            {
                "manifest_status": c01b["status"],
                "rows": len(combo),
                "random_forest_run": c01b["random_forest_run"],
            },
        ),
        abtest(
            "c02_replicates_complete",
            c02["status"] == "complete"
            and bootstrap_macro["n_valid_boot"].eq(2000).all()
            and bootstrap_macro["n_failed_boot"].eq(0).all(),
            {
                "manifest_status": c02["status"],
                "valid_replicates": sorted(bootstrap_macro["n_valid_boot"].unique().tolist()),
                "failed_replicates": sorted(bootstrap_macro["n_failed_boot"].unique().tolist()),
            },
        ),
        abtest(
            "all_vector_and_raster_outputs_exist",
            len(figure_files) == 24 and all(path.is_file() for path in figure_files),
            {
                "expected": 24,
                "found": sum(path.is_file() for path in figure_files),
            },
        ),
        abtest(
            "no_internal_rq_labels_in_figure_vectors",
            not any(
                re.search(r"\bRQ[1-6][A-Z]?\b", path.read_text(encoding="utf-8"))
                for path in figure_files
                if path.suffix == ".svg"
            ),
            "Generated SVG panels contain no RQ work labels.",
        ),
        abtest(
            "supplement_s17_and_s11",
            "### Table S17." in supp
            and "### Supplementary Figure S11." in supp
            and all(f"### Table S{index}." in supp for index in range(1, 18)),
            "S1-S16 order is retained and S17/S11 are present.",
        ),
        abtest(
            "pdf_qc_passed",
            pdf_qc["status"] == "passed",
            {
                "pages": pdf_qc["pages"],
                "main_figures": pdf_qc["main_figures"],
                "supplementary_figures": pdf_qc["supplementary_figures"],
                "supplementary_table_guides": pdf_qc["supplementary_table_guides"],
            },
        ),
        abtest(
            "pdf_metadata_title",
            reader.metadata.title
            == "Context-dependent calibration and experimental decision value of Boltz-2 ipTM in binder design",
            reader.metadata.title,
        ),
    ]
)

result = {
    "status": "passed" if all(item["passed"] for item in checks.values()) else "failed",
    "main": {"path": MAIN.relative_to(ROOT).as_posix(), "sha256": sha256(MAIN)},
    "supplement": {"path": SUPP.relative_to(ROOT).as_posix(), "sha256": sha256(SUPP)},
    "pdf": {
        "path": PDF.relative_to(ROOT).as_posix(),
        "sha256": sha256(PDF),
        "pages": len(reader.pages),
        "bytes": PDF.stat().st_size,
    },
    "checks": checks,
}
RESULT.write_text(json.dumps(result, ensure_ascii=False, indent=2) + "\n", encoding="utf-8")
print(json.dumps(result, ensure_ascii=False, indent=2))
if result["status"] != "passed":
    raise SystemExit(1)
