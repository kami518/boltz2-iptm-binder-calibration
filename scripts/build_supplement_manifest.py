"""Generate portable supplement provenance manifests."""
from __future__ import annotations

import hashlib
import json
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
OUT_MD = ROOT / "publication" / "manuscript" / "supplement_source_manifest.md"
OUT_JSON = ROOT / "publication" / "manuscript" / "supplement_build_manifest.json"


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        for chunk in iter(lambda: stream.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def main() -> None:
    dependency = json.loads((ROOT / "data" / "postreview_dependency_manifest.json").read_text(encoding="utf-8"))
    workbook = ROOT.parent / "assets" / "Supplementary_Tables_S1_S16_v1.1.0.xlsx"
    if not workbook.is_file():
        raise SystemExit(f"Missing cleaned workbook: {workbook}")

    supplementary = ROOT / "publication" / "figures" / "supplementary"
    figure_rows = []
    for index in range(1, 11):
        path = supplementary / f"Supplementary_Figure_S{index}.png"
        if not path.is_file():
            raise SystemExit(f"Missing supplementary figure: {path}")
        figure_rows.append((path, "PNG"))
    for suffix in ("png", "svg", "pdf"):
        path = supplementary / f"Supplementary_Figure_S11_ipSAE_postreview.{suffix}"
        if not path.is_file():
            raise SystemExit(f"Missing supplementary figure: {path}")
        figure_rows.append((path, suffix.upper()))

    lines = [
        "# Supplement source and output manifest",
        "",
        "Release: `v1.1.0`",
        "",
        "This manifest records the portable dependency closure and generated supplement outputs. Raw upstream data and third-party source code are excluded.",
        "",
        "## Pinned upstream sources",
        "",
        f"- Proteinbase: `{dependency['external_requirements']['proteinbase']}`",
        f"- Anthropic binder-design release: `{dependency['external_requirements']['anthropic_release']}`",
        f"- IPSAE: `{dependency['external_requirements']['ipsae']}`",
        f"- Pinc: `{dependency['external_requirements']['pinc']}`",
        "",
        "## Dependency closure",
        "",
        "The complete machine-readable dependency manifest is `data/postreview_dependency_manifest.json`.",
        "",
        "| Path | Source bytes | Source SHA-256 | Published bytes | Published SHA-256 |",
        "|---|---:|---|---:|---|",
    ]
    for item in dependency["dependencies"]:
        lines.append(
            f"| `{item['path']}` | {item['source_bytes']:,} | `{item['source_sha256']}` | "
            f"{item['published_bytes']:,} | `{item['published_sha256']}` |"
        )

    lines += [
        "",
        "## Release workbook",
        "",
        f"- File: `Supplementary_Tables_S1_S16_v1.1.0.xlsx`",
        f"- Bytes: {workbook.stat().st_size:,}",
        f"- SHA-256: `{sha256(workbook)}`",
        "",
        "## Supplementary figures",
        "",
        "| File | Format | Bytes | SHA-256 |",
        "|---|---|---:|---|",
    ]
    for path, format_name in figure_rows:
        lines.append(
            f"| `{path.relative_to(ROOT).as_posix()}` | {format_name} | "
            f"{path.stat().st_size:,} | `{sha256(path)}` |"
        )
    OUT_MD.write_text("\n".join(lines) + "\n", encoding="utf-8")

    build_manifest = {
        "status": "complete",
        "release": "v1.1.0",
        "dependencies": "data/postreview_dependency_manifest.json",
        "workbook": {
            "name": workbook.name,
            "bytes": workbook.stat().st_size,
            "sha256": sha256(workbook),
        },
        "figures": {
            path.name: {
                "path": path.relative_to(ROOT).as_posix(),
                "bytes": path.stat().st_size,
                "sha256": sha256(path),
            }
            for path, _ in figure_rows
        },
    }
    OUT_JSON.write_text(json.dumps(build_manifest, indent=2, ensure_ascii=False) + "\n", encoding="utf-8")
    print(f"wrote {OUT_MD.relative_to(ROOT)} and {OUT_JSON.relative_to(ROOT)}")


if __name__ == "__main__":
    main()
