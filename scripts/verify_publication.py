"""Verify a public reproducibility candidate without network access."""
from __future__ import annotations

import argparse
import ast
import csv
import hashlib
import json
import re
from collections import Counter
from pathlib import Path
from urllib.parse import unquote
from zipfile import ZipFile
from xml.etree import ElementTree as ET


ROOT = Path(__file__).resolve().parents[1]
MANIFEST = ROOT / "data" / "PUBLICATION_MANIFEST.csv"
TEXT_SUFFIXES = {".py", ".md", ".txt", ".json", ".csv", ".yaml", ".yml", ".toml", ".cff"}
LOCAL_PATH = re.compile(r"(?i)(?:\bfile:(?://|[/\\])|(?<![A-Za-z0-9])[A-Z]:[\\/]|/Users/|\\Users\\|/home/|\.codex)")
CREDENTIALS = re.compile(
    r"(?:BEGIN (?:RSA |OPENSSH )?PRIVATE KEY|gh[pousr]_[A-Za-z0-9]{20,}|"
    r"github_pat_[A-Za-z0-9_]{20,})"
)


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        for chunk in iter(lambda: stream.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def load_manifest() -> list[dict[str, str]]:
    if not MANIFEST.is_file():
        raise SystemExit(f"Missing manifest: {MANIFEST}")
    with MANIFEST.open(encoding="utf-8", newline="") as stream:
        return list(csv.DictReader(stream))


def verify_manifest(rows: list[dict[str, str]]) -> list[str]:
    errors: list[str] = []
    paths = [row["published_path"] for row in rows]
    duplicates = [path for path, count in Counter(paths).items() if count > 1]
    if duplicates:
        errors.append(f"duplicate manifest paths: {duplicates}")
    for row in rows:
        published = row["published_path"]
        if published.startswith(("/", "\\")) or re.match(r"^[A-Za-z]:", published):
            errors.append(f"absolute manifest path: {published}")
            continue
        if ".." in Path(published).parts:
            errors.append(f"path traversal in manifest: {published}")
            continue
        path = ROOT / published
        if not path.is_file():
            errors.append(f"missing: {published}")
        elif path.stat().st_size != int(row["published_bytes"]):
            errors.append(f"size mismatch: {published}")
        elif sha256(path) != row["published_sha256"]:
            errors.append(f"hash mismatch: {published}")
    actual = {
        path.relative_to(ROOT).as_posix()
        for path in ROOT.rglob("*")
        if path.is_file()
        and ".git" not in path.relative_to(ROOT).parts
        and "__pycache__" not in path.relative_to(ROOT).parts
        and path.relative_to(ROOT).as_posix() != "data/PUBLICATION_MANIFEST.csv"
    }
    expected = set(paths)
    if actual - expected:
        errors.append(f"files outside manifest: {sorted(actual - expected)}")
    if expected - actual:
        errors.append(f"manifest files absent: {sorted(expected - actual)}")
    return errors


def scan_python() -> list[str]:
    errors: list[str] = []
    for path in sorted(ROOT.rglob("*.py")):
        try:
            ast.parse(path.read_text(encoding="utf-8"), filename=str(path))
        except SyntaxError as error:
            errors.append(f"python syntax: {path.relative_to(ROOT)}: {error}")
    return errors


def scan_text() -> list[str]:
    errors: list[str] = []
    self_path = Path(__file__).resolve()
    tool_patterns = {
        (ROOT / "scripts" / "build_release.py").resolve(),
        (ROOT / "scripts" / "build_supplement_manifest.py").resolve(),
        (ROOT / "publication" / "manuscript" / "qc_review_pdf.py").resolve(),
    }
    for path in sorted(ROOT.rglob("*")):
        if not path.is_file() or path.suffix.lower() not in TEXT_SUFFIXES or path.resolve() in tool_patterns or path.resolve() == self_path:
            continue
        text = path.read_text(encoding="utf-8", errors="replace")
        relative = path.relative_to(ROOT).as_posix()
        if LOCAL_PATH.search(text):
            errors.append(f"local path or file URI: {relative}")
        if path.suffix.lower() != ".csv" and CREDENTIALS.search(text):
            errors.append(f"credential-like text: {relative}")
    return errors


def scan_markdown_links() -> list[str]:
    errors: list[str] = []
    for path in sorted(ROOT.rglob("*.md")):
        text = path.read_text(encoding="utf-8")
        for match in re.finditer(r"(?<!!)\[[^\]]+\]\(([^)]+)\)", text):
            target = match.group(1).split("#", 1)[0]
            if not target or re.match(r"^[A-Za-z][A-Za-z0-9+.-]*:", target) or target.startswith("#"):
                continue
            resolved = (path.parent / unquote(target)).resolve()
            if not resolved.exists():
                errors.append(f"broken Markdown link: {path.relative_to(ROOT)} -> {target}")
            try:
                resolved.relative_to(ROOT.resolve())
            except ValueError:
                errors.append(f"Markdown link escapes root: {path.relative_to(ROOT)} -> {target}")
    return errors


def scan_xlsx(path: Path) -> list[str]:
    errors: list[str] = []
    with ZipFile(path) as archive:
        shared: list[str] = []
        if "xl/sharedStrings.xml" in archive.namelist():
            ns = {"m": "http://schemas.openxmlformats.org/spreadsheetml/2006/main"}
            root = ET.fromstring(archive.read("xl/sharedStrings.xml"))
            for item in root.findall("m:si", ns):
                shared.append("".join(node.text or "" for node in item.iter("{http://schemas.openxmlformats.org/spreadsheetml/2006/main}t")))
        for name in archive.namelist():
            if not re.fullmatch(r"xl/worksheets/sheet\d+\.xml", name):
                continue
            root = ET.fromstring(archive.read(name))
            for cell in root.iter("{http://schemas.openxmlformats.org/spreadsheetml/2006/main}c"):
                value = cell.find("{http://schemas.openxmlformats.org/spreadsheetml/2006/main}v")
                text = ""
                if cell.attrib.get("t") == "s" and value is not None:
                    text = shared[int(value.text)]
                elif cell.attrib.get("t") == "inlineStr":
                    text = "".join(
                        node.text or ""
                        for node in cell.iter("{http://schemas.openxmlformats.org/spreadsheetml/2006/main}t")
                    )
                elif value is not None:
                    text = value.text or ""
                if LOCAL_PATH.search(text):
                    errors.append(f"{path.name}:{name}:{cell.attrib.get('r')}")
    return errors


def scan_pdf(path: Path) -> list[str]:
    try:
        from pypdf import PdfReader
    except ImportError as error:
        raise SystemExit("pypdf is required for PDF verification") from error
    errors: list[str] = []
    reader = PdfReader(str(path))
    for page_number, page in enumerate(reader.pages, start=1):
        text = page.extract_text() or ""
        if LOCAL_PATH.search(text):
            errors.append(f"PDF text local path: page {page_number}")
        for ref in page.get("/Annots", []):
            annotation = ref.get_object()
            if annotation.get("/Subtype") != "/Link":
                continue
            action = annotation.get("/A")
            action = action.get_object() if action else None
            uri = action.get("/URI") if action and action.get("/S") == "/URI" else None
            if uri and not str(uri).startswith(("http://", "https://")):
                errors.append(f"PDF non-HTTP(S) URI: page {page_number}: {str(uri)[:200]}")
    return errors


def verify_dependency_manifest() -> list[str]:
    path = ROOT / "data" / "postreview_dependency_manifest.json"
    if not path.is_file():
        return ["missing postreview dependency manifest"]
    try:
        data = json.loads(path.read_text(encoding="utf-8"))
    except json.JSONDecodeError as error:
        return [f"invalid dependency manifest: {error}"]
    errors: list[str] = []
    dependencies = data.get("dependencies", [])
    if len(dependencies) != 11:
        errors.append(f"dependency count mismatch: {len(dependencies)}")
    for item in dependencies:
        relative = item.get("path")
        if not relative:
            errors.append("dependency missing path")
            continue
        published = ROOT / relative
        if not published.is_file():
            errors.append(f"dependency missing: {relative}")
            continue
        if published.stat().st_size != item.get("published_bytes"):
            errors.append(f"dependency published size mismatch: {relative}")
        if sha256(published) != item.get("published_sha256"):
            errors.append(f"dependency published hash mismatch: {relative}")
        if not item.get("source_sha256") or not item.get("published_sha256"):
            errors.append(f"dependency missing dual hash: {relative}")
        if item.get("transform") == "path-redaction" and item.get("source_sha256") == item.get("published_sha256"):
            errors.append(f"dependency transformation not reflected in hash: {relative}")
    return errors


def negative_tests() -> list[str]:
    samples = [r"Q:\private\data.csv", "Q:/private/data.csv", "file:///Q:/private/data.csv"]
    return [] if all(LOCAL_PATH.search(sample) for sample in samples) else ["negative privacy tests failed"]


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--asset-dir", type=Path)
    parser.add_argument("--output-json", type=Path)
    args = parser.parse_args()

    rows = load_manifest()
    errors = []
    errors.extend(verify_manifest(rows))
    errors.extend(scan_python())
    errors.extend(scan_text())
    errors.extend(scan_markdown_links())
    errors.extend(verify_dependency_manifest())
    errors.extend(negative_tests())

    asset_dir = args.asset_dir.resolve() if args.asset_dir else None
    if asset_dir:
        xlsx = asset_dir / "Supplementary_Tables_S1_S16_v1.1.0.xlsx"
        pdf = asset_dir / "PEDS_Main_and_Supplement_Review_Draft_v1.1.0.pdf"
        if not xlsx.is_file():
            errors.append(f"missing release workbook: {xlsx.name}")
        else:
            errors.extend(scan_xlsx(xlsx))
        if not pdf.is_file():
            errors.append(f"missing release PDF: {pdf.name}")
        else:
            errors.extend(scan_pdf(pdf))

    result = {
        "status": "passed" if not errors else "failed",
        "root": str(ROOT),
        "manifest_files": len(rows),
        "errors": errors,
    }
    rendered = json.dumps(result, indent=2, ensure_ascii=False)
    print(rendered)
    if args.output_json:
        args.output_json.parent.mkdir(parents=True, exist_ok=True)
        args.output_json.write_text(rendered + "\n", encoding="utf-8")
    if errors:
        raise SystemExit(1)


if __name__ == "__main__":
    main()
