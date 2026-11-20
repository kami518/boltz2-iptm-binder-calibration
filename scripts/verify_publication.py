from __future__ import annotations

import csv
import hashlib
import py_compile
import re
from pathlib import Path

ROOT = next(path for path in Path(__file__).resolve().parents if (path / ".release-root").is_file())
MANIFEST = ROOT / "data/PUBLICATION_MANIFEST.csv"


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        for chunk in iter(lambda: stream.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def verify_manifest() -> None:
    if not MANIFEST.is_file():
        raise SystemExit(f"Missing manifest: {MANIFEST}")
    with MANIFEST.open(encoding="utf-8", newline="") as stream:
        rows = list(csv.DictReader(stream))
    errors: list[str] = []
    for row in rows:
        path = ROOT / row["published_path"]
        if not path.is_file():
            errors.append(f"missing: {row['published_path']}")
            continue
        if path.stat().st_size != int(row["published_bytes"]):
            errors.append(f"size mismatch: {row['published_path']}")
        if sha256(path) != row["published_sha256"]:
            errors.append(f"hash mismatch: {row['published_path']}")
    if errors:
        raise SystemExit("\n".join(errors))
    print(f"manifest: {len(rows)} files passed")


def compile_python() -> None:
    files = sorted(ROOT.glob("analysis/**/*.py")) + sorted(ROOT.glob("scripts/*.py"))
    for path in files:
        py_compile.compile(str(path), doraise=True)
    print(f"python: {len(files)} files compiled")


def scan_text() -> None:
    local_path = re.compile(r"(?:[A-Za-z]:\\(?:Users|diyipian|KimiData)\\|/home/|/Users/)")
    credentials = re.compile(r"(?:BEGIN (?:RSA |OPENSSH )?PRIVATE KEY|gh[pousr]_[A-Za-z0-9]{20,}|github_pat_[A-Za-z0-9_]{20,})")
    suffixes = {".py", ".md", ".txt", ".json", ".yaml", ".yml", ".toml", ".csv", ".cff"}
    findings: list[str] = []
    for path in ROOT.rglob("*"):
        if not path.is_file() or path.suffix.lower() not in suffixes:
            continue
        if path.resolve() == Path(__file__).resolve():
            continue
        text = path.read_text(encoding="utf-8", errors="replace")
        if local_path.search(text):
            findings.append(f"local path: {path.relative_to(ROOT)}")
        if credentials.search(text):
            findings.append(f"credential-like text: {path.relative_to(ROOT)}")
    if findings:
        raise SystemExit("\n".join(findings))
    print("privacy scan: passed")


if __name__ == "__main__":
    verify_manifest()
    compile_python()
    scan_text()
