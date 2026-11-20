"""Install the pinned MMseqs2 Windows binary inside the workspace and verify its SHA-256."""
from __future__ import annotations

import hashlib
import json
import urllib.request
import zipfile
from pathlib import Path

ROOT = next(_p for _p in Path(__file__).resolve().parents if (_p / ".release-root").is_file())
URL = "https://github.com/soedinglab/MMseqs2/releases/download/18-8cc5c/mmseqs-win64.zip"
EXPECTED_SHA256 = "9752087ca4e25e93b830775832a75044dc94ae348d66d4681f9ff0355eecfcbf"
OUT_DIR = ROOT / "tools" / "mmseqs2-18-8cc5c"
ARCHIVE = OUT_DIR / "mmseqs-win64.zip"
MANIFEST = ROOT / "analysis" / "audit_outputs" / "mmseqs_install_manifest_rq12.json"


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def main() -> None:
    OUT_DIR.mkdir(parents=True, exist_ok=True)
    if not ARCHIVE.exists() or sha256(ARCHIVE) != EXPECTED_SHA256:
        request = urllib.request.Request(URL, headers={"User-Agent": "rq12-reproducible-setup/1.0"})
        with urllib.request.urlopen(request, timeout=120) as response:
            payload = response.read()
        ARCHIVE.write_bytes(payload)
    actual_sha256 = sha256(ARCHIVE)
    if actual_sha256 != EXPECTED_SHA256:
        raise RuntimeError(f"SHA-256 mismatch: {actual_sha256}")
    extract_dir = OUT_DIR / "extracted"
    with zipfile.ZipFile(ARCHIVE) as archive:
        archive.extractall(extract_dir)
    executables = sorted(extract_dir.rglob("mmseqs.exe"))
    if len(executables) != 1:
        raise RuntimeError(f"Expected one mmseqs.exe, found: {executables}")
    manifest = {
        "url": URL,
        "release": "18-8cc5c",
        "archive": str(ARCHIVE),
        "bytes": ARCHIVE.stat().st_size,
        "sha256_expected": EXPECTED_SHA256,
        "sha256_observed": actual_sha256,
        "executable": str(executables[0]),
        "status": "download_and_checksum_passed",
    }
    MANIFEST.write_text(json.dumps(manifest, indent=2), encoding="utf-8")
    print(json.dumps(manifest, indent=2))


if __name__ == "__main__":
    main()
