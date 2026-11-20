from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path
from zipfile import ZipFile

ROOT = next(path for path in Path(__file__).resolve().parents if (path / ".release-root").is_file())
MANIFEST = ROOT / "data/release_assets_manifest.json"


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        for chunk in iter(lambda: stream.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def verify(path: Path, expected_bytes: int, expected_sha: str) -> None:
    if path.stat().st_size != expected_bytes:
        raise RuntimeError(f"Size mismatch: {path}")
    if sha256(path) != expected_sha:
        raise RuntimeError(f"SHA-256 mismatch: {path}")


def safe_extract(archive: Path, output_dir: Path) -> None:
    output_dir.mkdir(parents=True, exist_ok=True)
    with ZipFile(archive) as bundle:
        for member in bundle.infolist():
            candidate = (output_dir / member.filename).resolve()
            if output_dir.resolve() not in candidate.parents and candidate != output_dir.resolve():
                raise RuntimeError(f"Unsafe archive member: {member.filename}")
        bundle.extractall(output_dir)


def main() -> None:
    parser = argparse.ArgumentParser(description="Reassemble and verify GitHub Release assets.")
    parser.add_argument("asset_dir", type=Path, help="directory containing downloaded release assets")
    parser.add_argument("--output-dir", type=Path, default=ROOT / "release_assets")
    args = parser.parse_args()
    args.asset_dir = args.asset_dir.resolve()
    args.output_dir = args.output_dir.resolve()
    args.output_dir.mkdir(parents=True, exist_ok=True)

    manifest = json.loads(MANIFEST.read_text(encoding="utf-8"))
    for item in manifest["assets"]:
        archive = args.asset_dir / item["archive"]
        if not archive.exists():
            with archive.open("wb") as output:
                for piece in item["pieces"]:
                    path = args.asset_dir / piece["name"]
                    verify(path, piece["bytes"], piece["sha256"])
                    with path.open("rb") as stream:
                        while chunk := stream.read(1024 * 1024):
                            output.write(chunk)
        verify(archive, item["archive_bytes"], item["archive_sha256"])
        safe_extract(archive, args.output_dir)
        result = args.output_dir / item["output"]["name"]
        verify(result, item["output"]["bytes"], item["output"]["sha256"])
        print(f"verified: {result}")


if __name__ == "__main__":
    main()
