"""Create deterministic ZIPs, 1 MB pieces, and the release-asset manifest."""
from __future__ import annotations

import argparse
import hashlib
import json
import shutil
from pathlib import Path
from zipfile import ZIP_DEFLATED, ZipFile, ZipInfo


PIECE_BYTES = 1_000_000


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        for chunk in iter(lambda: stream.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def deterministic_zip(source: Path, archive: Path) -> None:
    info = ZipInfo(source.name, date_time=(1980, 1, 1, 0, 0, 0))
    info.compress_type = ZIP_DEFLATED
    info.external_attr = 0o644 << 16
    with ZipFile(archive, "w", compression=ZIP_DEFLATED, compresslevel=9) as bundle:
        bundle.writestr(info, source.read_bytes())


def split_archive(archive: Path, asset_dir: Path) -> list[dict[str, object]]:
    pieces: list[dict[str, object]] = []
    data = archive.read_bytes()
    for index, start in enumerate(range(0, len(data), PIECE_BYTES), start=1):
        piece = asset_dir / f"{archive.name}.piece{index:03d}"
        piece.write_bytes(data[start : start + PIECE_BYTES])
        pieces.append(
            {
                "name": piece.name,
                "bytes": piece.stat().st_size,
                "sha256": sha256(piece),
            }
        )
    return pieces


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--pdf", type=Path, required=True)
    parser.add_argument("--workbook", type=Path, required=True)
    parser.add_argument("--pilot", type=Path, required=True)
    parser.add_argument("--asset-dir", type=Path, required=True)
    parser.add_argument("--manifest", type=Path, required=True)
    args = parser.parse_args()
    asset_dir = args.asset_dir.resolve()
    if asset_dir.exists() and any(asset_dir.iterdir()):
        raise SystemExit(f"Asset directory must be absent or empty: {asset_dir}")
    asset_dir.mkdir(parents=True, exist_ok=True)

    outputs = [
        (
            args.pdf.resolve(),
            "PEDS_Main_and_Supplement_Review_Draft_v1.1.0.pdf",
        ),
        (
            args.workbook.resolve(),
            "Supplementary_Tables_S1_S16_v1.1.0.xlsx",
        ),
        (
            args.pilot.resolve(),
            "rq4c_pilot_parameters_v1.1.0.csv",
        ),
    ]
    assets = []
    for source, published_name in outputs:
        if not source.is_file():
            raise FileNotFoundError(source)
        published = asset_dir / published_name
        shutil.copy2(source, published)
        archive = asset_dir / f"{published_name}.zip"
        deterministic_zip(published, archive)
        assets.append(
            {
                "archive": archive.name,
                "archive_bytes": archive.stat().st_size,
                "archive_sha256": sha256(archive),
                "pieces": split_archive(archive, asset_dir),
                "output": {
                    "name": published_name,
                    "bytes": published.stat().st_size,
                    "sha256": sha256(published),
                },
            }
        )

    manifest = {
        "release": "v1.1.0",
        "piece_bytes": PIECE_BYTES,
        "assets": assets,
    }
    args.manifest.parent.mkdir(parents=True, exist_ok=True)
    args.manifest.write_text(json.dumps(manifest, indent=2, ensure_ascii=False) + "\n", encoding="utf-8")
    print(json.dumps(manifest, indent=2, ensure_ascii=False))


if __name__ == "__main__":
    main()
