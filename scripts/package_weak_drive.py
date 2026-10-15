"""Bundle measured results, immutable coordinates and the reproducible sources."""
from __future__ import annotations
import hashlib
import json
from pathlib import Path
import zipfile

ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / "outputs/weak-drive-20nm"


def main():
    screen = json.loads((OUT / "screen-report.json").read_text(encoding="utf-8"))
    extended = json.loads((OUT / "extended-report.json").read_text(encoding="utf-8"))
    validation = json.loads((OUT / "validation.json").read_text(encoding="utf-8"))
    assert all(report["status"] == "passed" for report in (screen, extended, validation))
    for name, expected in screen["source_sha256"].items():
        assert hashlib.sha256((ROOT / name).read_bytes()).hexdigest() == expected
    archive = OUT / "weak-drive-20nm-results.zip"
    allowed = {".json", ".csv", ".npz", ".png", ".svg", ".md", ".xyz", ".ILP"}
    with zipfile.ZipFile(archive, "w", compression=zipfile.ZIP_DEFLATED) as bundle:
        for path in sorted(OUT.rglob("*")):
            if path.is_file() and path != archive and "matplotlib-cache" not in path.parts and path.suffix in allowed:
                bundle.write(path, path.relative_to(OUT).as_posix())
        sources = set(screen["source_sha256"]) | {
            "scripts/validate_weak_drive.py", "scripts/plot_weak_drive.py", "scripts/package_weak_drive.py",
            "scripts/verify_20nm_game.py", "scripts/verify_20nm_game_browser.cjs",
            "scripts/capture_weak_drive_preview.cjs",
            "tests/test_track_gameplay.py", "tests/test_weak_drive_analysis.py", "requirements.txt",
        }
        for name in sorted(sources):
            bundle.write(ROOT / name, "implementation/" + name)
        for dirname in ("track-resize-check", "20nm-game-check"):
            for path in sorted((ROOT / "outputs" / dirname).glob("*")):
                if path.is_file() and path.suffix in {".json", ".png"}:
                    bundle.write(path, "game-acceptance/" + dirname + "/" + path.name)
    with zipfile.ZipFile(archive) as bundle:
        assert bundle.testzip() is None
        count = len(bundle.namelist())
    print(json.dumps({"archive": str(archive), "files": count, "size_MiB": archive.stat().st_size/2**20,
                      "sha256": hashlib.sha256(archive.read_bytes()).hexdigest()}, ensure_ascii=False))


if __name__ == "__main__":
    main()
