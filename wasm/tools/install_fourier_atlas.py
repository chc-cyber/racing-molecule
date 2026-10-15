"""Install verified, generated v2 assets alongside the reproducible v1 assets."""
from pathlib import Path
import gzip
import hashlib
import json
import shutil

ROOT = Path(__file__).resolve().parents[2]


def main():
    source = ROOT / "outputs/c60-proxy-study"
    target = ROOT / "wasm/dist/browser"
    metadata = json.loads((source / "assets/fourier-atlas.json").read_text())
    model = target / "assets/model.json"
    assert metadata["schema"] == "fourier-saip-browser-v2"
    assert metadata["base_model_sha256"] == hashlib.sha256(model.read_bytes()).hexdigest()
    for descriptor in metadata["arrays"]:
        relative = Path(descriptor["path"])
        assert relative.parts == ("assets", "fourier-atlas-0.bin.gz")
        data = gzip.decompress((source / relative).read_bytes())
        assert len(data) == descriptor["bytes"]
        assert hashlib.sha256(data).hexdigest() == descriptor["sha256"]
        shutil.copy2(source / relative, target / relative)
    shutil.copy2(source / "assets/fourier-atlas.json", target / "assets/fourier-atlas.json")
    print("Installed C-only three-chart atlas; v1 reference files retained.")


if __name__ == "__main__":
    main()
