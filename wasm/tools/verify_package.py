"""Verify the final static package and compiled-source provenance, without MD."""
from __future__ import annotations
import gzip
import hashlib
import json
from pathlib import Path
import zipfile

ROOT = Path(__file__).resolve().parents[2]


def main():
    site = ROOT / "wasm" / "dist"
    manifest = json.loads((site / "site-manifest.json").read_text(encoding="utf-8"))
    archive = ROOT / "outputs" / "wasm-preparation" / "racing-molecule-web-preview.zip"
    with zipfile.ZipFile(archive) as bundle:
        assert bundle.testzip() is None, "ZIP CRC failed"
        assert json.loads(bundle.read("site-manifest.json")) == manifest
        for item in manifest["files"]:
            data = (site / item["path"]).read_bytes()
            assert len(data) == item["bytes"] and hashlib.sha256(data).hexdigest() == item["sha256"], item["path"]
            assert bundle.read(item["path"]) == data, item["path"]
            assert not any(part in (".local", ".venv", "node_modules") for part in Path(item["path"]).parts)
    for metadata_name in ("model.json", "fourier-atlas.json"):
        meta = json.loads((site / "browser/assets" / metadata_name).read_text())
        descriptors = meta["surface"]["arrays"].values() if metadata_name == "model.json" else meta["arrays"]
        for descriptor in descriptors:
            if not descriptor.get("chunks"):
                continue
            pieces = []
            for part in descriptor["chunks"]:
                data = (site / "browser" / part["path"]).read_bytes()
                assert len(data) == part["bytes"] and hashlib.sha256(data).hexdigest() == part["sha256"]
                pieces.append(data)
            assert b"".join(pieces) == (site / "browser" / descriptor["path"]).read_bytes()
    build = json.loads((site / "browser" / "runtime" / "build.json").read_text(encoding="utf-8"))
    for source, digest in build["own_sources"].items():
        assert hashlib.sha256((ROOT / source).read_bytes()).hexdigest() == digest, f"WASM is stale: {source}"
    model = json.loads((site / "browser" / "assets" / "model.json").read_text(encoding="utf-8"))
    for descriptor in model["surface"]["arrays"].values():
        data = gzip.decompress((site / "browser" / descriptor["path"]).read_bytes())
        assert len(data) == descriptor["bytes"] and hashlib.sha256(data).hexdigest() == descriptor["sha256"]
    for filename in ("fourier.json", "fourier-atlas.json"):
        proxy_file = site / "browser/assets" / filename
        if not proxy_file.exists():
            continue
        proxy = json.loads(proxy_file.read_text())
        assert proxy['base_model_sha256'] == hashlib.sha256((site / 'browser/assets/model.json').read_bytes()).hexdigest()
        assert proxy['parameter_sha256'] == model['surface']['parameter_sha256']
        for descriptor in proxy['arrays']:
            data = gzip.decompress((site / 'browser' / descriptor['path']).read_bytes())
            assert len(data) == descriptor['bytes'] and hashlib.sha256(data).hexdigest() == descriptor['sha256']
    home = (site / 'index.html').read_text(encoding='utf-8')
    assert 'src="./ui/app.js"' in home and 'src="/app.js"' not in home
    settings_source=(site/'ui/app.js').read_text(encoding='utf-8')
    assert 'id="settings-race"' not in settings_source and 'id="mode-lock-note"' not in settings_source
    html = (site / "race.html").read_text(encoding="utf-8")
    assert 'data-engine="wasm"' in html and 'src="/app.js"' not in html and 'href="/style.css"' not in html
    assert '工坊待接入' not in html and 'data-app-route="garage"' in html
    workshop=(site/'editor/index.html').read_text(encoding='utf-8')
    assert 'data-engine="wasm"' in workshop and 'workshop-theme.css' in workshop
    assert 'id="competition-dialog"' not in workshop
    assert 'data-app-route="settings"' not in workshop
    ketcher=site/'editor/vendor/ketcher-c60'
    assert '/editor/vendor/ketcher-c60/' not in (ketcher/'index.html').read_text(encoding='utf-8')
    for asset in (ketcher/'assets').glob('*.js'):
        assert '/editor/vendor/ketcher-c60/' not in asset.read_text(encoding='utf-8'), f'Root-only editor URL: {asset.name}'
    for filename in ("app.js", "transport.js", "gold-surface.js", "style.css", "start-pose.js", "start-pose.css", "competition.mjs", "competition.css", "app-state.mjs", "design-system.css", "race-theme.css", "race-result.mjs", "race-result.css", "playback.mjs"):
        assert (site / filename).read_bytes() == (ROOT / "frontend" / filename).read_bytes(), f"Frontend is stale: {filename}"
    for filename in ("worker.mjs", "transport.mjs", "editor-worker.mjs", "embedded.mjs"):
        assert (site / "browser" / filename).read_bytes() == (ROOT / "frontend" / "browser" / filename).read_bytes()
    for filename in ('workshop.js','browser-service.mjs','ketcher-host.js','workshop-theme.css','ketcher-theme.css'):
        assert (site/'editor'/filename).read_bytes()==(ROOT/'frontend/editor'/filename).read_bytes()
    for filename in ('app.js','integration.css','minimal.css'):
        assert (site/'ui'/filename).read_bytes()==(ROOT/'web-ui/dist'/filename).read_bytes()
    output = {"passed": True, "manifest_files": len(manifest["files"]),
              "site_bytes": sum(item["bytes"] for item in manifest["files"]), "zip_bytes": archive.stat().st_size,
              "zip_sha256": hashlib.sha256(archive.read_bytes()).hexdigest(),
              "source_provenance": True, "potential_checksums": True, "zip_crc": True}
    (ROOT / "outputs" / "wasm-preparation" / "package-validation.json").write_text(json.dumps(output, indent=2)+"\n", encoding="utf-8")
    print(json.dumps(output, indent=2))


if __name__ == "__main__":
    main()
