"""Rebuild the auditable local ModelKit source patch and current file hashes.

The upstream zip metadata and pristine backups stay unchanged. Run this script
only after intentional changes to the five explicitly listed runtime sources.
"""
from pathlib import Path
import difflib
import hashlib
import json

source = Path(__file__).resolve().parent
vendor = source.parent / "vendor" / "jsmol"
paths = ["J/modelkit/ModelKitPopup.js", "J/modelkit/ModelKitPopupResourceBundle.js", "JS/JMenu.js", "core/coreswing.js", "core/coreswing.z.js"]
manifest_path = vendor / "manifest.json"
manifest = json.loads(manifest_path.read_text(encoding="utf-8"))
entries = {entry["path"]: entry for entry in manifest["files"]}
patch = []
modified = []
for name in paths:
    entry = entries["j2s/" + name]
    original = source / "upstream" / name
    current = vendor / "j2s" / name
    upstream_sha = hashlib.sha256(original.read_bytes()).hexdigest()
    assert upstream_sha == entry.get("upstreamSha256", entry["sha256"]), name
    entry["upstreamSha256"] = upstream_sha
    entry["sha256"] = hashlib.sha256(current.read_bytes()).hexdigest()
    entry["bytes"] = len(current.read_bytes())
    patch.extend(difflib.unified_diff(original.read_text(encoding="utf-8").splitlines(keepends=True), current.read_text(encoding="utf-8").splitlines(keepends=True), fromfile="a/j2s/" + name, tofile="b/j2s/" + name))
    modified.append({"path": "j2s/" + name, "upstreamSha256": upstream_sha, "sha256": entry["sha256"]})
patch_path = source / "modelkit-native-bridge.patch"
patch_path.write_text("".join(patch), encoding="utf-8", newline="\n")
manifest["included"] = "JSmol 16.3.31 local runtime with five documented ModelKit menu bridge source changes; all other upstream files unmodified."
manifest["extension"] = {
    "version": "molecular-car-modelkit-3",
    "date": "2026-10-13",
    "license": "LGPL-2.1-or-later (same as the modified upstream library)",
    "sourcePatch": "../../jsmol-source/modelkit-native-bridge.patch",
    "sourcePatchSha256": hashlib.sha256(patch_path.read_bytes()).hexdigest(),
    "originalSources": "../../jsmol-source/upstream/",
    "modified": modified,
}
manifest_path.write_text(json.dumps(manifest, indent=2, ensure_ascii=False) + "\n", encoding="utf-8", newline="\n")
print(json.dumps(manifest["extension"], indent=2))
