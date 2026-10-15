"""Install pinned browser assets from the official npm registry, without Node."""

import io
from pathlib import Path
import tarfile
import urllib.request

VERSION = "0.180.0"
DESTINATION = Path(__file__).resolve().parents[1] / "frontend" / "vendor"
FILES = {
    "package/build/three.module.js": "three.module.js",
    "package/build/three.core.js": "three.core.js",
    "package/examples/jsm/controls/OrbitControls.js": "OrbitControls.js",
    "package/LICENSE": "THREE-LICENSE.txt",
}


def main():
    url = f"https://registry.npmjs.org/three/-/three-{VERSION}.tgz"
    with urllib.request.urlopen(url, timeout=60) as response:
        archive = response.read()
    DESTINATION.mkdir(parents=True, exist_ok=True)
    with tarfile.open(fileobj=io.BytesIO(archive), mode="r:gz") as package:
        # Read exact known files; never extract an archive's arbitrary paths.
        for source, target in FILES.items():
            member = package.extractfile(source)
            if member is None:
                raise RuntimeError(f"Missing browser dependency: {source}")
            (DESTINATION / target).write_bytes(member.read())
    print(f"Three.js {VERSION}: browser assets ready")


if __name__ == "__main__":
    main()
