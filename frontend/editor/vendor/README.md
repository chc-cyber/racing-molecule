# Editor local dependencies

- `ketcher/`: unmodified official Ketcher 3.18.0 standalone, main 2D runtime. Apache-2.0, LICENSE / NOTICE retained. Sources and SHA-256 are recorded in `ketcher-manifest.json`; no npm build or CDN required.
- `jsmol/`: unmodified JSmol 16.3.31 runtime for the independent ModelKit probe. LGPL-2.1-or-later plus included jQuery MIT. Original README, both licenses, uncompressed bootstrap, `NOTICE.md` and file hashes are retained. It is not loaded by the main workshop.
- `lucide.min.js`: local Lucide 1.8.0 icons, with full ISC / Feather notices in `LUCIDE-LICENSE.txt`.

See [THIRD_PARTY.md](../THIRD_PARTY.md) for component sources and scope.

## Editor transform controls

- Component: Three.js `TransformControls`, release `r180` / `0.180.0`.
- Source: https://raw.githubusercontent.com/mrdoob/three.js/r180/examples/jsm/controls/TransformControls.js
- The source file is unmodified and uses the editor's existing `three` import map.
- SHA-256: `621fb8894c0f748c2be2a6ff8a6bdd72c18a22b24ddfe21d506514544df5c84e`.
- License: MIT; copyright 2010–2025 three.js authors. See `THREE-LICENSE.txt`.
- Purpose: real atom/fragment translation and rotation handles. The attached object is a temporary pivot; completed gestures apply to the canonical molecular model once.

Runtime loads this local file and does not use a CDN or add a package dependency.
