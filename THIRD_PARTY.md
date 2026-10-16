# License scope

Project-owned code is MIT (see LICENSE). Third-party files and modifications retain their original licenses. The linked browser WASM module contains LAMMPS and is distributed under GPL-2.0, not under MIT. Corresponding source archives and build inputs accompany the GitHub release; see RELEASING.md.

# Third-party components

The prototype uses these separately distributed open-source components:

| Component | Version used | License / source |
|---|---|---|
| LAMMPS | 30 Sep 2026 (20260930) | GPLv2; installed separately, no binaries copied into this project. [Official source](https://github.com/lammps/lammps) |
| RDKit | 2025.9.3 | BSD 3-Clause. [Official source and license](https://github.com/rdkit/rdkit) |
| NumPy | 1.26.4 | BSD 3-Clause. [Official source](https://github.com/numpy/numpy) |
| Numba | 0.61.2 | BSD 2-Clause. Compiles the continuous external fields and frozen-Au environment calculation on the CPU, installed in `.venv`. [Official source](https://github.com/numba/numba/tree/0.61.2) |
| Three.js | 0.180.0 | MIT. Pinned browser files in `frontend/vendor/`; copyright notice in `THREE-LICENSE.txt`. [Official source](https://github.com/mrdoob/three.js) |

The Python dependencies are installed inside `.venv`; their distributions retain
their own license notices. RDKit also installs Pillow as a transitive dependency.
The browser files are downloaded from the official npm registry by
`scripts/fetch_frontend.py`, and work locally without a CDN connection.

Example MOL files were generated locally from molecular graphs with RDKit.
The four-arm tetraphenylbenzene is a demonstration structure, not an experimental
Nanocar Race molecule. The molecular workshop uses Ketcher/Indigo and JSmol;
its exact versions, licenses, local patches and notices are recorded in
`frontend/editor/THIRD_PARTY.md` and `frontend/editor/vendor/`.

## Explicit gold reference experiment

`scripts/verify_au_surface.py` reads these potential files from the separately
installed LAMMPS `Potentials/` directory. They are not copied or redistributed by
this project. Each run records the filenames, SHA-256 hashes and source headers.

- `Au_u3.eam`: S. M. Foiles, M. I. Baskes and M. S. Daw, *Physical Review B*
  **33**, 7983–7991 (1986), [DOI](https://doi.org/10.1103/PhysRevB.33.7983).
  [NIST repository entry](https://www.ctcms.nist.gov/potentials/entry/1986--Foiles-S-M-Baskes-M-I-Daw-M-S--Au/).
- `CHAu.ILP`: W. Ouyang, O. Hod and R. Guerra, *Registry-Dependent Potential for
  Interfaces of Gold with Graphitic Systems*, *Journal of Chemical Theory and
  Computation* **17**, 7215–7223 (2021),
  [DOI](https://doi.org/10.1021/acs.jctc.1c00622).
  [LAMMPS SAIP documentation](https://docs.lammps.org/pair_saip_metal.html).

This benchmark tests a rigid planar benzene molecule. It does not establish a
force field for arbitrary player-designed molecules or replace the game's
continuum surface. Optional Matplotlib plots are analysis outputs, not game
runtime dependencies.

## Frozen periodic Au environment

`molecular_car/periodic_au.py` evaluates the same published SAIP mathematical
potential and local-normal derivatives against six frozen, periodic Au layers.
The definitions were checked against the separately installed LAMMPS baseline
implementation: [SAIP](https://github.com/lammps/lammps/blob/develop/src/INTERLAYER/pair_saip_metal.cpp)
and [local normals](https://github.com/lammps/lammps/blob/develop/src/INTERLAYER/pair_ilp_graphene_hbn.cpp).
The C/H parameters are read from the installed `CHAu.ILP`, rather than copied
into this repository. Numerical potential tables and compiled-function caches
are generated locally and ignored by Git. Numba installs llvmlite and its
LLVM runtime as separately licensed dependencies; their distribution notices
remain in the virtual environment. No LAMMPS or LLVM binaries are shipped here.

`molecular_car/fourier_au.py` and `fourier_au_kernel.py` approximate that same
classical frozen-SAIP potential using a bounded spatial Fourier expansion and
Hermite interpolation of its full local-normal response. Outside the bounded
domain they blend conservatively to the validated direction lattice sum. The
tables are generated from the installed parameters locally; this does not add
a new parameter set or redistribute potential files. Validation includes
native LAMMPS forces and complete Cartesian energy derivatives.

## Browser WASM development build

The browser build uses the same LAMMPS 30 Sep 2026 and RDKit 2025.09.3 source
versions as the native reference. `wasm/dependencies.json` pins the upstream
source objects and Emscripten SDK. `wasm/tools/build.py` links the selected
LAMMPS and RDKit C++ libraries with this project's conservative force kernels.
The resulting browser package redistributes compiled LAMMPS and RDKit code;
the earlier statements about separately installed binaries describe the native
version only. The compiled module is subject to LAMMPS's GPLv2 terms. A public
release must include the corresponding source and reproducible build inputs.

- Emscripten 6.0.12: MIT / University of Illinois-NCSA; runtime/library notices
  from the official SDK are included in the browser package.
- Boost 1.85.0: Boost Software License 1.0; selected headers compile into RDKit.
- zlib: the Emscripten SDK's pinned port, with its zlib license retained.
- better-enums 0.11.3: BSD 2-Clause; pinned to the commit in the dependency file.
- Catch2 3.4.0: build configuration dependency only, not linked into this game.

Float64 scalar Au potential tables are generated by `export_assets.py` from
the pinned LAMMPS `potentials/CHAu.ILP`, the bundled six-layer Au reference,
and the separately versioned hetero-Au game model. Browser assets include
their dimensions, model parameters and SHA-256 hashes. They retain the SAIP
paper attribution above. The current build includes the versioned Fourier proxy v2; it does not introduce a Python runtime.

Compiler downloads and upstream development trees live in ignored
`.local/wasm/`; generated websites live in ignored `wasm/dist/`. Neither tree
belongs in a normal source checkout. The static website package contains
only the required browser runtime, shared renderer, model tables and notices.

## Brotli resource decoder

- brotli-wasm 3.0.1, Apache-2.0: https://github.com/httptoolkit/brotli-wasm
- Browser wrapper and losslessly gzipped WASM: `frontend/browser/vendor/brotli/`. License and pinned npm source are preserved there.
- Used only to decode losslessly compressed simulation tables; it does not alter the numerical model.
