"""Provisional nonreactive X-Au game interaction, separate from CHAu SAIP.

Atom-type UFF self vdW parameters use geometric mixing with the bundled Au
vdW entry. No heteroatom is relabelled C or added to the SAIP parameter table.
The perfect frozen Au lattice is shared with SAIP. Two universal potential
basis tables allow adding/changing parameters without building ten 3D grids.
"""
from __future__ import annotations
import hashlib
import json
import os
from pathlib import Path
import time
import numpy as np

from .element_models import PROVISIONAL_AU_ELEMENTS

PARAMETER_FILE = Path(__file__).with_name("data")/"hetero_au_v1.json"
VERSION = "tapered-uff-vdw-lattice-basis-v1"


class HeteroAuPotential:
    def __init__(self, atoms, substrate, *, accelerator="numpy"):
        from .periodic_au import EV_TO_KCAL_MOL
        raw = PARAMETER_FILE.read_bytes()
        self.parameters = json.loads(raw)
        self.parameter_sha256 = hashlib.sha256(raw).hexdigest()
        if self.parameters["schema"] != "hetero-au-uff-vdw-game-v1":
            raise ValueError("Unknown provisional Au parameter schema.")
        gold = [self.parameters["gold_rmin_angstrom"], self.parameters["gold_well_depth_kcal_mol"]]
        if (not np.isfinite(gold).all() or min(gold) <= 0
                or self.parameters["cutoff_angstrom"] != 16.
                or set(self.parameters["supported_elements"]) != set(PROVISIONAL_AU_ELEMENTS)):
            raise ValueError("Invalid provisional Au parameter values, element registry or cutoff.")
        self.indices = np.array([i for i, a in enumerate(atoms) if a["element"] not in ("C", "H")], dtype=np.int64)
        self.coefficients = np.empty((len(self.indices), 2))
        self.pairs = []
        for k, i in enumerate(self.indices):
            atom = atoms[i]
            if atom["element"] not in PROVISIONAL_AU_ELEMENTS:
                raise ValueError("Missing provisional Au parameters for "+atom["element"])
            x, d = atom.get("uff_vdw_rmin_angstrom"), atom.get("uff_vdw_well_depth_kcal_mol")
            if x is None or d is None or not np.isfinite([x, d]).all() or min(x, d) <= 0:
                raise ValueError("Hetero-Au game surface requires validated atom-type UFF vdW parameters.")
            rmin = np.sqrt(x*self.parameters["gold_rmin_angstrom"])
            epsilon = np.sqrt(d*self.parameters["gold_well_depth_kcal_mol"])
            self.coefficients[k] = [epsilon*rmin**12/EV_TO_KCAL_MOL, -2*epsilon*rmin**6/EV_TO_KCAL_MOL]
            self.pairs.append({"atom_id": int(i)+1, "element": atom["element"],
                               "rmin_angstrom": float(rmin), "well_depth_kcal_mol": float(epsilon),
                               "self_rmin_angstrom": float(x), "self_well_depth_kcal_mol": float(d)})
        self.substrate = substrate
        self.used_atom_parameters_sha256 = hashlib.sha256(json.dumps(self.pairs, sort_keys=True).encode()).hexdigest()
        self._grid = None
        self.cache_loaded = False
        self.cache_path = None
        self.cache_build_seconds = 0.
        self.accelerator = accelerator
        if len(self.indices) and substrate.use_cache:
            self._prepare_cache()

    def _prepare_cache(self):
        from .periodic_au import CACHE_Z_MIN, CACHE_Z_MAX, CACHE_BLEND_WIDTH, CUTOFF, _taper
        s = self.substrate
        settings = {"version": VERSION, "period": s.period.tolist(), "heights": s.layer_heights.tolist(),
                    "nx": s.grid_nx, "ny": s.grid_ny, "dz": s.grid_dz,
                    "z_min": CACHE_Z_MIN, "z_max": CACHE_Z_MAX, "blend": CACHE_BLEND_WIDTH,
                    "cutoff": CUTOFF, "basis": ["Tap/r^12", "Tap/r^6"]}
        identity = hashlib.sha256(json.dumps(settings, sort_keys=True).encode()).hexdigest()
        self.cache_path = s.cache_dir/f"hetero-basis-{identity}.npz"
        self.z_start = CACHE_Z_MIN-3*s.grid_dz
        nz = int(np.ceil((CACHE_Z_MAX-CACHE_Z_MIN+6*s.grid_dz)/s.grid_dz))+1
        shape = (2, nz, s.grid_nx, s.grid_ny)
        if self.cache_path.is_file():
            try:
                with np.load(self.cache_path, allow_pickle=False) as saved:
                    grid = saved["potential"]
                    if (grid.shape == shape and np.isfinite(grid).all()
                            and str(saved["settings"].item()) == json.dumps(settings, sort_keys=True)):
                        self._grid = grid
                        self.cache_loaded = True
                        return
            except (OSError, ValueError, EOFError, KeyError):
                pass
        started = time.perf_counter()
        grid = np.zeros(shape)
        dx, dy = s.period/[s.grid_nx, s.grid_ny]
        qx, qy = np.meshgrid(np.arange(-int(np.ceil(CUTOFF/dx)), int(np.ceil(CUTOFF/dx))+1),
                             np.arange(-int(np.ceil(CUTOFF/dy)), int(np.ceil(CUTOFF/dy))+1), indexing="ij")
        lateral = (qx*dx)**2+(qy*dy)**2
        bins = (qx % s.grid_nx)*s.grid_ny+(qy % s.grid_ny)
        for iz in range(nz):
            z = self.z_start+iz*s.grid_dz
            for layer, height in enumerate(s.layer_heights):
                r2 = lateral+(z-height)**2
                keep = r2 < CUTOFF**2
                taper = _taper(np.sqrt(r2[keep]))[0]
                r6 = r2[keep]**-3
                for code, kernel in enumerate((taper*r6*r6, taper*r6)):
                    field = np.bincount(bins[keep], weights=kernel,
                                        minlength=s.grid_nx*s.grid_ny).reshape(s.grid_nx, s.grid_ny)
                    sx, sy = (layer % 3)*s.grid_nx//2, (layer % 3)*s.grid_ny//6
                    grid[code, iz] += np.roll(field, (sx, sy), axis=(0, 1))
                    grid[code, iz] += np.roll(field, (sx+s.grid_nx//2, sy+s.grid_ny//2), axis=(0, 1))
        self._grid = grid
        self.cache_build_seconds = time.perf_counter()-started
        s.cache_dir.mkdir(parents=True, exist_ok=True)
        temporary = self.cache_path.with_name(self.cache_path.name+f".{os.getpid()}.{time.time_ns()}.tmp")
        try:
            with temporary.open("wb") as output:
                np.savez_compressed(output, potential=grid, settings=json.dumps(settings, sort_keys=True))
            os.replace(temporary, self.cache_path)
        finally:
            if temporary.exists():
                temporary.unlink()

    def _interpolate(self, positions, active):
        from .periodic_au import _cubic_weights
        s = self.substrate
        spacing = np.array([*s.period/[s.grid_nx, s.grid_ny], s.grid_dz])
        values = positions[self.indices[active]]
        scaled = values/spacing
        scaled[:, 2] = (values[:, 2]-self.z_start)/s.grid_dz
        base = np.floor(scaled).astype(np.int64)
        fraction = scaled-base
        wx, dx = _cubic_weights(fraction[:, 0]); wy, dy = _cubic_weights(fraction[:, 1]); wz, dz = _cubic_weights(fraction[:, 2])
        offset = np.arange(-1, 3)
        x = (base[:, 0, None]+offset) % s.grid_nx
        y = (base[:, 1, None]+offset) % s.grid_ny
        z = base[:, 2, None]+offset
        # Coefficients multiply scalar potentials, never independently fitted forces.
        c = self._grid[:, z[:, :, None, None], x[:, None, :, None], y[:, None, None, :]]
        c = np.einsum("bnijk,nb->nijk", c, self.coefficients[active])
        e = np.einsum("nijk,ni,nj,nk->n", c, wz, wx, wy)
        g = np.column_stack([np.einsum("nijk,ni,nj,nk->n", c, wz, dx, wy),
                             np.einsum("nijk,ni,nj,nk->n", c, wz, wx, dy),
                             np.einsum("nijk,ni,nj,nk->n", c, dz, wx, wy)])
        return e, -g/spacing

    def energy_and_forces_eV(self, positions, center_weights=None, center_dz=None):
        from .periodic_au import _cache_weight, _taper
        s = self.substrate
        if self.accelerator == "numba":
            from .hetero_au_kernel import energy_and_forces
            return energy_and_forces(positions, self.indices, self.coefficients, s._tile,
                self._grid if self._grid is not None else np.zeros((2, 1, 4, 6)),
                float(s.period[0]), float(s.period[1]), s.grid_dz,
                self.z_start if self._grid is not None else 0., self._grid is not None, center_weights, center_dz)
        w, dw = _cache_weight(positions[self.indices, 2]) if self._grid is not None else (np.zeros(len(self.indices)), np.zeros(len(self.indices)))
        ce, cf = np.zeros(len(self.indices)), np.zeros((len(self.indices), 3))
        # Interpolation requires interior z nodes; never index exact-only atoms.
        active = np.flatnonzero(w > 0)
        if len(active):
            ce[active], cf[active] = self._interpolate(positions, active)
        ee, ef = np.zeros(len(self.indices)), np.zeros((len(self.indices), 3))
        exact = np.flatnonzero(w < 1)
        if len(exact):
            ids, delta, radius = s._pairs(positions, self.indices[exact])
            lookup = np.full(len(positions), -1, dtype=np.int64); lookup[self.indices] = np.arange(len(self.indices))
            rows = lookup[ids]
            a, b = self.coefficients[rows].T
            r6 = radius**-6
            raw = a*r6*r6+b*r6
            t, dt = _taper(radius)
            derivative = (-12*a*r6*r6-6*b*r6)/radius
            force = -((dt*raw+t*derivative)/radius)[:, None]*delta
            ee += np.bincount(rows, weights=t*raw, minlength=len(self.indices))
            for axis in range(3):
                ef[:, axis] += np.bincount(rows, weights=force[:, axis], minlength=len(self.indices))
        result = np.zeros_like(positions)
        result[self.indices] = w[:, None]*cf+(1-w[:, None])*ef
        result[self.indices, 2] -= dw*(ce-ee)
        atomic_energy = w*ce+(1-w)*ee
        if center_weights is not None:
            result[self.indices] *= center_weights[self.indices, None]
            result[self.indices, 2] -= center_dz[self.indices]*atomic_energy
            atomic_energy *= center_weights[self.indices]
        return float(np.sum(atomic_energy)), result

    def metadata(self):
        return {"model": self.parameters["schema"], "algorithm_version": VERSION,
                "status": self.parameters["status"], "parameter_file": str(PARAMETER_FILE),
                "parameter_sha256": self.parameter_sha256, "source": self.parameters["source"],
                "mixing": self.parameters["mixing"], "pair_potential": self.parameters["pair_potential"],
                "used_pairs": self.pairs, "cutoff_angstrom": 16.,
                "used_atom_parameters_sha256": self.used_atom_parameters_sha256,
                "cache_basis": ["Tap/r^12", "Tap/r^6"], "cache_loaded": self.cache_loaded,
                "cache_file": str(self.cache_path) if self.cache_path else None,
                "cache_build_seconds": self.cache_build_seconds,
                "cache_grid_shape": list(self._grid.shape) if self._grid is not None else None,
                "limitations": self.parameters["limitations"]}
