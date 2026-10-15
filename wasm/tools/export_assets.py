"""Build browser-readable immutable tables from the native reference model.

Only a developer runs this script. Browser races never require Python or an
HTTP calculation API. Double arrays remain little-endian float64, gzip is
lossless, and the existing full directional SAIP term is not tabulated away.
"""
from __future__ import annotations
import gzip
import hashlib
import json
from pathlib import Path
import shutil
import sys

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT))
import numpy as np
from rdkit import Chem
from molecular_car.molecule import prepare_molecule
from molecular_car.periodic_au import create_periodic_surface, BONDED_NORMAL_POLICY
from molecular_car.hetero_au import HeteroAuPotential
from molecular_car.element_models import EFFECTIVE_POLARIZABILITY
from molecular_car import fields, engine, relaxation
from molecular_car.track import TRACK
from molecular_car.start_pose import POSE_VERSION, DEFAULT_START_POSE, HEADING_RATIO, PLANE_RATIO


def write_json(path, data):
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(data, ensure_ascii=False, indent=2, allow_nan=False) + "\n", encoding="utf-8")


def main():
    target = ROOT / "wasm" / "dist" / "browser" / "assets"
    target.mkdir(parents=True, exist_ok=True)
    example = ROOT / "examples" / "tetraphenylbenzene.mol"
    source_parameters = ROOT / ".local" / "wasm" / "lammps" / "potentials" / "CHAu.ILP"
    surface_options = {"parameter_path": source_parameters} if source_parameters.is_file() else {}
    prepared = prepare_molecule(example.read_text(), example.name)
    surface = create_periodic_surface(prepared.atoms, prepared.positions,
        bonds=prepared.bonds, normal_policy=BONDED_NORMAL_POLICY, **surface_options)
    # Universal hetero basis, independent of player molecular parameters.
    hetero_atoms = [{**prepared.atoms[0], "element": "O"}]
    hetero = HeteroAuPotential(hetero_atoms, surface)
    arrays = {}
    for name, values in (("tile", surface._tile), ("saip-grid", surface._grid),
                         ("hetero-grid", hetero._grid)):
        payload = np.asarray(values, dtype="<f8", order="C").tobytes()
        compressed = gzip.compress(payload, compresslevel=9, mtime=0)
        (target / f"{name}.bin.gz").write_bytes(compressed)
        arrays[name] = {"path": f"assets/{name}.bin.gz", "shape": list(values.shape),
                        "dtype": "float64-le", "encoding": "gzip", "bytes": len(payload),
                        "download_bytes": len(compressed), "sha256": hashlib.sha256(payload).hexdigest()}
    metadata = {
        "schema": "molecular-car-browser-assets-v1", "rdkit_version": "2025.09.3",
        "lammps_version": 20260930, "track": TRACK,
        # Shared native defaults; browser/transport.mjs uses the same game preset.
        "settings": engine.DEFAULT_SETTINGS, "drive": engine.DEFAULT_DRIVE,
        "timestep_fs": engine.TIMESTEP_FS, "frame_steps": engine.FRAME_STEPS,
        "max_steps": engine.MAX_STEPS, "seed": engine.DEFAULT_RANDOM_SEED,
        "start_pose": {"protocol_version": POSE_VERSION, "default": DEFAULT_START_POSE,
                       "heading_ratio": HEADING_RATIO, "plane_ratio": PLANE_RATIO},
        "world_bounds": engine.WORLD_BOUNDS, "effective_polarizabilities": EFFECTIVE_POLARIZABILITY,
        "field_parameters": {k: getattr(fields, k) for k in (
            "ADSORPTION_DEPTH", "ADSORPTION_HEIGHT", "ADSORPTION_ALPHA", "SURFACE_SPACING",
            "CORRUGATION_DECAY", "TIP_SOFTENING", "FIELD_SCALE", "GRADIENT_BASE_FIELD",
            "GRADIENT_FIELD_SLOPE", "GRADIENT_POLARIZATION_GAIN", "BOUNDARY_WIDTH", "BOUNDARY_STIFFNESS",
            "TRACK_X_MIN", "TRACK_X_MAX", "LANE_HALF_WIDTH")},
        "surface": {"mode": "periodic-au", "model": "periodic-saip-2",
            "normal_policy": BONDED_NORMAL_POLICY, "normal_area_epsilon": 1e-6,
            "period": surface.period.tolist(), "layer_heights_angstrom": surface.layer_heights.tolist(),
            "surface_period_xy_angstrom": surface.period.tolist(),
            "lattice_spacing_angstrom": surface.spacing,
            "grid_dz": surface.grid_dz, "z_start": surface._z_start,
            "parameters": surface._kernel_parameters.tolist(), "parameter_sha256": surface.parameter_sha256,
            "hetero_au": hetero.parameters, "arrays": arrays},
        "preparation": {"protocol_version": relaxation.PROTOCOL_VERSION,
            "force_tolerance": relaxation.GLOBAL_FORCE_TOLERANCE,
            "wall_time_limit_seconds": relaxation.WALL_TIME_LIMIT_SECONDS,
            "evaluation_limit": relaxation.FORCE_EVALUATION_LIMIT,
            "stages": relaxation.STAGES,
            "bfgs_rescue_options": relaxation.BFGS_OPTIONS},
    }
    for filename in ("benzene.mol", "tetraphenylbenzene.mol"):
        shutil.copyfile(ROOT / "examples" / filename, target / filename)
    write_json(target / "model.json", metadata)
    # Identical fixed coordinates exercise all terms, both field modes, negative
    # lattice queries, multi-cell molecules, and both cache blending edges.
    fixtures = []
    samples = [(example.name, example.read_text()),
               ("benzene.mol", (ROOT / "examples" / "benzene.mol").read_text())]
    for name, smiles in (("ethanol", "CCO"), ("fluorobenzene", "Fc1ccccc1"),
                         ("linear-alkyne", "c1ccccc1C#Cc1ccccc1"),
                         ("phosphine-oxide", "CP(=O)(C)C"),
                         ("pyridine", "c1ccncc1"), ("boric-acid", "B(O)(O)O"),
                         ("tetramethylsilane", "C[Si](C)(C)C"),
                         ("dimethyl-sulfide", "CSC"), ("chlorobenzene", "Clc1ccccc1"),
                         ("bromobenzene", "Brc1ccccc1"), ("iodobenzene", "Ic1ccccc1"),
                         ("isotope-methane", "[13CH3][2H]")):
        samples.append((name + ".mol", Chem.MolToMolBlock(Chem.MolFromSmiles(smiles))))
    for filename, content in samples:
        p = prepare_molecule(content, filename)
        s = create_periodic_surface(p.atoms, p.positions, bonds=p.bonds,
                                   normal_policy=BONDED_NORMAL_POLICY, **surface_options)
        (target / filename).write_text(content, encoding="utf-8")
        for shift in ((0., 0., 0.), (-17.3, 6.1, 0.), (11., -8., 2.15-p.positions[:, 2].min()),
                      (0., 0., 11.85-p.positions[:, 2].max())):
            positions = p.positions + shift
            energy = float(p.forcefield.CalcEnergy(positions.ravel().tolist()))
            force = -np.asarray(p.forcefield.CalcGrad(positions.ravel().tolist())).reshape(-1, 3)
            se, sf = s.energy_and_forces(positions)
            from molecular_car.external_fields import ExternalFieldEvaluator
            for continuum in (False, True):
              evaluator = ExternalFieldEvaluator([a["charge"] for a in p.atoms],
                  [a["polarizability"] for a in p.atoms], surface_enabled=continuum)
              for mode in (0, 1):
                fe, ff = evaluator.energy_and_forces(positions, mode, .75, 6., 0., 7., .08)
                surface_energy, surface_forces = (0., np.zeros_like(sf)) if continuum else (se, sf)
                fixtures.append({"filename": filename, "continuum": continuum, "positions": positions.ravel().tolist(),
                    "charges": [a["charge"] for a in p.atoms], "masses": [a["mass"] for a in p.atoms],
                    "mode": mode, "uff_energy": energy, "uff_forces": force.ravel().tolist(),
                    "surface_energy": surface_energy, "surface_forces": surface_forces.ravel().tolist(),
                    "field_energy": fe, "field_forces": ff.ravel().tolist(),
                    "energy": energy+surface_energy+fe, "forces": (force+surface_forces+ff).ravel().tolist(),
                    "normal_counts": s._normal_counts.tolist(), "normal_neighbours": s._normal_neighbours.tolist()})
    write_json(ROOT / "outputs" / "wasm-preparation" / "force-fixtures.json", fixtures)
    rejected = [{"filename": "empty.mol", "content": ""},
                {"filename": "wrong.xyz", "content": samples[0][1]},
                {"filename": "multiple.sdf", "content": samples[0][1]+"$$$$\n"+samples[0][1]+"$$$$\n"}]
    for name, smiles in (("radical", "[CH3]"), ("fragments", "CC.O"),
                         ("unsupported-element", "[Na+]"), ("three-coordinate-P", "CP(C)C")):
        rejected.append({"filename": name+".mol", "content": Chem.MolToMolBlock(Chem.MolFromSmiles(smiles))})
    write_json(ROOT / "outputs" / "wasm-preparation" / "admission-fixtures.json", rejected)
    print(json.dumps({"assets": str(target), "cases": len(fixtures),
        "download_mib": sum(a["download_bytes"] for a in arrays.values()) / 1024**2}, indent=2))


if __name__ == "__main__":
    main()
