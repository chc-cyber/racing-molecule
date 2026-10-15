"""Finite release races for five generated/imported C/H test structures.

Use --prepare-only while production fixes are in progress. Actual MD must be
started after the release source is frozen. Analysis dependencies are never
added to this process. Every production race uses the real game Simulation.
"""
from __future__ import annotations

import os
for variable in ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS", "NUMEXPR_NUM_THREADS"):
    os.environ[variable] = "1"

import argparse
from concurrent.futures import ProcessPoolExecutor, as_completed
from contextlib import contextmanager
import csv
from datetime import datetime, timezone
import hashlib
import json
from pathlib import Path
import re
import sys
import time
import traceback

ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / "outputs/md-release-validation"
INPUTS = OUT / "inputs"
RACES = OUT / "races"
sys.path.insert(0, str(ROOT))
import numpy as np
from rdkit import Chem, rdBase
from rdkit.Chem import AllChem, rdMolDescriptors

SEEDS = (98245165, 2718281, 1234567)
EMBED_SEED = 20261012
FIXTURE_IDS = ("test-car-52", "benzene-12", "fullerene-c60", "ethynyl-rigid", "folded-eicosane")
PRODUCTION_SOURCES = (
    "molecular_car/engine.py", "molecular_car/molecule.py", "molecular_car/track.py",
    "molecular_car/periodic_au.py", "molecular_car/periodic_au_kernel.py",
    "molecular_car/fourier_au.py", "molecular_car/fourier_au_kernel.py",
    "molecular_car/data/au111_reference.json",
    "molecular_car/fields.py", "molecular_car/fields_kernel.py", "molecular_car/external_fields.py",
    "scripts/lammps_runtime.py", "scripts/validate_release_races.py",
)


def write_json(path, value):
    temporary = path.with_suffix(path.suffix + ".tmp")
    temporary.write_text(json.dumps(value, ensure_ascii=False, indent=2, allow_nan=False) + "\n", encoding="utf-8")
    temporary.replace(path)


def checksum(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def fullerene():
    """Truncate an icosahedron; the 30 surviving old edges are double bonds.

    Twelve degree-five pentagons have single bonds. The twenty resulting
    hexagons alternate single and double bonds. This is a C60 Kekule graph,
    with exactly sixty degree-three carbons and no hydrogen.
    """
    phi = (1. + np.sqrt(5.)) / 2.
    vertices = np.array([(0., a, b * phi) for a in (-1., 1.) for b in (-1., 1.)]
                        + [(a, b * phi, 0.) for a in (-1., 1.) for b in (-1., 1.)]
                        + [(a * phi, 0., b) for a in (-1., 1.) for b in (-1., 1.)])
    distances = np.linalg.norm(vertices[:, None] - vertices[None, :], axis=2)
    original_edges = [(i, j) for i in range(12) for j in range(i + 1, 12) if np.isclose(distances[i, j], 2.)]
    directed = sorted([(i, j) for i, j in original_edges] + [(j, i) for i, j in original_edges])
    positions = np.array([(2. * vertices[i] + vertices[j]) / 3. for i, j in directed])
    shortest = np.min(np.linalg.norm(positions[:, None] - positions[None, :], axis=2)
                      + np.eye(60) * 100.)
    graph = Chem.RWMol()
    for _ in range(60):
        graph.AddAtom(Chem.Atom("C"))
    for i in range(60):
        for j in range(i + 1, 60):
            if np.isclose(np.linalg.norm(positions[i] - positions[j]), shortest):
                old_edge = directed[i] == tuple(reversed(directed[j]))
                graph.AddBond(i, j, Chem.BondType.DOUBLE if old_edge else Chem.BondType.SINGLE)
    molecule = graph.GetMol()
    Chem.SanitizeMol(molecule)
    if (molecule.GetNumBonds() != 90 or any(a.GetDegree() != 3 for a in molecule.GetAtoms())
            or rdMolDescriptors.CalcMolFormula(molecule) != "C60"):
        raise AssertionError("The generated fullerene graph is not closed-shell C60")
    conformer = Chem.Conformer(60)
    conformer.Set3D(True)
    positions *= 1.44 / shortest
    for i, position in enumerate(positions):
        conformer.SetAtomPosition(i, tuple(position))
    molecule.AddConformer(conformer)
    return molecule


def embed_smiles(smiles, seed):
    molecule = Chem.AddHs(Chem.MolFromSmiles(smiles))
    parameters = AllChem.ETKDGv3()
    parameters.randomSeed = seed
    parameters.numThreads = 1
    with rdBase.BlockLogs():
        if AllChem.EmbedMolecule(molecule, parameters) != 0:
            raise ValueError(f"ETKDG could not embed {smiles}")
    return molecule


def optimize(molecule):
    with rdBase.BlockLogs():
        if not AllChem.UFFHasAllMoleculeParams(molecule):
            raise ValueError("A generated fixture is missing UFF parameters")
        field = AllChem.UFFGetMoleculeForceField(molecule, vdwThresh=100., ignoreInterfragInteractions=False)
        field.Initialize()
        status = field.Minimize(maxIts=3000, forceTol=.0005, energyTol=1e-9)
        energy = float(field.CalcEnergy())
    if status != 0 or not np.isfinite(energy):
        raise ValueError("Fixture UFF minimization did not converge")
    return energy


def orient_flat(molecule):
    """Fix a reproducible imported 3D pose: longest principal axis along x.

    The game's own importer preserves this orientation and still performs its
    normal UFF and surface preparation. No rotational constraint is applied
    during MD.
    """
    conformer = molecule.GetConformer()
    positions = np.asarray(conformer.GetPositions())
    masses = np.array([a.GetMass() for a in molecule.GetAtoms()])
    centered = positions - np.average(positions, axis=0, weights=masses)
    _, axes = np.linalg.eigh((centered * masses[:, None]).T @ centered)
    rotation = axes[:, ::-1]
    if np.linalg.det(rotation) < 0:
        rotation[:, -1] *= -1
    oriented = centered @ rotation
    for index, position in enumerate(oriented):
        conformer.SetAtomPosition(index, tuple(position))


def legacy_neighbour_diagnostic(molecule):
    positions = np.asarray(molecule.GetConformer().GetPositions())
    elements = [a.GetSymbol() for a in molecule.GetAtoms()]
    cutoffs = {("C", "C"): 2., ("C", "H"): 1.5, ("H", "C"): 1.5, ("H", "H"): 1.2}
    rcut = np.array([[cutoffs[a, b] for b in elements] for a in elements])
    square = np.sum((positions[:, None] - positions[None, :])**2, axis=2)
    mask = (square > 0.) & (square < rcut**2)
    counts = mask.sum(axis=1)
    return {"legacy_distance_normal_max_neighbours": int(counts.max()),
            "legacy_distance_normal_atoms_over_3": (np.flatnonzero(counts > 3) + 1).tolist(),
            "chemical_degree_max": max(a.GetDegree() for a in molecule.GetAtoms())}


def prepare_fixtures():
    INPUTS.mkdir(parents=True, exist_ok=True)
    entries = []
    for fixture, original, label in (
        ("test-car-52", "tetraphenylbenzene.mol", "Four-arm test molecule"),
        ("benzene-12", "benzene.mol", "Benzene"),
    ):
        target = INPUTS / (fixture + ".mol")
        target.write_bytes((ROOT / "examples" / original).read_bytes())
        molecule = Chem.MolFromMolFile(str(target), removeHs=False)
        entries.append((fixture, label, molecule, {"source": "existing locally generated game fixture", "file": target.name}))
    for fixture, label, molecule, detail in (
        ("fullerene-c60", "C60 cage", fullerene(), {"source": "procedural truncated icosahedron; valid C60 Kekule bond graph"}),
        ("ethynyl-rigid", "Ethynyl-linked three phenyl rings",
         embed_smiles("c1ccccc1C#Cc1ccc(C#Cc2ccccc2)cc1", EMBED_SEED),
         {"source": "RDKit ETKDGv3", "smiles": "c1ccccc1C#Cc1ccc(C#Cc2ccccc2)cc1", "embedding_seed": EMBED_SEED}),
    ):
        detail["UFF_optimized_energy_kcal_mol"] = optimize(molecule)
        orient_flat(molecule)
        target = INPUTS / (fixture + ".mol")
        molecule.SetProp("_Name", label)
        target.write_text(Chem.MolToMolBlock(molecule), encoding="utf-8")
        detail["file"] = target.name
        entries.append((fixture, label, molecule, detail))
    # Select a compact, bent UFF minimum from a small deterministic conformer
    # set. Even ordinary sp3 C already has four short bonded neighbours under
    # the old distance rule; no artificial overlaps are introduced.
    candidates = []
    for seed in (EMBED_SEED, EMBED_SEED + 1, EMBED_SEED + 2, EMBED_SEED + 3):
        molecule = embed_smiles("C" * 20, seed)
        energy = optimize(molecule)
        carbon_positions = molecule.GetConformer().GetPositions()[:20]
        end_distance = float(np.linalg.norm(carbon_positions[-1] - carbon_positions[0]))
        candidates.append((end_distance, seed, energy, molecule))
    end_distance, seed, energy, molecule = min(candidates, key=lambda item: item[0])
    orient_flat(molecule)
    molecule.SetProp("_Name", "Folded n-eicosane C20H42")
    target = INPUTS / "folded-eicosane.mol"
    target.write_text(Chem.MolToMolBlock(molecule), encoding="utf-8")
    entries.append(("folded-eicosane", "Folded n-eicosane", molecule,
                    {"source": "RDKit ETKDGv3; most compact of four relaxed seeds", "smiles": "C" * 20,
                     "embedding_seed": seed, "candidate_seeds": [c[1] for c in candidates],
                     "carbon_end_to_end_angstrom": end_distance, "UFF_optimized_energy_kcal_mol": energy,
                     "file": target.name}))
    manifest = {"status": "passed", "fixture_generator": "scripts/validate_release_races.py",
                "rdkit_version": rdBase.rdkitVersion, "fixtures": []}
    for fixture, label, molecule, detail in entries:
        path = INPUTS / detail["file"]
        positions = molecule.GetConformer().GetPositions()
        record = {"id": fixture, "label": label, "formula": rdMolDescriptors.CalcMolFormula(molecule),
                  "atoms": molecule.GetNumAtoms(), "file_sha256": checksum(path),
                  "elements": [a.GetSymbol() for a in molecule.GetAtoms()],
                  "masses_u": [a.GetMass() for a in molecule.GetAtoms()],
                  "bonds_zero_based": [[b.GetBeginAtomIdx(), b.GetEndAtomIdx()] for b in molecule.GetBonds()],
                  "span_angstrom": np.ptp(positions, axis=0).tolist(), **detail,
                  **legacy_neighbour_diagnostic(molecule)}
        if max(record["span_angstrom"]) >= 45 or record["atoms"] > 256:
            raise AssertionError("Generated representative violates current import limits")
        manifest["fixtures"].append(record)
    flexible = next(row for row in manifest["fixtures"] if row["id"] == "folded-eicosane")
    if flexible["legacy_distance_normal_max_neighbours"] <= 3:
        raise AssertionError("The flexible diagnostic does not reproduce the old neighbour limit")
    write_json(INPUTS / "manifest.json", manifest)
    return manifest


def source_hashes(manifest):
    files = list(PRODUCTION_SOURCES) + ["outputs/md-release-validation/inputs/manifest.json"]
    files += ["outputs/md-release-validation/inputs/" + f["file"] for f in manifest["fixtures"]]
    return {name: checksum(ROOT / name) for name in files}


def preparation_log(log_path):
    content = log_path.read_text(encoding="utf-8", errors="replace") if log_path.exists() else ""
    criterion = re.search(r"Stopping criterion\s*=\s*([^\r\n]+)", content)
    norms = re.search(r"Force two-norm initial, final\s*=\s*\r?\n\s*(\S+)\s+(\S+)", content)
    return {"joint_relaxation_stopping_criterion": criterion.group(1).strip() if criterion else None,
            "joint_relaxation_force_norm_initial_final": [float(norms.group(1)), float(norms.group(2))] if norms else None,
            "joint_relaxation_protocol": "Unmodified game's FIRE minimize and timestep reset; detailed native log retained"}


def principal_and_rotation(positions, masses, reference):
    center = np.average(positions, axis=0, weights=masses)
    centered = positions - center
    eigenvalues, axes = np.linalg.eigh((centered * masses[:, None]).T @ centered / masses.sum())
    heading = float(np.degrees(np.arccos(np.clip(abs(axes[0, -1]), 0., 1.)))) if eigenvalues[-1] > 1.15 * eigenvalues[-2] else None
    tilt = float(np.degrees(np.arccos(np.clip(abs(axes[2, 0]), 0., 1.)))) if eigenvalues[1] > 1.2 * max(eigenvalues[0], 1e-10) else None
    u, _, vh = np.linalg.svd((reference * masses[:, None]).T @ centered)
    if np.linalg.det(u @ vh) < 0:
        u[:, -1] *= -1
    rotation = u @ vh
    fitted_angle = float(np.degrees(np.arccos(np.clip((np.trace(rotation) - 1.) / 2., -1., 1.))))
    residual = reference @ rotation - centered
    fitted_rmsd = float(np.sqrt(np.average(np.sum(residual**2, axis=1), weights=masses)))
    return {"long_axis_track_angle_degrees": heading, "plane_normal_tilt_degrees": tilt,
            "reference_fit_rotation_angle_degrees": fitted_angle, "reference_fit_RMSD_angstrom": fitted_rmsd}


def integrated_masses(simulation):
    n = simulation.n
    rmass = simulation.lmp.extract_atom("rmass")
    if rmass:
        # The game keeps atom sorting disabled; still reorder by stable IDs.
        tags = np.ctypeslib.as_array(simulation.lmp.extract_atom("id"), shape=(n,))
        result = np.empty(n)
        result[np.asarray(tags, dtype=int) - 1] = np.ctypeslib.as_array(rmass, shape=(n,))
        return result
    types = np.ctypeslib.as_array(simulation.lmp.gather_atoms("type", 0, 1), shape=(n,))
    table = simulation.lmp.extract_atom("mass")
    return np.array([table[int(kind)] for kind in types])


def run_trial(config, hashes, fixture):
    import molecular_car.engine as engine
    from scripts.lammps_runtime import open_lammps
    key = f"{config['phase']}_{fixture['id']}_seed{config['seed']}_{config['duration_ps']:g}ps"
    folder = RACES / "trials"
    folder.mkdir(parents=True, exist_ok=True)
    checkpoint, coordinates_path, log_path = (folder / (key + suffix) for suffix in (".json", ".npz", ".log"))
    if checkpoint.exists() and coordinates_path.exists():
        prior = json.loads(checkpoint.read_text(encoding="utf-8"))
        if (prior.get("status") == "passed" and prior["source_sha256"] == hashes and prior["config"] == config
                and checksum(coordinates_path) == prior.get("coordinates_sha256")):
            return prior
    simulation = None
    original_context = engine.open_lammps

    @contextmanager
    def logged_context():
        with open_lammps() as (lmp, runtime):
            lmp.command('log "' + log_path.as_posix() + '"')
            yield lmp, runtime

    engine.open_lammps = logged_context
    began = time.perf_counter()
    samples, coordinates, velocity_samples = [], [], []
    result = {"key": key, "config": config, "fixture": fixture["id"], "status": "running", "source_sha256": hashes}
    try:
        content = (INPUTS / fixture["file"]).read_text(encoding="utf-8")
        simulation = engine.Simulation(content, fixture["file"],
            {"surface": "periodic-au", "temperature": 600., "damping_fs": 350., "corrugation": .08},
            random_seed=config["seed"], drive={"mode": "gradient", "strength": 0.})
        result["runtime"] = simulation.runtime
        prepared_seconds = time.perf_counter() - began
        masses = integrated_masses(simulation)
        mass_error = float(np.max(np.abs(masses - simulation.masses)))
        if mass_error > 1e-8:
            raise AssertionError(f"LAMMPS integrates masses inconsistent with the molecular input: {mass_error}")
        initial_pose = simulation.coordinates()
        sample_steps = round(config["sample_ps"] * 1000 / engine.TIMESTEP_FS)
        warmup_steps = round(config["warmup_ps"] * 1000 / engine.TIMESTEP_FS)
        for offset in range(0, warmup_steps, sample_steps):
            simulation.advance(min(sample_steps, warmup_steps - offset))
        simulation.update({"drive": {"mode": "gradient", "strength": config["strength"]}})
        # One explicit setup refresh after *every* phase's strength assignment:
        # includes diffusion (zero -> zero), preserving a common RNG cadence.
        simulation.lmp.command("run 0")
        simulation._check_callback()
        frame = simulation.frame()
        initial_step = frame["step"]
        reference_positions = simulation.coordinates()
        reference = reference_positions - np.average(reference_positions, axis=0, weights=simulation.masses)
        bond_ids = np.array([[b["a"] - 1, b["b"] - 1] for b in simulation.prepared.bonds])
        bond_low, bond_high, force_max, conservative_max = float("inf"), 0., 0., 0.

        def append(state):
            nonlocal bond_low, bond_high, force_max, conservative_max
            positions = np.asarray(state["positions"]).reshape(-1, 3)
            forces = np.ctypeslib.as_array(simulation.lmp.gather_atoms("f", 1, 3), shape=(simulation.n * 3,)).reshape(-1, 3)
            velocities = np.ctypeslib.as_array(simulation.lmp.gather_atoms("v", 1, 3), shape=(simulation.n * 3,)).reshape(-1, 3).copy()
            conservative = simulation._molecular_forces
            if not np.isfinite(positions).all() or not np.isfinite(forces).all() or not np.isfinite(conservative).all() or not np.isfinite(state["energy"]):
                raise AssertionError("Non-finite sampled coordinates, energy or forces")
            independent_com = np.average(positions, axis=0, weights=masses)
            if not np.allclose(independent_com, state["center"], rtol=0., atol=1e-10):
                raise AssertionError("Reported COM does not match native integrated masses")
            lengths = np.linalg.norm(positions[bond_ids[:, 0]] - positions[bond_ids[:, 1]], axis=1)
            bond_low, bond_high = min(bond_low, float(lengths.min())), max(bond_high, float(lengths.max()))
            if bond_low < .6 or bond_high > 2.5:
                raise AssertionError("Sampled CH bond geometry left [0.6, 2.5] angstrom diagnostic limits")
            force_max = max(force_max, float(np.linalg.norm(forces, axis=1).max()))
            conservative_max = max(conservative_max, float(np.linalg.norm(conservative, axis=1).max()))
            coordinates.append(positions.copy())
            velocity_samples.append(velocities)
            samples.append({"time_ps": (state["step"] - initial_step) * engine.TIMESTEP_FS / 1000.,
                "step": state["step"], "center_nm": (independent_com / 10.).tolist(), "energy_kcal_mol": state["energy"],
                "internal_temperature_K": state["internal_temperature"], "lowest_atom_z_nm": float(positions[:, 2].min()) / 10.,
                "highest_atom_z_nm": float(positions[:, 2].max()) / 10., **principal_and_rotation(positions, masses, reference)})

        append(frame)
        propagation_seconds = 0.
        propagation_steps = 0
        production_steps = round(config["duration_ps"] * 1000 / engine.TIMESTEP_FS)
        for offset in range(0, production_steps, sample_steps):
            count = min(sample_steps, production_steps - offset)
            started = time.perf_counter()
            frame = simulation.advance(count)
            propagation_seconds += time.perf_counter() - started
            propagation_steps += count
            append(frame)
            if config["phase"] == "race" and frame["center"][0] >= 200.:
                break
        time_array = np.array([s["time_ps"] for s in samples])
        centers = np.array([s["center_nm"] for s in samples])
        relative = centers - centers[0]
        reached = bool(config["phase"] == "race" and centers[-1, 0] >= 20.)
        lifted_start = None
        longest_lift = 0.
        for sample in samples:
            if sample["lowest_atom_z_nm"] > .8:
                if lifted_start is None:
                    lifted_start = sample["time_ps"]
                longest_lift = max(longest_lift, sample["time_ps"] - lifted_start)
            else:
                lifted_start = None
        result.update(status="passed", runtime=simulation.runtime, atoms=simulation.n, formula=fixture["formula"],
            prepared_atom_notes=simulation.prepared.notes, native_mass_max_error_u=mass_error,
            initial_pose_angstrom=initial_pose.tolist(), prepared_seconds=prepared_seconds,
            preparation={"game_metadata": simulation.runtime.get("preparation"), **preparation_log(log_path)}, production_steps=propagation_steps,
            production_wall_seconds=propagation_seconds, production_us_per_step=propagation_seconds / propagation_steps * 1e6,
            bond_length_range_angstrom=[bond_low, bond_high], sampled_total_force_max_kcal_mol_angstrom=force_max,
            sampled_conservative_force_max_kcal_mol_angstrom=conservative_max,
            outcome="finish" if reached else ("diffusion_window_complete" if config["phase"] == "diffusion" else "time_limit"),
            finish_sample_time_ps=float(time_array[-1]) if reached else None,
            observed_time_ps=float(time_array[-1]), displacement_nm=relative[-1].tolist(),
            final_reference_fit_rotation_degrees=samples[-1]["reference_fit_rotation_angle_degrees"],
            maximum_lateral_excursion_nm=float(np.abs(relative[:, 1]).max()),
            longest_all_atoms_above_0_8_nm_ps=longest_lift, sustained_lift=bool(longest_lift >= 5. - 1e-8),
            sampled_atom_height_range_nm=[min(s["lowest_atom_z_nm"] for s in samples), max(s["highest_atom_z_nm"] for s in samples)],
            samples=samples)
        np.savez_compressed(coordinates_path, time_ps=time_array, positions_angstrom=np.stack(coordinates),
                            velocities_angstrom_fs=np.stack(velocity_samples), masses_u=masses, bonds_zero_based=bond_ids)
        result["coordinates_file"] = coordinates_path.relative_to(OUT).as_posix()
        result["coordinates_sha256"] = checksum(coordinates_path)
    except BaseException as error:
        result.update(status="failed", error=str(error), traceback=traceback.format_exc(), samples=samples)
    finally:
        if simulation is not None:
            simulation.close()
        engine.open_lammps = original_context
        if log_path.exists():
            result["preparation"] = {"game_metadata": result.get("runtime", {}).get("preparation"), **preparation_log(log_path)}
        result["native_log"] = log_path.relative_to(OUT).as_posix()
        result["total_wall_seconds"] = time.perf_counter() - began
        write_json(checkpoint, result)
    return result


def summary(rows):
    groups = {}
    for row in rows:
        groups.setdefault((row["config"]["phase"], row["fixture"]), []).append(row)
    result = []
    for (phase, fixture), trials in sorted(groups.items()):
        passed = [r for r in trials if r["status"] == "passed"]
        times = [r["finish_sample_time_ps"] for r in passed if r["finish_sample_time_ps"] is not None]
        result.append({"phase": phase, "fixture": fixture, "trials": len(trials), "passed": len(passed),
            "finishers": len(times), "observed_finish_times_ps": times,
            "mean_forward_displacement_nm": float(np.mean([r["displacement_nm"][0] for r in passed])) if passed else None,
            "final_transverse_RMS_nm": float(np.sqrt(np.mean([r["displacement_nm"][1]**2 for r in passed]))) if passed else None,
            "sustained_lift_trials": sum(r["sustained_lift"] for r in passed), "failed_trial_keys": [r["key"] for r in trials if r["status"] != "passed"]})
    return result


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--prepare-only", action="store_true")
    parser.add_argument("--workers", type=int, default=2)
    parser.add_argument("--duration-ps", type=float, default=100.)
    parser.add_argument("--diffusion-ps", type=float, default=10.)
    parser.add_argument("--fixtures", nargs="+", choices=FIXTURE_IDS)
    args = parser.parse_args()
    if args.prepare_only:
        manifest = prepare_fixtures()
        print(json.dumps({"manifest": str(INPUTS / "manifest.json"), "fixtures": [{k: row[k] for k in ("id", "formula", "atoms", "legacy_distance_normal_max_neighbours")} for row in manifest["fixtures"]]}))
        return 0
    if not 1 <= args.workers <= 3 or not 1 <= args.duration_ps <= 100 or not 0 <= args.diffusion_ps <= 20:
        parser.error("Use 1–3 workers, 1–100 ps race and 0–20 ps diffusion")
    from molecular_car.engine import TIMESTEP_FS
    from molecular_car.track import TRACK
    from molecular_car.periodic_au import parameter_file
    if TIMESTEP_FS != .5 or TRACK["finish_x"] != 200. or TRACK["half_width"] != 100.:
        raise ValueError("Release validation requires the current 0.5 fs and 20 nm square track")
    manifest = json.loads((INPUTS / "manifest.json").read_text(encoding="utf-8"))
    fixtures = [f for f in manifest["fixtures"] if args.fixtures is None or f["id"] in args.fixtures]
    hashes = source_hashes(manifest)
    phases = [("race", .75, args.duration_ps)] + ([("diffusion", 0., args.diffusion_ps)] if args.diffusion_ps else [])
    tasks = [(fixture, {"phase": phase, "strength": strength, "seed": seed, "duration_ps": duration,
                       "warmup_ps": 5., "sample_ps": .1}) for phase, strength, duration in phases for fixture in fixtures for seed in SEEDS]
    RACES.mkdir(parents=True, exist_ok=True)
    report = {"status": "running", "created_utc": datetime.now(timezone.utc).isoformat(), "source_sha256": hashes,
        "parameter_file_sha256": checksum(parameter_file()), "fixtures": fixtures, "planned_trials": len(tasks),
        "protocol": {"surface": "complete existing frozen periodic Au SAIP with the release normal-policy metadata recorded per trial",
            "temperature_K": 600., "damping_fs": 350., "timestep_fs": .5, "race_strength": .75,
            "track_nm": [20., 20.], "race_limit_ps": args.duration_ps, "diffusion_window_ps": args.diffusion_ps,
            "seeds": SEEDS, "warmup_ps": 5., "sample_ps": .1,
            "random_number_protocol": "Original game constructor and bath; 50 runs of 200 steps for zero-drive warmup; strength assignment; exactly one run 0 refresh for both phases; runs of 200 steps until finish/window. No bath reset or extra force setup.",
            "preparation": "Original imported geometry -> game UFF relaxation -> unchanged native joint FIRE surface relaxation -> zero-drive thermal preparation",
            "preparation_convergence": "Native log stopping criterion distinguishes energy tolerance, force tolerance and limits. force_tolerance_met=false alone does not mean no convergence criterion was met.",
            "energy_force_checks": "Game checks finite conservative energy/force every callback; independent sampled total LAMMPS forces include Langevin. No energy-conservation inference under Langevin/driving.",
            "bond_height_checks": "CH bonds and mass-weighted COM audited every 0.1 ps. Overall lift: every atom above 0.8 nm for at least 5 ps.",
            "orientation": "Mass-weighted principal-axis heading only when largest/second eigenvalue >1.15; plane normal when middle/smallest >1.2. Kabsch reference-fit angle is bounded 0..180 degrees and flexible deformation can confound it.",
            "qualification": "Finite smoke-validation trajectories of five shapes and three seeds; absence of observed failure is not a proof for all player molecules or converged diffusion."},
        "trials": [], "summary": []}
    path = RACES / "report.json"
    write_json(path, report)
    began = time.perf_counter()
    with ProcessPoolExecutor(max_workers=args.workers) as pool:
        futures = {pool.submit(run_trial, config, hashes, fixture): (fixture, config) for fixture, config in tasks}
        for future in as_completed(futures):
            row = future.result()
            report["trials"].append(row)
            report["summary"] = summary(report["trials"])
            report["elapsed_wall_seconds"] = time.perf_counter() - began
            write_json(path, report)
            print(f"{len(report['trials'])}/{len(tasks)} {row['key']} {row['status']} {row.get('outcome', '')} {row.get('error', '')}", flush=True)
    report["status"] = "passed" if all(r["status"] == "passed" for r in report["trials"]) else "failed"
    if source_hashes(manifest) != hashes:
        report.update(status="failed", error="Source changed during the validation")
    report["trials"].sort(key=lambda r: r["key"])
    write_json(path, report)
    with (RACES / "summary.csv").open("w", encoding="utf-8", newline="") as handle:
        columns = list(report["summary"][0])
        writer = csv.DictWriter(handle, fieldnames=columns)
        writer.writeheader()
        writer.writerows(report["summary"])
    print(f"{report['status']}: {path}", flush=True)
    return 0 if report["status"] == "passed" else 1


if __name__ == "__main__":
    raise SystemExit(main())
