"""Isolated atom-count scaling measurements of the optimized game engine.

Only benchmark fixtures and reports are written. No service, integration step,
force-field parameter or production source is changed. Four equal alkyl arms
around a central carbon give a connected chemically comparable C/H family.
Input conformers are relaxed offline before normal game preparation.
"""
from __future__ import annotations

import argparse
import csv
import hashlib
import itertools
import json
import os
from pathlib import Path
import sys
import time

for name in ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS", "NUMEXPR_NUM_THREADS"):
    os.environ[name] = "1"
ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / "outputs/atom-count-scaling"
sys.path.insert(0, str(ROOT))

import numpy as np
from rdkit import Chem
from rdkit.Chem import AllChem
from molecular_car.engine import Simulation, TIMESTEP_FS

SEEDS = (98245165, 1234567, 2718281)
ARMS = (1, 2, 4, 6, 8, 12, 16, 20)


def write_json(path, value):
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(value, ensure_ascii=False, indent=2) + "\n", encoding="utf-8")


def source_hashes():
    names = ("molecular_car/engine.py", "molecular_car/molecule.py", "molecular_car/fields.py",
             "molecular_car/fields_kernel.py", "molecular_car/external_fields.py", "scripts/benchmark_atom_scaling.py")
    return {name: hashlib.sha256((ROOT / name).read_bytes()).hexdigest() for name in names}


def make_fixture(arm):
    graph = Chem.RWMol()
    root = graph.AddAtom(Chem.Atom("C"))
    for branch in range(4):
        previous = root
        for _ in range(arm):
            index = graph.AddAtom(Chem.Atom("C"))
            graph.AddBond(previous, index, Chem.BondType.SINGLE)
            previous = index
    heavy = graph.GetMol()
    Chem.SanitizeMol(heavy)
    smiles = Chem.MolToSmiles(heavy)
    mol = Chem.AddHs(heavy)
    parameters = AllChem.ETKDGv3()
    parameters.randomSeed = 20261011
    parameters.numThreads = 1
    status = AllChem.EmbedMolecule(mol, parameters)
    if status:
        parameters.useRandomCoords = True
        status = AllChem.EmbedMolecule(mol, parameters)
    if status:
        raise ValueError("ETKDG input embedding failed")
    ff = AllChem.UFFGetMoleculeForceField(mol, vdwThresh=100., ignoreInterfragInteractions=False)
    ff.Initialize()
    iterations = 0
    status = 1
    while status and iterations < 9000:
        status = ff.Minimize(maxIts=300, forceTol=.0005, energyTol=1e-9)
        iterations += 300
    if status:
        raise ValueError("Offline input UFF relaxation did not converge within 9000 iterations")
    name = f"four_arm_alkane_L{arm}"
    mol.SetProp("_Name", name)
    path = OUT / "molecules" / f"{name}.mol"
    path.parent.mkdir(parents=True, exist_ok=True)
    text = Chem.MolToMolBlock(mol)
    path.write_text(text, encoding="utf-8")
    positions = mol.GetConformer().GetPositions()
    return {"name": name, "family": "four-arm saturated alkane", "arm_carbons": arm,
            "carbons": 1 + 4*arm, "atoms": mol.GetNumAtoms(), "smiles": smiles,
            "formula": f"C{1+4*arm}H{2*(1+4*arm)+2}",
            "input_file": str(path.relative_to(ROOT)), "input_sha256": hashlib.sha256(text.encode()).hexdigest(),
            "offline_relaxation_iteration_budget_used": iterations,
            "input_cartesian_span_A": float(np.ptp(positions, axis=0).max()), "text": text}


def topology_counts(mol):
    distances = Chem.GetDistanceMatrix(mol)
    candidates = int(np.count_nonzero(np.triu(distances >= 3., k=1)))
    degrees = [atom.GetDegree() for atom in mol.GetAtoms()]
    return {"nonbonded_candidate_pairs_graph_distance_ge_3": candidates,
            "all_unordered_pairs": mol.GetNumAtoms()*(mol.GetNumAtoms()-1)//2,
            "bonds": mol.GetNumBonds(), "angle_candidates": sum(d*(d-1)//2 for d in degrees),
            "torsion_candidates": sum((degrees[b.GetBeginAtomIdx()]-1)*(degrees[b.GetEndAtomIdx()]-1)
                                       for b in mol.GetBonds())}


def median_loop(method, iterations, batches=3):
    for _ in range(100):
        method()
    times = []
    for _ in range(batches):
        started = time.perf_counter()
        for i in range(iterations):
            method()
        times.append((time.perf_counter()-started)/iterations*1e6)
    return {"median_us": float(np.median(times)), "batches_us": times}


def measure(fixture, seed, steps, warmup):
    sim = Simulation(fixture["text"], Path(fixture["input_file"]).name,
        settings={"surface": "continuum", "temperature": 600., "corrugation": .08, "damping_fs": 350.},
        drive={"mode": "gradient", "strength": .75}, random_seed=seed)
    try:
        if sim.field_evaluator.accelerator != "numba":
            raise AssertionError("Expected compiled external field")
        positions = sim.coordinates()
        if not np.isfinite(positions).all():
            raise AssertionError("Nonfinite initial state")
        sim.lmp.command(f"run {warmup} pre no post no")
        sim._check_callback()
        start_step = int(sim.lmp.get_thermo("step"))
        times = []
        for _ in range(3):
            started = time.perf_counter()
            sim.lmp.command(f"run {steps} pre no post no")
            sim._check_callback()
            times.append((time.perf_counter()-started)/steps*1e6)
        end_step = int(sim.lmp.get_thermo("step"))
        end_positions = sim.coordinates()
        energy = float(sim.lmp.get_thermo("etotal"))
        if end_step-start_step != 3*steps or not np.isfinite(end_positions).all() or not np.isfinite(energy):
            raise AssertionError("Actual propagation/state validation failed")
        flat = end_positions.ravel().tolist()
        ff = sim.prepared.forcefield
        def uff():
            ff.CalcEnergy(flat)
            ff.CalcGrad(flat)
        def field():
            sim.field_evaluator.energy_and_forces(end_positions, 1, .75, 0., 0., 0., .08)
        static_iterations = max(300, min(2000, 60000//sim.n))
        timings = {"UFF_energy_plus_gradient_static": median_loop(uff, static_iterations),
                   "compiled_external_static": median_loop(field, static_iterations),
                   "coordinate_list_conversion_static": median_loop(lambda: end_positions.ravel().tolist(), static_iterations)}
        gradient = np.asarray(ff.CalcGrad(flat))
        external_energy, external_force = sim.field_evaluator.energy_and_forces(end_positions, 1, .75, 0., 0., 0., .08)
        if not np.isfinite(gradient).all() or not np.isfinite(external_force).all() or not np.isfinite(external_energy):
            raise AssertionError("Nonfinite fixed-coordinate force")
        return {"fixture": fixture["name"], "family": fixture["family"], "atoms": sim.n,
                "seed": seed, "batches_us_per_step": times, "median_us_per_step": float(np.median(times)),
                "start_step": start_step, "end_step": end_step, "initial_span_A": float(np.ptp(positions, axis=0).max()),
                "final_span_A": float(np.ptp(end_positions, axis=0).max()), "final_energy_kcal_mol": energy,
                "static_components": timings, "topology": topology_counts(sim.prepared._rdkit_molecule),
                "runtime": sim.runtime, "passed": True}
    finally:
        sim.close()


def nonnegative_fit(columns, target, labels, scales):
    x = np.column_stack(columns)/np.asarray(scales)
    y = np.asarray(target)
    best = None
    for count in range(1, x.shape[1]+1):
        for indices in itertools.combinations(range(x.shape[1]), count):
            coefficients = np.linalg.lstsq(x[:, indices], y, rcond=None)[0]
            if (coefficients < 0).any():
                continue
            all_coefficients = np.zeros(x.shape[1])
            all_coefficients[list(indices)] = coefficients
            predicted = x@all_coefficients
            sse = float(np.sum((predicted-y)**2))
            if best is None or sse < best[0]:
                best = (sse, all_coefficients, predicted)
    sse, coefficients, predicted = best
    return {"coefficients": dict(zip(labels, (coefficients/np.asarray(scales)).tolist())),
            "r_squared": float(1-sse/np.sum((y-y.mean())**2)),
            "rmse_us": float(np.sqrt(sse/len(y))), "predicted_us": predicted.tolist(),
            "maximum_relative_error": float(np.max(np.abs(predicted-y)/y))}


def analysis(summary):
    family = [row for row in summary if row["family"] == "four-arm saturated alkane"]
    n = np.asarray([row["atoms"] for row in family], dtype=float)
    full = np.asarray([row["full_MD_us_per_step"] for row in family])
    uff = np.asarray([row["UFF_static_us"] for row in family])
    pairs = np.asarray([row["topology"]["nonbonded_candidate_pairs_graph_distance_ge_3"] for row in family])
    full_fit = nonnegative_fit([np.ones(len(n)), n, n*n], full,
        ["constant_us", "linear_us_per_atom", "quadratic_us_per_atom_squared"], [1., 100., 10000.])
    uff_fit = nonnegative_fit([np.ones(len(n)), n, pairs], uff,
        ["constant_us", "linear_us_per_atom", "nonbonded_us_per_candidate_pair"], [1., 100., 10000.])
    fitted = np.polyfit(np.log(n), np.log(full), 1)
    local = [{"lower_atoms": int(a), "upper_atoms": int(b), "effective_power": float(np.log(y/x)/np.log(b/a))}
             for a, b, x, y in zip(n[:-1], n[1:], full[:-1], full[1:])]
    return {"family_atoms": n.astype(int).tolist(),
            "descriptive_full_cost_model": "T_us = a + b*N + c*N^2; nonnegative constrained least squares",
            "full_cost_fit": full_fit,
            "descriptive_UFF_cost_model": "T_UFF_us = a + b*N + c*P; P is graph-distance >=3 candidate nonbonded pairs",
            "UFF_cost_fit": uff_fit,
            "descriptive_single_power_fit": {"effective_exponent": float(fitted[0]), "prefactor_us": float(np.exp(fitted[1]))},
            "adjacent_effective_powers": local,
            "limitations": ["This is a finite-size descriptive fit, not a universal scaling exponent.",
                "Explicit hydrogen counts, molecular topology and conformation affect cost.",
                "Static component timings include Python binding cost and may not sum to trajectory MD costs.",
                "The current UFF vdwThresh=100 builds effectively all eligible compact intramolecular pairs; no distance neighbor-list benchmark is claimed.",
                "LAMMPS setup/minimization, JIT initialization, frames, rendering and fixture preparation are excluded."]}


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--steps", type=int, default=1500)
    parser.add_argument("--warmup", type=int, default=500)
    parser.add_argument("--seeds", type=int, default=3)
    args = parser.parse_args()
    if not 1 <= args.seeds <= len(SEEDS) or args.steps <= 0 or args.warmup <= 0:
        raise ValueError("Invalid positive timing protocol")
    OUT.mkdir(parents=True, exist_ok=True)
    started = time.perf_counter()
    hashes = source_hashes()
    fixtures, omitted = [], []
    for arm in ARMS:
        try:
            fixture = make_fixture(arm)
            if fixture["input_cartesian_span_A"] > 45.:
                raise ValueError("Offline input span exceeds 45 A")
            fixtures.append(fixture)
            print(f"INPUT {fixture['name']}: {fixture['atoms']} atoms, {fixture['input_cartesian_span_A']:.2f} A span", flush=True)
        except Exception as exc:
            omitted.append({"arm_carbons": arm, "reason": str(exc)})
            print(f"OMIT arm {arm}: {exc}", flush=True)
    for filename in ("benzene.mol", "tetraphenylbenzene.mol"):
        path = ROOT / "examples" / filename
        fixtures.append({"name": path.stem, "family": "existing game fixture", "input_file": str(path.relative_to(ROOT)),
                         "input_sha256": hashlib.sha256(path.read_bytes()).hexdigest(), "text": path.read_text(encoding="utf-8")})
    write_json(OUT / "molecule-provenance.json", {"fixtures": [{k:v for k,v in f.items() if k != "text"} for f in fixtures], "omitted": omitted})
    trials = []
    for index, seed in enumerate(SEEDS[:args.seeds]):
        order = list(fixtures)
        if index % 2:
            order.reverse()
        for fixture in order:
            try:
                row = measure(fixture, seed, args.steps, args.warmup)
                trials.append(row)
                print(f"TIMING N={row['atoms']} {row['fixture']} seed={seed}: MD={row['median_us_per_step']:.2f} us, UFF={row['static_components']['UFF_energy_plus_gradient_static']['median_us']:.2f} us", flush=True)
            except Exception as exc:
                omitted.append({"name": fixture["name"], "seed": seed, "reason": str(exc)})
                print(f"OMIT {fixture['name']} seed {seed}: {exc}", flush=True)
            write_json(OUT / "trials-progress.json", {"trials": trials, "omitted": omitted})
    summary = []
    for fixture in fixtures:
        rows = [row for row in trials if row["fixture"] == fixture["name"]]
        if len(rows) != args.seeds:
            continue
        costs = [row["median_us_per_step"] for row in rows]
        static = {key: float(np.median([row["static_components"][key]["median_us"] for row in rows]))
                  for key in rows[0]["static_components"]}
        summary.append({"fixture": fixture["name"], "family": fixture["family"], "atoms": rows[0]["atoms"],
            "full_MD_us_per_step": float(np.median(costs)), "trial_min_MD_us": min(costs), "trial_max_MD_us": max(costs),
            "UFF_static_us": static["UFF_energy_plus_gradient_static"], "external_static_us": static["compiled_external_static"],
            "coordinate_conversion_static_us": static["coordinate_list_conversion_static"], "topology": rows[0]["topology"]})
    summary.sort(key=lambda row: (row["family"] != "four-arm saturated alkane", row["atoms"]))
    if len([r for r in summary if r["family"] == "four-arm saturated alkane"]) < 5:
        raise AssertionError("Too few complete chemically comparable family sizes")
    if source_hashes() != hashes:
        raise AssertionError("Production or benchmark source changed during measurement")
    report = {"status": "passed", "protocol": {"surface": "continuum", "external_backend": "numba",
        "temperature_K": 600., "drive_mode": "gradient", "strength": .75, "corrugation": .08, "damping_fs": 350.,
        "time_step_fs": TIMESTEP_FS, "threads": 1, "warmup_steps": args.warmup, "timed_steps_per_batch": args.steps,
        "batches_per_seed": 3, "seeds": list(SEEDS[:args.seeds]), "size_order": "ascending, descending, ascending per seed; existing fixtures at opposite ends",
        "aggregation": "median of 3 timed batches per seed, then median across seeds", "no_product_model_change": True},
        "fixtures": [{k:v for k,v in f.items() if k != "text"} for f in fixtures], "omitted": omitted,
        "trials": trials, "summary": summary, "analysis": analysis(summary),
        "source_sha256": hashes, "wall_seconds": time.perf_counter()-started}
    write_json(OUT / "report.json", report)
    keys = ("fixture", "family", "atoms", "full_MD_us_per_step", "trial_min_MD_us", "trial_max_MD_us", "UFF_static_us", "external_static_us", "coordinate_conversion_static_us")
    with (OUT / "timing-summary.csv").open("w", newline="", encoding="utf-8") as stream:
        writer = csv.DictWriter(stream, fieldnames=list(keys)+["candidate_nonbonded_pairs"])
        writer.writeheader()
        for row in summary:
            writer.writerow({**{key:row[key] for key in keys}, "candidate_nonbonded_pairs": row["topology"]["nonbonded_candidate_pairs_graph_distance_ge_3"]})
    print(f"PASS {len(summary)} sizes/fixtures, {len(trials)} valid propagated trials in {report['wall_seconds']:.1f}s", flush=True)
    print(json.dumps(report["analysis"], ensure_ascii=False), flush=True)


if __name__ == "__main__":
    main()
