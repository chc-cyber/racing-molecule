"""Measure a no-force integration floor, never an alternative racing model.

This diagnostic deliberately omits molecular, surface, field and wall forces.
Its trajectory has no chemical meaning and is neither saved nor compared with
the game.  It estimates how much work remains in LAMMPS and its Python callback
boundary, so that moving the *same* forces to native code can be assessed.

Free H atoms would eventually leave the original fixed game box.  Restore
initial positions and velocities outside timing every 250 steps, keeping the
world, bath, pair-zero neighbor setup and timestep identical to the game.
All requested warmup/timing steps are still executed.  Segmented command
overhead is included; rendering, JSON, initialization and restoration are not.
"""

from __future__ import annotations

import argparse
import hashlib
import json
import os
from pathlib import Path
import statistics
import sys
import time

for _name in ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS", "NUMBA_NUM_THREADS"):
    os.environ[_name] = "1"

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))

import numpy as np

from molecular_car.engine import TIMESTEP_FS, WORLD_BOUNDS
from molecular_car.molecule import prepare_molecule
from scripts.lammps_runtime import open_lammps


def sha256(relative):
    return hashlib.sha256((ROOT / relative).read_bytes()).hexdigest()


def uff_api_probe(content, filename, iterations=2000):
    """Same static UFF, compare supplied-list and conformer-coordinate API.

    Both routes still execute RDKit's C++ UFF and create a Python gradient
    tuple.  Their difference only suggests the supplied-coordinate conversion
    cost; the no-argument path scatters conformer positions internally and
    cannot directly be used for the evolving game trajectory.
    """
    prepared = prepare_molecule(content, filename)
    flat = prepared.positions.ravel().tolist()
    forcefield = prepared.forcefield
    list_energy = forcefield.CalcEnergy(flat)
    list_gradient = np.asarray(forcefield.CalcGrad(flat))
    conformer_energy = forcefield.CalcEnergy()
    conformer_gradient = np.asarray(forcefield.CalcGrad())
    de = abs(list_energy - conformer_energy)
    dg = float(np.max(np.abs(list_gradient - conformer_gradient)))
    if de > 1e-9 or dg > 1e-9:
        raise RuntimeError("Static UFF API probe routes did not evaluate identical coordinates.")

    def calls(supplied):
        if supplied:
            forcefield.CalcEnergy(flat)
            forcefield.CalcGrad(flat)
        else:
            forcefield.CalcEnergy()
            forcefield.CalcGrad()

    for _ in range(250):
        calls(True)
        calls(False)
    timings = {True: [], False: []}
    # Alternate order rather than comparing two long runs at different loads.
    for round_index in range(3):
        order = (True, False) if round_index % 2 == 0 else (False, True)
        for supplied in order:
            began = time.perf_counter_ns()
            for _ in range(iterations):
                calls(supplied)
            timings[supplied].append((time.perf_counter_ns() - began) / iterations / 1000)
    listed = statistics.median(timings[True])
    conformer = statistics.median(timings[False])
    return {
        "file": filename, "atoms": len(prepared.atoms),
        "static_coordinates_only": True,
        "iterations_per_batch": iterations,
        "supplied_list_uff_energy_gradient_batches_us": timings[True],
        "conformer_noargs_uff_energy_gradient_batches_us": timings[False],
        "supplied_list_uff_energy_gradient_median_us": listed,
        "conformer_noargs_uff_energy_gradient_median_us": conformer,
        "supplied_coordinate_api_extra_us": listed - conformer,
        "energy_difference_kcal_mol": de,
        "max_gradient_difference_kcal_mol_angstrom": dg,
        "interpretation": "C++ UFF numerical work remains in both. No-argument API avoids extracting supplied Python coordinates but scatters conformer positions inside C++; it still returns a Python gradient tuple. This static proxy is not a replacement for evolving MD coordinates or a native force kernel timing.",
    }


def measure(content, filename, mode, warmup=1000, steps=2000, segment=250):
    prepared = prepare_molecule(content, filename)
    n = len(prepared.atoms)
    types = list(dict.fromkeys(atom["element"] for atom in prepared.atoms))
    callbacks = 0
    with open_lammps() as (lmp, runtime):
        lmp.commands_list([
            "units real", "atom_style atomic", "boundary f f f",
            "atom_modify map array sort 0 0.0",
            "region world block " + " ".join(str(value) for bounds in WORLD_BOUNDS for value in bounds) + " units box",
            f"create_box {len(types)} world",
        ])
        for index, element in enumerate(types, 1):
            mass = next(atom["mass"] for atom in prepared.atoms if atom["element"] == element)
            lmp.command(f"mass {index} {mass:.12g}")
        count = lmp.create_atoms(n, list(range(1, n + 1)),
                                 [types.index(atom["element"]) + 1 for atom in prepared.atoms],
                                 prepared.positions.ravel().tolist())
        if count != n:
            raise RuntimeError("Native floor could not create the full molecule.")
        lmp.commands_list([
            "pair_style zero 10.0", "pair_coeff * *",
            "compute internal_temperature all temp/com",
            "thermo_style custom step temp pe ke etotal", "thermo 1000000",
            f"timestep {TIMESTEP_FS}",
        ])
        if mode == "empty_python_callback":
            lmp.command("fix molecular all external pf/callback 1 1")

            def callback(caller, step, nlocal, tags, positions, output):
                nonlocal callbacks
                callbacks += 1
                output[:] = 0.0
                lmp.fix_external_set_energy_global("molecular", 0.0)

            lmp.set_fix_external_callback("molecular", callback)
        elif mode != "no_python_callback":
            raise ValueError(mode)
        lmp.commands_list([
            "velocity all create 600.0 17329 mom yes rot yes dist gaussian",
            "fix integrate all nve",
            "fix bath all langevin 600.0 600.0 350.0 98245165 tally yes",
            "run 0",
        ])
        initial_x = lmp.numpy.extract_atom("x", nelem=n).copy()
        initial_v = lmp.numpy.extract_atom("v", nelem=n).copy()

        def execute(total_steps):
            measured = 0.0
            remaining = total_steps
            while remaining:
                count = min(segment, remaining)
                lmp.numpy.extract_atom("x", nelem=n)[:] = initial_x
                lmp.numpy.extract_atom("v", nelem=n)[:] = initial_v
                lmp.command("run 0")
                began = time.perf_counter_ns()
                lmp.command(f"run {count} pre no post no")
                measured += time.perf_counter_ns() - began
                if lmp.get_natoms() != n:
                    raise RuntimeError("Native floor lost an atom.")
                remaining -= count
            return measured / total_steps / 1000.0

        execute(warmup)
        calls_before = callbacks
        batches = [execute(steps) for _ in range(3)]
        return {
            "file": filename, "atoms": n, "mode": mode,
            "batches_us_per_step": batches,
            "median_us_per_step": statistics.median(batches),
            "callback_calls_including_untimed_run0": callbacks - calls_before,
            "runtime": runtime,
        }


def add_model(report, optimized_path, benchmark_path=None):
    """Optionally attach a model using fresh optimized cost decomposition.

    The decomposition report must contain results with file names and
    components_us_per_step; the keys are those of the historical diagnostic.
    All measured RDKit UFF and compiled-field costs are retained.  In reality
    a native path could speed their conversions but still needs their physics.
    This is a cost model, not a native implementation benchmark or prediction.
    """
    optimized = json.loads(optimized_path.read_text(encoding="utf-8"))
    benchmark = None if benchmark_path is None else json.loads(benchmark_path.read_text(encoding="utf-8"))
    profile_sources = optimized.get("source_sha256", {})
    benchmark_sources = {} if benchmark is None else benchmark.get("source_sha256", {})
    if benchmark is not None and profile_sources != benchmark_sources:
        raise ValueError("The final optimized profile and benchmark must share identical recorded source hashes.")
    report["provenance"] = {
        "historical_floor_and_static_uff_measurement": {
            "source_sha256": dict(report.get("source_sha256", {})),
            "interpretation": "The top-level source_sha256 belongs to the earlier no-force floor and static UFF API probe, not the final optimized MD cost model. The floor did not instantiate or execute engine.Simulation; it imported only TIMESTEP_FS and WORLD_BOUNDS from the historical engine. This historical engine hash is intentionally retained rather than relabeled as the final validated engine.",
            "imported_engine_values": {
                "timestep_fs": 0.5,
                "world_bounds_angstrom": [[-110., 1110.], [-350., 350.], [-40., 120.]],
            },
        },
        "final_optimized_native_cost_model_inputs": {
            "profile_source_sha256": profile_sources,
            "benchmark_source_sha256": benchmark_sources,
            "validated_engine_sha256": profile_sources.get("molecular_car/engine.py"),
            "profile_and_benchmark_sources_identical": benchmark is not None and profile_sources == benchmark_sources,
            "interpretation": "Current full-MD cost, retained UFF cost and retained compiled external-field cost come from the final validated input reports, whose individual file hashes and full source maps are recorded. The historical no-force floor is a separately measured diagnostic used only as a cost-model input; it is not a final native racing benchmark.",
        },
    }
    modeled = []
    for row in optimized["results"]:
        if row.get("implementation", "optimized") != "optimized":
            continue
        native = next(result for result in report["results"]
                      if result["file"] == row["file"] and result["mode"] == "no_python_callback")
        empty = next(result for result in report["results"]
                     if result["file"] == row["file"] and result["mode"] == "empty_python_callback")
        components = row["components_us_per_step"]
        uff = components["UFF_energy"] + components["UFF_gradient"]
        field = components["surface_drive_walls"]
        current = row.get("baseline_median_us_per_step", row["instrumented_us_per_step"])
        current_kind = "profile_instrumented"
        if benchmark is not None:
            current = statistics.median(trial["median_us_per_step"]["optimized"]
                                        for trial in benchmark["trials"] if trial["file"] == row["file"])
            current_kind = "median_of_three_seed_uninstrumented_batch_medians"
        cost = native["median_us_per_step"] + uff + field
        boundary = empty["median_us_per_step"] - native["median_us_per_step"]
        static_probe = next(probe for probe in report["static_uff_api_probe"] if probe["file"] == row["file"])
        conversion_proxy = max(0., static_probe["supplied_coordinate_api_extra_us"])
        modeled.append({
            "file": row["file"], "atoms": row["atoms"],
            "current_optimized_full_md_us_per_step": current,
            "current_full_md_cost_kind": current_kind,
            "profile_instrumented_full_md_us_per_step": row["instrumented_us_per_step"],
            "retained_rdkit_uff_energy_gradient_us_per_step": uff,
            "retained_compiled_external_fields_us_per_step": field,
            "native_no_force_floor_us_per_step": native["median_us_per_step"],
            "measured_empty_python_boundary_extra_us_per_step": boundary,
            "modeled_native_with_retained_force_cost_us_per_step": cost,
            "modeled_speedup_if_retaining_measured_force_cost": current / cost,
            "static_uff_coordinate_conversion_proxy_us_per_step": conversion_proxy,
            "modeled_native_if_static_conversion_proxy_also_removed_us_per_step": cost - conversion_proxy,
            "modeled_speedup_if_static_conversion_proxy_also_removed": current / (cost - conversion_proxy),
            "all_force_cost_zero_theoretical_speedup_ceiling": current / native["median_us_per_step"],
            "interpretation": "Two bookkeeping scenarios, not an expected performance interval: retain measured force costs, or also subtract a static supplied-coordinate API conversion proxy. The latter compares different internal coordinate routes and is particularly uncertain. Actual native topology work and force implementation may cost more or less. Neither speedup is measured for a working native UFF racer.",
        })
    report["native_path_cost_model"] = modeled
    report["optimized_cost_report"] = str(optimized_path.resolve())
    report["optimized_cost_report_sha256"] = hashlib.sha256(optimized_path.read_bytes()).hexdigest()
    if benchmark_path is not None:
        report["optimized_md_benchmark"] = str(benchmark_path.resolve())
        report["optimized_md_benchmark_sha256"] = hashlib.sha256(benchmark_path.read_bytes()).hexdigest()


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--output", type=Path, default=ROOT / "outputs" / "native-path-estimate" / "report.json")
    parser.add_argument("--optimized-cost-report", type=Path)
    parser.add_argument("--optimized-md-benchmark", type=Path)
    parser.add_argument("--segment", type=int, default=250)
    parser.add_argument("--update-model-only", action="store_true",
                        help="Attach a fresh optimized profile to the saved floor; do not repeat timings.")
    args = parser.parse_args()
    if args.update_model_only:
        if args.optimized_cost_report is None or not args.output.is_file():
            parser.error("--update-model-only needs an existing --output and --optimized-cost-report")
        report = json.loads(args.output.read_text(encoding="utf-8"))
        add_model(report, args.optimized_cost_report, args.optimized_md_benchmark)
        report["native_cost_model_source_sha256"] = sha256("scripts/estimate_native_path.py")
        args.output.write_text(json.dumps(report, ensure_ascii=False, indent=2), encoding="utf-8")
        print(f"Updated modeled estimate without repeating timing: {args.output}")
        return
    if not 1 <= args.segment <= 250:
        parser.error("segment must lie within 1..250 to avoid free-atom escape")
    report = {
        "status": "diagnostic_only",
        "protocol": "Single-thread LAMMPS, original game world and pair_style zero 10; 0.5 fs; 600 K Langevin; damping 350 fs; velocity seed 17329; Langevin seed 98245165; 1000 warmup steps; three 2000-step batches; initial x/v restored outside timing every 250 steps; no frames, JSON or network in timing.",
        "interpretation": "No-force floor deliberately omits all UFF and external forces. Its trajectory is nonphysical and is NOT a native racing benchmark. Removing the Python callback does not remove the cost of native energy/force evaluation. The optional native-path model retains all measured optimized UFF energy+gradient and compiled field cost and adds the native floor; the ratio is an estimate, not measured speedup of a migrated force field.",
        "limitations": [
            "No native UFF migration was implemented or validated in this diagnostic.",
            "The floor does not include future native molecular topology or neighbor-list work.",
            "Segmented run commands add small driver overhead; repeated restoration excludes free-atom changes in neighbor load and lost atoms.",
            "One seed and two molecules do not establish general scaling or a guaranteed native-path gain.",
        ],
        "official_sources": {
            "RDKit_python_forcefield": "https://www.rdkit.org/docs/source/rdkit.ForceField.rdForceField.html",
            "RDKit_cpp_forcefield": "https://www.rdkit.org/docs/cppapi/classForceFields_1_1ForceField.html",
            "LAMMPS_external_callback": "https://docs.lammps.org/fix_external.html",
        },
        "results": [],
        "static_uff_api_probe": [],
        "source_sha256": {relative: sha256(relative) for relative in (
            "scripts/estimate_native_path.py", "scripts/lammps_runtime.py", "molecular_car/engine.py", "molecular_car/molecule.py")},
    }
    if args.segment != 250:
        report["protocol"] = report["protocol"].replace("every 250 steps", f"every {args.segment} steps")
    for filename in ("benzene.mol", "tetraphenylbenzene.mol"):
        content = (ROOT / "examples" / filename).read_text(encoding="utf-8")
        for mode in ("no_python_callback", "empty_python_callback"):
            row = measure(content, filename, mode, segment=args.segment)
            report["results"].append(row)
            print(f'{filename} {mode}: {row["median_us_per_step"]:.3f} us/step', flush=True)
        report["static_uff_api_probe"].append(uff_api_probe(content, filename))
    if args.optimized_cost_report:
        add_model(report, args.optimized_cost_report, args.optimized_md_benchmark)
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(report, ensure_ascii=False, indent=2), encoding="utf-8")
    print(f"Saved diagnostic: {args.output}")


if __name__ == "__main__":
    main()
