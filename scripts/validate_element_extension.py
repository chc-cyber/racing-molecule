"""Probe the actual game preparation path with locally generated element fixtures.

This first phase never imports or starts Simulation/LAMMPS. Passing UFF and
Gasteiger preparation is a software applicability check, not evidence that UFF,
effective polarizabilities or heteroatom/Au interactions are experimentally
accurate for these molecular environments.
"""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
from pathlib import Path
import platform
import sys
import time

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))

import numpy as np
from rdkit import Chem, rdBase
from rdkit.Chem import AllChem, rdMolDescriptors

from molecular_car.molecule import EFFECTIVE_POLARIZABILITY, EMBED_SEED, STARTING_MIN_Z, prepare_molecule


CASES = [
    ("H", "methane", "C", "representative", "passed"),
    ("C", "benzene", "c1ccccc1", "representative", "passed"),
    ("B", "phenylboronic-acid", "OB(O)c1ccccc1", "representative", "passed"),
    ("N", "pyridine", "c1ccncc1", "representative", "passed"),
    ("O", "ethanol", "CCO", "representative", "passed"),
    ("F", "fluorobenzene", "Fc1ccccc1", "representative", "passed"),
    ("Si", "tetramethylsilane", "C[Si](C)(C)C", "representative", "passed"),
    ("P", "trimethylphosphine-oxide", "CP(C)(C)=O", "representative", "passed"),
    ("S", "dimethyl-sulfide", "CSC", "representative", "passed"),
    ("Cl", "chlorobenzene", "Clc1ccccc1", "representative", "passed"),
    ("Br", "bromobenzene", "Brc1ccccc1", "representative", "passed"),
    ("I", "iodobenzene", "Ic1ccccc1", "representative", "passed"),
    ("B", "trimethylborane", "B(C)(C)C", "environment", "passed"),
    ("B", "phenylboronic-dimethyl-ester", "COB(OC)c1ccccc1", "environment", "passed"),
    ("N", "trimethylamine", "CN(C)C", "environment", "passed"),
    ("N", "acetonitrile", "CC#N", "environment", "passed"),
    ("N", "acetamide", "CC(=O)N", "environment", "passed"),
    ("N", "nitromethane", "C[N+](=O)[O-]", "environment", "passed"),
    ("O", "dimethyl-ether", "COC", "environment", "passed"),
    ("O", "acetone", "CC(=O)C", "environment", "passed"),
    ("O", "dimethyl-peroxide", "COOC", "environment", "passed"),
    ("Si", "hexamethyldisiloxane", "C[Si](C)(C)O[Si](C)(C)C", "environment", "passed"),
    ("P", "trimethylphosphine", "CP(C)C", "parameter-boundary", "rejected-other"),
    ("P", "trimethyl-phosphate", "COP(=O)(OC)OC", "environment", "passed"),
    ("S", "ethanethiol", "CCS", "environment", "passed"),
    ("S", "thiophene", "c1ccsc1", "environment", "passed"),
    ("S", "dimethyl-sulfoxide", "CS(C)=O", "environment", "passed"),
    ("S", "dimethyl-sulfone", "CS(C)(=O)=O", "environment", "passed"),
    # A valid neutral, closed-shell graph that the current UFF type checker
    # rejects. Keep this real boundary case rather than claiming all valences
    # are supported once an element has one passing example.
    ("S", "methyl-sulfur-pentafluoride", "CS(F)(F)(F)(F)F", "parameter-boundary", "rejected-uff"),
    # Newly exposed by the full UFF energy/force gate: the old finite-only
    # path accepted I(III), despite a 4.75746 kcal/mol/A FD mismatch.
    ("I", "iodosylbenzene", "O=Ic1ccccc1", "parameter-boundary", "rejected-other"),
    ("I", "iodoxybenzene", "O=I(=O)c1ccccc1", "parameter-boundary", "rejected-other"),
]


def sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def charge_probe(molecule: Chem.Mol) -> dict:
    """Check charges independently, even if UFF rejects this environment."""
    copied = Chem.Mol(molecule)
    try:
        with rdBase.BlockLogs():
            AllChem.ComputeGasteigerCharges(copied, nIter=12, throwOnParamFailure=True)
        charges = np.array([float(a.GetProp("_GasteigerCharge")) for a in copied.GetAtoms()])
        total_formal_charge = sum(a.GetFormalCharge() for a in copied.GetAtoms())
        finite = bool(np.isfinite(charges).all())
        balanced = bool(finite and np.isclose(charges.sum(), total_formal_charge, atol=1e-5))
        return {
            "status": "passed" if finite and balanced else "failed",
            "finite": finite,
            "formal_charge_e": total_formal_charge,
            "partial_charge_sum_e": float(charges.sum()) if finite else None,
            "charge_range_e": [float(charges.min()), float(charges.max())] if finite else None,
            "n_iter": 12,
            "throw_on_parameter_failure": True,
        }
    except Exception as exc:
        return {"status": "failed", "exception_type": type(exc).__name__, "error": str(exc)}


def atom_environments(molecule: Chem.Mol, element: str) -> list[dict]:
    return [
        {
            "atom_index_zero_based": atom.GetIdx(),
            "element": atom.GetSymbol(),
            "degree_with_explicit_H": atom.GetDegree(),
            "hybridization": str(atom.GetHybridization()),
            "aromatic": atom.GetIsAromatic(),
            "formal_charge_e": atom.GetFormalCharge(),
            "neighbours": [a.GetSymbol() for a in atom.GetNeighbors()],
        }
        for atom in molecule.GetAtoms() if atom.GetSymbol() == element
    ]


def standalone_uff_fd(explicit: Chem.Mol) -> dict:
    """Diagnose UFF itself even if the production gate rejects this molecule.

    Repeat the documented fixed-seed geometry procedure on a separate RDKit
    molecule, then differentiate all Cartesian components of its scalar energy.
    This never substitutes a different force field or modifies game admission.
    """
    molecule = Chem.Mol(explicit)
    try:
        parameters = AllChem.ETKDGv3()
        parameters.randomSeed = EMBED_SEED
        parameters.numThreads = 1
        with rdBase.BlockLogs():
            embedded = AllChem.EmbedMolecule(molecule, parameters)
            if embedded != 0:
                parameters.useRandomCoords = True
                embedded = AllChem.EmbedMolecule(molecule, parameters)
            if embedded != 0:
                raise ValueError("Independent ETKDG embedding failed")
            field = AllChem.UFFGetMoleculeForceField(molecule, vdwThresh=100., ignoreInterfragInteractions=False)
            field.Initialize()
            minimizer_status = field.Minimize(maxIts=300, forceTol=.001, energyTol=1e-7)
        x = np.array(molecule.GetConformer().GetPositions(), copy=True)
        masses = np.array([atom.GetMass() for atom in molecule.GetAtoms()])
        center = np.average(x, axis=0, weights=masses)
        x[:, :2] -= center[:2]
        x[:, 2] += STARTING_MIN_Z - x[:, 2].min()
        for index, point in enumerate(x):
            molecule.GetConformer().SetAtomPosition(index, tuple(point))
        with rdBase.BlockLogs():
            field = AllChem.UFFGetMoleculeForceField(molecule, vdwThresh=100., ignoreInterfragInteractions=False)
            field.Initialize()
        flat = x.ravel().copy()
        energy = float(field.CalcEnergy(tuple(flat)))
        force = -np.asarray(field.CalcGrad(tuple(flat)))
        if not np.isfinite(energy) or not np.isfinite(force).all():
            raise ValueError("Independent UFF energy or forces are nonfinite")
        epsilon = 1e-5
        errors, normalized = [], []
        for component in range(flat.size):
            plus, minus = flat.copy(), flat.copy()
            plus[component] += epsilon
            minus[component] -= epsilon
            fd = -(field.CalcEnergy(tuple(plus)) - field.CalcEnergy(tuple(minus))) / (2 * epsilon)
            error = abs(force[component] - fd)
            errors.append(float(error))
            normalized.append(float(error / (.001 + 2e-6 * abs(fd))))
        return {
            "status": "passed" if max(normalized) <= 1. else "force-energy-inconsistent",
            "geometry_protocol": "Independent copy: production ETKDG seed, 300 UFF steps, game COM/min-height placement, rebuild UFF",
            "minimizer_return_code": int(minimizer_status),
            "energy_kcal_mol": energy,
            "maximum_atom_force_kcal_mol_A": float(np.linalg.norm(force.reshape(-1, 3), axis=1).max()),
            "global_force_norm_kcal_mol_A": float(np.linalg.norm(force)),
            "checked_cartesian_components": flat.size,
            "epsilon_A": epsilon,
            "absolute_tolerance_kcal_mol_A": .001,
            "relative_tolerance": 2e-6,
            "maximum_force_fd_error_kcal_mol_A": max(errors),
            "maximum_error_over_tolerance": max(normalized),
            "worst_flat_component_zero_based": int(np.argmax(normalized)),
            "force_energy_consistent": bool(max(normalized) <= 1.),
        }
    except Exception as exc:
        return {"status": "diagnostic-failed", "error": str(exc), "exception_type": type(exc).__name__}


def probe_case(case: tuple[str, str, str, str, str], output_dir: Path) -> dict:
    element, name, smiles, role, expected = case
    item = {
        "target_element": element, "name": name, "smiles": smiles,
        "role": role, "expected_outcome": expected,
    }
    begin = time.perf_counter()
    with rdBase.BlockLogs():
        molecule = Chem.MolFromSmiles(smiles)
    if molecule is None:
        raise ValueError(f"Fixture SMILES is invalid: {name}")
    formal_charge = sum(a.GetFormalCharge() for a in molecule.GetAtoms())
    radicals = sum(a.GetNumRadicalElectrons() for a in molecule.GetAtoms())
    if formal_charge or radicals or len(Chem.GetMolFrags(molecule)) != 1:
        raise ValueError(f"Fixture must be one neutral closed-shell molecule: {name}")
    if not any(a.GetSymbol() == "C" for a in molecule.GetAtoms()):
        raise ValueError(f"Fixture must have a carbon skeleton: {name}")
    molecule.SetProp("_Name", name)
    AllChem.Compute2DCoords(molecule)
    text = Chem.MolToMolBlock(molecule)
    input_file = output_dir / "inputs" / f"{name}.mol"
    input_file.write_text(text, encoding="utf-8", newline="\n")
    explicit = Chem.AddHs(molecule)
    item.update({
        "input_file": str(input_file.relative_to(ROOT)).replace("\\", "/"),
        "input_sha256": sha256(input_file),
        "input_coordinates": "2D local RDKit depiction; actual preparation embeds ETKDGv3",
        "formula": rdMolDescriptors.CalcMolFormula(explicit),
        "atom_count_with_H": explicit.GetNumAtoms(),
        "elements": sorted({a.GetSymbol() for a in explicit.GetAtoms()}),
        "formal_charge_e": formal_charge,
        "radical_electrons": radicals,
        "target_atom_environments": atom_environments(explicit, element),
    })
    with rdBase.BlockLogs():
        uff_available = bool(AllChem.UFFHasAllMoleculeParams(explicit))
    item["uff_all_parameters_available"] = uff_available
    item["gasteiger_independent"] = charge_probe(explicit)
    item["uff_independent_full_cartesian_fd"] = (standalone_uff_fd(explicit) if uff_available
        else {"status": "not-applicable", "reason": "UFF atom typing has missing parameters"})
    try:
        # This is the unmodified production preparation function: parameter
        # checks, ETKDG, UFF minimization and Gasteiger charge assignment.
        prepared = prepare_molecule(text, input_file.name)
        x = tuple(float(v) for v in prepared.positions.ravel())
        energy = float(prepared.forcefield.CalcEnergy(x))
        forces = -np.asarray(prepared.forcefield.CalcGrad(x)).reshape((-1, 3))
        masses = np.array([a["mass"] for a in prepared.atoms])
        charges = np.array([a["charge"] for a in prepared.atoms])
        polarizability = np.array([a["polarizability"] for a in prepared.atoms])
        if not (np.isfinite(energy) and np.isfinite(forces).all() and
                np.isfinite(prepared.positions).all() and np.isfinite(charges).all() and
                np.isfinite(masses).all() and (masses > 0).all() and
                np.isfinite(polarizability).all()):
            raise ValueError("Prepared UFF energy/forces, positions, charges or masses are invalid")
        if not np.isclose(charges.sum(), formal_charge, atol=1e-5):
            raise ValueError("Prepared charges do not sum to the formal charge")
        fixture_file = output_dir / "prepared" / f"{name}.mol"
        fixture_file.write_text(Chem.MolToMolBlock(prepared._rdkit_molecule), encoding="utf-8", newline="\n")
        item.update({
            "preparation_status": "passed", "outcome": "passed",
            "prepared_file": str(fixture_file.relative_to(ROOT)).replace("\\", "/"),
            "prepared_sha256": sha256(fixture_file),
            "prepared_atoms": prepared.atoms,
            "bonds_zero_based": [[b["a"] - 1, b["b"] - 1, b["order"]] for b in prepared.bonds],
            "notes": prepared.notes,
            "uff_energy_kcal_mol": energy,
            "uff_force_max_atom_kcal_mol_A": float(np.linalg.norm(forces, axis=1).max()),
            "uff_force_global_norm_kcal_mol_A": float(np.linalg.norm(forces)),
            "partial_charge_sum_e": float(charges.sum()),
            "partial_charge_range_e": [float(charges.min()), float(charges.max())],
            "mass_range_u": [float(masses.min()), float(masses.max())],
            "prepared_position_span_A": np.ptp(prepared.positions, axis=0).tolist(),
            "prepared_min_atom_height_A": float(prepared.positions[:, 2].min()),
            # These are deliberately explicit software screening thresholds,
            # not claimed physical stability criteria or a change to game
            # acceptance. The extreme-I(V) boundary is finite yet unsuitable
            # as an unchecked starting state at the game's fixed timestep.
            "preparation_force_screen": {
                "residual_force_alert": bool(np.linalg.norm(forces, axis=1).max() > 1.0),
                "residual_force_alert_threshold_kcal_mol_A": 1.0,
                "extreme_energy_or_force": bool(abs(energy) > 1e6 or np.linalg.norm(forces, axis=1).max() > 1e4),
                "extreme_energy_threshold_kcal_mol": 1e6,
                "extreme_atom_force_threshold_kcal_mol_A": 1e4,
                "interpretation": "Diagnostic screening only. Actual preparation acceptance is recorded separately; joint force convergence and MD readiness are not established here.",
            },
        })
    except Exception as exc:
        (output_dir / "prepared" / f"{name}.mol").unlink(missing_ok=True)
        item.update({
            "preparation_status": "rejected",
            "outcome": "rejected-uff" if not uff_available else "rejected-other",
            "exception_type": type(exc).__name__, "error": str(exc),
        })
    item["matches_expected_outcome"] = item["outcome"] == expected
    item["elapsed_seconds"] = time.perf_counter() - begin
    return item


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--probe-only", action="store_true", help="Documented first phase; no MD is ever started by this version")
    parser.add_argument("--output", type=Path, default=ROOT / "outputs/md-preparation-elements/elements")
    args = parser.parse_args()
    output_dir = args.output.resolve()
    if not output_dir.is_relative_to(ROOT / "outputs"):
        raise ValueError("Validation output must stay inside this workspace's outputs directory")
    (output_dir / "inputs").mkdir(parents=True, exist_ok=True)
    (output_dir / "prepared").mkdir(parents=True, exist_ok=True)
    source_files = [ROOT / "molecular_car/molecule.py", ROOT / "molecular_car/element_models.py", Path(__file__).resolve()]
    before = {str(p.relative_to(ROOT)).replace("\\", "/"): sha256(p) for p in source_files}
    start = time.perf_counter()
    cases = [probe_case(case, output_dir) for case in CASES]
    after = {str(p.relative_to(ROOT)).replace("\\", "/"): sha256(p) for p in source_files}
    per_element = []
    for element in EFFECTIVE_POLARIZABILITY:
        selected = [c for c in cases if c["target_element"] == element]
        primary = [c for c in selected if c["role"] == "representative"]
        per_element.append({
            "element": element,
            "effective_game_polarizability": EFFECTIVE_POLARIZABILITY[element],
            "representative_passed": bool(primary and all(c["preparation_status"] == "passed" for c in primary)),
            "tested_cases": len(selected),
            "prepared_cases": [c["name"] for c in selected if c["preparation_status"] == "passed"],
            "rejected_cases": [c["name"] for c in selected if c["preparation_status"] == "rejected"],
        })
    success = before == after and all(c["matches_expected_outcome"] for c in cases)
    report = {
        "status": "passed" if success else "failed",
        "phase": "production preparation probe only; no Simulation or LAMMPS",
        "created_utc": time.strftime("%Y-%m-%dT%H:%M:%SZ", time.gmtime()),
        "python_version": platform.python_version(), "rdkit_version": rdBase.rdkitVersion,
        "source_sha256_before": before, "source_sha256_after": after,
        "source_unchanged": before == after,
        "effective_game_polarizability": dict(EFFECTIVE_POLARIZABILITY),
        "md_readiness": "Not assessed; no joint relaxation, surface forces, or time integration in this phase",
        "protocol": {
            "input_origin": "Locally authored SMILES, sanitized RDKit graph and 2D MOL; no downloaded structures",
            "molecule_scope": "One carbon-containing molecule, net neutral and closed-shell; nitro has cancelling local formal charges",
            "preparation": "Actual prepare_molecule: UFF parameter guard, fixed-seed ETKDGv3, at most 300 UFF iterations, finite UFF energy/gradient, strict 12-iteration Gasteiger",
            "independent_parameter_checks": "UFFHasAllMoleculeParams, copied-graph Gasteiger and all-3N UFF scalar-energy finite differences on independent fixed-seed geometry",
            "expected_boundaries": "Methyl sulfur pentafluoride is rejected for missing UFF parameters. Tricoordinate P and iodoxybenzene are checked for actual strict preparation rejection and their UFF energy/force mismatch is recorded independently.",
            "no_md": True,
        },
        "limitations": [
            "Passing parameter availability and finite preparation does not establish chemical accuracy, stability, force-field transferability or STM behaviour.",
            "One supported environment does not imply all valences, bond orders, charges or radicals of that element are supported.",
            "UFF internal forces contain no Gasteiger pair electrostatics; these fixed charges feed the separate game external field.",
            "Atomwise polarizabilities are fixed effective game coefficients, not calibrated molecular polarizabilities.",
            "Heteroatom/Au forces are not tested here. Future isotropic LJ extensions must be labelled as game approximations, distinct from the CHAu SAIP parameter set.",
            "Hypervalent iodine examples passing preparation only means the implemented atom typing/charge routines accept them; it is no evidence of model accuracy.",
            "The legacy finite-only preparation accepted iodoxybenzene with extreme UFF values and tricoordinate phosphine with inconsistent UFF energy/gradient. These are explicit boundaries rather than supported racing environments.",
            "The imported UFF optimization and later joint UFF/surface relaxation are different stages; this probe only checks the former.",
        ],
        "per_element": per_element, "cases": cases,
        "summary": {
            "representative_elements": len(per_element),
            "representative_elements_passed": sum(e["representative_passed"] for e in per_element),
            "cases": len(cases),
            "prepared": sum(c["preparation_status"] == "passed" for c in cases),
            "rejected": sum(c["preparation_status"] == "rejected" for c in cases),
            "unexpected_outcomes": [c["name"] for c in cases if not c["matches_expected_outcome"]],
            "accepted_with_residual_force_alert": [c["name"] for c in cases if c.get("preparation_force_screen", {}).get("residual_force_alert")],
            "accepted_with_extreme_energy_or_force": [c["name"] for c in cases if c.get("preparation_force_screen", {}).get("extreme_energy_or_force")],
            "independent_uff_force_energy_inconsistent": [c["name"] for c in cases if c["uff_independent_full_cartesian_fd"]["status"] == "force-energy-inconsistent"],
        },
        "elapsed_seconds": time.perf_counter() - start,
    }
    report_file = output_dir / "probe-report.json"
    report_file.write_text(json.dumps(report, ensure_ascii=False, indent=2) + "\n", encoding="utf-8")
    fields = ["target_element", "name", "smiles", "role", "formula", "atom_count_with_H", "uff_all_parameters_available", "preparation_status", "outcome", "uff_energy_kcal_mol", "uff_force_global_norm_kcal_mol_A", "uff_force_max_atom_kcal_mol_A", "error"]
    with (output_dir / "probe-summary.csv").open("w", encoding="utf-8-sig", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, extrasaction="ignore")
        writer.writeheader()
        writer.writerows(cases)
    print(json.dumps({"status": report["status"], **report["summary"], "report": str(report_file)}, indent=2))
    if not success:
        raise SystemExit(1)


if __name__ == "__main__":
    main()
