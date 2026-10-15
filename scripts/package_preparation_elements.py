"""Verify completed acceptance data and save a compact reproducible archive."""
import hashlib
import json
from pathlib import Path
import zipfile

ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT/"outputs/md-preparation-elements"


def read(name):
    return json.loads((OUT/name).read_text(encoding="utf-8"))


def digest(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def main():
    convergence = read("preparation/convergence-report.json")
    hetero = read("preparation/hetero-acceptance-report.json")
    probe = read("elements/probe-report.json")
    md = read("elements/md/report-all.json")
    audit = read("elements/md/coordinate-audit.json")
    numerics = read("numerics/report-all.json")
    review = read("review/hetero-static-review.json")
    api = read("api/port-8782/report.json")
    assert convergence["status"] == md["status"] == probe["status"] == numerics["status"] == audit["status"] == "passed"
    assert hetero["passed"] and review["passed"] and api["passed"]
    assert len(convergence["cases"]) == 15 and all(c["strict"]["checks_passed"] for c in convergence["cases"])
    accepted = [c for c in hetero["cases"] if not c["expected_rejection"]]
    assert len(accepted) == 36 and all(c["passed"] for c in hetero["cases"])
    assert len(md["driven"]) == 15 and len(md["nve"]) == 8 and len(md["rejections"]) == 4
    assert all(c["passed"] for c in md["driven"]+md["nve"]+md["rejections"])
    sources = {}
    for report, key in ((convergence, "source_sha256_before"), (hetero, "source_sha256_before"),
                        (probe, "source_sha256_before"), (md, "source_sha256"),
                        (numerics, "source_sha256"), (review, "fingerprint_start"), (api, "source_sha256")):
        for name, value in report[key].items():
            name = name.replace("\\", "/")
            if name in sources:
                assert sources[name] == value, "Inconsistent source in reports: "+name
            assert digest(ROOT/name) == value, "Source changed after acceptance: "+name
            sources[name] = value
    for case in md["driven"]+md["nve"]:
        assert digest(ROOT/case["trajectory_file"]) == case["trajectory_sha256"]
    for case in read("numerics/nve-report.json")["cases"]:
        assert digest(Path(case["initial_state_and_trajectories_npz"])) == case["trajectory_sha256"]
    strict = [c["strict"]["preparation"] for c in convergence["cases"]]
    summary = {
        "completed_requested_checks": True,
        "settings": {"surface": "periodic-au", "temperature_K": 600., "gradient_strength": .75,
                     "timestep_fs": .5, "langevin_damping_fs": 350.},
        "preparation": {"protocol_version": convergence["protocol_version"], "CH_paired_cases": 15,
                        "element_surface_cases": len(accepted), "accepted_cases_total": 51,
                        "force_norm_range_kcal_mol_A": [min(p["global_force_norm_kcal_mol_A"] for p in strict),
                                                        max(p["global_force_norm_kcal_mol_A"] for p in strict)],
                        "criterion": "Actual model full 3N conservative force two-norm <=1e-4; final-coordinate refresh mandatory; exhausted budgets reject."},
        "elements": {"supported": ["H", "B", "C", "N", "O", "F", "Si", "P", "S", "Cl", "Br", "I"],
                     "environment_probe_cases": len(probe["cases"]), "probe_summary": probe["summary"],
                     "known_rejected_environments": ["six-coordinate S missing UFF parameters", "three-coordinate P UFF energy/gradient inconsistency",
                                                     "tested I(III) UFF gradient inconsistency", "tested I(V) pathological optimized energy/force"],
                     "Au_interface": "CH keeps existing SAIP; X uses versioned provisional nonreactive UFF-vdW lattice LJ."},
        "MD": {"driven_5ps_cases": 15, "mixed_NVE_2ps_pairs": 8, "CH_NVE_5ps_pairs": 9,
               "mixed_full_coordinate_samples_audited": audit["full_coordinate_samples"],
               "fine_to_coarse_mixed_energy_error_ratio": [min(c["fine_to_coarse_energy_error_ratio"] for c in md["nve"]),
                                                         max(c["fine_to_coarse_energy_error_ratio"] for c in md["nve"])]},
        "tests": {"Python_passed": 142, "API_passed": len(api["checks"])},
        "preview": {"url": api["base_url"], "server_pid": api["server_pid"], "state": api["preview"]},
        "remaining_limits": ["Geometric local convergence, not global minimum or 600 K equilibrium sampling.",
                             "Hetero-Au adsorption parameters are not calibrated; no chemistry, charge transfer or image response.",
                             "Known rejected environments require separate UFF repair/validation; arbitrary player states remain outside these fixtures.",
                             "Contact-dependent friction, tensor polarizability, Au motion and defects remain deferred."],
        "source_sha256": sources,
    }
    (OUT/"summary.json").write_text(json.dumps(summary, ensure_ascii=False, allow_nan=False, indent=2)+"\n", encoding="utf-8")
    paths = {OUT/"summary.json", OUT/"python-tests.log", OUT/"elements/probe-report.json", OUT/"elements/probe-summary.csv"}
    for directory in ("elements/md", "elements/inputs", "elements/prepared", "numerics", "review", "api/port-8782"):
        paths.update(p for p in (OUT/directory).rglob("*") if p.is_file())
    for pattern in ("*-initial.npz", "preparation-before-after.*", "convergence-report.json",
                    "hetero-acceptance-report.json", "native-failure-cleanup-report.json", "source-snapshot.zip"):
        paths.update((OUT/"preparation").glob(pattern))
    archive_path = ROOT/"outputs/md-preparation-elements-results.zip"
    with zipfile.ZipFile(archive_path, "w", compression=zipfile.ZIP_DEFLATED) as archive:
        for path in sorted(paths):
            archive.write(path, "results/"+path.relative_to(OUT).as_posix())
        for name in sorted(sources):
            archive.write(ROOT/name, "source/"+name)
        for name in ("README.md", "scripts/validate_preparation_convergence.py", "scripts/validate_preparation_api.py",
                     "scripts/validate_element_extension.py", "scripts/validate_element_md.py", "scripts/plot_element_md.py",
                     "scripts/validate_release_numerics.py", "scripts/package_preparation_elements.py"):
            if name not in sources:
                archive.write(ROOT/name, "source/"+name)
    with zipfile.ZipFile(archive_path) as archive:
        assert archive.testzip() is None
        entries = len(archive.namelist())
    result = {"path": str(archive_path), "bytes": archive_path.stat().st_size,
              "sha256": digest(archive_path), "entries": entries, "CRC_passed": True,
              "verified_source_files": len(sources)}
    (OUT/"archive.json").write_text(json.dumps(result, indent=2)+"\n", encoding="utf-8")
    print(json.dumps(result, indent=2))


if __name__ == "__main__":
    main()
