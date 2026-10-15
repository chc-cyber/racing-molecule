"""Audit and package the completed MD release checks without rerunning MD."""
from datetime import datetime, timezone
import hashlib
import json
from pathlib import Path
import zipfile

ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / "outputs/md-release-validation"


def digest(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def main():
    reports = {"numerics": OUT / "numerics/report-all.json",
               "races": OUT / "races/report.json",
               "runtime": OUT / "runtime/port-8781/report.json",
               "normal": OUT / "bonded-normal-check.json"}
    data = {name: json.loads(path.read_text(encoding="utf-8")) for name, path in reports.items()}
    assert data["numerics"]["status"] == data["races"]["status"] == "passed"
    assert data["runtime"]["passed"]
    assert len(data["races"]["trials"]) == data["races"]["planned_trials"] == 30
    assert all(trial["status"] == "passed" for trial in data["races"]["trials"])
    sources = {}
    for report in data.values():
        for name, expected in report.get("source_sha256", {}).items():
            assert name not in sources or sources[name] == expected, f"Reports used different sources: {name}"
            assert digest(ROOT / name) == expected, f"Source changed after validation: {name}"
            sources[name] = expected
    for trial in data["races"]["trials"]:
        coordinates = Path(trial["coordinates_file"])
        if not coordinates.is_absolute():
            coordinates = OUT / coordinates
        assert digest(coordinates) == trial["coordinates_sha256"]

    nve = data["numerics"]["reports"]["nve"]["cases"]
    pairs = [case["dt_comparisons"] for case in nve]
    preparation = {trial["fixture"]: trial["preparation"] for trial in data["races"]["trials"] if trial["config"]["phase"] == "race"}
    summary = {"schema": "md-release-validation-v1", "completed_requested_checks": True,
        "created_utc": datetime.now(timezone.utc).isoformat(),
        "default_configuration": data["runtime"]["defaults"], "default_drive": data["runtime"]["drive"],
        "game_timestep_fs": .5, "thermostat": "Unchanged per-atom Langevin, damping 350 fs",
        "mass_policy": "Ignore isotope labels only in private race copy; use identical standard-element masses throughout",
        "surface_normal_policy": "bonded-saip-game-v1; original distance-selected reference path retained",
        "evidence": {name: {"file": str(path.relative_to(OUT)).replace('\\', '/'), "sha256": digest(path)} for name, path in reports.items()},
        "source_sha256": sources,
        "paired_NVE": {"pairs": len(nve), "duration_ps_each": 5., "initial_temperature_K": 600.,
            "same_initial_positions_and_velocities": True,
            "coarse_max_energy_error_percent_initial_ke_range": [min(p[0]["max_drift_fraction_initial_ke"] for p in pairs)*100, max(p[0]["max_drift_fraction_initial_ke"] for p in pairs)*100],
            "fine_max_energy_error_percent_initial_ke_range": [min(p[1]["max_drift_fraction_initial_ke"] for p in pairs)*100, max(p[1]["max_drift_fraction_initial_ke"] for p in pairs)*100],
            "fine_to_coarse_error_ratio_range": [min(c["fine_to_coarse_energy_drift_ratio"] for c in nve), max(c["fine_to_coarse_energy_drift_ratio"] for c in nve)]},
        "race_and_zero_drive_summary": data["races"]["summary"],
        "preparation_diagnostics": preparation,
        "remaining_work": [
            "Improve joint UFF/surface preparation convergence. Most representatives hit 350 FIRE iterations; C60 met energy tolerance but not the strict force tolerance. Numerical MD stability is not proof of adequate minimization or equilibration.",
            "Near-collinear regularization is conservative but can still yield large forces. These samples do not establish stability for every user-designed geometry at 0.5 fs.",
            "The sp/sp3 surface-axis rule is an explicit game adaptation, not experimental validation of SAIP parameter transferability.",
            "The validation covers five neutral C/H shapes and three seeds, not charged cars, all supported sizes, experimental race rankings or converged diffusion coefficients. Contact-dependent friction and anisotropic polarization remain deferred."
        ],
        "preview": {"url": data["runtime"]["base_url"], "server_pid": data["runtime"]["server_pid"], **data["runtime"]["preview"]}}
    (OUT / "summary.json").write_text(json.dumps(summary, indent=2, ensure_ascii=False)+"\n", encoding="utf-8")

    archive = ROOT / "outputs/md-release-validation-results.zip"
    files = []
    for folder in ("inputs", "races", "numerics", "figures", "runtime/port-8781"):
        files.extend(path for path in (OUT/folder).rglob("*") if path.is_file() and path.suffix != ".pyc")
    files.extend(OUT/name for name in ("summary.json", "bonded-normal-check.json", "python-tests.log"))
    with zipfile.ZipFile(archive, "w", zipfile.ZIP_DEFLATED, compresslevel=6) as package:
        for path in sorted(set(files)):
            package.write(path, str(path.relative_to(OUT)).replace('\\', '/'))
        for name in ("scripts/validate_release_runtime.py", "scripts/validate_release_numerics.py",
                     "scripts/validate_release_races.py", "scripts/plot_release_races.py",
                     "scripts/package_md_release_validation.py"):
            path = ROOT/name
            if path.is_file():
                package.write(path, "validation-source/"+name)
        for name in sources:
            if not name.startswith("molecular_car/"):
                continue
            package.write(ROOT/name, "validation-source/"+name)
    with zipfile.ZipFile(archive) as package:
        assert package.testzip() is None
        count = len(package.namelist())
    artifact = {"path": str(archive), "bytes": archive.stat().st_size, "sha256": digest(archive), "entries": count, "crc_passed": True}
    (OUT / "archive.json").write_text(json.dumps(artifact, indent=2)+"\n", encoding="utf-8")
    print(json.dumps({"checks_completed": True, "artifact": artifact}, ensure_ascii=False))


if __name__ == "__main__":
    main()
