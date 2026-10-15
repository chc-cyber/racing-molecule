"""Independent saved-trajectory audit; no MD, plotting, or timing work."""
from pathlib import Path
import gzip
import hashlib
import json
import math
import statistics

ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / "outputs/parallel-md-validation"
NEW = ROOT / "outputs/c60-proxy-study/races"
BASES = {
    "periodic": ROOT / "outputs/crystal-contrast-guarded/browser-trials",
    "v1": ROOT / "outputs/wasm-fourier-validation/races/browser-trials",
}


def read(path):
    return json.loads(path.read_text())


def metrics(row):
    samples = row["samples"]
    y0 = samples[0]["center_nm"][1]
    dy = [r["center_nm"][1] - y0 for r in samples]
    slow = [samples[k]["center_nm"][0] - samples[k - 20]["center_nm"][0] < .02
            for k in range(20, len(samples))]
    return {"finish_ps": row["finish_time_ps"],
            "lateral_rms_nm": math.sqrt(sum(y * y for y in dy) / len(dy)),
            "max_lateral_nm": max(map(abs, dy)),
            "low_progress_fraction": sum(slow) / len(slow)}


report = read(NEW / "browser-confirm-report.json")
assert report["status"] == "complete" and len(report["trials"]) == 24
audit = {"status": "running", "scope": "24 saved atlas races; historical matched periodic/v1 references, n=3 per condition", "trials": [], "groups": []}
for row in report["trials"]:
    key = row["key"]
    assert row["status"] == "passed" and row["outcome"] == "finish", key
    assert row["competitive_valid"] and not row["sustained_lift"] and not row["right_censored"], key
    assert row["runtime"]["timestep_fs"] == .5
    assert row["runtime"]["preparation"]["accepted"]
    assert "v2" in row["runtime"]["surface"]["implementation"]
    samples = row["samples"]
    assert samples[0]["step"] == 0 and samples[-2]["center_nm"][0] < 20 <= samples[-1]["center_nm"][0]
    for i, s in enumerate(samples):
        assert s["step"] == i * 100 and abs(s["time_ps"] - i * .05) < 1e-10
        assert all(math.isfinite(x) for x in [*s["center_nm"], s["energy_kcal_mol"], s["internal_temperature_K"]])
    assert abs(row["finish_time_ps"] - samples[-1]["time_ps"]) < 1e-10
    assert max(abs(a - b) for a, b in zip(row["finish_crossing_interval_ps"], [samples[-2]["time_ps"], samples[-1]["time_ps"]])) < 1e-10
    run = peak = 0
    for s in samples[1:]:
        run = run + 1 if s["min_z_A"] > 8 else 0
        peak = max(peak, run)
    assert abs(peak * .05 - row["maximum_lift_duration_ps"]) < 1e-10
    assert peak < 100
    assert row["bond_ratio_range"][0] >= .65 and row["bond_ratio_range"][1] <= 1.45
    coordinate_path = NEW / "browser-trials" / row["coordinates_file"]
    coordinates = json.loads(gzip.decompress(coordinate_path.read_bytes()))
    masses = coordinates["masses"]
    total_mass = sum(masses)
    maximum_com_error = 0.
    observed_samples = {s["step"]: s for s in samples}
    for step, positions in zip(coordinates["steps"], coordinates["positions"]):
        assert len(positions) == 3 * len(masses) and step in observed_samples
        com = [sum(m * positions[3 * i + d] for i, m in enumerate(masses)) / total_mass / 10 for d in range(3)]
        maximum_com_error = max(maximum_com_error, max(abs(a - b) for a, b in zip(com, observed_samples[step]["center_nm"])))
    assert maximum_com_error < 1e-10
    item = {"key": key, "car": row["config"]["molecule"], "scale": row["config"]["lattice_scale"],
            "seed": row["config"]["seed"], "max_com_error_nm": maximum_com_error,
            "atlas": metrics(row), "bond_range": row["bond_ratio_range"], "height_range_A": row["atom_height_range_A"],
            "maximum_lift_duration_ps": peak * .05, "comparisons": {}}
    item["proxy_counts"] = [sum(s["proxy_counts"][i] for s in samples) for i in range(3)]
    for label, folder in BASES.items():
        old = read(folder / (key + ".json"))
        assert old["config"] == row["config"] and old["status"] == "passed" and old["outcome"] == "finish"
        item[label] = metrics(old)
        item["comparisons"][label] = {"finish_delta_pct": 100 * (row["finish_time_ps"] / old["finish_time_ps"] - 1),
            "initial_com_difference_nm": max(abs(a-b) for a,b in zip(old["initial_center_nm"], row["initial_center_nm"]))}
    audit["trials"].append(item)

for car in ["benzene-12", "pyridine-11", "test-car-52", "fullerene-c60"]:
    for scale in [1.5, 2]:
        a = [r for r in audit["trials"] if r["car"] == car and r["scale"] == scale]
        assert len(a) == 3
        g = {"car": car, "scale": scale, "n": len(a)}
        counts = [sum(r["proxy_counts"][i] for r in a) for i in range(3)]
        g["proxy_coverage"] = [n / sum(counts) for n in counts]
        for model in ["periodic", "v1", "atlas"]:
            g[model] = {}
            for metric in a[0][model]:
                values = [r[model][metric] for r in a]
                g[model][metric] = {"mean": statistics.mean(values), "sd": statistics.stdev(values)}
        for ref in ["periodic", "v1"]:
            g["delta_vs_" + ref + "_pct"] = 100 * (g["atlas"]["finish_ps"]["mean"] / g[ref]["finish_ps"]["mean"] - 1)
        audit["groups"].append(g)
audit["rankings"] = {str(scale): {m: [g["car"] for g in sorted([g for g in audit["groups"] if g["scale"] == scale], key=lambda g:g[m]["finish_ps"]["mean"])] for m in ["periodic", "v1", "atlas"]} for scale in [1.5, 2]}
audit["max_com_error_nm"] = max(r["max_com_error_nm"] for r in audit["trials"])
audit["status"] = "passed"
audit["interpretation"] = "Passed saved-data integrity/geometry checks, not statistical equivalence. Historical preparation and floating-point arithmetic are not identical; each group has only three seeds. No 600 K acceptance is inferred."
audit["source_hashes"] = {str(p.relative_to(ROOT)): hashlib.sha256(p.read_bytes()).hexdigest() for p in [Path(__file__), NEW / "browser-confirm-report.json"]}
OUT.mkdir(exist_ok=True)
(OUT / "race-review.json").write_text(json.dumps(audit, indent=2) + "\n")
print(json.dumps({"status": audit["status"], "trials":len(audit["trials"]), "max_com_error_nm":audit["max_com_error_nm"], "groups":audit["groups"], "rankings":audit["rankings"]}))
