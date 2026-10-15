"""Audit persisted element-MD trajectories and plot finite measured diagnostics."""
from __future__ import annotations
import hashlib
import json
import os
from pathlib import Path
import sys

ROOT = Path(__file__).resolve().parents[1]
BASE = ROOT / "outputs/md-preparation-elements/elements"
OUT = BASE / "md"
sys.path.insert(0, str(ROOT / ".local/analysis-libs"))
os.environ["MPLCONFIGDIR"] = str(OUT / "matplotlib-cache")
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D

BLUE, ORANGE, GREEN = "#0072B2", "#D55E00", "#009E73"


def digest(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def audit(report, probe):
    if report["status"] != "passed" or probe["status"] != "passed":
        raise ValueError("Only completed passed reports may be plotted")
    for name, expected in report["source_sha256"].items():
        if digest(ROOT / name) != expected:
            raise ValueError("MD source/input fingerprint changed: " + name)
    prepared = {c["name"]: c for c in probe["cases"] if c["preparation_status"] == "passed"}
    rows = []
    for trial in [*report["driven"], *report["nve"]]:
        path = ROOT / trial["trajectory_file"]
        if digest(path) != trial["trajectory_sha256"]:
            raise ValueError("Coordinate fingerprint changed: " + trial["key"])
        with np.load(path, allow_pickle=False) as data:
            expected_mass = np.array([a["mass"] for a in prepared[trial["fixture"]]["prepared_atoms"]])
            masses = data["masses_u"]
            mass_error = float(np.abs(masses - expected_mass).max())
            if mass_error > 1e-8:
                raise ValueError("Stored actual integration masses disagree with prepared masses")
            pairs = data["bonds_zero_based"]
            methods = trial.get("methods", [trial])
            summaries = []
            for method in methods:
                suffix = str(method.get("dt_fs", "")).replace(".", "_")
                x = data["positions_dt_" + suffix] if "dt_fs" in method else data["positions_A"]
                v = data["velocities_dt_" + suffix] if "dt_fs" in method else data["velocities_A_fs"]
                if not np.isfinite(x).all() or not np.isfinite(v).all() or x.shape != v.shape:
                    raise ValueError("Saved trajectory has invalid positions/velocities")
                if x.shape != (len(method["samples"]), len(masses), 3):
                    raise ValueError("Saved trajectory dimensions do not match samples")
                time_error = float(np.abs(data["times_ps"] - np.array([s["time_ps"] for s in method["samples"]])).max())
                com = (x * masses[None, :, None]).sum(axis=1) / masses.sum()
                com_error = float(np.abs(com - np.array([s["COM_A"] for s in method["samples"]])).max())
                lengths = np.linalg.norm(x[:, pairs[:, 0]] - x[:, pairs[:, 1]], axis=2)
                ratios = lengths / lengths[0]
                ratio_range = [float(ratios.min()), float(ratios.max())]
                ratio_error = float(np.abs(np.array(ratio_range) - method["bond_ratio_range"]).max())
                min_heights = x[:, :, 2].min(axis=1)
                max_heights = x[:, :, 2].max(axis=1)
                height_error = max(float(np.abs(min_heights - np.array([s["min_atom_height_A"] for s in method["samples"]])).max()),
                                   float(np.abs(max_heights - np.array([s["max_atom_height_A"] for s in method["samples"]])).max()))
                if max(time_error, com_error, ratio_error, height_error) > 1e-10:
                    raise ValueError("Persisted full-coordinate audit disagrees with scalar report")
                if "dt_fs" in method:
                    if not np.array_equal(x[0], data["initial_positions_A"]) or not np.array_equal(v[0], data["initial_velocities_A_fs"]):
                        raise ValueError("Saved paired NVE starts do not match bit for bit")
                summaries.append({"dt_fs": method.get("dt_fs", .5), "sample_count": len(x),
                                  "mass_weighted_COM_max_error_A": com_error,
                                  "bond_ratio_report_max_error": ratio_error,
                                  "atom_height_report_max_error_A": height_error,
                                  "time_report_max_error_ps": time_error,
                                  "finite_full_positions_and_velocities": True,
                                  "paired_initial_xv_bit_identical": True if "dt_fs" in method else None})
            rows.append({"key": trial["key"], "coordinate_sha256": digest(path), "mass_max_error_u": mass_error,
                         "methods": summaries, "passed": True})
    return {"status": "passed", "source_and_coordinate_hashes_passed": True,
            "trajectories": len(rows), "drive_trajectories": len(report["driven"]), "paired_nve_trajectories": len(report["nve"]),
            "full_coordinate_samples": sum(m["sample_count"] for r in rows for m in r["methods"]), "cases": rows}


def plot(report):
    plt.rcParams.update({"font.family": "DejaVu Sans", "font.size": 9, "axes.spines.top": False,
                         "axes.spines.right": False, "svg.fonttype": "none"})
    driven = report["driven"]
    labels = [("P(V)" if r["target_element"] == "P" else r["target_element"]) + ("*" if r["surface"] == "fourier-au" else "") for r in driven]
    fig, axes = plt.subplots(2, 2, figsize=(14.5, 10.5))
    fig.subplots_adjust(left=.075, right=.975, top=.86, bottom=.19, wspace=.28, hspace=.46)
    fig.suptitle("Mixed-element game MD: preparation and finite numerical acceptance", fontsize=16, y=.982)
    fig.text(.5, .943, "600 K | Fixed gradient S = 0.75 | 0.5 fs game timestep | 15 driven 5 ps trials (one seed each) | 8 paired NVE checks", ha="center", fontsize=10)
    x = np.arange(len(driven))
    ax = axes[0, 0]
    initial = [r["preparation"]["game_metadata"]["initial"]["global_force_norm_kcal_mol_A"] for r in driven]
    final = [r["preparation"]["independent"]["global_force_two_norm_kcal_mol_A"] for r in driven]
    for index, (before, after) in enumerate(zip(initial, final)):
        ax.plot([index, index], [before, after], color="#aaaaaa", lw=.8)
    ax.scatter(x, initial, color=ORANGE, marker="x", label="Before joint minimization", s=30)
    ax.scatter(x, final, color=BLUE, label="Independent final force norm", s=25)
    ax.axhline(1e-4, color="#333333", ls="--", lw=.8, label="Strict force gate")
    ax.set(yscale="log", title="A  Deterministic joint force convergence", ylabel="Global 3N force two-norm (kcal mol$^{-1}$ Å$^{-1}$)")
    ax.set_xticks(x, labels)
    ax.legend(fontsize=8, frameon=True, facecolor="white", edgecolor="white", framealpha=.96, loc="center right")
    ax.grid(axis="y", alpha=.17)
    ax = axes[0, 1]
    for snapshot_index, (color, label) in enumerate(zip((BLUE, GREEN, ORANGE), ("Initial", "1 ps", "5 ps"))):
        errors = [r["samples"][int(round((0., 1., 5.)[snapshot_index] / .05))]["full_conservative_force_fd"]["maximum_force_error_kcal_mol_A"] for r in driven]
        ax.scatter(x + (snapshot_index - 1) * .16, errors, s=22, color=color, label=label)
    ax.axhline(.001, color="#333333", ls="--", lw=.8, label="Declared absolute tolerance")
    ax.set(yscale="log", title="B  Complete conservative energy/force finite differences", ylabel="Maximum Cartesian force mismatch (kcal mol$^{-1}$ Å$^{-1}$)")
    ax.set_xticks(x, labels)
    ax.legend(fontsize=8, frameon=False, loc="upper right")
    ax.grid(axis="y", alpha=.17)
    ax = axes[1, 0]
    nve_names = list(dict.fromkeys(r["fixture"] for r in report["nve"]))
    for index, name in enumerate(nve_names):
        pairs = [r for r in report["nve"] if r["fixture"] == name]
        for seed_index, row in enumerate(pairs):
            dx = (seed_index - .5) * .12
            values = [100 * m["maximum_energy_error_fraction_initial_ke"] for m in row["methods"]]
            ax.plot([index + dx - .16, index + dx + .16], values, color="#aaaaaa", lw=.8)
            for dt_index, (color, value) in enumerate(zip((BLUE, ORANGE), values)):
                ax.scatter(index + dx + (dt_index - .5) * .32, value, color=color, s=35,
                           marker=("o", "s")[seed_index])
    ax.set(title="C  2 ps NVE: same coordinates and velocities for each pair", ylabel="Maximum |E(t) − E(0)| / initial kinetic energy (%)", ylim=(0., None))
    ax.set_xticks(np.arange(len(nve_names)), ["Pyridine\nN", "SiMe$_4$\nSi", "Me$_2$S\nS", "Iodobenzene\nI"])
    ax.legend(handles=[Line2D([0], [0], marker="o", ls="", color=BLUE, label="0.5 fs"),
                       Line2D([0], [0], marker="o", ls="", color=ORANGE, label="0.25 fs (diagnostic)")], frameon=False, fontsize=8)
    ax.grid(axis="y", alpha=.17)
    ax = axes[1, 1]
    for color, row in zip(plt.cm.tab20(np.linspace(0., 1., 12)), [r for r in driven if r["surface"] == "periodic-au"]):
        samples = row["samples"]
        t = [s["time_ps"] for s in samples]
        centers = np.array([s["COM_A"] for s in samples])
        ax.plot(t, (centers[:, 0] - centers[0, 0]) / 10., lw=1.1, color=color,
                label="P(V)" if row["target_element"] == "P" else row["target_element"])
    ax.set(title="D  Actual short driven trajectories on periodic Au", xlabel="Driven time (ps)", ylabel="COM forward displacement (nm)", xlim=(0., 5.))
    ax.legend(ncol=4, fontsize=8, frameon=False, loc="upper left")
    ax.grid(alpha=.17)
    fig.text(.075, .112, "All 23 constructors met the independently checked 1e−4 joint-force gate before integration. * Fourier Au; other labels periodic Au. Labels identify whole-molecule fixtures.", fontsize=9)
    fig.text(.075, .080, "NVE uses 600 K initial velocities, zero drive and no bath. Circles/squares are two individual seeds; no confidence intervals or trajectory extrapolation.", fontsize=9)
    fig.text(.075, .048, "C/H retains CHAu SAIP game normals; other elements use provisional nonreactive UFF/Au lattice LJ. Numerical passes do not calibrate adsorption or STM physics.", fontsize=9)
    fig.text(.075, .018, "Supported examples include P(V). Tricoordinate P, I(III), I(V) boundary fixtures and six-coordinate S are explicitly rejected before MD.", fontsize=9)
    for extension in ("png", "svg"):
        fig.savefig(OUT / ("element-md-validation." + extension), dpi=180, facecolor="white")
    plt.close(fig)


def main():
    report_path = OUT / "report-all.json"
    probe_path = BASE / "probe-report.json"
    report = json.loads(report_path.read_text(encoding="utf-8"))
    probe = json.loads(probe_path.read_text(encoding="utf-8"))
    result = audit(report, probe)
    result.update(report_sha256=digest(report_path), probe_sha256=digest(probe_path))
    (OUT / "coordinate-audit.json").write_text(json.dumps(result, indent=2) + "\n", encoding="utf-8")
    plot(report)
    print(json.dumps({"status": result["status"], "trajectories": result["trajectories"],
                      "full_coordinate_samples": result["full_coordinate_samples"],
                      "plot": str(OUT / "element-md-validation.png")}, indent=2))


if __name__ == "__main__":
    main()
