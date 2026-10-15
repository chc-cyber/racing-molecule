"""Scientific figures from accepted, isolated pre/post game benchmarks."""

from __future__ import annotations

import csv
import json
import os
from pathlib import Path
import sys
import zipfile

ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / "outputs/continuous-optimization-check"
FIGURES = OUT / "figures"
sys.path.insert(0, str(ROOT / ".local/analysis-libs"))
os.environ["MPLCONFIGDIR"] = str(OUT / "matplotlib-cache")
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

COLORS = {"reference": "#0072B2", "optimized": "#D55E00"}
LABELS = {"reference": "Before", "optimized": "After"}
FILES = ("benzene.mol", "tetraphenylbenzene.mol")
MOLECULE_LABELS = {FILES[0]: "12-atom benzene", FILES[1]: "52-atom test racer"}


def read(name):
    return json.loads((OUT / name).read_text(encoding="utf-8"))


def save(fig, name):
    fig.savefig(FIGURES / f"{name}.png", dpi=180, facecolor="white")
    fig.savefig(FIGURES / f"{name}.svg", facecolor="white")
    plt.close(fig)


def finish_axes(axes):
    for axis in np.asarray(axes).ravel():
        axis.grid(alpha=.18)
        axis.spines[["top", "right"]].set_visible(False)


def main():
    reports = {name: read(f"{name}-report.json") for name in ("static", "benchmark", "motion", "nve", "profile")}
    if not all(value["passed"] for value in reports.values()):
        raise ValueError("All acceptance reports must pass before preparing final figures")
    # Different stages may run separately, but must refer to the same source.
    hashes = reports["benchmark"]["source_sha256"]
    for name, report in reports.items():
        if report["source_sha256"] != hashes:
            raise ValueError(f"Source mismatch in {name}; repeat that stage")
    FIGURES.mkdir(parents=True, exist_ok=True)
    performance = reports["benchmark"]["summary"]
    summary = {"performance": performance, "static": reports["static"]["maxima"],
               "motion": [{key: row[key] for key in ("file", "atoms", "seed", "max_position_difference_A", "first_ps_max_position_difference_A", "final_center_difference_A")}
                          for row in reports["motion"]["trials"]],
               "nve": [{key: row[key] for key in ("file", "atoms", "max_energy_drift_kcal_mol", "max_paired_energy_difference_kcal_mol")}
                       for row in reports["nve"]["trials"]],
               "third_step_estimates": [{"file": row["file"], "atoms": row["atoms"],
                   "UFF_fraction": row["component_fraction"]["UFF_energy"] + row["component_fraction"]["UFF_gradient"],
                   **row["third_step_estimates"]} for row in reports["profile"]["results"] if row["implementation"] == "optimized"],
               "benchmark_protocol": reports["benchmark"]["protocol"], "source_sha256": hashes}
    (OUT / "summary.json").write_text(json.dumps(summary, ensure_ascii=False, indent=2) + "\n", encoding="utf-8")
    with (OUT / "performance-summary.csv").open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=("file", "atoms", "before_us_per_step", "after_us_per_step", "speedup", "cost_reduction_percent"))
        writer.writeheader()
        for row in performance:
            writer.writerow({"file": row["file"], "atoms": row["atoms"],
                "before_us_per_step": row["median_us_per_step"]["reference"],
                "after_us_per_step": row["median_us_per_step"]["optimized"],
                "speedup": row["speedup"], "cost_reduction_percent": row["cost_reduction_percent"]})
    with (OUT / "motion-samples.csv").open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=("file", "seed", "implementation", "time_ps", "center_x_A", "center_y_A", "center_z_A", "minimum_z_A", "carbon_plane_tilt_deg", "total_energy_kcal_mol"))
        writer.writeheader()
        for trial in reports["motion"]["trials"]:
            for label, samples in trial["samples"].items():
                for sample in samples:
                    writer.writerow({"file": trial["file"], "seed": trial["seed"], "implementation": label,
                        "time_ps": sample["time_ps"], "center_x_A": sample["center_A"][0],
                        "center_y_A": sample["center_A"][1], "center_z_A": sample["center_A"][2],
                        "minimum_z_A": sample["minimum_z_A"], "carbon_plane_tilt_deg": sample["carbon_plane_tilt_deg"],
                        "total_energy_kcal_mol": sample["total_energy_kcal_mol"]})
    with (OUT / "static-energy-force.csv").open("w", encoding="utf-8", newline="") as handle:
        keys = ("file", "atoms", "case", "drive", "strength", "pulse", "reference_energy_kcal_mol", "energy_error_kcal_mol", "max_force_error_kcal_mol_A", "permuted_force_error_kcal_mol_A", "restored_identity_force_error_kcal_mol_A")
        writer = csv.DictWriter(handle, fieldnames=keys, extrasaction="ignore")
        writer.writeheader()
        writer.writerows(reports["static"]["cases"])

    fig, axes = plt.subplots(1, 2, figsize=(11.6, 5.0))
    fig.subplots_adjust(left=.075, right=.985, top=.78, bottom=.22, wspace=.30)
    fig.suptitle("Continuous surface optimization: actual game MD cost", fontsize=15, y=.97)
    fig.text(.5, .885, "Unchanged UFF, Langevin dynamics and 0.5 fs step | Three paired seeds | 600 K, fixed gradient 0.75", ha="center", fontsize=9)
    for index, row in enumerate(performance):
        filename = row["file"]
        subset = [trial for trial in reports["benchmark"]["trials"] if trial["file"] == filename]
        for offset, label in ((-.18, "reference"), (.18, "optimized")):
            x = index + offset
            value = row["median_us_per_step"][label]
            axes[0].bar(x, value, width=.32, color=COLORS[label], label=LABELS[label] if index == 0 else None)
            samples = [trial["median_us_per_step"][label] for trial in subset]
            axes[0].scatter(np.full(len(samples), x), samples, color="#222222", s=20, zorder=4)
            axes[0].text(x, max(samples) * 1.045, f"{value:.1f}", ha="center", fontsize=9)
        axes[1].bar(index, row["speedup"], color="#009E73", width=.55)
        axes[1].scatter(np.full(len(subset), index), [trial["speedup"] for trial in subset], color="#222222", s=20, zorder=4)
        axes[1].text(index, max(trial["speedup"] for trial in subset) * 1.04, f"{row['speedup']:.2f}x", ha="center", fontsize=11)
    axes[0].set(title="A   Full propagation cost", ylabel="Wall-clock cost (µs / step)")
    axes[1].set(title="B   Measured full-MD speedup", ylabel="Before / after (higher is faster)")
    for axis in axes:
        axis.set_xticks(np.arange(2), ["12-atom\nbenzene", "52-atom\ntest racer"])
        axis.set_ylim(0, axis.get_ylim()[1] * 1.18)
    axes[0].legend(frameon=False, fontsize=9)
    fig.text(.075, .055, "1000-step warmup + three alternating 2000-step intervals per seed. Preparation, JIT, frames and rendering excluded.", fontsize=8.4)
    finish_axes(axes)
    save(fig, "continuous-performance")

    fig, axes = plt.subplots(1, 3, figsize=(14., 4.9))
    fig.subplots_adjust(left=.055, right=.985, top=.77, bottom=.22, wspace=.35)
    fig.suptitle("Energy and Cartesian forces: same-coordinate verification", fontsize=15, y=.97)
    fig.text(.5, .875, "Actual callback including UFF, surface, drivers and walls | Gradient, tip, pulse, boundaries and atom-order changes", ha="center", fontsize=9)
    cases = reports["static"]["cases"]
    for filename, color, marker in ((FILES[0], "#0072B2", "o"), (FILES[1], "#D55E00", "s")):
        subset = [row for row in cases if row["file"] == filename]
        indices = np.array([row["case"] for row in subset]) + 1
        axes[0].scatter(indices, [max(abs(row["energy_error_kcal_mol"]), 1e-16) for row in subset], color=color, marker=marker, label=MOLECULE_LABELS[filename], s=32)
        axes[1].scatter(indices, [max(row["max_force_error_kcal_mol_A"], 1e-16) for row in subset], color=color, marker=marker, s=32)
    for axis in axes[:2]:
        axis.set(yscale="log", xlabel="Configuration index")
        axis.set_xticks(range(1, 7))
        axis.axhline(1e-9, color="#777777", lw=1, ls="--")
    axes[0].set(title="A   Before/after energy difference", ylabel="Absolute difference (kcal / mol)")
    axes[1].set(title="B   Maximum force difference", ylabel="Difference (kcal / mol / Å)")
    axes[0].legend(frameon=False, fontsize=7.5)
    derivatives = reports["static"]["cartesian_derivative_checks"]
    for offset, label in ((-.12, "reference"), (.12, "optimized")):
        subset = [row for row in derivatives if row["implementation"] == label]
        axes[2].scatter(np.arange(len(subset)) + 1 + offset,
            [row["max_force_derivative_error_kcal_mol_A"] for row in subset],
            color=COLORS[label], marker="o", label=LABELS[label], s=30)
    axes[2].set(title="C   Complete finite-difference check", xlabel="Derivative fixture index", ylabel="Force vs −dE/dr (kcal / mol / Å)", yscale="log")
    axes[2].axhline(1e-4, color="#777777", lw=1, ls="--")
    axes[2].legend(frameon=False, fontsize=8)
    axes[2].set_xticks(range(1, 6))
    fig.text(.055, .055, "Dashed lines: declared acceptance thresholds. Values below 1e−16 are displayed at 1e−16. Finite-difference displacement: 1e−5 Å.", fontsize=8.4)
    finish_axes(axes)
    save(fig, "continuous-energy-force-check")

    fig, axes = plt.subplots(2, 4, figsize=(14., 7.9))
    fig.subplots_adjust(left=.055, right=.985, top=.82, bottom=.12, hspace=.48, wspace=.34)
    fig.suptitle("Same model, paired racing trajectories", fontsize=15, y=.97)
    fig.text(.5, .9, "10 ps | Three paired seeds per molecule | Identical starting positions and velocities | 600 K, fixed gradient 0.75", ha="center", fontsize=9)
    for index, filename in enumerate(FILES):
        trials = [row for row in reports["motion"]["trials"] if row["file"] == filename]
        for label in ("reference", "optimized"):
            times = np.array([row["time_ps"] for row in trials[0]["samples"][label]])
            centers = np.array([[row["center_A"] for row in trial["samples"][label]] for trial in trials])
            displacement = (centers - centers[:, :1, :]) / 10
            height = centers[:, :, 2]
            tilt = np.array([[row["carbon_plane_tilt_deg"] for row in trial["samples"][label]] for trial in trials])
            for axis, value, title, ylabel in zip(axes[index], (displacement[:, :, 0], displacement[:, :, 1], height, tilt),
                ("Forward displacement", "Lateral displacement", "Center height", "Carbon-plane tilt"),
                ("Δx (nm)", "Δy (nm)", "z (Å)", "Tilt (degree)")):
                mean, sd = value.mean(axis=0), value.std(axis=0, ddof=1)
                axis.fill_between(times, mean - sd, mean + sd, color=COLORS[label], alpha=.10)
                axis.plot(times, mean, color=COLORS[label], lw=1.6, ls="-" if label == "reference" else "--", label=LABELS[label])
                axis.set(title=f"{chr(65 + index * 4 + list(axes[index]).index(axis))}   {title}", xlabel="Simulation time (ps)", ylabel=ylabel)
        axes[index, 0].text(.04, .95, MOLECULE_LABELS[filename], transform=axes[index, 0].transAxes, va="top", fontsize=8.4)
        axes[index, 3].legend(frameon=False, fontsize=8)
    fig.text(.055, .035, "Lines: means of three seeds; shading: ±1 sample SD. Overlap is expected for unchanged equations. A nonplanar carbon fit is a descriptive angle.", fontsize=8.3)
    finish_axes(axes)
    save(fig, "continuous-racing-comparison")

    fig, axes = plt.subplots(1, 3, figsize=(14., 5.0))
    fig.subplots_adjust(left=.055, right=.985, top=.76, bottom=.24, wspace=.37)
    fig.suptitle("Residual cost and the potential of native UFF integration", fontsize=15, y=.97)
    fig.text(.5, .87, "Separate instrumented diagnostic | Migration not implemented | Performance scenarios are Amdahl estimates", ha="center", fontsize=9)
    profiles = reports["profile"]["results"]
    keys = ("UFF_energy", "UFF_gradient", "surface_drive_walls", "other_inside_callback", "outside_callback")
    colors = ("#0072B2", "#56B4E9", "#E69F00", "#CC79A7", "#777777")
    names = ("UFF energy", "UFF gradient", "External field", "Other callback", "Outside callback")
    for index, filename in enumerate(FILES):
        subset = [next(row for row in profiles if row["file"] == filename and row["implementation"] == label) for label in ("reference", "optimized")]
        bottom = np.zeros(2)
        for key, color, name in zip(keys, colors, names):
            values = np.array([row["components_us_per_step"][key] for row in subset])
            axes[index].bar([0, 1], values, bottom=bottom, width=.60, color=color, label=name)
            bottom += values
        axes[index].set(title=f"{chr(65 + index)}   {MOLECULE_LABELS[filename]}", ylabel="Instrumented cost (µs / step)")
        axes[index].set_xticks([0, 1], ["Before", "After"])
    for filename, color in ((FILES[0], "#0072B2"), (FILES[1], "#D55E00")):
        row = next(row for row in profiles if row["file"] == filename and row["implementation"] == "optimized")
        scenarios = row["third_step_estimates"]["UFF_only_scenarios"]
        axes[2].plot([value["assumed_UFF_speedup"] for value in scenarios], [value["whole_MD_speedup"] for value in scenarios], "o-", color=color, label=MOLECULE_LABELS[filename])
        axes[2].axhline(row["third_step_estimates"]["UFF_only_zero_cost_upper_bound_speedup"], color=color, ls=":", lw=1)
    axes[2].set(title="C   Hypothetical UFF-only acceleration", xlabel="Assumed UFF speedup", ylabel="Predicted additional full-MD speedup")
    axes[2].set_xticks([2, 4, 10], ["2x", "4x", "10x"])
    axes[2].legend(frameon=False, fontsize=7.5)
    axes[0].legend(frameon=False, fontsize=7.8, loc="upper right")
    fig.text(.055, .05, "Diagnostic includes timer overhead and one seed. Dotted ceilings remove all UFF cost; transfer/callback savings are not assumed here.", fontsize=8.5)
    finish_axes(axes)
    save(fig, "continuous-residual-cost")

    fig, axes = plt.subplots(1, 2, figsize=(11.6, 4.9))
    fig.subplots_adjust(left=.075, right=.985, top=.78, bottom=.22, wspace=.32)
    fig.suptitle("NVE conservation and short-trajectory numerical agreement", fontsize=15, y=.97)
    fig.text(.5, .88, "Time step unchanged at 0.5 fs | NVE: 10 K, no thermostat or drive | Racing: identical paired stochastic seeds", ha="center", fontsize=9)
    for filename, color in ((FILES[0], "#0072B2"), (FILES[1], "#D55E00")):
        nve = next(row for row in reports["nve"]["trials"] if row["file"] == filename)
        for label in ("reference", "optimized"):
            values = nve["samples"][label]
            energies = np.array([row["total_energy_kcal_mol"] for row in values])
            axes[0].plot([row["time_ps"] for row in values], energies - energies[0], color=color,
                ls="-" if label == "reference" else "--", label=f"{nve['atoms']} atoms, {LABELS[label].lower()}")
        trials = [row for row in reports["motion"]["trials"] if row["file"] == filename]
        maximum = np.max([[sample["max_position_difference_A"] for sample in row["comparison"]] for row in trials], axis=0)
        axes[1].plot([row["time_ps"] for row in trials[0]["comparison"]], maximum, color=color, label=f"{nve['atoms']} atoms")
    axes[0].set(title="A   NVE total-energy drift", xlabel="Simulation time (ps)", ylabel="E(t) − E(0) (kcal / mol)")
    axes[1].set(title="B   Maximum paired position difference", xlabel="Simulation time (ps)", ylabel="Maximum Cartesian difference (Å)")
    if all(row["max_position_difference_A"] == 0 for row in reports["motion"]["trials"]):
        axes[1].set_yticks([0.])
        axes[1].text(.5, .25, "Exactly 0 at every sampled step", transform=axes[1].transAxes, ha="center", fontsize=10)
    for axis in axes:
        axis.legend(frameon=False, fontsize=8)
    fig.text(.075, .055, "NVE drift is the pre-existing integrator error. Racing position difference is the maximum over all atoms, coordinates and three paired seeds.", fontsize=8.3)
    finish_axes(axes)
    save(fig, "continuous-conservation-agreement")

    archive = OUT / "continuous-optimization-results.zip"
    with zipfile.ZipFile(archive, "w", compression=zipfile.ZIP_DEFLATED) as handle:
        for path in sorted(OUT.rglob("*")):
            if (path.is_file() and path != archive
                    and not {"matplotlib-cache", "__pycache__", "diagnostics"}.intersection(path.parts)
                    and not (path.name.startswith("report-") and path.name != "report-all.json")
                    and path.suffix.lower() in (".json", ".csv", ".png", ".svg", ".py", ".md")):
                handle.write(path, path.relative_to(OUT).as_posix())
        for name in ("README.md", "THIRD_PARTY.md", "scripts/lammps_runtime.py",
                     "scripts/verify_continuum_optimization.py", "scripts/plot_continuum_optimization.py",
                     "scripts/verify_continuum_game.py", "scripts/verify_continuum_game_browser.cjs",
                     "scripts/estimate_native_path.py", "tests/test_fields_kernel.py",
                     "tests/test_force_workflow.py", "examples/benzene.mol", "examples/tetraphenylbenzene.mol",
                     "molecular_car/engine.py", "molecular_car/external_fields.py",
                     "molecular_car/fields_kernel.py", "molecular_car/fields.py", "molecular_car/molecule.py"):
            path = ROOT / name
            if path.exists():
                handle.write(path, "implementation/" + name)
        for directory in ("continuous-game-check", "native-path-estimate"):
            supplement = ROOT / "outputs" / directory
            for path in sorted(supplement.glob("*")):
                if path.is_file() and path.suffix.lower() in (".json", ".png", ".md"):
                    handle.write(path, f"supplement/{directory}/{path.name}")
    with zipfile.ZipFile(archive) as handle:
        if handle.testzip() is not None:
            raise ValueError("Results archive CRC verification failed")
    print(json.dumps(summary, ensure_ascii=False, indent=2))
    print(f"Figures and verified archive saved: {archive}")


if __name__ == "__main__":
    main()
