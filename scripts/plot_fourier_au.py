"""Scientific summary of accepted Fourier Au proxy benchmarks."""

from __future__ import annotations

import csv
import json
import os
from pathlib import Path
import sys
import zipfile

ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / "outputs/fourier-au-check"
FIGURES = OUT / "figures"
sys.path.insert(0, str(ROOT / ".local/analysis-libs"))
os.environ["MPLCONFIGDIR"] = str(OUT / "matplotlib-cache")
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np


def read(name):
    return json.loads((OUT / name).read_text(encoding="utf-8"))


def save(fig, name):
    fig.savefig(FIGURES / f"{name}.png", dpi=180, facecolor="white")
    fig.savefig(FIGURES / f"{name}.svg", facecolor="white")
    plt.close(fig)


def main():
    static, holdout, paired = read("static-report.json"), read("holdout-report.json"), read("paired-md-report.json")
    if not static["passed"] or not holdout["passed"]:
        raise ValueError("Accept the numerical proxy before plotting its performance.")
    FIGURES.mkdir(parents=True, exist_ok=True)
    labels = {"continuum": "Original\ncontinuum", "periodic": "Periodic cache\n+ exact direction", "fourier": "Fourier proxy\n+ smooth fallback"}
    colors = {"continuum": "#0072B2", "periodic": "#E69F00", "fourier": "#009E73"}
    rows = []
    for name in labels:
        report = paired["reports"][name]
        samples = np.array([trial["median_wall_seconds_per_step"] * 1e6 for trial in report["trials"]])
        rows.append({"mode": name, "median_us_per_step": float(np.median(samples)),
                     "mean_us_per_step": float(np.mean(samples)), "sample_sd_us_per_step": float(np.std(samples, ddof=1)),
                     "seed_medians_us": samples.tolist(),
                     "NVE_maximum_energy_drift_eV": report["nve"]["maximum_energy_drift_eV"]})
    summary = {"performance": rows, "protocol": paired["protocol"],
               "fourier_speedup_vs_periodic": paired["fourier_speedup_vs_periodic"],
               "fourier_cost_vs_continuum": paired["fourier_cost_vs_continuum"],
               "static_max_force_error_eV_A": static["max_force_error_eV_A"],
               "holdout_max_force_error_eV_A": holdout["max_force_error_eV_A"],
               "holdout_max_force_derivative_error_eV_A": holdout["max_force_derivative_error_eV_A"],
               "surface": static["surface"]}
    (OUT / "summary.json").write_text(json.dumps(summary, ensure_ascii=False, indent=2) + "\n", encoding="utf-8")
    with (OUT / "performance-summary.csv").open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=[key for key in rows[0] if key != "seed_medians_us"], extrasaction="ignore")
        writer.writeheader(); writer.writerows(rows)
    fig, axes = plt.subplots(1, 2, figsize=(11.6, 4.9))
    fig.subplots_adjust(left=.075, right=.985, top=.79, bottom=.23, wspace=.30)
    fig.suptitle("Static Au Fourier proxy: measured MD performance", fontsize=15, y=.97)
    fig.text(.5, .89, "Same MD box and callback | 12-atom benzene | 600 K | 0.5 fs | Three paired seeds", ha="center", fontsize=9)
    indices = np.arange(len(rows))
    costs = np.array([row["median_us_per_step"] for row in rows])
    axes[0].bar(indices, costs, color=[colors[row["mode"]] for row in rows], width=.6)
    axes[0].set(title="A   Full propagation cost", ylabel="Wall-clock cost (µs / step)")
    ratios = costs / costs[0]
    axes[1].bar(indices, ratios, color=[colors[row["mode"]] for row in rows], width=.6)
    axes[1].set(title="B   Cost relative to original continuum", ylabel="Relative cost (lower is faster)")
    for index, row in enumerate(rows):
        axes[0].scatter(np.full(3, index), row["seed_medians_us"], color="black", s=20, zorder=4)
        axes[0].text(index, max(row["seed_medians_us"])*1.035, f"{costs[index]:.1f}", ha="center", fontsize=10)
        axes[1].text(index, ratios[index]*1.035, f"{ratios[index]:.2f}×", ha="center", fontsize=10)
    axes[0].set_ylim(0, max(max(row["seed_medians_us"]) for row in rows)*1.2)
    axes[1].set_ylim(0, max(ratios)*1.2)
    for axis in axes:
        axis.set_xticks(indices, [labels[row["mode"]] for row in rows], fontsize=9)
        axis.grid(axis="y", alpha=.2)
    fig.text(.075, .045, "1000-step warmup + three 2000-step intervals per seed. Preparation, compilation, sampling and rendering excluded.", fontsize=8.6)
    save(fig, "fourier-performance")

    cases = static["cases"]
    fig, axes = plt.subplots(1, 2, figsize=(11.6, 4.9))
    fig.subplots_adjust(left=.075, right=.985, top=.78, bottom=.21, wspace=.30)
    fig.suptitle("Fourier proxy compared with native frozen-Au SAIP", fontsize=15, y=.97)
    fig.text(.5, .89, "Includes conservative reference fallback outside the approximation domain", ha="center", fontsize=9)
    for tilt, color in ((0., "#0072B2"), (30., "#009E73"), (60., "#E69F00"), (85., "#CC79A7")):
        subset = [case for case in cases if case["atoms"] == 12 and case["tilt_deg"] == tilt]
        axes[0].scatter([case["min_z_A"] for case in subset], [case["max_force_error_eV_A"] for case in subset],
                        s=28, color=color, label=f"Benzene {tilt:g} deg")
    for atoms, color, marker in ((52, "#222222", "D"), (212, "#D55E00", "X")):
        subset = [case for case in cases if case["atoms"] == atoms]
        axes[0].scatter([case["min_z_A"] for case in subset], [case["max_force_error_eV_A"] for case in subset],
                        s=52, color=color, marker=marker, label=f"{atoms}-atom molecule")
    independent = holdout["cases"]
    axes[1].scatter([case["min_z_A"] for case in independent], [case["max_force_error_eV_A"] for case in independent],
                    s=26, color="#0072B2", label="Independent geometries")
    for axis, title in zip(axes, ("A   Established 36 configurations", "B   Random orientations and H-normal holdouts")):
        axis.set(title=title, xlabel="Minimum atom height (Å)", ylabel="Maximum Cartesian force error (eV / Å)", yscale="log")
        axis.axhline(.002, color="#777777", lw=1, ls="--", label="Acceptance 0.002 eV / Å")
        axis.grid(alpha=.2)
        axis.legend(frameon=False, fontsize=7.8, ncol=2)
    fig.text(.075, .045, "Validation uses identical molecular coordinates; reference SAIP/OPT is disabled. Error bounds apply to sampled configurations.", fontsize=8.4)
    save(fig, "fourier-native-accuracy")

    motion = read("paired-motion-report.json")
    if not motion["passed"]:
        raise ValueError("Motion checks have not passed.")
    fig, axes = plt.subplots(2, 2, figsize=(11.6, 7.7))
    fig.subplots_adjust(left=.075, right=.985, top=.82, bottom=.13, hspace=.42, wspace=.30)
    fig.suptitle("Static Au proxies: short MD motion and contact", fontsize=15, y=.97)
    fig.text(.5, .915, "12-atom benzene | 600 K | 0.75 fixed gradient | 10 ps | Three paired initial velocity seeds", ha="center", fontsize=9)
    handles, motion_rows = [], []
    for name, trials in motion["reports"].items():
        time_ps = np.array([sample["time_ps"] for sample in trials[0]["samples"]])
        centers = np.array([[sample["center_angstrom"] for sample in trial["samples"]] for trial in trials])
        displacement = (centers - centers[:, :1, :]) / 10
        height = centers[:, :, 2]
        tilt = np.array([[sample["carbon_plane_tilt_deg"] for sample in trial["samples"]] for trial in trials])
        values = (displacement[:, :, 0], displacement[:, :, 1], height, tilt)
        for index, (axis, samples) in enumerate(zip(axes.ravel(), values)):
            if index == 1:
                # RMS across seeds is not a sample standard deviation band.
                mean = np.sqrt(np.mean(samples**2, axis=0))
                line, = axis.plot(time_ps, mean, color=colors[name], lw=1.7)
            else:
                mean, sd = np.mean(samples, axis=0), np.std(samples, axis=0, ddof=1)
                line, = axis.plot(time_ps, mean, color=colors[name], lw=1.7)
                axis.fill_between(time_ps, mean-sd, mean+sd, color=colors[name], alpha=.13)
            if index == 0:
                handles.append(line)
        for trial_index, trial in enumerate(trials):
            for sample_index, sample in enumerate(trial["samples"]):
                motion_rows.append({"mode": name, "seed": trial["description"]["seed"], "time_ps": sample["time_ps"],
                                    "forward_displacement_nm": float(displacement[trial_index, sample_index, 0]),
                                    "lateral_displacement_nm": float(displacement[trial_index, sample_index, 1]),
                                    "COM_height_A": float(height[trial_index, sample_index]),
                                    "carbon_tilt_deg": float(tilt[trial_index, sample_index])})
    panels = (("A   Forward displacement", "Displacement (nm)"),
              ("B   Lateral wandering", "RMS lateral displacement (nm)"),
              ("C   Surface contact", "Center-of-mass height (Å)"),
              ("D   Molecular orientation", "Carbon plane tilt (degrees)"))
    for axis, (title, ylabel) in zip(axes.ravel(), panels):
        axis.set(title=title, xlabel="Time (ps)", ylabel=ylabel, xlim=(0, 10))
        axis.grid(alpha=.2)
    fig.legend(handles, [labels[name].replace("\n", " ") for name in motion["reports"]],
               loc="upper center", bbox_to_anchor=(.5, .885), ncol=3, frameon=False, fontsize=9)
    fig.text(.075, .037, "A/C/D: mean ± sample SD (n=3). B: RMS across seeds. A short dynamic check, not converged race or diffusion statistics.", fontsize=8.4)
    save(fig, "fourier-motion-comparison")
    with (OUT / "motion-summary.csv").open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(motion_rows[0]))
        writer.writeheader(); writer.writerows(motion_rows)

    destination = OUT / "fourier-results.zip"
    with zipfile.ZipFile(destination, "w", compression=zipfile.ZIP_DEFLATED) as bundle:
        for folder in (OUT, ROOT / "outputs/fourier-game-check", ROOT / "outputs/fourier-prototypes"):
            for path in sorted(folder.rglob("*")):
                if path.is_file() and path.suffix in (".json", ".csv", ".png", ".svg") and "matplotlib-cache" not in path.parts:
                    bundle.write(path, path.relative_to(ROOT / "outputs"))
    with zipfile.ZipFile(destination) as bundle:
        if bundle.testzip() is not None:
            raise AssertionError("Result bundle CRC failed.")
    print(json.dumps({"figures": str(FIGURES), "bundle": str(destination),
                      "performance": rows}, indent=2))


if __name__ == "__main__":
    main()
