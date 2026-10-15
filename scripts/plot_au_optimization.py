"""Scientific figures from measured Au optimization and accuracy reports."""

from __future__ import annotations

import csv
import json
import os
from pathlib import Path
import sys
import zipfile

ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / "outputs" / "au-optimization-check"
FAST = ROOT / "outputs" / "periodic-au-check"
FIGURES = OUT / "figures"
sys.path.insert(0, str(ROOT / ".local" / "analysis-libs"))
os.environ["MPLCONFIGDIR"] = str(OUT / "matplotlib-cache")
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np


def save(fig, name):
    fig.savefig(FIGURES / f"{name}.png", dpi=180, facecolor="white")
    fig.savefig(FIGURES / f"{name}.svg", facecolor="white")
    plt.close(fig)


def main():
    FIGURES.mkdir(parents=True, exist_ok=True)
    explicit = json.loads((OUT / "report.json").read_text(encoding="utf-8"))
    cached = json.loads((FAST / "md-report.json").read_text(encoding="utf-8"))
    numpy_cached = json.loads((FAST / "md-report-numpy.json").read_text(encoding="utf-8"))
    accuracy = json.loads((FAST / "static-report.json").read_text(encoding="utf-8"))
    if explicit["status"] != "passed" or not explicit["timings"] or not accuracy["passed"]:
        raise ValueError("Complete native optimization benchmarks before plotting.")
    rows = []
    for trial in explicit["timings"]:
        config = trial["configuration"]
        if config["model"] != "frozen_au":
            continue
        cost = np.array([batch["seconds_per_step"] * 1e6 for batch in trial["batches"]])
        rows.append({"model": config["name"], "style": config["style"],
                     "Au_atoms_in_LAMMPS": trial["Au_atoms"], "median_us_per_step": float(np.median(cost)),
                     "sample_sd_us_per_step": float(np.std(cost, ddof=1)), "cost_samples_us": cost.tolist(),
                     "preparation_seconds": trial["preparation_wall_seconds"],
                     "approved": config["style"] == "base"})
    for label, report in (("numpy", numpy_cached), ("numba", cached)):
        trial = report["trials"][0]
        cost = np.array([batch["seconds_per_step"] * 1e6 for batch in trial["batches"]])
        rows.append({"model": f"cached-periodic-au-{label}", "style": f"cached-{label}",
                     "Au_atoms_in_LAMMPS": 0, "median_us_per_step": float(np.median(cost)),
                     "sample_sd_us_per_step": float(np.std(cost, ddof=1)), "cost_samples_us": cost.tolist(),
                     "preparation_seconds": trial["description"]["preparation_wall_seconds"], "approved": True})
    fast_trial = cached["trials"][0]
    baseline = next(row["median_us_per_step"] for row in rows
                    if row["model"] == "frozen_au-base-big")
    for row in rows:
        row["speedup_vs_baseline"] = baseline / row["median_us_per_step"]
    summary = {"performance": rows, "static_accuracy": {
        key: accuracy[key] for key in ("max_energy_error_eV", "max_force_error_eV_A",
                                      "finite_difference_force_error_eV_A")},
        "cache": fast_trial["description"]["environment"],
        "protocol": "1 thread; 12 atom benzene; 600 K; 0.5 fs; Langevin 350 fs; strength 0.75; three contiguous 2000 step batches after 1000 step warmup; preparation and sampling excluded",
        "NVE_maximum_energy_drift_eV": cached["nve"]["maximum_energy_drift_eV"],
        "compiled_three_seeds_mean_us_per_step": [trial["microseconds_per_step"] for trial in cached["trials"]],
        "source_sha256": json.loads((FAST / "report-all.json").read_text(encoding="utf-8"))["source_sha256"],
        "OPT_approved": explicit["static_checks"].get("OPT_approved", False)}
    (OUT / "optimization-summary.json").write_text(json.dumps(summary, indent=2) + "\n", encoding="utf-8")
    with (OUT / "optimization-summary.csv").open("w", newline="", encoding="utf-8") as handle:
        fields = [key for key in rows[0] if key != "cost_samples_us"]
        writer = csv.DictWriter(handle, fieldnames=fields, extrasaction="ignore")
        writer.writeheader(); writer.writerows(rows)

    colors = ["#757575" if not row["approved"] else
              "#009E73" if row["style"] == "cached-numba" else "#E69F00" if row["style"] == "cached-numpy" else "#0072B2" for row in rows]
    labels = [f"{row['Au_atoms_in_LAMMPS']} Au\n" +
              ("SAIP/OPT*" if row["style"] == "opt" else "Cached\nNumPy" if row["style"] == "cached-numpy" else "Cached\nNumba" if row["style"] == "cached-numba" else "SAIP")
              for row in rows]
    fig, axes = plt.subplots(1, 2, figsize=(11.8, 4.7))
    fig.subplots_adjust(left=.075, right=.985, top=.78, bottom=.24, wspace=.3)
    fig.suptitle("Measured speedup of frozen Au implementations", fontsize=15, y=.97)
    fig.text(.5, .89, "12-atom benzene | 600 K | 0.5 fs | Three timed MD batches, each 2000 steps", ha="center", fontsize=9)
    indices = np.arange(len(rows))
    bars = axes[0].bar(indices, [row["median_us_per_step"] for row in rows], color=colors, width=.65)
    axes[0].set(yscale="log", ylabel="Propagation cost (us / step)", title="A   MD cost (lower is faster)")
    speed_bars = axes[1].bar(indices, [row["speedup_vs_baseline"] for row in rows], color=colors, width=.65)
    axes[1].set(ylabel="Speedup relative to 2880-Au SAIP", title="B   Relative speed")
    for axis in axes:
        axis.set_xticks(indices, labels, fontsize=9)
        axis.grid(axis="y", alpha=.2)
    for index, row in enumerate(rows):
        if not row["approved"]:
            bars[index].set_hatch("///")
            speed_bars[index].set_hatch("///")
        axes[0].scatter(np.full(3, index), row["cost_samples_us"], s=20, color="black", zorder=3)
        axes[0].text(index, row["median_us_per_step"] * 1.15, f"{row['median_us_per_step']:.0f}", ha="center", fontsize=9)
        axes[1].text(index, row["speedup_vs_baseline"] * 1.03, f"{row['speedup_vs_baseline']:.2f}x", ha="center", fontsize=9)
    axes[0].set_ylim(top=max(row["median_us_per_step"] for row in rows) * 1.7)
    axes[1].set_ylim(0, max(row["speedup_vs_baseline"] for row in rows) * 1.18)
    fig.text(.07, .035, "*OPT failed the Cartesian force check in the installed build; its timings are experimental and it is disabled by default.", fontsize=8.8)
    save(fig, "optimization-performance")

    cases = accuracy["cases"]
    fig, axes = plt.subplots(1, 2, figsize=(11.8, 4.7))
    fig.subplots_adjust(left=.075, right=.985, top=.80, bottom=.20, wspace=.29)
    fig.suptitle("Cached periodic surface compared with native SAIP", fontsize=15, y=.97)
    fig.text(.5, .89, "Identical coordinates | Six relaxed Au layers | 16 Å tapered potential | Local-normal forces retained", ha="center", fontsize=9)
    tilt_colors = {0.: "#0072B2", 30.: "#009E73", 60.: "#E69F00", 85.: "#CC79A7", None: "#222222"}
    for tilt, color in tilt_colors.items():
        subset = [case for case in cases if case["tilt_deg"] == tilt and case["atoms"] == 12]
        if not subset:
            continue
        axes[0].scatter([case["min_z_A"] for case in subset],
                        [case["energy_error_eV"] * 1000 for case in subset],
                        label=f"Benzene {tilt:g} deg", s=28, color=color, alpha=.8)
        axes[1].scatter([case["min_z_A"] for case in subset],
                        [max(case["max_force_error_eV_A"], 1e-13) for case in subset], s=28, color=color, alpha=.8)
    for atoms, marker, color in ((52, "D", "#222222"), (212, "X", "#D55E00")):
        subset = [case for case in cases if case["atoms"] == atoms]
        axes[0].scatter([case["min_z_A"] for case in subset], [case["energy_error_eV"] * 1000 for case in subset],
                        label=f"{atoms}-atom connected molecule", s=55, marker=marker, color=color, zorder=4)
        axes[1].scatter([case["min_z_A"] for case in subset], [max(case["max_force_error_eV_A"], 1e-13) for case in subset],
                        s=55, marker=marker, color=color, zorder=4)
    axes[0].set(title="A   Adsorption energy error", xlabel="Minimum molecular atom height (Å)", ylabel="Fast - native energy (meV)")
    axes[0].axhline(0, color="#777777", lw=.7)
    axes[1].set(title="B   Largest Cartesian force error", xlabel="Minimum molecular atom height (Å)",
                ylabel="Maximum atom-component error (eV / Å)", yscale="log")
    for axis in axes:
        axis.grid(alpha=.2)
    axes[0].legend(frameon=False, fontsize=8, ncol=2)
    fig.text(.075, .035, "Includes nonplanar carbon forces, cell-corner geometry and a connected molecule spanning multiple primitive cells.", fontsize=8.8)
    save(fig, "periodic-potential-accuracy")

    destination = OUT / "optimization-results.zip"
    files = []
    for folder in (OUT, FAST, ROOT / "outputs" / "surface-mode-check"):
        files.extend(path for path in folder.rglob("*")
                     if path.is_file() and path.suffix in (".json", ".csv", ".png", ".svg")
                     and "matplotlib-cache" not in path.parts)
    with zipfile.ZipFile(destination, "w", compression=zipfile.ZIP_DEFLATED) as bundle:
        for path in sorted(set(files)):
            bundle.write(path, path.relative_to(ROOT / "outputs"))
    with zipfile.ZipFile(destination) as bundle:
        if bundle.testzip() is not None:
            raise AssertionError("Result bundle did not pass CRC verification.")
    print(json.dumps({"figures": str(FIGURES), "bundle": str(destination),
                      "performance": [{key: row[key] for key in ("model", "median_us_per_step", "speedup_vs_baseline")} for row in rows]}, indent=2))


if __name__ == "__main__":
    main()
