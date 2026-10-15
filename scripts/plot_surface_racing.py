"""Publication-style figures from real LAMMPS surface-racing trials.

Matplotlib is an optional, isolated analysis dependency. No synthetic or
extrapolated motion is drawn. Main race figures require all 18 matched trials.
"""

from __future__ import annotations

import csv
import json
import os
from pathlib import Path
import sys
import zipfile

for variable in ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS", "NUMEXPR_NUM_THREADS"):
    os.environ[variable] = "1"
ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / "outputs/surface-racing-check"
FIGURES = OUT / "figures"
sys.path.insert(0, str(ROOT / ".local/analysis-libs"))
os.environ["MPLCONFIGDIR"] = str(OUT / "matplotlib-cache")

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
from matplotlib.lines import Line2D

MODELS = ("continuous", "frozen_au", "mobile_au")
STRENGTHS = (.75, .25)
COLORS = {"continuous": "#0072B2", "frozen_au": "#D55E00", "mobile_au": "#009E73"}
LABELS = {"continuous": "Continuum", "frozen_au": "Frozen Au", "mobile_au": "Mobile Au"}
SEED_STYLES = ("-", "--", ":")


def save_figure(fig, name):
    fig.savefig(FIGURES / f"{name}.png", dpi=180, facecolor="white")
    fig.savefig(FIGURES / f"{name}.svg", facecolor="white")
    plt.close(fig)


def figure_heading(fig, title, subtitle, *, subtitle_y=.932):
    fig.suptitle(title, fontsize=15, y=.98)
    fig.text(.5, subtitle_y, subtitle, ha="center", fontsize=9, color="#535353")


def model_legend(fig):
    handles = [Line2D([0], [0], color=COLORS[m], lw=2, label=LABELS[m]) for m in MODELS]
    fig.legend(handles=handles, loc="upper center", bbox_to_anchor=(.5, .902), ncol=3,
               frameon=False, handlelength=2.7)


def style_axes(ax, panel, title, ylabel, *, zero=False):
    ax.text(-.12, 1.035, panel, transform=ax.transAxes, fontweight="bold", fontsize=12)
    ax.set(title=title, xlabel="Simulation time (ps)", ylabel=ylabel)
    ax.grid(alpha=.18, lw=.6)
    if zero:
        ax.axhline(0, color="#777777", lw=.7, zorder=0)


def mean_sd(ax, times, values, model, individual=False):
    color = COLORS[model]
    if individual:
        for data in values:
            ax.plot(times, data, color=color, lw=.7, alpha=.22)
    mean = np.mean(values, axis=0)
    sd = np.std(values, axis=0, ddof=1)
    ax.fill_between(times, mean - sd, mean + sd, color=color, alpha=.14, lw=0)
    ax.plot(times, mean, color=color, lw=1.8)


report = json.loads((OUT / "report.json").read_text(encoding="utf-8"))
trials = [t for t in report["trials"] if t["configuration"]["stage"] == "races"]
if len(trials) != 18 or any(t["status"] != "passed" for t in trials):
    raise ValueError("Complete all 18 surface-racing trials before drawing the main comparison.")
if len({t["duration_ps"] for t in trials}) != 1:
    raise ValueError("Race comparisons require the same physical window for all three models.")
duration = trials[0]["duration_ps"]
groups = {}
for strength in STRENGTHS:
    for model in MODELS:
        group = sorted((t for t in trials if t["surface"] == model and np.isclose(t["strength"], strength)),
                       key=lambda t: t["seed"])
        if len(group) != 3:
            raise ValueError(f"Expected three seeds for {model}, S={strength}.")
        data = [np.genfromtxt(t["trajectory_csv"], delimiter=",", names=True, encoding="utf-8") for t in group]
        times = data[0]["time_ps"]
        if any(not np.array_equal(times, row["time_ps"]) for row in data):
            raise ValueError("Sampling clocks differ across seeds.")
        groups[(strength, model)] = {"trials": group, "data": data, "times": times}

FIGURES.mkdir(parents=True, exist_ok=True)
plt.rcParams.update({"font.family": "DejaVu Sans", "font.size": 10,
    "axes.spines.top": False, "axes.spines.right": False, "axes.titlesize": 11,
    "svg.fonttype": "none", "savefig.bbox": "tight"})

# Motion: distinguish displacement, unsigned transverse excursion and velocity.
fig, axes = plt.subplots(2, 3, figsize=(13.6, 8.0), sharex=True)
figure_heading(fig, "Surface model changes in early molecular racing",
               f"C6H6 | 600 K | dt = 0.5 fs | damping = 350 fs | {duration:g} ps | 3 seeds per condition")
model_legend(fig)
for row_index, strength in enumerate(STRENGTHS):
    for model in MODELS:
        group = groups[(strength, model)]
        data, times = group["data"], group["times"]
        mean_sd(axes[row_index, 0], times, np.array([d["dx_A"] / 10 for d in data]), model, individual=True)
        dy = np.array([d["dy_A"] / 10 for d in data])
        for y in dy:
            axes[row_index, 1].plot(times, np.abs(y), color=COLORS[model], lw=.65, alpha=.2)
        axes[row_index, 1].plot(times, np.sqrt(np.mean(dy**2, axis=0)), color=COLORS[model], lw=1.8)
        interval = float(times[1] - times[0])
        width = int(round(2.0 / interval))
        width += (width + 1) % 2  # odd, centered approximately 2 ps window
        velocity = np.array([np.convolve(d["vx_A_ps"] / 10, np.ones(width) / width, mode="valid") for d in data])
        half = width // 2
        mean_sd(axes[row_index, 2], times[half:-half], velocity, model)
    style_axes(axes[row_index, 0], "ad"[row_index], f"Forward displacement | S = {strength:g}", "Delta x (nm)", zero=True)
    style_axes(axes[row_index, 1], "be"[row_index], f"Lateral wandering | S = {strength:g}", "RMS Delta y (nm)")
    style_axes(axes[row_index, 2], "cf"[row_index], f"Forward velocity | S = {strength:g}", "v_x (nm/ps)", zero=True)
    for ax in axes[row_index]:
        ax.set_xlim(0, duration)
fig.text(.065, .035, "Shading: mean +/- sample SD (n = 3), not a confidence interval. Thin lines: individual seeds.", fontsize=9)
fig.text(.065, .012, "Velocity uses a centered 2 ps moving average; lateral RMS is a finite-window statistic, not a diffusion coefficient.", fontsize=9)
fig.subplots_adjust(left=.07, right=.985, top=.825, bottom=.13, wspace=.31, hspace=.48)
save_figure(fig, "race-dynamics")

# Paired contrasts make close trajectories visible without subtracting ensemble
# averages before pairing. Seed correspondence matches input velocities only.
fig, axes = plt.subplots(2, 2, figsize=(11.2, 8.0), sharex=True)
figure_heading(fig, "Displacement differences relative to the continuum model",
               f"Matched initial velocities by seed | 600 K | {duration:g} ps | 3 seed pairs")
handles = [Line2D([0], [0], color=COLORS[m], lw=2, label=LABELS[m]) for m in MODELS[1:]]
handles.append(Line2D([0], [0], color="#777777", lw=.8, label="Continuum reference = 0"))
fig.legend(handles=handles, loc="upper center", bbox_to_anchor=(.5, .902), ncol=3, frameon=False)
for row_index, strength in enumerate(STRENGTHS):
    reference = groups[(strength, "continuous")]
    for model in MODELS[1:]:
        group = groups[(strength, model)]
        if [t["seed"] for t in reference["trials"]] != [t["seed"] for t in group["trials"]]:
            raise ValueError("Seed pairing differs between models.")
        if not np.array_equal(reference["times"], group["times"]):
            raise ValueError("Paired model sampling clocks differ.")
        for col, field in enumerate(("dx_A", "dy_A")):
            differences = np.array([(d[field] - r[field]) / 10
                                    for d, r in zip(group["data"], reference["data"])])
            mean_sd(axes[row_index, col], group["times"], differences, model, individual=True)
    for col, component in enumerate(("x", "y")):
        style_axes(axes[row_index, col], "abcd"[row_index*2+col],
                   f"{'Forward' if col == 0 else 'Lateral'} contrast | S = {strength:g}",
                   f"Delta {component} difference (nm)", zero=True)
        axes[row_index, col].set_xlim(0, duration)
fig.text(.075, .035, "Shading: paired difference mean +/- sample SD (n = 3). Thin lines: individual seed pairs.", fontsize=9)
fig.text(.075, .012, "Matching inputs does not imply identical stepwise Langevin noise; these are finite-window model contrasts.", fontsize=9)
fig.subplots_adjust(left=.1, right=.985, top=.825, bottom=.13, wspace=.32, hspace=.48)
save_figure(fig, "race-model-differences")

# True plan-view paths with equal scales; do not visually amplify transverse motion.
fig, axes = plt.subplots(2, 3, figsize=(13.6, 8.0))
figure_heading(fig, "Unwrapped molecular center-of-mass trajectories",
               f"Equal x/y scale within each panel | 600 K | {duration:g} ps | 3 independent seeds")
handles = [Line2D([0], [0], color="#333333", linestyle=s, lw=1.3, label=f"Seed {i+1}")
           for i, s in enumerate(SEED_STYLES)]
handles.extend([Line2D([0], [0], marker="*", linestyle="none", color="#333333", label="Start"),
                Line2D([0], [0], marker="o", linestyle="none", color="#333333", label="End")])
fig.legend(handles=handles, loc="upper center", bbox_to_anchor=(.5, .9), ncol=5, frameon=False)
for row_index, strength in enumerate(STRENGTHS):
    all_dx = np.concatenate([d["dx_A"] / 10 for m in MODELS for d in groups[(strength, m)]["data"]])
    all_dy = np.concatenate([d["dy_A"] / 10 for m in MODELS for d in groups[(strength, m)]["data"]])
    xspan = [min(0, all_dx.min()) - .25, max(0, all_dx.max()) + .25]
    yspan = [min(0, all_dy.min()) - .25, max(0, all_dy.max()) + .25]
    for col, model in enumerate(MODELS):
        ax = axes[row_index, col]
        for index, data in enumerate(groups[(strength, model)]["data"]):
            x, y = data["dx_A"] / 10, data["dy_A"] / 10
            ax.plot(x, y, color=COLORS[model], linestyle=SEED_STYLES[index], lw=1.2, alpha=.8)
            ax.scatter(x[-1], y[-1], color=COLORS[model], s=20, zorder=4)
        ax.scatter([0], [0], marker="*", s=70, color="black", zorder=5)
        ax.set(xlabel="Delta x (nm)", ylabel="Delta y (nm)", xlim=xspan, ylim=yspan,
               aspect="equal", title=f"{LABELS[model]} | S = {strength:g}")
        ax.grid(alpha=.18, lw=.6)
        ax.text(-.12, 1.035, "abcdef"[row_index*3+col], transform=ax.transAxes, fontweight="bold", fontsize=12)
fig.text(.065, .025, "Small Au surface repeats laterally. Image flags retain physical displacement; no reset of the drive at cell crossings.", fontsize=9)
fig.subplots_adjust(left=.07, right=.985, top=.825, bottom=.11, wspace=.28, hspace=.36)
save_figure(fig, "race-trajectories")

# Contact and conformation, required when dynamic benzene tilts away from fit reference.
fig, axes = plt.subplots(2, 3, figsize=(13.6, 8.0), sharex=True)
figure_heading(fig, "Molecular contact and conformation diagnostics",
               "Reference: continuum z = 0 / instantaneous mean top-layer Au z | UFF fixed topology")
model_legend(fig)
for row_index, strength in enumerate(STRENGTHS):
    for model in MODELS:
        group = groups[(strength, model)]
        for col, field in enumerate(("COM_height_above_top_Au_A", "carbon_tilt_deg", "carbon_plane_rms_A")):
            mean_sd(axes[row_index, col], group["times"], np.array([d[field] for d in group["data"]]), model, individual=True)
    for col, (title, ylabel) in enumerate((("COM height above substrate", "Height (angstrom)"), ("Carbon-ring tilt", "Tilt (degrees)"),
                                         ("Carbon-ring planarity", "Plane RMS (angstrom)"))):
        style_axes(axes[row_index, col], "abcdef"[row_index*3+col], f"{title} | S = {strength:g}", ylabel)
        axes[row_index, col].set_xlim(0, duration)
fig.text(.065, .025, "Mean +/- sample SD; thin lines: individual seeds. Planarity: best-fit C plane. Tilt: nominal Au normal (global z).", fontsize=9)
fig.subplots_adjust(left=.07, right=.985, top=.825, bottom=.12, wspace=.31, hspace=.48)
save_figure(fig, "surface-contact")

# Performance: identical simulated windows and callback monitoring; setup is separate.
fig, axes = plt.subplots(1, 3, figsize=(13.6, 5.4))
figure_heading(fig, "Wall-clock performance of the three surface implementations",
               "Intel Core i9-13900HX | one LAMMPS / BLAS thread | 6 runs/model (2 strengths x 3 seeds)",
               subtitle_y=.9)
metrics = (("integration_seconds_per_step", 1e6, "MD propagation cost", "Wall time per step (microseconds)"),
           ("ps_per_wall_second", 1, "MD throughput", "Simulated ps / wall second"),
           ("preparation_wall_seconds", 1, "Preparation cost", "Wall time (s)"))
performance_rows = []
performance_maxima = np.zeros(len(metrics))
for model in MODELS:
    model_trials = [t for t in trials if t["surface"] == model]
    for col, (key, factor, title, ylabel) in enumerate(metrics):
        values = np.array([(t[key] if key == "preparation_wall_seconds" else t["performance"][key]) * factor for t in model_trials])
        performance_maxima[col] = max(performance_maxima[col], values.max(), values.mean() + values.std(ddof=1))
        index = MODELS.index(model)
        axes[col].bar(index, values.mean(), color=COLORS[model], width=.58, alpha=.82)
        axes[col].errorbar(index, values.mean(), yerr=values.std(ddof=1), color="#202020", capsize=4, fmt="none", lw=1)
        axes[col].scatter(index + np.linspace(-.13, .13, len(values)), values, color="#202020", s=17, zorder=4)
        axes[col].text(index, values.max() * 1.12, f"{values.mean():.2f}", ha="center", va="bottom", fontsize=10)
    wall = np.array([t["performance"]["integration_wall_seconds"] for t in model_trials])
    cost = np.array([t["performance"]["integration_seconds_per_step"] * 1e6 for t in model_trials])
    throughput = np.array([t["performance"]["ps_per_wall_second"] for t in model_trials])
    performance_rows.append({"model": model, "n_runs": len(model_trials),
        "propagation_wall_seconds_mean": float(wall.mean()), "propagation_wall_seconds_SD": float(wall.std(ddof=1)),
        "microseconds_per_step_mean": float(cost.mean()), "microseconds_per_step_SD": float(cost.std(ddof=1)),
        "ps_per_wall_second_mean": float(throughput.mean()), "ps_per_wall_second_SD": float(throughput.std(ddof=1)),
        "total_atoms": model_trials[0]["total_atoms"], "integrated_atoms": model_trials[0]["integrated_atoms"]})
for col, (_, _, title, ylabel) in enumerate(metrics):
    axes[col].set(title=title, ylabel=ylabel, xticks=range(3), xticklabels=[LABELS[m] for m in MODELS], ylim=(0, performance_maxima[col] * 1.32))
    axes[col].grid(axis="y", alpha=.18, lw=.6)
    axes[col].set_axisbelow(True)
    axes[col].text(-.12, 1.035, "abc"[col], transform=axes[col].transAxes, fontweight="bold", fontsize=12)
fig.text(.065, .09, "Bars: mean +/- sample SD of six runs; points: individual runs. Production timer includes UFF callbacks and neighbour work.", fontsize=9)
fig.text(.065, .048, "MD timer excludes sampling and CSV I/O. Preparation includes input loading and mobile-Au thermalization; cached slab relaxation excluded.", fontsize=9)
fig.text(.065, .007, "Atoms (total / integrated): continuum 12 / 12; frozen Au 2892 / 12; mobile Au 2892 / 1932.", fontsize=9)
fig.subplots_adjust(left=.075, right=.985, top=.805, bottom=.23, wspace=.32)
save_figure(fig, "compute-performance")

# Machine-readable figure summaries are derived from the raw trials, not from graph pixels.
summary_rows = []
for strength in STRENGTHS:
    for model in MODELS:
        group = groups[(strength, model)]
        dx = np.array([d["dx_A"][-1] / 10 for d in group["data"]])
        dy = np.array([d["dy_A"][-1] / 10 for d in group["data"]])
        baseline_dx = np.array([d["dx_A"][-1] / 10 for d in groups[(strength, "continuous")]["data"]])
        paired_dx = dx - baseline_dx
        speed = dx / duration
        summary_rows.append({"model": model, "strength": strength, "temperature_K": 600, "duration_ps": duration,
            "n_seeds": 3, "forward_nm_mean": float(dx.mean()), "forward_nm_SD": float(dx.std(ddof=1)),
            "lateral_nm_RMS": float(np.sqrt(np.mean(dy**2))), "forward_speed_nm_ps_mean": float(speed.mean()),
            "forward_speed_nm_ps_SD": float(speed.std(ddof=1)),
            "paired_forward_difference_nm_mean": float(paired_dx.mean()),
            "paired_forward_difference_nm_SD": float(paired_dx.std(ddof=1)),
            "sampled_maximum_tilt_deg": max(float(d["carbon_tilt_deg"].max()) for d in group["data"]),
            "sampled_maximum_carbon_plane_RMS_A": max(float(d["carbon_plane_rms_A"].max()) for d in group["data"]),
            "minimum_bond_A": min(t["production_monitor"]["bond_min_angstrom"] for t in group["trials"]),
            "maximum_bond_A": max(t["production_monitor"]["bond_max_angstrom"] for t in group["trials"])})
baseline = performance_rows[0]["microseconds_per_step_mean"]
for row in performance_rows:
    row["cost_relative_to_continuum"] = row["microseconds_per_step_mean"] / baseline
for name, rows in (("motion-summary.csv", summary_rows), ("performance-summary.csv", performance_rows)):
    with (OUT / name).open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)
(OUT / "figure-summary.json").write_text(json.dumps({"motion": summary_rows, "performance": performance_rows,
    "source_report": str(OUT / "report.json"), "statistics": "Sample SD; no significance test or asymptotic diffusion fit."},
    indent=2, allow_nan=False) + "\n", encoding="utf-8")
archive_files = [OUT / name for name in ("report.json", "figure-summary.json", "motion-summary.csv", "performance-summary.csv")]
archive_files.extend(sorted(FIGURES.glob("*.png")))
archive_files.extend(sorted(FIGURES.glob("*.svg")))
for trial in report["trials"]:
    directory = OUT / "trials" / trial["trial_id"]
    archive_files.extend(directory / name for name in ("trajectory.csv", "summary.json", "configuration.json")
                         if (directory / name).exists())
with zipfile.ZipFile(OUT / "comparison-results.zip", "w", compression=zipfile.ZIP_DEFLATED) as archive:
    for path in archive_files:
        archive.write(path, path.relative_to(OUT))
print(FIGURES)
