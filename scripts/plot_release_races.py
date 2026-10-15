"""Plot only measured, accepted finite release trajectories and diagnostics."""
from __future__ import annotations

import hashlib
import json
import os
from pathlib import Path
import sys

for variable in ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS", "NUMEXPR_NUM_THREADS"):
    os.environ[variable] = "1"
ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / "outputs/md-release-validation"
FIGURES = OUT / "races/figures"
sys.path.insert(0, str(ROOT / ".local/analysis-libs"))
os.environ["MPLCONFIGDIR"] = str(OUT / "races/matplotlib-cache")
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
from matplotlib.lines import Line2D

COLORS = ("#0072B2", "#D55E00", "#009E73")
LABELS = {"test-car-52": "Four-arm car", "benzene-12": "Benzene", "fullerene-c60": "C60 cage",
          "ethynyl-rigid": "Ethynyl-linked rings", "folded-eicosane": "Folded eicosane"}


def group(report, fixture, phase):
    rows = sorted((r for r in report["trials"] if r["fixture"] == fixture and r["config"]["phase"] == phase),
                  key=lambda r: r["config"]["seed"])
    if len(rows) != 3:
        raise ValueError(f"Expected three seeds for {fixture}/{phase}")
    return rows


def heading(fig, title, subtitle):
    fig.suptitle(title, fontsize=15, y=.985)
    fig.text(.5, .920, subtitle, ha="center", fontsize=9, color="#555555")


def save(fig, name):
    for extension in ("png", "svg"):
        fig.savefig(FIGURES / f"{name}.{extension}", dpi=180, facecolor="white")
    plt.close(fig)


def seed_legend(fig, seeds, y=.88):
    fig.legend(handles=[Line2D([0], [0], color=color, label=f"Seed {seed}") for seed, color in zip(sorted(seeds), COLORS)],
               loc="upper center", bbox_to_anchor=(.5, y), ncol=3, frameon=False, fontsize=8.5)


def overview(report):
    fixtures = report["fixtures"]
    fig, axes = plt.subplots(1, 2, figsize=(12.5, 5.8))
    fig.subplots_adjust(left=.075, right=.985, top=.80, bottom=.31, wspace=.26)
    heading(fig, "Release validation: five different C/H test structures",
            "Full periodic Au | Fixed-gradient S = 0.75 | 600 K | 350 fs Langevin damping | 0.5 fs timestep | Three seeds per structure")
    seed_legend(fig, report["protocol"]["seeds"])
    limit = report["protocol"]["race_limit_ps"]
    for j, fixture in enumerate(fixtures):
        rows = group(report, fixture["id"], "race")
        for i, row in enumerate(rows):
            x = j + (i - 1) * .10
            complete = row["outcome"] == "finish"
            axes[0].scatter(x, row["observed_time_ps"], color=COLORS[i], s=45,
                            marker="o" if complete else "^", facecolors=COLORS[i] if complete else "none")
            axes[1].scatter(x, row["maximum_lateral_excursion_nm"], color=COLORS[i], s=42)
        axes[0].text(j, limit * 1.09, f"{sum(r['outcome'] == 'finish' for r in rows)}/3 finished", ha="center", fontsize=8)
    ticks = [f"{LABELS[f['id']]}\n{f['formula']} ({f['atoms']} atoms)" for f in fixtures]
    for ax in axes:
        ax.set_xticks(range(len(fixtures)), ticks, rotation=18, ha="right")
        ax.grid(axis="y", alpha=.18)
    axes[0].set(title="A   Observed race completion or time limit", ylabel="Observed driven time (ps)", ylim=(0., limit * 1.22))
    axes[1].set(title="B   Transverse excursion during the observed race", ylabel="Maximum |COM y − starting y| (nm)", ylim=(0., None))
    fig.text(.075, .105, "Filled circle: first sampled finish at absolute x ≥20 nm. Open triangle: still unfinished at the stated time limit; no extrapolated finish time.", fontsize=8.5)
    fig.text(.075, .055, "Generated/imported C/H test structures, not experimental nanocars. Native integration masses, finite energy/forces, bonds and height were audited.", fontsize=8.5)
    save(fig, "01-five-shape-race-overview")


def race_trajectories(report):
    fixtures = report["fixtures"]
    fig, axes = plt.subplots(2, len(fixtures), figsize=(3.25 * len(fixtures), 7.8), squeeze=False)
    fig.subplots_adjust(left=.065, right=.985, top=.78, bottom=.19, wspace=.38, hspace=.38)
    heading(fig, "Actual release-race forward and transverse trajectories",
            "5 ps zero-drive thermal preparation | Fixed S = 0.75 | Same 0.1 ps sampling and random-number setup cadence for every structure")
    seed_legend(fig, report["protocol"]["seeds"])
    limit = report["protocol"]["race_limit_ps"]
    for j, fixture in enumerate(fixtures):
        rows = group(report, fixture["id"], "race")
        for i, row in enumerate(rows):
            time = np.array([s["time_ps"] for s in row["samples"]])
            center = np.array([s["center_nm"] for s in row["samples"]])
            relative = center - center[0]
            axes[0, j].plot(time, center[:, 0], color=COLORS[i], lw=1.1)
            axes[0, j].scatter(time[-1], center[-1, 0], color=COLORS[i], marker="*" if row["outcome"] == "finish" else "^", s=35)
            axes[1, j].plot(relative[:, 0], relative[:, 1], color=COLORS[i], lw=1.0)
            axes[1, j].scatter(relative[-1, 0], relative[-1, 1], color=COLORS[i], marker="*" if row["outcome"] == "finish" else "^", s=30)
        axes[0, j].set(title=f"{LABELS[fixture['id']]}\n{fixture['formula']} | {fixture['atoms']} atoms",
                       xlabel="Driven time (ps)", ylabel="Absolute COM x (nm)", xlim=(0., limit))
        axes[0, j].axhline(20., color="#555555", lw=.7, ls=":")
        axes[1, j].set(xlabel="COM Δx (nm)", ylabel="COM Δy (nm)")
        axes[1, j].axhline(0., color="#777777", lw=.7, ls=":")
        axes[1, j].set_title(f"Finish {sum(r['outcome'] == 'finish' for r in rows)}/3; lift {sum(r['sustained_lift'] for r in rows)}/3", fontsize=8.5)
        for ax in axes[:, j]:
            ax.grid(alpha=.17, lw=.6)
    fig.text(.065, .112, "Top: absolute COM x versus time. Bottom: XY displacement; x/y scales differ to resolve small lateral deviations. Star = finish; triangle = time limit.", fontsize=8.5)
    fig.text(.065, .065, "Every seed remains a separate measured curve. Lift requires every atom above 0.8 nm for ≥5 ps; atomic and COM heights are retained in the report.", fontsize=8.5)
    fig.text(.065, .022, "High atoms in an upright cage or tilted ring do not by themselves imply desorption. No energy-conservation claim applies to Langevin and driven races.", fontsize=8.5)
    save(fig, "01-race-displacement-and-xy")


def no_drive(report):
    fixtures = report["fixtures"]
    fig, axes = plt.subplots(2, len(fixtures), figsize=(3.25 * len(fixtures), 7.8), squeeze=False)
    fig.subplots_adjust(left=.065, right=.985, top=.78, bottom=.19, wspace=.38, hspace=.40)
    duration = report["protocol"]["diffusion_window_ps"]
    heading(fig, "Short no-drive displacement and orientation diagnostics",
            f"S = 0 | 600 K | Full periodic Au | {duration:g} ps after the same 5 ps thermal preparation | Three finite trajectories per structure")
    seed_legend(fig, report["protocol"]["seeds"])
    for j, fixture in enumerate(fixtures):
        rows = group(report, fixture["id"], "diffusion")
        for i, row in enumerate(rows):
            time = np.array([s["time_ps"] for s in row["samples"]])
            center = np.array([s["center_nm"] for s in row["samples"]])
            relative = center - center[0]
            angle = np.array([s["reference_fit_rotation_angle_degrees"] for s in row["samples"]])
            axes[0, j].plot(relative[:, 0], relative[:, 1], color=COLORS[i], lw=1.0)
            axes[0, j].scatter(relative[-1, 0], relative[-1, 1], color=COLORS[i], marker="^", s=28)
            axes[1, j].plot(time, angle, color=COLORS[i], lw=1.0)
        axes[0, j].set(title=LABELS[fixture["id"]], xlabel="COM Δx (nm)", ylabel="COM Δy (nm)")
        axes[0, j].set_aspect("equal", adjustable="datalim")
        axes[1, j].set(xlabel="No-drive time (ps)", ylabel="Reference-fit angle (degrees)", xlim=(0., duration), ylim=(0., 180.))
        for ax in axes[:, j]:
            ax.grid(alpha=.17, lw=.6)
    fig.text(.065, .112, "XY paths use equal x/y scales within each panel. Rotation is a mass-weighted, atom-labelled best rigid fit to the initial structure, bounded 0–180°.", fontsize=8.5)
    fig.text(.065, .065, "For flexible structures, shape changes also influence the fitted angle; fitted RMSD and conditional principal-axis statistics are retained in the report.", fontsize=8.5)
    fig.text(.065, .022, "The finite track has soft walls near its start and edges. These short trajectories are diagnostics, not converged translational or rotational diffusion coefficients.", fontsize=8.5)
    save(fig, "02-no-drive-motion-and-orientation")


def nve_energy():
    path = OUT / "numerics/nve-report.json"
    if not path.exists():
        return None
    report = json.loads(path.read_text(encoding="utf-8"))
    if not report.get("passed") or any(not row["passed"] for row in report["cases"]):
        raise ValueError("Only accepted paired conservative NVE checks can be plotted")
    for name, expected in report["source_sha256"].items():
        if hashlib.sha256((ROOT / name).read_bytes()).hexdigest() != expected:
            raise ValueError(f"Source differs from the conservative NVE measurement: {name}")
    fixtures = list(dict.fromkeys(row["fixture"] for row in report["cases"]))
    dt_colors = {"0.5": "#D55E00", "0.25": "#0072B2"}
    styles = ("-", "--", ":")
    fig, axes = plt.subplots(1, len(fixtures), figsize=(4.25 * len(fixtures), 5.8), squeeze=False)
    fig.subplots_adjust(left=.075, right=.985, top=.79, bottom=.23, wspace=.30)
    heading(fig, "Conservative NVE: energy error decreases when timestep is halved",
            "Same initial positions and velocities for each pair | 600 K initial kinetic condition | Zero drive | Langevin bath disabled")
    handles = [Line2D([0], [0], color=dt_colors[str(dt)], label=f"dt = {dt:g} fs") for dt in (.5, .25)]
    handles += [Line2D([0], [0], color="#555555", ls=style, label=f"Seed {seed}")
                for seed, style in zip(sorted({row["seed"] for row in report["cases"]}), styles)]
    fig.legend(handles=handles, loc="upper center", bbox_to_anchor=(.5, .885), ncol=5, frameon=False, fontsize=8)
    ratios = []
    for j, fixture in enumerate(fixtures):
        rows = sorted((row for row in report["cases"] if row["fixture"] == fixture), key=lambda row: row["seed"])
        ratios.extend(row["fine_to_coarse_energy_drift_ratio"] for row in rows)
        for i, row in enumerate(rows):
            trajectory = Path(row["initial_state_and_trajectories_npz"])
            if hashlib.sha256(trajectory.read_bytes()).hexdigest() != row["trajectory_sha256"]:
                raise ValueError("A paired NVE coordinate file differs from its checksum")
            for method in row["dt_comparisons"]:
                time = np.array([s["time_ps"] for s in method["samples"]])
                energy = np.array([s["total_energy_kcal_mol"] for s in method["samples"]])
                relative_percent = 100. * (energy - energy[0]) / method["initial_ke_kcal_mol"]
                axes[0, j].plot(time, relative_percent, color=dt_colors[str(method["dt_fs"])], ls=styles[i], lw=1.)
        short_label = {"benzene.mol": "Benzene", "tetraphenylbenzene.mol": "Four-arm test car", "folded-eicosane.mol": "Folded eicosane"}.get(fixture, Path(fixture).stem)
        axes[0, j].set(title=f"{short_label} | {rows[0]['atoms']} atoms", xlabel="Conservative simulation time (ps)",
                       ylabel="[E(t) − E(0)] / initial KE (%)")
        axes[0, j].axhline(0., color="#777777", lw=.7)
        axes[0, j].grid(alpha=.18)
    fig.text(.075, .115, f"Fine/coarse maximum |energy error| ratio across the nine pairs: {min(ratios):.3f}–{max(ratios):.3f}. Individual seed errors remain visible.", fontsize=8.5)
    fig.text(.075, .068, "Only the validation process changes dt to 0.25 fs. The game's 0.5 fs timestep and dynamics remain unchanged.", fontsize=8.5)
    fig.text(.075, .027, "Finite conservative numerical checks test force/integrator consistency; they do not calibrate UFF/SAIP or guarantee long-time physical accuracy.", fontsize=8.5)
    save(fig, "03-conservative-energy-timestep-check")
    return {"report_sha256": hashlib.sha256(path.read_bytes()).hexdigest(), "pairs": len(report["cases"]),
            "fine_to_coarse_energy_error_ratio_range": [min(ratios), max(ratios)]}


def main():
    path = OUT / "races/report.json"
    report = json.loads(path.read_text(encoding="utf-8"))
    if report.get("status") != "passed" or len(report["trials"]) != report["planned_trials"]:
        raise ValueError("Only a complete, numerically accepted release report can be plotted")
    for source, expected in report["source_sha256"].items():
        if hashlib.sha256((ROOT / source).read_bytes()).hexdigest() != expected:
            raise ValueError(f"Measured source or input changed: {source}")
    for row in report["trials"]:
        if row["status"] != "passed" or row["source_sha256"] != report["source_sha256"]:
            raise ValueError(f"Failed/inconsistent measurement: {row['key']}")
        if hashlib.sha256((OUT / row["coordinates_file"]).read_bytes()).hexdigest() != row["coordinates_sha256"]:
            raise ValueError(f"Coordinate checksum mismatch: {row['key']}")
    FIGURES.mkdir(parents=True, exist_ok=True)
    plt.rcParams.update({"font.family": "DejaVu Sans", "font.size": 9,
                         "axes.spines.top": False, "axes.spines.right": False,
                         "svg.fonttype": "none", "axes.titlesize": 10})
    race_trajectories(report)
    if report["protocol"]["diffusion_window_ps"]:
        no_drive(report)
    nve = nve_energy()
    metadata = {"report_sha256": hashlib.sha256(path.read_bytes()).hexdigest(),
                "figures": sorted(p.name for p in FIGURES.glob("*.png")),
                "sustained_lift_trials": sum(row["sustained_lift"] for row in report["trials"]),
                "failed_trials": 0, "qualification": report["protocol"]["qualification"], "paired_nve": nve}
    (OUT / "races/plot-metadata.json").write_text(json.dumps(metadata, indent=2) + "\n", encoding="utf-8")
    print(json.dumps(metadata))


if __name__ == "__main__":
    main()
