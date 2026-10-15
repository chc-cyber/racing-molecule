"""Figures from completed, actual-game 20 nm weak-drive scans.

This plotting process alone uses optional isolated analysis packages. It does
not import the simulation or alter its MD dependencies. No synthetic motion,
extrapolated finish times or confidence intervals are plotted.
"""
from __future__ import annotations

import argparse
import hashlib
import json
import os
from pathlib import Path
import sys

for variable in ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS", "NUMEXPR_NUM_THREADS"):
    os.environ[variable] = "1"
ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / "outputs/weak-drive-20nm"
sys.path.insert(0, str(ROOT / ".local/analysis-libs"))
os.environ["MPLCONFIGDIR"] = str(OUT / "matplotlib-cache")

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.colors import TwoSlopeNorm
from matplotlib.lines import Line2D
import numpy as np

SURFACES = ("periodic-au", "continuum")
SURFACE_LABELS = {"periodic-au": "Full periodic Au", "continuum": "Corrugated continuum"}
SURFACE_COLORS = {"periodic-au": "#D55E00", "continuum": "#0072B2"}
SEED_COLORS = ("#0072B2", "#D55E00", "#009E73")
SEED_STYLES = ("-", "--", ":")
TEMPERATURES = (150, 300, 600)
STRENGTHS = (.05, .10, .15, .25)
FIGURES = OUT / "figures"


def load_report(path):
    report = json.loads(path.read_text(encoding="utf-8"))
    if report.get("status") != "passed":
        raise ValueError(f"Only a complete passed report can be plotted: {path}")
    rows = report["trials"]
    if len(rows) != report["planned_trials"] or len({r["key"] for r in rows}) != len(rows):
        raise ValueError("The report is incomplete or contains duplicate trajectories")
    for name, expected in report["source_sha256"].items():
        if hashlib.sha256((ROOT / name).read_bytes()).hexdigest() != expected:
            raise ValueError(f"Source differs from the measured version: {name}")
    for row in rows:
        if row["status"] != "passed" or row["source_sha256"] != report["source_sha256"]:
            raise ValueError(f"Failed or inconsistent trajectory: {row['key']}")
        coordinates = OUT / row["coordinates_file"]
        if hashlib.sha256(coordinates.read_bytes()).hexdigest() != row["coordinates_sha256"]:
            raise ValueError(f"Coordinates no longer match their recorded checksum: {coordinates}")
    return report


def trials_for(report, surface, temperature, strength, corrugation=.08):
    rows = sorted((r for r in report["trials"]
                   if r["config"]["surface"] == surface
                   and np.isclose(r["config"]["temperature_K"], temperature)
                   and np.isclose(r["config"]["strength"], strength)
                   and np.isclose(r["config"]["corrugation"], corrugation)),
                  key=lambda r: r["config"]["seed"])
    if len(rows) != 3:
        raise ValueError(f"Expected three seeds: {surface}, {temperature} K, S={strength}")
    return rows


def series(row):
    time = np.array([s["time_ps"] for s in row["samples"]])
    center = np.array([s["center_nm"] for s in row["samples"]])
    return time, center, center - center[0]


def lift_diagnostic(row):
    """Explicit geometric warning, not a thermodynamic desorption criterion.

    A sustained lift flag requires *every atom* to be above 0.8 nm
    for at least 5 ps in the sampled trajectory. Rotation of one phenyl ring
    or a large single-atom z excursion therefore does not trigger this flag.
    """
    with np.load(OUT / row["coordinates_file"], allow_pickle=False) as data:
        time = data["time_ps"]
        positions = data["positions_angstrom"]
        minimum_atom_height_nm = positions[:, :, 2].min(axis=1) / 10.
        all_atoms_above = minimum_atom_height_nm > .8
        first = None
        first_sustained = None
        longest = 0.
        for index, above in enumerate(all_atoms_above):
            if above:
                if first is None:
                    first = index
                duration = float(time[index] - time[first])
                longest = max(longest, duration)
                if duration >= 5. - 1e-8 and first_sustained is None:
                    first_sustained = float(time[first])
            else:
                first = None
    _, centers, _ = series(row)
    return {"sustained_lift": first_sustained is not None,
            "first_lift_ps": first_sustained, "longest_all_atoms_above_0_8_nm_ps": longest,
            "lowest_atom_height_range_nm": [float(minimum_atom_height_nm.min()), float(minimum_atom_height_nm.max())],
            "com_height_range_nm": [float(centers[:, 2].min()), float(centers[:, 2].max())]}


def heading(fig, title, subtitle):
    fig.suptitle(title, fontsize=15, y=.982)
    fig.text(.5, .938, subtitle, ha="center", fontsize=9, color="#525252")


def save(fig, name):
    for extension in ("png", "svg"):
        fig.savefig(FIGURES / f"{name}.{extension}", dpi=190, facecolor="white")
    plt.close(fig)


def group_stats(rows):
    dx = np.array([r["metrics"]["displacement_nm"][0] for r in rows])
    dy = np.array([r["metrics"]["displacement_nm"][1] for r in rows])
    return (float(dx.mean()), float(np.sqrt(np.mean(dy**2))),
            float(np.mean([r["metrics"]["backward_fraction"] for r in rows])) * 100,
            float(np.mean([r["metrics"]["pause_fraction"] for r in rows])) * 100)


def parameter_matrix(report, diagnostics):
    values = np.empty((4, len(TEMPERATURES), len(STRENGTHS)))
    lifted = np.zeros((len(TEMPERATURES), len(STRENGTHS)), dtype=int)
    for i, temperature in enumerate(TEMPERATURES):
        for j, strength in enumerate(STRENGTHS):
            rows = trials_for(report, "periodic-au", temperature, strength)
            values[:, i, j] = group_stats(rows)
            lifted[i, j] = sum(diagnostics[r["key"]]["sustained_lift"] for r in rows)
    fig, axes = plt.subplots(2, 2, figsize=(11.2, 8.2))
    fig.subplots_adjust(left=.09, right=.95, top=.855, bottom=.18, hspace=.44, wspace=.32)
    heading(fig, "Weak forward bias on the complete existing Au lattice",
            "20 nm × 20 nm track | C30H22, 52 atoms | Up to 50 ps after 5 ps zero-drive equilibration | Three seeds per cell")
    titles = ("A   Mean end-point forward displacement (nm)",
              "B   End-point transverse RMS across seeds (nm)",
              "C   Backward longitudinal windows (%)",
              "D   Near-zero longitudinal progress windows (%)")
    for k, ax in enumerate(axes.flat):
        if k == 0:
            magnitude = max(float(np.abs(values[k]).max()), 1e-4)
            norm = TwoSlopeNorm(vmin=-magnitude, vcenter=0., vmax=magnitude)
            cmap = "RdBu"
        else:
            norm = matplotlib.colors.Normalize(0., 100. if k >= 2 else max(float(values[k].max()), 1e-4))
            cmap = "cividis"
        image = ax.imshow(values[k], cmap=cmap, norm=norm, aspect="auto")
        ax.set(title=titles[k], xlabel="Fixed-gradient strength S", ylabel="Bath temperature (K)",
               xticks=np.arange(len(STRENGTHS)), xticklabels=[f"{s:.2f}" for s in STRENGTHS],
               yticks=np.arange(len(TEMPERATURES)), yticklabels=TEMPERATURES)
        ax.tick_params(length=0)
        for i in range(len(TEMPERATURES)):
            for j in range(len(STRENGTHS)):
                color = "white" if np.mean(matplotlib.colormaps[cmap](norm(values[k, i, j]))[:3]) < .50 else "#111111"
                label = f"{values[k, i, j]:.2f}" if k < 2 else f"{values[k, i, j]:.1f}"
                if lifted[i, j]:
                    label += " †"
                ax.text(j, i, label, ha="center", va="center", color=color, fontsize=10)
        fig.colorbar(image, ax=ax, fraction=.043, pad=.035)
    fig.text(.09, .105, "Backward: Δx < −0.02 nm per 1 ps. Near-zero: |Δx| ≤ 0.02 nm per 1 ps; sideways motion is allowed.", fontsize=8.5)
    fig.text(.09, .066, "End-point RMS = √mean[(yend − ystart)²], not a diffusion coefficient. No confidence intervals or finish-time extrapolations.", fontsize=8.5)
    if lifted.any():
        fig.text(.09, .027, "† At least one seed has sustained lift: all atoms above 0.8 nm for ≥5 ps. Flagged data retained; see height diagnostics.", fontsize=8)
    else:
        fig.text(.09, .027, "No sustained-lift flags: no trial had all atoms above 0.8 nm for ≥5 ps. This is a geometric diagnostic.", fontsize=8)
    save(fig, "01-periodic-au-parameter-map")


def surface_comparison(report, diagnostics):
    fig, axes = plt.subplots(2, 3, figsize=(13.2, 8.0), sharex=True)
    fig.subplots_adjust(left=.075, right=.985, top=.80, bottom=.19, hspace=.34, wspace=.30)
    heading(fig, "Full periodic Au versus the existing corrugated continuum",
            "Matched temperature, drive, timestep and seeds | Complete molecular MD | Up to 50 ps | Three seeds per condition")
    handles = [Line2D([0], [0], color=SURFACE_COLORS[s], marker="o", label=SURFACE_LABELS[s]) for s in SURFACES]
    fig.legend(handles=handles, loc="upper center", bbox_to_anchor=(.5, .89), ncol=2, frameon=False)
    x = np.array(STRENGTHS)
    for j, temperature in enumerate(TEMPERATURES):
        for surface, offset in zip(SURFACES, (-.0025, .0025)):
            dxs, lateral, flags = [], [], []
            for strength in STRENGTHS:
                rows = trials_for(report, surface, temperature, strength)
                dxs.append([r["metrics"]["displacement_nm"][0] for r in rows])
                lateral.append([abs(r["metrics"]["displacement_nm"][1]) for r in rows])
                flags.append(any(diagnostics[r["key"]]["sustained_lift"] for r in rows))
            for k, individual in enumerate((np.array(dxs), np.array(lateral))):
                aggregate = individual.mean(axis=1) if k == 0 else np.sqrt(np.mean(individual**2, axis=1))
                axes[k, j].errorbar(x + offset, aggregate,
                                   yerr=[aggregate - individual.min(axis=1), individual.max(axis=1) - aggregate],
                                   fmt="o-", color=SURFACE_COLORS[surface], ms=4, lw=1.5, capsize=3)
                for seed_index in range(3):
                    axes[k, j].scatter(x + offset + (seed_index - 1) * .001, individual[:, seed_index],
                                       color=SURFACE_COLORS[surface], s=15, alpha=.45, zorder=3)
                for index, flag in enumerate(flags):
                    if flag:
                        axes[k, j].annotate("†", (x[index] + offset, individual[index].max()),
                                            xytext=(0, 5), textcoords="offset points", ha="center", color=SURFACE_COLORS[surface])
        axes[0, j].set_title(f"{chr(65 + j)}   {temperature} K")
        axes[1, j].set_title(f"{chr(68 + j)}   {temperature} K")
        axes[0, j].set_ylabel("End-point forward displacement (nm)")
        axes[1, j].set_ylabel("End-point transverse displacement (nm)")
        axes[1, j].set_xlabel("Fixed-gradient strength S")
        axes[0, j].axhline(0., lw=.7, color="#777777")
        for k in (0, 1):
            axes[k, j].set_xticks(x, [f"{s:.2f}" for s in STRENGTHS])
            axes[k, j].set_xlim(.03, .27)
            axes[k, j].grid(alpha=.18, lw=.6)
    fig.text(.075, .110, "Top line: mean Δx. Bottom line: RMS Δy. Dots: individual Δx or |Δy|. Bars: seed minimum–maximum, not confidence intervals.", fontsize=8.5)
    fig.text(.075, .069, "Periodic Au retains its complete existing frozen-lattice SAIP potential. Continuum retains its separate A = 0.08 kcal/mol corrugation.", fontsize=8.5)
    fig.text(.075, .029, "Strong- and zero-drive reference trajectories appear in the next figure. † marks sustained lift when present; flagged values are retained.", fontsize=8.5)
    save(fig, "02-surface-model-comparison")


def trajectory_panels(report, conditions, diagnostics, *, name, compare_surfaces=False):
    fig, axes = plt.subplots(2, len(conditions), figsize=(max(10., 3.3 * len(conditions)), 8.0), squeeze=False)
    fig.subplots_adjust(left=.075, right=.985, top=.795, bottom=.19, hspace=.39, wspace=.35)
    horizon = max(r["config"]["duration_ps"] for r in report["trials"])
    heading(fig, "Actual molecular trajectories: wandering, pauses and forward progress",
            f"{'Matched surface comparison' if compare_surfaces else 'Full existing periodic Au'} | Up to {horizon:g} ps | Sampled every 0.1 ps | Three separate seed trajectories")
    seed_ids = sorted(report["protocol"]["seeds"])
    if compare_surfaces:
        handles = [Line2D([0], [0], color=SURFACE_COLORS[s], lw=2, label=SURFACE_LABELS[s]) for s in SURFACES]
        handles += [Line2D([0], [0], color="#555555", ls=ls, label=f"Seed {seed}") for seed, ls in zip(seed_ids, SEED_STYLES)]
    else:
        handles = [Line2D([0], [0], color=color, lw=1.8, label=f"Seed {seed}") for seed, color in zip(seed_ids, SEED_COLORS)]
    fig.legend(handles=handles, loc="upper center", bbox_to_anchor=(.5, .891), ncol=len(handles), frameon=False, fontsize=8)
    for j, condition in enumerate(conditions):
        temperature, strength = condition["temperature_K"], condition["strength"]
        surfaces = SURFACES if compare_surfaces else ("periodic-au",)
        status = []
        for surface in surfaces:
            rows = trials_for(report, surface, temperature, strength, condition.get("corrugation", .08))
            flags = sum(diagnostics[r["key"]]["sustained_lift"] for r in rows)
            finishers = sum(r["metrics"]["finish_observed"] for r in rows)
            status.append(f"{'Au' if surface == 'periodic-au' else 'Continuum'}: finish {finishers}/3; lift {flags}/3")
            for i, row in enumerate(rows):
                time, _, relative = series(row)
                color = SURFACE_COLORS[surface] if compare_surfaces else SEED_COLORS[i]
                style = SEED_STYLES[i] if compare_surfaces else "-"
                axes[0, j].plot(relative[:, 0], relative[:, 1], color=color, ls=style, lw=1.0, alpha=.86)
                axes[0, j].scatter(relative[0, 0], relative[0, 1], color=color, marker="o", s=20)
                axes[0, j].scatter(relative[-1, 0], relative[-1, 1], color=color,
                                   marker="*" if row["metrics"]["finish_observed"] else "^", s=38)
                axes[1, j].plot(time, relative[:, 0], color=color, ls=style, lw=1.1, alpha=.86)
                axes[1, j].scatter(time[-1], relative[-1, 0], color=color,
                                   marker="*" if row["metrics"]["finish_observed"] else "^", s=28)
                first_lift = diagnostics[row["key"]]["first_lift_ps"]
                if first_lift is not None:
                    position = np.searchsorted(time, first_lift)
                    axes[0, j].scatter(relative[position, 0], relative[position, 1], marker="x", color="#111111", s=40, zorder=5)
                    axes[1, j].scatter(time[position], relative[position, 0], marker="x", color="#111111", s=40, zorder=5)
        axes[0, j].set(title=f"{temperature:g} K | S = {strength:g}", xlabel="COM Δx (nm)", ylabel="COM Δy (nm)")
        axes[0, j].set_aspect("equal", adjustable="datalim")
        axes[0, j].axhline(0., lw=.6, color="#999999", zorder=0)
        axes[0, j].axvline(0., lw=.6, color="#999999", zorder=0)
        axes[1, j].set(xlabel="Driven simulation time (ps)", ylabel="COM Δx (nm)", xlim=(0., horizon))
        axes[1, j].axhline(0., lw=.6, color="#999999", zorder=0)
        axes[1, j].set_title("\n".join(status), fontsize=8., pad=5.)
        for ax in axes[:, j]:
            ax.grid(alpha=.17, lw=.6)
            ax.ticklabel_format(axis="both", style="plain", useOffset=False)
    fig.text(.075, .111, "Each XY panel uses equal x/y scales within that panel; numeric axis ranges differ between panels. Start: circle. End: triangle; finish: star.", fontsize=8.5)
    fig.text(.075, .069, "Curves stop at the first sampled COM crossing of absolute x = 20 nm. No crossing within the window is an observed non-finish, not an inferred time.", fontsize=8.5)
    fig.text(.075, .028, "Black ×: first sustained-lift flag (all atoms above 0.8 nm for ≥5 ps). Static bias sets drift; Langevin supplies random forces.", fontsize=8.5)
    save(fig, name)


def height_figure(reports, diagnostics):
    flagged = [r for report in reports for r in report["trials"] if diagnostics[r["key"]]["sustained_lift"]]
    if not flagged:
        return
    fig, axes = plt.subplots(1, 2, figsize=(12., 5.6))
    fig.subplots_adjust(left=.08, right=.98, top=.81, bottom=.24, wspace=.30)
    heading(fig, "Sustained lift must be distinguished from surface wandering",
            "Geometric diagnostic: every atom above 0.8 nm for at least 5 ps | Only flagged actual trajectories shown")
    for row in flagged:
        time, center, relative = series(row)
        config = row["config"]
        color = SURFACE_COLORS[config["surface"]]
        axes[0].plot(time, center[:, 2], color=color, alpha=.45, lw=.9)
        axes[1].plot(time, relative[:, 1], color=color, alpha=.45, lw=.9)
    axes[0].set(xlabel="Driven simulation time (ps)", ylabel="COM height above Au plane (nm)")
    axes[1].set(xlabel="Driven simulation time (ps)", ylabel="COM transverse displacement (nm)")
    for ax in axes:
        ax.grid(alpha=.18)
    fig.legend(handles=[Line2D([0], [0], color=SURFACE_COLORS[s], label=SURFACE_LABELS[s]) for s in SURFACES],
               loc="upper center", bbox_to_anchor=(.5, .89), frameon=False, ncol=2)
    fig.text(.08, .125, f"Flagged trajectories: {len(flagged)}. Their displacement metrics remain in the main figures and are marked by † or ×.", fontsize=8.5)
    fig.text(.08, .074, "A large single-atom height or molecular rotation alone does not trigger this rule. This threshold is not a thermodynamic desorption measurement.", fontsize=8.5)
    save(fig, "05-sustained-lift-diagnostics")


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--selected", help="Optional JSON list of temperature_K/strength conditions to show in detailed panels")
    args = parser.parse_args()
    report = load_report(OUT / "screen-report.json")
    extended_path = OUT / "extended-report.json"
    extended = None
    if extended_path.exists() and json.loads(extended_path.read_text(encoding="utf-8")).get("status") == "passed":
        extended = load_report(extended_path)
    reports = [report] + ([extended] if extended is not None else [])
    diagnostics = {r["key"]: lift_diagnostic(r) for item in reports for r in item["trials"]}
    FIGURES.mkdir(parents=True, exist_ok=True)
    plt.rcParams.update({"font.family": "DejaVu Sans", "font.size": 9,
                         "axes.spines.top": False, "axes.spines.right": False,
                         "axes.titlesize": 10, "svg.fonttype": "none"})
    parameter_matrix(report, diagnostics)
    surface_comparison(report, diagnostics)
    conditions = (json.loads(Path(args.selected).read_text(encoding="utf-8")) if args.selected else
                  [{"temperature_K": 300, "strength": 0.},
                   {"temperature_K": 600, "strength": .15},
                   {"temperature_K": 600, "strength": .25},
                   {"temperature_K": 300, "strength": .75}])
    trajectory_panels(report, conditions, diagnostics, name="03-selected-actual-trajectories")
    if extended is not None:
        selected = []
        for row in extended["summary"]:
            condition = {key: row[key] for key in ("temperature_K", "strength", "corrugation")}
            if condition not in selected:
                selected.append(condition)
        trajectory_panels(extended, selected, diagnostics, name="04-extended-surface-trajectories", compare_surfaces=True)
    height_figure(reports, diagnostics)
    metadata = {
        "reports": {item["phase"]: hashlib.sha256((OUT / (item["phase"] + "-report.json")).read_bytes()).hexdigest() for item in reports},
        "selected_screen_conditions": conditions,
        "lift_rule": "Every atom above z=0.8 nm for >=5 ps; geometric flag, not a thermodynamic desorption criterion",
        "diagnostics": diagnostics,
        "figures": sorted(path.name for path in FIGURES.glob("*.png")),
        "qualification": "Actual finite trajectories of C30H22 and three seeds; no confidence intervals, converged diffusion estimates or extrapolated finish times",
    }
    (OUT / "plot-metadata.json").write_text(json.dumps(metadata, ensure_ascii=False, indent=2) + "\n", encoding="utf-8")
    print(json.dumps({"figures": metadata["figures"], "sustained_lift_trials": sum(d["sustained_lift"] for d in diagnostics.values()),
                      "extended_included": extended is not None}, ensure_ascii=False))


if __name__ == "__main__":
    main()
