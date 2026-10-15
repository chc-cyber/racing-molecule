"""Publication-style plot from the measured paired preparation report only."""
from pathlib import Path
import json
import os
import sys
ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / "outputs/md-preparation-elements/preparation"
sys.path.insert(0, str(ROOT / ".local/analysis-libs"))
os.environ["MPLCONFIGDIR"] = str(OUT / "matplotlib-cache")
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import Patch

report = json.loads((OUT / "convergence-report.json").read_text(encoding="utf-8"))
if report["status"] != "passed" or not report["sources_unchanged_during_run"]:
    raise ValueError("Only a complete frozen-source acceptance report may be plotted.")
names = ("benzene-12", "test-car-52", "fullerene-c60", "ethynyl-rigid", "folded-eicosane")
labels = ("Benzene\n12 atoms", "Four-arm\n52 atoms", "C60\n60 atoms", "Ethynyl\n36 atoms", "Chain\n62 atoms")
surfaces = (("periodic-au", "Periodic Au"), ("fourier-au", "Fourier Au"), ("continuum", "Continuum"))
plt.rcParams.update({"font.family": "DejaVu Sans", "font.size": 10,
                     "axes.spines.top": False, "axes.spines.right": False,
                     "svg.fonttype": "none"})
fig, axes = plt.subplots(2, 3, figsize=(13.4, 7.2), sharex="col")
before_color, after_color = "#91A4B8", "#178568"
x = np.arange(len(names)); width = .34
for column, (surface, title) in enumerate(surfaces):
    rows = [next(row for row in report["cases"] if row["surface"] == surface and row["fixture"] == name) for name in names]
    for protocol, color, offset in (("former", before_color, -width/2), ("strict", after_color, width/2)):
        force = [row[protocol]["preparation"]["global_force_norm_kcal_mol_A"] for row in rows]
        seconds = [row[protocol]["preparation"]["wall_seconds"] for row in rows]
        axes[0, column].bar(x+offset, force, width, color=color)
        axes[1, column].bar(x+offset, seconds, width, color=color)
    axes[0, column].axhline(1e-4, color="#BA3D34", lw=1.2, ls="--")
    axes[0, column].set_yscale("log")
    axes[0, column].set_ylim(1e-6, 12.)
    axes[0, column].set_title(title, fontweight="bold")
    axes[1, column].set_yscale("log")
    axes[1, column].set_ylim(.004, 4.)
    axes[1, column].set_xticks(x, labels)
    for ax in axes[:, column]:
        ax.grid(axis="y", alpha=.17)
        ax.set_axisbelow(True)
axes[0, 0].set_ylabel("Global conservative force two-norm\n(kcal mol⁻¹ Å⁻¹)")
axes[1, 0].set_ylabel("Native joint preparation time (s)")
axes[0, 2].text(.97, .27, "Required: ≤ 10⁻⁴", transform=axes[0, 2].transAxes, ha="right", color="#BA3D34", fontsize=9)
fig.suptitle("Joint geometry relaxation — force convergence before race thermalization", fontsize=14, fontweight="bold", y=.99)
fig.legend(handles=[Patch(color=before_color, label="Former FIRE: 350 max, energy tolerance 10⁻⁶"),
                    Patch(color=after_color, label="Strict force gate: FIRE + HFTN, optional CG")],
           loc="upper center", bbox_to_anchor=(.5, .955), ncol=2, frameon=False, fontsize=10)
fig.text(.04, .025, "Each bar: one measured native preparation on the same input; fixed 0.5 fs race timestep and Langevin settings.\n"
         "Times exclude molecule/surface/cache construction and race MD. Both protocols were run separately; all strict cases passed independent force and native cleanup checks.", fontsize=9, color="#475569")
fig.tight_layout(rect=(0., .10, 1., .89))
for extension in ("png", "svg"):
    fig.savefig(OUT / f"preparation-before-after.{extension}", dpi=190)
plt.close(fig)
