"""Scientific finite-size scaling figure from isolated actual-game timings."""
from __future__ import annotations

import hashlib
import json
import os
from pathlib import Path
import sys
import zipfile

ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / "outputs/atom-count-scaling"
sys.path.insert(0, str(ROOT / ".local/analysis-libs"))
os.environ["MPLCONFIGDIR"] = str(OUT / "matplotlib-cache")
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np


def main():
    report = json.loads((OUT / "report.json").read_text(encoding="utf-8"))
    if report["status"] != "passed":
        raise ValueError("Only accepted finite-size measurements can be plotted")
    for name, expected in report["source_sha256"].items():
        if hashlib.sha256((ROOT / name).read_bytes()).hexdigest() != expected:
            raise ValueError(f"Source changed after measurement: {name}")
    rows = [r for r in report["summary"] if r["family"] == "four-arm saturated alkane"]
    n = np.array([r["atoms"] for r in rows])
    full = np.array([r["full_MD_us_per_step"] for r in rows])
    uff = np.array([r["UFF_static_us"] for r in rows])
    external = np.array([r["external_static_us"] for r in rows])
    fit = report["analysis"]["full_cost_fit"]
    coeff = fit["coefficients"]
    a, b, c = (coeff[name] for name in ("constant_us", "linear_us_per_atom", "quadratic_us_per_atom_squared"))
    line = np.linspace(n.min(), n.max(), 250)
    predicted = a + b * line + c * line**2
    fig, axes = plt.subplots(1, 2, figsize=(12.5, 5.2))
    fig.subplots_adjust(left=.07, right=.985, top=.78, bottom=.23, wspace=.30)
    fig.suptitle("Actual-game MD cost versus explicit atom count", fontsize=15, y=.97)
    fig.text(.5, .87, "Same four-arm alkane family | One thread | 0.5 fs | 600 K | Three seeds, three timed batches per seed", ha="center", fontsize=9)
    low = np.array([r["trial_min_MD_us"] for r in rows])
    high = np.array([r["trial_max_MD_us"] for r in rows])
    axes[0].errorbar(n, full, yerr=[full-low, high-full], fmt="o", color="#0072B2", capsize=3,
                     label="Actual full MD: median; seed range", zorder=3)
    axes[0].plot(line, predicted, color="#D55E00", label="Descriptive a + bN + cN² fit")
    axes[0].text(.04, .96, f"T = {a:.2f} + {b:.3f}N + {c:.5f}N² µs\nR² = {fit['r_squared']:.4f}",
                 transform=axes[0].transAxes, va="top", fontsize=9)
    controls = [r for r in report["summary"] if r["family"] == "existing game fixture"]
    axes[0].scatter([r["atoms"] for r in controls], [r["full_MD_us_per_step"] for r in controls],
                    marker="x", s=50, color="#777777", label="Game fixtures: excluded from fit")
    axes[0].set(title="A   Complete propagation cost", xlabel="Explicit atoms, including H", ylabel="Time (µs / MD step)")
    axes[0].legend(frameon=False, fontsize=8, loc="lower right")
    axes[1].plot(n, uff, "o-", color="#0072B2", label="UFF energy + gradient, static diagnostic")
    axes[1].plot(n, external, "s-", color="#009E73", label="Compiled external fields, static diagnostic")
    axes[1].set(title="B   Force-component diagnostic", xlabel="Explicit atoms, including H", ylabel="Time (µs / evaluation)", yscale="log")
    axes[1].legend(frameon=False, fontsize=8)
    for ax in axes:
        ax.grid(alpha=.18)
        ax.spines[["top", "right"]].set_visible(False)
        ax.set_xticks(n[::2].tolist() + [int(n[-1])])
    fig.text(.07, .095, "Preparation, minimization, compilation, frames and rendering excluded. Fit describes only this measured family and size range.", fontsize=8.5)
    fig.text(.07, .047, "Static component timings include bindings and do not add exactly to full MD. All eligible intramolecular vdW pairs retained; no model changes.", fontsize=8.5)
    fig.savefig(OUT / "atom-count-scaling.png", dpi=180, facecolor="white")
    fig.savefig(OUT / "atom-count-scaling.svg", facecolor="white")
    plt.close(fig)
    archive = OUT / "atom-count-scaling-results.zip"
    with zipfile.ZipFile(archive, "w", compression=zipfile.ZIP_DEFLATED) as handle:
        for path in sorted(OUT.rglob("*")):
            if (path.is_file() and path != archive and "matplotlib-cache" not in path.parts
                    and path.name != "trials-progress.json"
                    and path.suffix.lower() in (".json", ".csv", ".mol", ".png", ".svg", ".md")):
                handle.write(path, path.relative_to(OUT).as_posix())
        for name in ("scripts/benchmark_atom_scaling.py", "scripts/plot_atom_scaling.py"):
            handle.write(ROOT / name, "implementation/" + name)
    with zipfile.ZipFile(archive) as handle:
        if handle.testzip() is not None:
            raise ValueError("Scaling archive CRC verification failed")
    print(json.dumps({"figure": str(OUT / "atom-count-scaling.png"), "archive": str(archive), "fit": fit}, ensure_ascii=False))


if __name__ == "__main__":
    main()
