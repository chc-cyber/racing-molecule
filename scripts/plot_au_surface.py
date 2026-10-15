"""Plot the completed explicit Au/benzene reference using optional Matplotlib.

The analysis environment is separate from LAMMPS and the game dependencies.
CSV energies are actual LAMMPS SAIP evaluations or the unchanged game surface.
"""

import json
import os
from pathlib import Path
import sys

ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / "outputs/au-surface-check"
sys.path.insert(0, str(ROOT / ".local/analysis-libs"))
os.environ["MPLCONFIGDIR"] = str(OUT / "matplotlib-cache")

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.colors import Normalize
import numpy as np


def rows(filename):
    return np.genfromtxt(OUT / filename, delimiter=",", names=True, encoding="utf-8")


def grid(data):
    """Arrange samples by recorded fractional coordinates, not inferred order."""
    fx = np.unique(data["fx"])
    fy = np.unique(data["fy"])
    arrays = {key: np.empty((len(fy), len(fx))) for key in ("x_A", "y_A", "energy_eV")}
    for sample in data:
        ix = int(np.searchsorted(fx, sample["fx"]))
        iy = int(np.searchsorted(fy, sample["fy"]))
        for key in arrays:
            arrays[key][iy, ix] = sample[key]
    x, y, energy = (arrays[key] for key in ("x_A", "y_A", "energy_eV"))
    return x - x[0, 0], y - y[0, 0], (energy - energy.min()) * 1000


report = json.loads((OUT / "report.json").read_text(encoding="utf-8"))
if report["status"] != "passed" or report.get("slab_only"):
    raise ValueError("Complete verify_au_surface.py including the benzene scans first.")
benzene = report["benzene"]
comparison = benzene["continuous_game_surface_comparison"]
height_saip = rows("benzene-saip-height.csv")
height_continuous = rows("benzene-continuous-height.csv")
saip_data = rows("benzene-saip-lateral.csv")
saip_data = saip_data[np.isclose(saip_data["angle_deg"], 0)]
continuous_data = rows("benzene-continuous-lateral.csv")
continuous_data = continuous_data[np.isclose(
    continuous_data["height_A"], comparison["at_SAIP_angle0_atop_reference_height_angstrom"], atol=1e-8)]
saip_grid, continuous_grid = grid(saip_data), grid(continuous_data)
normalization = Normalize(0, max(saip_grid[2].max(), continuous_grid[2].max()))

plt.rcParams.update({"font.family": "Microsoft YaHei", "font.size": 10,
                     "axes.spines.top": False, "axes.spines.right": False})
fig = plt.figure(figsize=(14, 5.4), facecolor="#f7f7f0")
layout = fig.add_gridspec(1, 3, width_ratios=(1.15, 1, 1))
axes = [fig.add_subplot(layout[i]) for i in range(3)]
for ax in axes:
    ax.set_facecolor("white")

axes[0].plot(height_saip["height_A"], height_saip["energy_eV"], color="#ad8427", lw=2,
             label="显式 Au + SAIP（苯中心位于顶位）")
axes[0].plot(height_continuous["height_A"], height_continuous["energy_eV"], color="#537765", lw=2,
             label="现有连续势（固定横向位置）")
for z, energy, color in (
    (benzene["atop_vertical_minimum_height_angstrom"], benzene["atop_vertical_minimum_energy_eV"], "#ad8427"),
    (comparison["anchor_vertical_minimum_height_angstrom"], comparison["anchor_vertical_minimum_energy_eV"], "#537765"),
):
    axes[0].scatter([z], [energy], s=32, color=color, zorder=4)
    axes[0].annotate(f"{z:.3f} Å / {energy:.3f} eV", (z, energy),
                     xytext=(7, -22 if color == "#ad8427" else 16),
                     textcoords="offset points", color=color, fontsize=9)
axes[0].axhline(0, color="#c7c8bf", lw=.8)
axes[0].set(xlabel="苯平面到顶层 Au 中心平面的高度（Å）", ylabel="分子–表面势能（eV）",
            title="吸附高度与势阱深度", xlim=(2.4, 5.4), ylim=(-.9, 1.0))
axes[0].legend(frameon=False, fontsize=8, loc="upper right")
axes[0].grid(alpha=.15)

for ax, values, title in zip(axes[1:], (saip_grid, continuous_grid),
                              ("显式 Au + SAIP", "现有连续表面")):
    x, y, relative_energy = values
    mesh = ax.pcolormesh(x, y, relative_energy, cmap="YlOrBr", norm=normalization, shading="gouraud")
    ax.plot([x[0, 0], x[0, -1], x[-1, -1], x[-1, 0], x[0, 0]],
            [y[0, 0], y[0, -1], y[-1, -1], y[-1, 0], y[0, 0]], color="#675c41", lw=.8)
    ax.set(xlabel="横向平移 x（Å）", ylabel="横向平移 y（Å）", aspect="equal",
           title=f"{title}\n网格最大值−最小值：{relative_energy.max():.2f} meV")
    ax.grid(alpha=.1)

fig.subplots_adjust(left=.065, right=.9, top=.8, bottom=.22, wspace=.35)
colorbar = fig.colorbar(mesh, ax=axes[1:], fraction=.028, pad=.04)
colorbar.set_label("相对各自网格最低点的势能（meV）")
fig.suptitle("小范围 Au(111) 参考计算 · 刚性平面苯 · 外场关闭", fontsize=15, y=.96)
height = benzene["atop_vertical_minimum_height_angstrom"]
fig.text(.49, .835, f"横向图：两模型均固定高度 {height:.3f} Å、同一苯朝向；共用色阶", fontsize=10)
fig.text(.065, .1, "2880 个 Au / 六层 / 底两层固定 / 横向周期边界。图中平移原胞各自从零计，未拟合两模型的晶格相位。", fontsize=9, color="#5c655c")
fig.text(.065, .055, "横向起伏是固定高度、固定朝向的网格范围，不是最小扩散路径的势垒；本结果未改写游戏参数。", fontsize=9, color="#5c655c")
fig.savefig(OUT / "surface-comparison.png", dpi=160, facecolor=fig.get_facecolor())
fig.savefig(OUT / "surface-comparison.svg", facecolor=fig.get_facecolor())
print(OUT / "surface-comparison.png")
