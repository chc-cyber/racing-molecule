"""Plot measured trajectories from verify_thermal_drive's report.

Matplotlib is an optional analysis dependency, independent of the game runtime.
"""
import json
import os
from pathlib import Path
import sys

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / ".local/analysis-libs"))
OUT = ROOT / "outputs/thermal-drive-check"
os.environ["MPLCONFIGDIR"] = str(OUT / "matplotlib-cache")
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

report = json.loads((OUT / "report.json").read_text(encoding="utf-8"))
if not report["passed"]:
    raise ValueError("Run the trajectory verification successfully first")
plt.rcParams.update({"font.family": "Microsoft YaHei", "font.size": 11,
                     "axes.spines.top": False, "axes.spines.right": False})
fig, axes = plt.subplots(1, 2, figsize=(13.2, 5.4), gridspec_kw={"width_ratios": [1, 1.05]})
fig.patch.set_facecolor("#f7f7f0")
colors = ("#64937b", "#b68c2e", "#bd6857", "#5368a4")
for group, color in zip(report["groups"], colors):
    label = f"{group['temperature_K']} K / 强度 {group['strength']}"
    for index, row in enumerate(group["trials"]):
        final = np.array(row["samples"][-1]["center"])
        initial = final-np.array(row["displacement_nm"])*10
        coords = np.array([initial] + [sample["center"] for sample in row["samples"]])
        times = [0] + [sample["time_fs"]/1000 for sample in row["samples"]]
        delta = (coords-initial)/10
        for ax in axes:
            ax.set_facecolor("#ffffff")
        axes[0].plot(times, delta[:, 1], color=color, lw=1.7, alpha=.75,
                     linestyle="--" if group["strength"] < .75 else "-",
                     label=label if index == 0 else None)
        axes[1].plot(delta[:, 0], delta[:, 1], color=color, lw=1.8, alpha=.8,
                     linestyle="--" if group["strength"] < .75 else "-")
        axes[1].plot(delta[-1, 0], delta[-1, 1], "o", color=color, ms=3)
axes[0].axhline(0, color="#d5d6cd", lw=.9)
axes[0].set(xlabel="仿真时间（ps）", ylabel="横向位移（nm）", title="升温增加横向游走：局部放大")
axes[0].legend(frameon=False, fontsize=9, loc="upper left")
axes[0].grid(alpha=.14)
axes[1].set(xlabel="向终点方向位移（nm）", ylabel="横向位移（nm）", title="俯视轨迹：两轴比例相同")
axes[1].set_xlim(-.8, 20)
axes[1].set_ylim(-5, 5)
axes[1].set_aspect("equal", adjustable="box")
axes[1].grid(alpha=.14)
axes[1].text(.02, .05, "赛道宽 50 nm；此图仅显示中央 10 nm\n蓝色虚线：弱驱动，前进慢，偏航占比更大",
             transform=axes[1].transAxes, fontsize=9, color="#596156")
fig.suptitle("真实 LAMMPS 轨迹 · 每组 3 个种子 · 50 ps · 固定广域场", fontsize=15, y=.99)
fig.text(.05, .025, "短轨迹对比，不能据此给出迷路概率；高温固定拓扑 UFF 不代表真实化学稳定性。",
         fontsize=9, color="#68705f")
fig.tight_layout(rect=[.015, .055, .99, .945], w_pad=2.5)
fig.savefig(OUT / "trajectory-comparison.png", dpi=150, facecolor=fig.get_facecolor())
fig.savefig(OUT / "trajectory-comparison.svg", facecolor=fig.get_facecolor())
print(OUT / "trajectory-comparison.png")
