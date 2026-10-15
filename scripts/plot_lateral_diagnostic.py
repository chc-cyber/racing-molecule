"""Plot an observed race COM trajectory without modifying the running game."""
import csv
import json
import os
from pathlib import Path
import sys

ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / "outputs/lateral-diagnostic"
sys.path.insert(0, str(ROOT / ".local/analysis-libs"))
os.environ["MPLCONFIGDIR"] = str(OUT / "matplotlib-cache")
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

with (OUT / "current-com.csv").open(encoding="utf-8-sig", newline="") as source:
    rows = list(csv.DictReader(source))
summary = json.loads((OUT / "current-summary.json").read_text(encoding="utf-8-sig"))
t = np.array([float(row["Time_ps"]) for row in rows])
x = np.array([float(row["X_nm"]) for row in rows])
y = np.array([float(row["Dy_nm"]) for row in rows])
plt.rcParams.update({"font.family": "Microsoft YaHei", "font.size": 11,
                     "axes.spines.top": False, "axes.spines.right": False})
fig, axes = plt.subplots(1, 2, figsize=(13.1, 5.3), gridspec_kw={"width_ratios": [1.05, 1]})
fig.patch.set_facecolor("#f7f7f0")
axes[0].plot(t, y, color="#598373", linewidth=1.45)
axes[0].axhline(0, color="#aaad9f", lw=.8, linestyle="--")
axes[0].set(xlabel="仿真时间（ps）", ylabel="横向位移（nm）", title="真实横向轨迹：偏离后又回到起点附近",
            xlim=(0, t[-1] * 1.04), ylim=(-.8, 1.5))
peak = int(y.argmax())
axes[0].scatter([t[peak], t[-1]], [y[peak], y[-1]], color="#b68a38", s=35, zorder=3)
axes[0].annotate(f"{t[peak]:.2f} ps：+{y[peak]:.2f} nm", (t[peak], y[peak]),
                 xytext=(30, 1.29), arrowprops={"arrowstyle":"-", "color":"#b68a38"}, fontsize=10)
axes[0].annotate(f"结束：{y[-1]:+.2f} nm", (t[-1], y[-1]),
                 xytext=(160, .53), arrowprops={"arrowstyle":"-", "color":"#b68a38"}, fontsize=10)
axes[0].grid(alpha=.13)
axes[1].axhspan(-25, 25, color="#efe4bd", alpha=.6)
axes[1].axhline(0, color="#d8c783", linestyle="--", lw=.8)
axes[1].plot(x, y, color="#598373", linewidth=1.5)
axes[1].set(xlabel="沿主干道的位置（nm）", ylabel="横向位置变化（nm）",
            title="放回 100 × 50 nm 赛道：偏移显得很小", xlim=(0, 101), ylim=(-25, 25))
axes[1].set_aspect("equal", adjustable="box")
axes[1].grid(alpha=.12)
fig.suptitle("你刚跑完的比赛 · 真实质心坐标 · 相机不会改变此数据", fontsize=15, y=.98)
fig.text(.04, .02, "单条轨迹不是扩散系数测量；比赛中驱动强度由 0.25 调到了 0.70。", fontsize=9, color="#626b5c")
fig.tight_layout(rect=[.01, .055, .99, .945], w_pad=2.3)
fig.savefig(OUT / "current-lateral-trajectory.png", dpi=150, facecolor=fig.get_facecolor())
fig.savefig(OUT / "current-lateral-trajectory.svg", facecolor=fig.get_facecolor())
print(OUT / "current-lateral-trajectory.png")
