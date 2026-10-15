"""Illustrative fixed-electrode field; no game or MD parameters are changed."""
from pathlib import Path
import os
import sys
ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / "outputs/static-field-concept"
OUT.mkdir(parents=True, exist_ok=True)
sys.path.insert(0, str(ROOT / ".local/analysis-libs"))
os.environ["MPLCONFIGDIR"] = str(OUT / "matplotlib-cache")
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib import font_manager, patches
import numpy as np

font_manager.fontManager.addfont("C:/Windows/Fonts/msyh.ttc")
plt.rcParams.update({"font.family":"Microsoft YaHei", "axes.unicode_minus":False})
h, wavelength, amplitude, base, slope, car_z = 1.5, 6., .45, 1., .03, .3
k = 2*np.pi/wavelength
x = np.linspace(-6,6,501)
z = np.linspace(0,h,151)
X,Z = np.meshgrid(x,z)

def values(X,Z):
    ratio = np.sinh(k*Z)/np.sinh(k*h)
    phi = -(base+slope*X)*Z+amplitude*ratio*np.cos(k*X)
    ex = slope*Z+amplitude*k*ratio*np.sin(k*X)
    ez = base+slope*X-amplitude*k*np.cosh(k*Z)/np.sinh(k*h)*np.cos(k*X)
    return phi,ex,ez

phi,ex,ez=values(X,Z)
fig,axes=plt.subplots(1,2,figsize=(12.4,5.3),gridspec_kw={"width_ratios":[1.45,1]})
fig.subplots_adjust(left=.065,right=.975,top=.75,bottom=.24,wspace=.28)
fig.suptitle("固定电极阵列在完整金面上方形成空间场纹理",fontsize=16,y=.97)
fig.text(.5,.87,"截面示意：条纹电极（ky=0），高度 h=1.5 nm，波长 6 nm；幅度为示意值，未经赛车标定",ha="center",fontsize=9)
ax=axes[0]
ax.contourf(X,Z,phi,levels=25,cmap="Blues_r",alpha=.5)
ax.streamplot(x,z,ex,ez,color="#315c72",density=(1.5,.65),linewidth=.7,arrowsize=.8)
for left in np.linspace(-6,5.5,24):
    center=left+.25
    color=plt.colormaps["coolwarm"]((np.cos(k*center)+1)/2)
    ax.add_patch(patches.Rectangle((left,h+.015),.46,.08,facecolor=color,edgecolor="white",lw=.5,clip_on=False))
ax.axhline(0,color="#b89939",lw=9)
ax.axhline(car_z,color="#cd6b26",ls="--",lw=1.2)
ax.text(0,h+.145,"上方：固定分段电极，分别保持不同电压",ha="center",fontsize=10)
ax.text(0,.055,"完整接地 Au 面，电势 φ=0",ha="center",fontsize=10,color="#7f6927",bbox={"facecolor":"white","edgecolor":"none","alpha":.9,"pad":2})
ax.text(-5.8,car_z+.06,"赛车所在高度",fontsize=9,color="#b34c0e")
ax.set(xlim=(-6,6),ylim=(-.03,h+.05),xlabel="沿赛道位置 x / nm",ylabel="离金面高度 z / nm")
_,ex_car,ez_car=values(x,np.full_like(x,car_z))
intensity=ex_car**2+ez_car**2
u=-.5*intensity
u-=np.mean(u)
axes[1].plot(x,u,color="#b34c0e",lw=2)
axes[1].axhline(0,color="#888888",lw=.6)
axes[1].set(xlabel="沿赛道位置 x / nm",ylabel="诱导能量（示意单位，扣除均值）",title="同一高度：强场区域对应低诱导能量")
axes[1].grid(alpha=.2)
axes[1].text(.04,.05,"Uind ∝ −|E|²\n热扰动可使赛车越过局部势垒",transform=axes[1].transAxes,va="bottom",fontsize=10)
for ax in axes:
    ax.spines[["top","right"]].set_visible(False)
fig.text(.065,.105,"左图线条是电场线，并非赛车轨迹。右图为标量可极化点探针，未加入整车构象、表面势或 Langevin 热浴。",fontsize=9)
fig.text(.065,.05,"电压图案与所有电极位置均固定。增加 y 方向变化后，条纹扩展为二维场纹理；此候选场尚未接入游戏。",fontsize=9)
for extension in ("png","svg"):
    fig.savefig(OUT/f"fixed-electrode-field.{extension}",dpi=180,facecolor="white")
print(OUT/"fixed-electrode-field.png")
