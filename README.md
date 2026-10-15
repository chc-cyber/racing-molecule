# Racing Molecule · 分子赛车

设计自己的分子，让它在金表面完成一场 20 nm 比赛。LAMMPS 负责分子动力学，RDKit 提供分子处理与 UFF 求力，Three.js 呈现真实轨迹。

**[在线游玩](https://chc-cyber.github.io/racing-molecule/) · [下载发行包](https://github.com/chc-cyber/racing-molecule/releases) · [构建与发布](RELEASING.md)**

## 怎么玩

1. 进入工坊，用二维或三维工具设计分子，也可导入 MOL、单记录 SDF 或工程文件。
2. 选择标准赛事或自定义条件，保存并进入赛道。
3. 准备完成后开始比赛，支持暂停、继续、1× / 2× 播放，以及球棍 / 空间填充显示。
4. 沿赛道净前进 20 nm 后显示成绩，平均速度以 m/s 表示，可保存原始轨迹。

网页的计算在玩家设备上完成，无需安装 Python 或 LAMMPS，也不需要计算服务器。初次进入工坊及准备赛车需要下载较大的运行库和势表，网络较慢时请耐心等待。暂不承诺离线运行。

## 当前版本

这是开源预览版。默认条件为 **150 K、驱动 0.6、Au 代理势 v2、晶格起伏 1.5 倍、1 fs 时间步**，Langevin 阻尼时间 350 fs。标准赛事固定规则；自定义赛事可在准备前改变部分条件。

默认 1× 下，100 ps 轨迹播放约 50 秒，2× 约 25 秒。计算或下载不足时会缓冲，倍速不改变 MD 步长和比赛成绩。

支持 H、B、C、N、O、F、Si、P、S、Cl、Br、I；比赛分子最多 256 原子（含氢），跨度不超过 45 Å。具体化学环境还须通过力场、价态和几何松弛检查，元素在名单内不代表任何结构均适用。

这是以可玩性为目标的经典 MD 游戏模型，采用固定拓扑及有效电场参数；C/H–Au 为 SAIP 游戏适配，其他元素为非反应范德华混合近似。成绩不代表实验纳米车速度或合成可行性。

Chrome、Edge、Firefox 的本地完整流程已验收；公开 HTTPS 站点完成 23 项流程与缓存检查，详见 [发布验收](docs/RELEASE_ACCEPTANCE.md)。手机真机与 Mac/Safari 仍待验证。刷新会丢失未下载的比赛轨迹，页面会提示确认；工坊工程建议下载备份。

## 开发与本地运行

- [纯浏览器构建入口与发行流程](RELEASING.md)
- [WASM 实现、固定版本与验证命令](wasm/README.md)
- [二维 Ketcher 扩展源码与构建](frontend/editor/ketcher-source/README.md)
- [三维 JSmol 修改与对应源码](frontend/editor/jsmol-source/README.md)
- [开发历史与 MD 数值验证记录](DEVELOPMENT_NOTES.md)

本地 Python 版本需要 Windows、兼容的 64 位 CPython 3.11 与 LAMMPS GUI 发行版：

```powershell
.\scripts\setup.ps1 -PythonExe 'C:\path\to\python.exe' -LammpsHome 'C:\path\to\LAMMPS'
.\start-game.ps1
```

打开启动器显示的本地地址。本地与网页入口共用同一项目，网页玩家不需要执行上述安装。

## 许可证

自有代码采用 [MIT](LICENSE)。第三方组件及其修改保留原许可证；含 LAMMPS 的链接 WASM 模块按 GPLv2 分发，MIT 不覆盖整个发行包。发行页提供对应源码、构建输入及校验清单。详见 [第三方说明](THIRD_PARTY.md) 和 [工坊第三方说明](frontend/editor/THIRD_PARTY.md)。
