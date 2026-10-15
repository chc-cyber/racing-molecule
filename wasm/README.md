# 纯浏览器 WASM 试玩版

当前已实现真实的浏览器 MD：LAMMPS 和 RDKit C++ UFF 编译成一个 WASM 模块，外场及冻结 Au 势在同一模块内求力，单独的 Web Worker 积分，现有 Three.js 界面播放真实采样坐标。玩家不安装 Python、LAMMPS、MCP 或本地计算服务；网站只发送静态文件。

这是核心试玩预览；在线入口与发布状态见根目录 README。保留现有本地版作为模型与数值验收基准，GitHub 使用同一仓库、两个运行入口，避免独立维护两个游戏。

当前规则版本为 `race-v1`，支持 20 nm 赛程。首次打开先选择赛事与分子，再准备；回放记录 `competition`（含赛事类型、规则版本、参数、起跑姿态及固定随机种子），不同类型成绩可分开处理。网页首页、工坊与赛道现已整合。工坊的导入／导出、补氢、ETKDG 构象、二维排版及 UFF 已接入浏览器 Worker，使用同一个编译模块中的独立 RDKit 工坊入口，不运行比赛求力。下文早期验收中的针尖／脉冲属于历史记录。

## 已支持

- 默认 20 nm × 20 nm，代理势 v2 / 150 K / 广域场 0.6 / 起伏 1.5 倍 / 350 fs 阻尼。
- 默认比赛时间步 1 fs、准备松弛 0.5 fs，每 50 fs 采样；LAMMPS NVE 积分加 Langevin 热浴；所有力与积分使用 float64，不启用 fast-math。
- 默认傅里叶代理 v2，可切换完整冻结周期 Au 模型或原连续表面；固定共价法向、法向链式力和保守缓存混合均保留。傅里叶支持域外回退完整方向项。
- 标准赛事与自定义赛事共用本地版界面。仅固定广域梯度场；针尖、手动目标与脉冲已移除。比赛支持开始、暂停、继续与重置。
- 标准参数固定；自定义在准备前选择温度、驱动、表面与起伏。比赛期间参数固定，修改后需要重新准备；步长 1 fs 与阻尼 350 fs 暂不开放编辑。
- 金晶格起伏倍率与广域场沿程增长分别调节，默认分别为 1.5 和 1；界面调整会重新准备并从第0步开始。
- 金原子跑道、球棍/空间填充切换、跟随/全跑道视角。
- 1×／2× 播放按钮位于模型切换下方，默认 100 ps 播放 50 秒，2× 为 25 秒；不足时缓冲，不外推坐标、不累计追赶时间。引擎按需提前计算约 10 ps，按顺序读取未播放采样。暂停／继续也适用于计算已经完成但动画尚未冲线的阶段。
- 倍速不修改 1 fs 时间步、求力、随机种子和比赛成绩；结果弹窗等最终原始坐标显示后出现。验证入口：`node --test tests/playback.test.mjs wasm/tools/transport.test.mjs`、`node scripts/verify_playback.cjs wasm`（需本机 Chrome 和 Playwright）；本机验收结果位于 `outputs/playback-check/`。
- 自动/手动起跑摆放、保留文件朝向，以及点选车尾/车头参考；摆放在原有零驱动联合松弛之前完成，比赛方向自由。
- 3D MOL/单记录 SDF 导入；2D 文件在 Worker 内补氢、ETKDGv3 构象生成、UFF 优化及联合松弛。上限为 256 原子（含氢）、45 Å 跨度、1 MB 文件，保持既有具体化学环境准入检查。
- H、B、C、N、O、F、Si、P、S、Cl、Br、I；C/H SAIP 游戏适配与其他元素范德华混合近似沿用本地模型，同位素质量不参与比赛。
- `joint-force-relaxation-v2`：原生阶段未收敛时增加曲率初始化 BFGS，使用同一保守目标及原准备预算；所有最终坐标重新验力，比赛积分与求力模型不变。
- 本机 JSON 轨迹导出；切到后台自动暂停计算，离开网页终止 Worker。

长时间运行的回放缓存限制约为 24 MB 原始坐标与基础帧数据；JavaScript 对象和导出 JSON 另有开销。超出后保留末尾轨迹，导出明确记录 `replay_complete=false`、丢弃帧数及第一保留步。普通 20 nm 比赛的完整轨迹可保留。

## 构建与本机预览

以下步骤已在当前 Windows / CPython 3.11 环境执行。开发者需要原项目 Python 依赖；玩家不需要这些工具。SDK、上游源码和编译缓存只安装在忽略目录 `.local/wasm/`，不改变全局环境。首次下载及编译较大，需预留数 GB；它们不属于网站载荷。

```powershell
.\.venv\Scripts\python.exe -E -B wasm/tools/bootstrap.py
.\.venv\Scripts\python.exe -E -B wasm/tools/build.py --jobs 4
.\.venv\Scripts\python.exe -E -B wasm/tools/export_assets.py
.\.venv\Scripts\python.exe -E -B wasm/tools/export_fourier.py
.\.venv\Scripts\python.exe -E -B wasm/tools/export_fourier_atlas.py
.\.venv\Scripts\python.exe -E -B wasm/tools/install_fourier_atlas.py
.\.venv\Scripts\python.exe -E -B wasm/tools/native_reference.py
node wasm/tools/validate.mjs
node --test wasm/tools/transport.test.mjs
.\.venv\Scripts\python.exe -E -B wasm/tools/package_web.py
.\.venv\Scripts\python.exe -E -B wasm/tools/verify_package.py
.\.venv\Scripts\python.exe -E -B wasm/tools/serve.py --port 8783
```

`native_reference.py` 需要现有 LAMMPS Python 配置，只用于开发验收。打开 `http://127.0.0.1:8783/`，从首页进入工坊，完成分子后准备比赛；这里的 Python 只发送文件，不执行 MD。也可把 `wasm/dist/` 交给任何静态服务器；不能直接双击 HTML 用 `file://` 游玩。选择未被现有服务占用的端口。

修改自有 C++ 后可用 `build.py --link-only` 复用依赖库。构建固定 Emscripten 6.0.12、LAMMPS stable_30Sep2026、RDKit Release_2025_09_3；版本见 `dependencies.json`，下载源 SHA-256 见 `source-checksums.json`。势表无损 gzip 压缩，解压后核对大小与 SHA-256。

生成结果：

- `wasm/dist/`：当前约 196 MB 静态网站（约 2100 个清单文件，包含完整工坊），主要为WASM、压缩势表与Three.js，含两代傅里叶参考表。
- `outputs/wasm-preparation/racing-molecule-web-preview.zip`：当前约 128 MB 的可运行静态网站包；默认 v2 模式首次准备时下载新碳表约30.58 MB。
- `outputs/wasm-preparation/wasm-validation.json`：实际 WASM 能量、受力、准入、松弛、传播与成对 NVE 验收结果。
- `outputs/wasm-preparation/browser-race.json`、`browser-preview.png`：实际浏览器完赛轨迹与画面。

当前构建路径以 Windows 为验收范围；尚未运行 GitHub Actions 或 Linux 构建，不宣称跨平台构建已验收。

## 晶格起伏与沿程增长验收，2026-10-14

环境设置新增`lattice_scale`（界面0–3）和`drive.growth`（0.1–2），默认均为1。周期势使用六相位平均和逐原子短程保护，从同一势能求力；非1倍率需要额外表面查询。增长调节保留x=0诱导前进力，同时改变沿程增幅及部分竖直/电荷作用。具体表达式、论文基准边界与147场扫描结果（含42场1.5倍起伏补测）见根目录README。现有0.5 fs步长、UFF及Langevin模型保持原值。

原生/WASM的28组新静态对照与12组2 ps配对NVE通过。真实Chrome已验证调整两项各自重建并归零、短程积分、暂停与模式切换；没有页面错误。55文件静态包的清单、源码、势表、ZIP内容和CRC均通过检查。开发比较推荐候选为150 K / 驱动0.60 / 起伏2 / 增长1；该轮验收时默认仍为600 K / 0.75 / 起伏1 / 增长1；现行默认见本页开头。

## 起跑摆放验收，2026-10-13

浏览器 Worker 调用 `prepareWithPose(content, filename, settings, startPose)`，与本地采用相同主轴阈值、符号和刚性旋转。旧三参数 `prepare` 及 `prepareChemistry` 继续保留文件朝向供既有数值参考使用。自动摆放只是几何起点；联合松弛可能改变朝向，不能保证最快起步。没有新增带场整备或比赛中的方向约束。

16 组原生/真实 WASM 几何对照最大坐标差 3.55e-15 Å；无效姿态明确拒绝。竖直自动、原始手动 90°、连接 C60 自动三个样例通过原 1e-4 kcal/mol/Å 准备门槛，从第 0 步完成 2000 步 / 1 ps。原有 224 组固定坐标力/能量及 5 ps 成对 NVE 回归通过。轨迹、数值范围和有限样本限制见根目录 README 与 `outputs/start-pose-validation/`。

```powershell
.\.venv\Scripts\python.exe scripts/validate_start_pose.py --output outputs/start-pose-validation/final
node wasm/tools/validate_start_pose.mjs
```

新导入默认自动；重置与更换表面保留已应用设置。重新应用姿态创建新 Worker、清空旧轨迹，再完成原准备流程。回放包含请求设置和实际摆放元数据。

## 连接 C60 的准备补救验收，2026-10-13

已将本地准备补救同步到浏览器 C++，并重新生成试玩包。72 原子的 C60—苯环样例在周期 Au 和连续表面上均通过原 1e-4 kcal/mol/Å 全局受力门槛，并完成各 10000 步 / 5 ps 试跑。预算耗尽正确拒绝，在同一个引擎对象上重新准备苯分子可继续积分。原有 224 组固定坐标力/能量及 5 ps 成对 NVE 对照重新通过。

开发验收可在已有构建后运行：

```powershell
.\.venv\Scripts\python.exe scripts/validate_preparation_rescue.py
.\.venv\Scripts\python.exe scripts/diagnose_joint_fullerene.py --output outputs/joint-fullerene-diagnosis/native-5ps --finite-difference --smoke-steps 10000
node wasm/tools/validate_fullerene.mjs
```

`validate_preparation_rescue.py` 读取前期元素验收在 `outputs/md-preparation-elements/elements/inputs/` 生成的代表分子；新环境可先运行 `scripts/validate_element_extension.py` 生成它们。浏览器补救报告写入 `outputs/joint-fullerene-diagnosis/wasm-report.json`，包括同坐标本地受力对照与该样例的 5 ps 质心对照。此处的 Node 验收执行真实 WASM，不代替多浏览器/移动设备验收。

## 数值与浏览器验收，2026-10-12

14 个分子 × 4 个固定构型位置 × 2 种驱动 × 2 种表面，共 224 组实际 WASM 对照，覆盖全部 12 种元素、同位素忽略、跨晶胞坐标、负坐标及两个缓存混合区。总能量最大绝对差为 2.85e-14 kcal/mol，总力最大分量差为 4.27e-14 kcal/mol/Å，均通过 1e-7 阈值。七类无效输入明确拒绝，并返回可读原因。

52 原子默认车通过完整联合松弛，残余全局力 8.67e-5 kcal/mol/Å，低于既有 1e-4 门槛。5 ps 成对 NVE 使用完全相同的初始坐标、速度、0.5 fs 步长及 100 步批次；最大逐原子位置差 1.09e-11 Å，最大总能量差 1.44e-11 kcal/mol。另一次关闭驱动与热浴的 5 ps 检查，最大能量偏差 0.140 kcal/mol，与有限时间步的非零误差相符；这不替代原有步长验证，也不声称任意分子已验证。

本机 Node 实际执行 10,000 步约 3.5 s，是单次计时，含帧采样，不含网页下载、GPU 绘制与播放节流；不作为跨浏览器、手机或不同 CPU 的性能承诺。长时间的热浴轨迹会对浮点差异敏感，不能要求不同平台逐帧位级相同。

Codex Chromium 浏览器里默认车实际以 113,100 步 / 56.55 ps 到达 x=200.218 Å，完整导出 1,132 帧。已核对重置、运行/暂停、两种驱动、手动针尖、两种表面和显示模型，二维乙醇 MOL 实际载入为 9 原子赛车。现有本地服务入口只读检查仍正常；此前 142 项 Python 回归检查通过。

## 非单位起伏求力优化，2026-10-14

`potentials.cpp` 的1.5/2倍等非单位起伏路径现在复用法向工作区、预平均六相位对齐网格，并用保守空间候选表减少 Au 位点遍历。候选表仅排除不可能位于16 Å截断内的位点，逐原子逐相位查询，保留原位点顺序和精确截断；长分子跨多个晶胞不依赖质心最小镜像。碳原子接触路径单独编译。1倍起伏继续原路径。

诊断接口`surfaceOptimization(0..3)`分别选择原算法、合并相位、增加预平均网格、再增加候选表；正式 Worker 默认3，未增加玩家设置。预平均表在首次非单位起伏求力时生成，新增约25 MB势表内存，另有候选索引；下载势表不变，本次WASM仅增加约13 KB。`forceTiming`是开发用固定坐标内核计时，不代表整个积分器耗时。

该轮本机串行MD耗时降低约13%–19%，未达到所有分子半分钟内计算完成。当时Worker的40 ms/100步节流尚未改变，后续傅里叶轮次已取消，见下节。完整记录见[优化报告](../outputs/surface-kernel-optimization/REPORT.md)。新环境需先按前文生成参考fixtures；本工作目录中可复现：

```powershell
node wasm/tools/validate_surface_optimization.mjs
node wasm/tools/validate_surface_elements.mjs
node wasm/tools/benchmark_surface_optimization.mjs
node wasm/tools/scan_surface_optimization.mjs --jobs scripts/data/surface-optimization-jobs.json
node scripts/verify_surface_optimization_browser.cjs
.\.venv\Scripts\python.exe scripts/analyze_surface_optimization.py
```

先单独运行串行性能测试，再运行多车比赛复核；两者不可同时用于单玩家耗时比较。优化前参考模块保存在`outputs/surface-kernel-optimization/reference/`，这些本地验收档案不属于普通源码发布内容。正式发布应另行分发需要复现的证据包。

## 首轮浏览器傅里叶实验模式，2026-10-14（历史）

本节更新前一轮优化的状态：浏览器已支持 `fourier-au`，且取消 Worker 每100步主动等待40 ms的节流。仍每100步采样并让出事件循环，0.5 fs时间步、UFF、热浴及驱动不变；浏览器调度仍有成本。

12个复谐波采用Float64及高度/方向Hermite插值，六相位直接合并到系数，保留完整位置/法向/混合权重导数。域外回退周期参考。只在选择时下载32.55 MB碳方向表（解压37.59 MB），氢方向项仍完整计算；完整开发包保留可选氢表。该轮发布时默认仍为周期势；现行默认已改为 v2。

当前机器1.5/2倍起伏、150 K、驱动0.6的串行2 ps样本中，苯与吡啶约6.6–7.5倍加速，52原子车约3.9–4.4倍，C60约1.4倍。真实Chrome Worker吡啶20 nm两场分别约18.4秒和20.6秒，不含下载、准备与Three.js绘制，不能推广成所有分子半分钟完成。完整结果及统计限制见[报告](../outputs/wasm-fourier-validation/REPORT.md)。

本工作目录的验收命令如下；独立新环境还需生成旧比赛参考输入/档案。先串行性能测试，再多任务比赛，不将并发耗时当成单玩家性能。

```powershell
.\.venv\Scripts\python.exe scripts/prepare_browser_fourier_validation.py
node wasm/tools/validate_browser_fourier.mjs
node wasm/tools/benchmark_browser_fourier.mjs
node wasm/tools/scan_browser_fourier.mjs
node scripts/verify_fourier_wasm_browser.cjs
node scripts/verify_fourier_ui.cjs
.\.venv\Scripts\python.exe scripts/analyze_browser_fourier.py
.\.venv\Scripts\python.exe wasm/tools/package_web.py
.\.venv\Scripts\python.exe wasm/tools/verify_package.py
```

## 三方向表与时间步研究，2026-10-14

浏览器选择傅里叶时现读取`fourier-atlas.json`，使用三张归一化平滑方向图，支持2.4–11.2 Å高度；C表解压36.99 MB、下载30.58 MB。v1系数仍保留给原数值测试，新`atlas_support.mjs`及验证工具直接读取独立生成目录。两代元数据和C++白名单均受检，包内保留两代资产，因此完整ZIP不是玩家首载量。

84组原生/WASM对照、18组全Cartesian有限差分、三图重叠、负法向、周期平移和退化法向均有验证。126个既有C60快照全部使用代理，最大力差0.03747 kcal/mol/Å；低高度/近退化压力构型误差较大，不是全域误差界。24场新代理比赛与旧周期/v1结果分开比较，详见[并行优化报告](../outputs/parallel-md-validation/REPORT.md)。

新增诊断`setTimestep`及准备设置`timestep_fs`支持0.25/0.5/1/2 fs，普通游戏UI固定默认1 fs，运行途中不允许用该诊断接口换步长。准备松弛仍固定0.5 fs，比赛采样保持50 fs、阻尼350 fs；成绩与超时以物理时间处理。此前128条2 ps NVE及逐步峰值补检作为精度与敏感性记录保留。现已按游戏需求采用1 fs，具体组合验收见[默认参数验收](../outputs/default-candidate-acceptance/report.html)；2 fs不作为默认。步长研究使用隔离的旧代理构建，最终主构建另做新代理与步长组合检查，勿混用不同扫描的并发耗时。

本工作目录保留可复现输入/参考档案，新环境须先生成它们。主要验证入口：

```powershell
node wasm/tools/validate_parallel_regression.mjs
node wasm/tools/validate_browser_atlas.mjs
node wasm/tools/scan_atlas_races.mjs
node wasm/tools/benchmark_parallel_md.mjs
node scripts/audit_parallel_md_worker.mjs
node scripts/verify_parallel_md_browser.cjs
node scripts/verify_atlas_ui.cjs
.\.venv\Scripts\python.exe scripts/analyze_parallel_md.py
```

串行benchmark和真实浏览器性能测试必须单独运行，不能与比赛扫描并行。时间步独立构建、运行命令及统计限制见`outputs/timestep-study/REPORT.md`。构建脚本现会检查编译期间源文件是否变动，避免并行开发生成错误源码哈希。

## GitHub 与网站发布方式

同一 GitHub 仓库保留本地 Python 入口、`wasm/src/` 与构建脚本，以及共用的 `molecular_car/` 参数资产和 `frontend/` 显示层。本地 Python/Numba 与浏览器 C++ 求力是两套实现，以参数版本、参考算例及数值验收约束一致性；目前不是同一份求力代码的两种绑定。

GitHub Releases 分别提供本地版安装/源码包、浏览器静态网站包及对应源代码材料。`wasm/dist/` 是构建产物，不放入普通源码提交；`.venv/`、`.local/`、测试轨迹和上游开发树同样不提交。

GitHub Pages 或其他静态托管只发布静态网站构建。资源 URL 已使用相对路径，可放在项目子路径。托管须为 `.js`/`.mjs` 设置 JavaScript MIME、`.wasm` 设置 `application/wasm`，提供 HTTPS；势表 `.bin.gz` 按原始 gzip 文件发送，浏览器自行解压，不要误加会导致再次解压的 `Content-Encoding`。本版单线程 Worker 不要求 SharedArrayBuffer 或 COOP/COEP。

正式发布前仍需：

1. 完成 Chrome/Edge、Firefox、Safari 与移动设备的加载、内存、发热及等待时间验收；当前未做全浏览器矩阵或弱设备性能基线。
2. 浏览器工坊已接入文件读写、二维／三维、构象准备与 UFF；仍需更广泛的车形、弱设备首载与长时操作验收。傅里叶采用现有默认代理。
3. 决定项目自身许可证，补齐公开分发 WASM 的对应源代码及构建材料入口。已保留编译依赖许可证与精确源版本，但当前预览 ZIP 不是完整的公开发行源码交付。
4. 在最终托管域名和项目子路径进行部署验收，包括文件类型、缓存版本一致性及真实网络首载；本轮只准备本地包，没有创建仓库或公开站点。
