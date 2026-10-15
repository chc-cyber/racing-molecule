# 保存并进入赛道

顶部保留“保存工程”和一个主按钮“保存并进入赛道”。主按钮顺序执行：保存工程 → 检查结构 → 准备赛车 → 切换赛车 → 同页进入赛道。进入后赛车处于准备就绪状态，点击“开始试跑”才开始模拟。“返回工坊”恢复浏览器保存的设计。

保存工程会下载 `.mcar`，同时保存同一个设计快照到浏览器。浏览器负责文件下载，网页不能确认用户是否取消下载或磁盘是否写入完成；本阶段没有新增文件选择确认。进入赛道依赖成功的浏览器暂存，准备、松弛和起跑摆放仅作用于赛车副本。未应用的几何预览继续使用已有的应用/取消按钮，解决后主按钮才可用。

准备中显示当前阶段和“取消”。准备失败显示实际原因和“返回编辑”；连接中断时加一个“重试”。没有通过检查后的二次确认。准备期间编辑区暂时冻结，取消后立即恢复，包括服务暂时无法访问的情况。取消和返回编辑不会重新打开文件或重置工坊的撤销记录。历史记录不跨页面持久化，返回工坊恢复保存的工程内容。

## 后端协议

新增模块 `molecular_car/editor/race_handoff.py`；现有 `/api/load` 和模拟准备算法保持原语义。新接口：

| 接口 | 内容 |
|---|---|
| POST `/api/editor/race/prepare` | `{id, filename, content}`；id 是客户端生成的 UUID |
| GET `/api/editor/race/status?id=…` | 独立任务状态：checking / preparing / ready / error / cancelled / committed |
| POST `/api/editor/race/cancel` | `{id}`；取消候选，已提交时返回提交状态 |
| POST `/api/editor/race/commit` | `{id}`；仅准备就绪候选可切换，返回新 generation |
| GET `/api/replay/previous` | 下载最近一次被候选替换的比赛轨迹 |

候选进程直接使用原 `engine.worker_main`。在独立子进程中围绕原 `prepare_molecule` 增加阶段通知，原返回对象直接交给 `Simulation` 使用，没有额外的预检查优化或二次构象生成。checking 包括当前化学/参数检查、补氢、必要的三维生成和真空 UFF；preparing 包括原表面、LAMMPS 和联合松弛。状态标签不是独立化学操作。

原比赛在候选准备期间继续运行或保持暂停，`/api/state` 和 `/api/frames` 始终代表它。候选失败和取消只释放候选进程。提交使用生命周期锁和状态锁：核对旧比赛 generation、环境与驱动设置、候选进程存活；保存旧回放；接管已经准备好的同一进程和第 0 帧；增加 generation。旧进程通过捕获的独立引用停止，不会停止新进程。提交响应丢失后再次提交同一 ID 不会重建或再次替换。

准备请求使用 ID 和文件指纹去重，网络重试先恢复状态；同一 ID 的不同文件被拒绝。取消可以先于迟到的准备请求到达。工坊用 sessionStorage 记录任务和保存的快照，刷新后继续恢复；暂存不可用时当前页面仍能重试，刷新恢复受浏览器存储可用性限制。断网取消先在本页屏蔽迟到结果并解锁编辑，保留取消记录；后台释放失败时，在刷新或下一次准备前重试释放，服务端也会回收过期任务。提交前再次比较工坊快照，防止旧结果覆盖新设计。更换赛道或准备期间修改环境/驱动会使候选提交失败，需要重新准备。

## 资源和轨迹

最多一辆未提交候选，同时最多两个模拟进程。候选准备最长 180 秒；就绪候选最多保留 10 分钟；最近 8 个任务保留状态供重试。取消、失败、过期和服务退出释放候选进程。

只保留最近一场被替换比赛的轨迹，在赛道“保存轨迹”菜单选择当前或上一场。旧轨迹仍沿用原 5010 帧缓冲限制；超出时回放标记 truncated。该记录存于服务内存，服务重启后丢失。原 `/api/load`、重置、切换表面等入口仍沿用原流程，此次事务保护适用于工坊候选流程。

## 运行和验证

沿用项目启动方式打开 `/editor/index.html`。独立验收实例使用 8786，默认表面验收实例使用 8787；不要对正在玩的服务运行验收，因为测试会主动比赛和替换赛车。

```powershell
.\.venv\Scripts\python.exe -X utf8 -E -B -m molecular_car.server --port 8786
node frontend/editor/tests/workshop-completion.cjs http://127.0.0.1:8786
node scripts/verify_editor_seed_browser.cjs
node --disable-warning=MODULE_TYPELESS_PACKAGE_JSON --test tests/editor_race_completion.test.mjs
.\.venv\Scripts\python.exe -X utf8 -E -B -m unittest tests.test_editor_race_handoff tests.test_surface_modes.SurfaceLifecycleTests tests.test_track_gameplay tests.test_start_pose
```

默认表面验收另开服务，然后运行：

```powershell
.\.venv\Scripts\python.exe -X utf8 -E -B -m molecular_car.server --port 8787 --molecule examples/benzene.mol
.\.venv\Scripts\python.exe -X utf8 -E -B scripts/verify_editor_handoff.py
```

新流程未增加第三方依赖，组件版本和许可证沿用 [THIRD_PARTY.md](THIRD_PARTY.md)。实际结果见 [VALIDATION.md](VALIDATION.md) 及 `outputs/editor-completion-check/report.json`、`outputs/editor-handoff-default/report.json`。
