# 实际验证记录

## 2026-10-14 当前版本：保存并进入赛道

主按钮串联同一设计快照的工程下载、检查、候选赛车准备、原子切换和同页进入赛道。仅在实际准备完成后切换到已准备好的同一 LAMMPS 进程，第 0 帧处于暂停状态。现有 `/api/load` 语义与模拟准备算法不变；工坊改用独立候选接口。

| 验证 | 结果 | 证据 |
|---|---:|---|
| 实际下载、预览阻拦、保存/构象/ID/撤销同步、失败保留比赛、断网取消、重试和刷新恢复、同页进入、旧轨迹下载、真实试跑、窄屏 | **10/10** | `outputs/editor-completion-check/report.json` |
| 最终页面简洁文案及具体错误指引 | **14/14** | `outputs/editor-copy-check/report.json` |
| 原生二维/三维绘制、C60/甲基、二面角、稳定 ID、共同撤销、过期异步结果 | **14/14** | `outputs/editor-seed-check/browser-report.json` |
| 分子模型/同步/几何与新完成状态机 | **96/96** | `outputs/editor-completion-unit.log` |
| 候选生命周期/原生命周期/起跑摆放/试跑控制与独立几何服务 | **82/82** | `outputs/editor-completion-backend.log` |
| 默认 fourier-au 实际联合松弛失败保留旧车/轨迹、成功提交去重、实际 LAMMPS 积分 | **3/3** | `outputs/editor-handoff-default/report.json` |

浏览器完成验收使用任务自管的 8786 实例和既有 continuum 表面，真实执行 LAMMPS。单独 8787 实例使用当前默认 fourier-au，未改势函数、准备阈值、温度或驱动。默认表面下四苯基苯候选真实经过 checking → preparing → error：30 秒软预算耗尽，全局残余力 2.68029 kcal/mol/Å，原比赛 generation 不变，轨迹 seq 从 0 继续到 17，已有帧前缀完全一致。这是实际模拟失败，未伪造错误响应。之后苯候选真实通过 checking → preparing → ready，提交两次只切换一次，真实积分到第 50 步；残余力和最终复核记录保存在报告内。

浏览器网络恢复验收有意丢弃一条已经被服务器接受的准备响应，并实际点击重试、刷新页面；恢复同一任务且没有重复准备。取消验收有意中断取消请求，页面立即恢复编辑，后续操作先完成后台释放。未调用旧 `/api/load`，没有外部请求或浏览器异常。已保存工程与回到工坊的原子、三维坐标逐项一致。最近一场轨迹通过真实下载菜单导出并逐项核对。

最终报告列出源码 SHA-256，验证期间源码无变化。文案验证的未收敛展示仍是对真实 UFF 响应附加的受控提示，不当作物理未收敛验收；真正的物理准备失败以上述 fourier-au 报告为准。

较早的全量 `tests.test_surface_modes` 运行有两个已过时的测试失败：一个仍要求默认 surface 为 periodic-au，另一个模拟旧 FourierAuPotential 入口，与当前 fourier-au/atlas 默认实现不符。本次涉及的 SurfaceLifecycleTests 单独通过，未为通过旧断言而改模拟模型。这两项旧测试不计入上表通过项。

没有新增第三方依赖。工程保存使用浏览器下载及暂存，不承诺可检测磁盘写入完成；撤销历史不跨页面持久化。上一场轨迹最多保留一场、仍遵守原帧缓冲限制，服务重启后丢失。复杂分子的物理准备仍可能不收敛；此次改善失败后的交互和原比赛保留，没有降低准入门槛。复现、路由、资源回收与限制见 [COMPLETION.md](COMPLETION.md)。

## 2026-10-13 历史版本：统一简写界面文案

本轮只更新工坊文字、独立几何接口的消息组织及 C60 模板提示的显示。图标、绘图、构象算法、原子身份、撤销、赛车准备和公共路由保持原逻辑。任务自管的 8786 服务已重载，其他服务没有停止。

- 图标提示统一为动作名称及快捷键；二维/三维底栏只显示模式、选择状态和原子/键数量。
- 重建、补氢、UFF、排版、导入、保存、替换等成功反馈只报告结果。几何接口新增兼容字段 `message`，正常完成文字从 `warnings` 中分离。
- 删除“可撤销、可保存、真实氢、原生、坐标/ID 保留、未执行其它操作”等重复保证；帮助窗口改为操作指引。
- 价态、立体冲突、未收敛、键长异常、文件兼容性和结束当前比赛的实际影响继续提示。UFF 的多个片段说明改为仅在有多个未连接片段时显示。
- C60 原生放置指引精简为“请先将 C₆₀ 放在空白处，再用单键连接。”实际测试发现 Ketcher 默认信息区域隐藏该提示，现通过公开消息订阅将这一条指引显示到工坊状态栏，其它原生调试消息不转接。
- 新/旧 C60 完成说明兼容过滤，真实手性提示不被过滤。

| 验证 | 结果 | 证据 |
|---|---:|---|
| 精简文案、实际提示、正常结果、错误指引及旧工程 | **14/14** | `outputs/editor-copy-check/report.json`；`frontend/editor/tests/workshop-copy.cjs` |
| 连接、C60、共同撤销、文件往返及异步优化 | **21/21** | `outputs/editor-simplified-check/report.json` |
| 二维绘制、补氢、C60–CH3 及三维种子 | **14/14** | `outputs/editor-seed-check/browser-report.json` |
| 共同模型及适配器 | **92/92** | `outputs/editor-copy-unit.log`；含旧 C60 提示和真实手性提示同时存在的兼容检查 |
| 独立后端几何与文件 | **59/59** | 原有后端命令；文案匹配改用当前问题/解决办法，构象、ID 与数据不变的断言继续通过 |

浏览器共 **49 项**。三份报告的源码 SHA-256 均与当前源码一致；页面错误、外部请求和 `/api/load` 调用均为零，文案专项的 generation 为 **1 → 1**。汇总在 `outputs/editor-copy-summary.json`。真实 UFF 收敛成功使用实际服务；未收敛提示的显示顺序使用带真实计算结果的受控警告返回值验证，不声称该算例实际未收敛。本轮没有新的 LAMMPS 试跑验收。

Ketcher 原生扩展从已有固定源码构建为 `3.18.0+molecular-car.c60.3`，仅修改模板放置提示；模板数据、来源及编号映射仍为版本 2。169 个发行文件哈希核对通过，144 个依赖及许可记录完整保留；没有改压缩 bundle、增加依赖或 CDN。构建输出在 `outputs/editor-copy-ketcher-build.log`。JSmol 库与图标资产没有改动。

已查看 `outputs/editor-copy-check/concise-2d.png` 与窄屏截图；专项还实际下载 `saved.mcar` 并核对保存的完整原子记录。其它此前限制见历史记录及 README。

## 2026-10-13 历史版本：三键直线排版与侧栏操作细节

本轮只修改独立工坊、独立几何模块、相关测试及说明。任务自管的 8786 服务已重新加载；其他服务、赛道加载和公共启动配置没有变动。最终源码冻结后核对全部浏览器报告的源码 SHA-256，六组报告均匹配当前版本；汇总在 `outputs/editor-polish-summary.json`。下方旧版本的功能描述是历史记录，不能代表当前菜单行为。

| 验证 | 结果 | 证据 |
|---|---:|---|
| 本轮排版、文字图标、右键关闭、构象和共同撤销 | **11/11** | `tests/workshop-polish.cjs`；`outputs/editor-polish-check/report.json` |
| 独立构象操作、连接、C60、文件往返及旧优化丢弃 | **21/21** | `outputs/editor-simplified-check/report.json` |
| 主工坊原生手势、侧栏工具及共同历史 | **17/17** | `outputs/editor-toolbar-check/report.json` |
| 周期表替换、元素栏、菜单关闭与窄屏 | **14/14** | `outputs/editor-ketcher-style-check/report.json` |
| 氢标签、新增区域、C60–CH3 与构象保留 | **14/14** | `outputs/editor-seed-check/browser-report.json` |
| 独立研究页原生菜单与原子身份回归 | **19/19** | `outputs/editor-jsmol-native-check/report.json`；菜单仅在研究页开启 |
| 共同模型与适配器 | **92/92** | `outputs/editor-polish-unit.log` |
| 独立后端几何、文件与三键排版 | **59/59** | 原有 57 项及 `tests/test_editor_linear_layout.py` 2 项；见汇总中的执行命令 |

浏览器共 **96 项**，页面错误与外部请求均为零。主工坊脚本未调用 `/api/load`；本轮专项检查的赛道 generation 为 **1 → 1**。独立研究页通过静态拦截运行，不接触任何 API 或正在运行的 8765 服务。本轮没有新的 LAMMPS 试跑验收。

本轮实际结果：

1. 原问题是未消毒的编辑图没有识别杂化状态，RDKit 将普通炔碳角排成 120°。现在只在私有绘图副本上识别杂化，不做完整化学消毒，也不改变 canonical 图。2-丁炔在真实 Ketcher 页面中两端均为 **180°**；普通炔链、二炔、芳基炔及包含未完成价态的草稿有后端覆盖。
2. 二维排版、撤销、重做和再返回三维均精确保留所有旧原子 ID、三维坐标和键级。2D 标签、C60 球壳、新基团尺度、手动二面角、显式补氢和旧异步结果丢弃的既有回归同时通过。
3. 原子工具显示字母 **Atom**，有单独的文字宽度以避免裁切；+H 使用普通 Arial 粗体 H。右侧所有常用元素字号保留、字重为 700；原生 Ketcher 元素栏字重为 600，两边均为粗体视觉强调。右下复位与顶部适配视图统一四角框。
4. 绕键旋转单环箭头头尾有明显缺口，键穿过环；UFF 箭头缩小，与曲线有间距。桌面、390px 窄屏和放大图标截图已逐一查看：`polish-3d.png`、`polish-narrow.png`、`icons-enlarged.png`、`linear-alkyne-2d.png`，均位于 `outputs/editor-polish-check/`。
5. 主工坊设置 `nativeMenus:false`，关闭三维右键、原生角落和 logo 菜单入口，并去掉失效入口的绘制。ModelKit 仍保持原生编辑状态；真实替换、添加、连接、删除、移动、基团旋转、绕键旋转和共同撤销通过。二维 SUP 的右键展开/折叠保留。独立研究页默认开启原生菜单，19 项回归也通过。
6. 图标清单为 `.3`；47 个本地图标、内联模块及许可证哈希核对通过，其中 28 个官方 Ketcher SVG 与保留的官方源码逐字节一致。本轮未改 JSmol 的 LGPL 库文件、添加依赖或使用 CDN。

限制：炔链直线排版不意味着能强制任意环内三键满足 180°，环约束仍由成熟绘图工具处理。手机触摸、复杂 C60 多点外接等此前限制仍在；本轮窄屏验收使用桌面鼠标。初始构象构建仍不自动执行 UFF。

## 2026-10-13 历史版本：独立构象操作、合并键工具与精简界面

本轮使用任务自己管理的 **8786** 独立实例，在加载新后端模块后，以最终冻结源码完成实际 Windows Chrome 验收。其他任务服务未停止，模拟与加载语义未修改。

| 验证 | 结果 | 可复现入口及证据 |
|---|---:|---|
| 共同模型、Ketcher 显式映射、几何、周期表与提示 | **92/92** | `outputs/editor-simplified-unit.log`；新增 4 个操作提示/C60/交错片段编号检查 |
| 独立后端文件、补氢、初始构建与几何 | **57/57** | `outputs/editor-simplified-backend.log`；新增 `tests/test_editor_separate_operations.py` 6 项 |
| 新界面、合并键工具、操作分拆、文件与真实优化 | **21/21** | `frontend/editor/tests/workshop-simplified.cjs`；`outputs/editor-simplified-check/report.json` |
| 原生二维绘制、氢标签、C60 及新增区域 | **14/14** | `scripts/verify_editor_seed_browser.cjs`；`outputs/editor-seed-check/browser-report.json` |
| 原生三维手势、共同历史与视图往返 | **17/17** | `frontend/editor/tests/jsmol-workshop-toolbar.cjs --roundtrip`；`outputs/editor-toolbar-check/report.json` |
| 周期表、元素替换、原生选择与菜单同步 | **14/14** | `frontend/editor/tests/ketcher-style-workshop.cjs`；`outputs/editor-ketcher-style-check/report.json` |
| 原生 ModelKit 独立菜单、变换与身份回归 | **19/19** | `frontend/editor/tests/jsmol-native.cjs`；`outputs/editor-jsmol-native-check/report.json` |

共 **85 项实际浏览器检查**通过。五份报告的源码 SHA-256 均与最终文件一致，`sourcesChangedDuringRun=[]`、`pageErrors=[]`、`externalRequests=[]`。验收强制阻止 `/api/load`，新独立实例的比赛 generation 始终 **1 → 1**；本轮没有通过工坊送入赛道，也未进行 LAMMPS 试跑验收；独立服务启动时的默认初始化不计为验收成功。此前真实试跑结果保留在历史节，不计入本轮结果。

实际验收重点：

1. `+H` 图标为加号在前的同基线符号，原子工具始终为 Atom 轨道图标。适配视图为四角框；绕键旋转为化学键穿过单箭头椭圆环。UFF 的能量曲线与二维排版的平面环图标明确不同。3D 只显示 UFF，2D 只显示二维排版；390 px 与桌面均无水平溢出，截图已目视核对。
2. Ketcher `hiddenControls` 原生配置整组去掉 R 基团、反应与图形注释的下拉入口，修正只隐藏内部按钮而留下重叠小三角的问题；选择、键型、电荷、周期表、环模板及原生 C60 SUP 继续有效。
3. 常驻橙色提示条已移除。C60 加成完成后白色状态栏只报告当前操作；旧工程里错存为长期结构警告的那条说明过滤。真实未解决的价态/手性问题显示可点击警告计数，实际弹窗保留原因。新验收加载混合旧 C60 说明与手性问题，确认只出现 1 个真实警告，不丢失手性提示。
4. 一个化学键按钮分别执行：点击已有键改为所选单/双/三键、拖到已有原子连接、拖到空白新增当前元素并连接；普通原子点击不换元素或增加历史。真实用同一个工具完成单双三键，以及 N 双键延伸。每次仅一份历史，原有原子/键 ID、坐标与其余拓扑保持，撤销/重做精确恢复。
5. 满价且含手动 H 的连接被明确拒绝，项目与历史不变。C60 外接单键按明确加成降一条笼双键级并检查价态，连接和降级整体撤销。补氢后真实 65 原子、95 键，60 笼原子与所有旧坐标逐值保留；原生二维折叠和返回三维球壳有效。
6. 独立 `rebuild` 只写三维坐标：六碳苯重建后仍六原子、六键，原 ID 与二维排版保留，无显式 H、无 energy、`optimized:false`。独立 `hydrogenate` 在二维仍保持二维视图，加入六个真实 H；三维只添加 H，旧坐标保持。独立 UFF 从当前构象优化，固定所选原子精确不动，图与 ID 不变，临时 H 不写入工程。三项分别可撤销。
7. 后端新测试禁止 UFF 调用，坐标重建仍通过；二维补氢测试禁止嵌入/优化仍通过。UFF 测试禁止嵌入，仍从当前坐标降低能量。输入失败不会修改原项目。当前页面与原生菜单不再调用旧 `generate` 组合操作；它只为已有独立接口/历史测试保留兼容。
8. TPB 通过真实 JSmol 绕键手势改变构象，保留全部键长、拓扑和 ID；撤销/重做逐值恢复。再真实拖动二维原子并整理二维，全部三维坐标保持；返回三维仍为手工编辑后的构象。
9. 延迟真实 UFF HTTP 响应，在等待时继续原生三维移动，旧响应被丢弃，新设计保持。旧 prepare 的同类边界继续通过。真实手性分子整体旋转、3D MOL 下载、文件选择器重导入，保留元素、电荷、同位素、四面体标记及坐标（V2000 舍入范围内）。MOL 导入获得新工程 UUID，不宣称普通 MOL 保存永久编辑器 ID。
10. 验收另发现交错排列的断连片段会被 Indigo 重排，原有引用表因顺序不同拒绝绑定。KET 传输现在按**已知 canonical 边及原记录顺序**排列连通部分，并同步明确原子/键引用表；源项目顺序、ID 和坐标不变，不匹配输出图、不猜对称 C60 身份。实际含四个手动 H 与孤立碳的原生项目加载及后续连接验证通过，并有独立映射单元测试。

新截图、逐项完整项目快照与实际导出的 `chiral-edited.mol` 在 `outputs/editor-simplified-check/`，包括 `simplified-2d.png`、`simplified-3d.png`、`simplified-3d-narrow.png`。原生菜单扩展为 `molecular-car-modelkit-3`，保留 LGPL 原始源码、补丁和哈希；47 份本地图标为 28 份原始 Ketcher、1 份 C60 原创扩展、18 份工坊原创。无新增公共依赖或 CDN。

已知限制：三维键工具直接提供单/双/三键，芳香与立体标记仍通过 Ketcher 或属性面板编辑。重新生成是初始构象构建，不保证无原子碰撞或能量最小；C60 复用球壳来源，复杂多点外接尚未完整验收。JSmol 左上原生入口仍露出细颜色条。手机触摸与完整剪贴板连续链未进行本轮专项验收。未修改赛车动画、力场、赛道准备、共享路由或公共启动脚本。

## 2026-10-13 上一轮：Ketcher 图标布局、统一周期表与原生选择修复

本节保留上一轮界面统一的过程。同名报告已由顶部新一轮验收覆盖，当前源码哈希和结果以顶部为准。上一轮仅修改前端，无服务重启。

| 验证 | 结果 | 入口及证据 |
|---|---:|---|
| 新图标布局、统一周期表、原生选择与元素替换 | **14/14** | `frontend/editor/tests/ketcher-style-workshop.cjs`；`outputs/editor-ketcher-style-check/report.json` |
| 主页面三维工具、原生手势、共同历史及往返 | **17/17** | `frontend/editor/tests/jsmol-workshop-toolbar.cjs --roundtrip`；`outputs/editor-toolbar-check/report.json` |
| 原生二维绘制、氢标签、C60 球壳及新增结构构建 | **14/14** | `scripts/verify_editor_seed_browser.cjs`；`outputs/editor-seed-check/browser-report.json` |
| 独立原生 ModelKit 菜单、手势与身份回归 | **19/19** | `frontend/editor/tests/jsmol-native.cjs`；`outputs/editor-jsmol-native-check/report.json` |

共 64 项实际浏览器检查通过。四份报告的源码 SHA-256 均与最终文件一致，`sourcesChangedDuringRun=[]`、`pageErrors=[]`、`externalRequests=[]`。主页面脚本阻止所有 `/api/load`，比赛 generation 为 **2 → 2**。本轮没有执行 LAMMPS 试跑，既有试跑结果保留在历史节，不计入本轮统计。

实际核对：

1. 顶部打开、保存、撤销、重做及三维通用工具直接使用固定 Ketcher 3.18.0 的原始 SVG。左侧分组工具、右侧常用元素与周期表、白色工具组及青绿色激活状态与二维保持一致；三维专属动作使用同规格原创图标，仍调用 JSmol 原生 picking。42 份图标及本地 Apache LICENSE / NOTICE 来源与哈希已核对，主工坊不再加载 Lucide。
2. 原生“元素周期表…”与右侧周期表入口打开同一个 118 元素弹窗。实际选择 Si / Cl，在空白添加真实原子、点击原子替换元素，稳定 ID、两套坐标和其余键逐项保持；共用撤销 / 重做精确恢复完整项目。
3. 选中一个 H 后直接点击弹窗里的“替换当前选中原子为 Cl”，无需再点击画布；只增加一份历史，可撤销。三个原生菜单类别都不再显示 `??` 空元素槽，键级、连接、元素及移动模式的菜单选择反向同步常驻图标。
4. 验收发现并复现了 JSmol 的进入画布假点击：`JSmol.min.js` 的 mouseenter 发出 504 后紧接 502，`J/awtjs2d/Mouse.js` 的 502 在返回上次按下位置时无新按下动作也会调用 clicked。证据在 `outputs/editor-ketcher-style-selection-failure/enter-diagnostic.json`。`jsmol-host.js` 对本 applet 实例过滤未按下的 502 点击，只保留原生 released 恢复，实际 501/506/502 继续原生处理；没有改库文件或选择语义。完整元素替换序列之后，真实移回同一 H 不再改变原生或共同选择、修订、历史及项目；随后真实点击能正确选中。
5. 原生移动、分子 Alt 旋转、TPB 绕键旋转、C60 整组移动、删除原子 / 键、键级和拖动连接继续通过。三维创建并手工移动 Si 后，真实二维 / 三维往返保留原有 13 个原子的 ID、pos2、pos3，只补入 4 个 H；真实 prepare 响应 `optimized:false`，无 UFF 能量。
6. 二维氢标签、旧苯与新增环尺度、TPB 手工二面角、C60 60C/90 键与球壳距离、C60–CH3 两条构建路径及延迟结果丢弃继续通过。当前本机 C60–CH3 首次转换 248 ms、已有球壳新增甲基转换 260 ms；这是一轮组件加载后的测量，不代表所有设备性能。自动构建均无 UFF。
7. 390 px 页面没有横向溢出。左右工具组、测量面板和周期表可实际点击。桌面和窄屏截图已目视检查：`outputs/editor-ketcher-style-check/ketcher-style-desktop.png`、`ketcher-style-narrow.png`。

原始 SVG、版权声明、来源清单与原生菜单可读补丁保存在 `vendor/ketcher-icons/`、`jsmol-source/`，运行不依赖 CDN。没有修改公共依赖、启动脚本、比赛加载语义或模拟文件。

已知外观限制：JSmol 原生左上菜单入口仍保留，在自定义工具组背后露出约 10 px 的灰色 / 黄色细条；原生右键菜单和操作有效。本轮验收针对桌面鼠标及窄屏布局，未重新验收手机触摸或完整剪贴板连续链。三维专属动作在 Ketcher 中没有直接对应图标，其图标来源明确标为工坊原创。

## 2026-10-13 上一轮：无自动 UFF 的初始构建与常驻原生工具栏

本节保留上一轮实现与验证过程。同名主页面输出已由顶部新一轮验收更新，当前源码哈希、布局及结果以顶部记录为准。下面“局部三维准备、原生 ModelKit 与折叠 C₆₀”属于更早的自动 UFF 版本历史验收。

| 验证 | 结果 | 入口及证据 |
|---|---:|---|
| 共同模型、Ketcher 映射、周期表、构象确认与撤销 | 88/88 | `outputs/editor-seed-unit.log`；新增 4 个坐标确认边界测试 |
| Python 文件、SUP、初始构建与几何边界 | 51/51 | `outputs/editor-seed-backend.log`；`tests/test_editor_initial_construction.py` 11 项新增检查 |
| 实际主页面二维绘制及初始三维构建 | 14/14 | `scripts/verify_editor_seed_browser.cjs`；`outputs/editor-seed-check/browser-report.json` |
| 实际主页面三维图标、原生菜单及构象往返 | 17/17 | `frontend/editor/tests/jsmol-workshop-toolbar.cjs --roundtrip`；`outputs/editor-toolbar-check/report.json` |
| 原生 JSmol 宿主 / 原生菜单回归 | 18/18、19/19 | `outputs/editor-jsmol-host-check/`、`outputs/editor-jsmol-native-check/` |

上一轮运行时两份主页面报告的 source SHA-256 与当时源码逐项一致，`sourcesChangedDuringRun=[]`、`pageErrors=[]`、`externalRequests=[]`。构象脚本强制阻止 `/api/load`，比赛 generation 始终为 1；工具栏脚本仅 1 次真实 `/api/editor/geometry`，没有赛道加载。该轮没有重新进行 LAMMPS 试跑。

实际验证内容：

1. 在原生 Ketcher 中绘制苯、进入三维补入 6 个真实 H、返回二维核对每个氢文字位置；原生新增并连接苯的尺度合理，旧原子坐标和稳定 ID 保持。手工 TPB 二面角 71.23456789° 在新增结构和往返后保持。
2. 原生底栏插入折叠 C₆₀，真实图为 60C/90 键；球壳所有 1770 对内部距离保持。实际用 Ketcher 键工具从缩写连接点拖出一个甲基，连接按明确加成降笼内双键级，再进入三维，得到 65 个原子、95 条键、4 个真实 H。
3. C60–CH3 首次转换实测 **257 ms**，已有球壳后新增甲基实测 **227 ms**，均包含页面转换与同步。已有 60 个笼坐标逐值保持，外接 C–C 为 1.52 Å，4 个 C–H 为 1.10 Å，甲基位于球壳外。组件加载后的本机单次测量不代表所有设备性能。
4. 每个自动转换的真实 HTTP 响应都有 `construction.optimized:false`，没有 `energy` 字段。后端测试将 UFF/MMFF 调用设为错误，所有自动路径仍通过；C60–CH3 同时禁止 `EmbedMolecule`，仍成功。普通新苯/链/手性片段只嵌入新片段及一个 cap，不嵌入旧笼。无 UFF 参数的 Fe²⁺、He 也可生成初始构象。
5. C60–CH3 两条路径的撤销、重做恢复完整精确快照，往返不重复加氢。延迟真实 prepare 响应并继续原生二维绘制后，旧结果被丢弃，当前设计保持。
6. 9 个常驻图标实际执行选择、元素添加/替换、键级、拖动连接、删除、原子移动、分子/基团移动、绕键旋转及测量。周期表、精确键长面板、Alt 分子旋转和原生右键菜单反向激活都经过真实鼠标/键盘操作。原生 TPB 绕键旋转改变二面角，同时保留各键长度、拓扑与 UUID。
7. 原生三维添加一个碳并手工移动，再经二维返回三维：所有原有 13 个原子的 ID、pos2、pos3 精确保持，真实后端只添加 4 个氢，最终 17 原子。两待构建新原子只拖动其中一个时，另一仍保留 pending；同一步撤销恢复坐标与两个标记。未知构象确认 ID 会原子性拒绝。
8. 未释放的原生拖动可通过键盘切换工具安全取消，晚到的鼠标释放不会覆盖设计。390 px 窄屏图标栏横排可滚动，实际宽 372 px、高 47 px；提示框宽 271.375 px，右边界 385.375 px，不溢出页面。截图已目视核对。

源码新增独立 `molecular_car/editor/construction.py`，未改共享模拟、加载语义、赛车画面或公共依赖。几何整理仍为显式 UFF；已有元素、电荷或键级变更也不自动移动旧原子。几何构建不是能量最小化，多个原有固定点约束的新桥接会明确拒绝。三维原生键图标支持单、双、三键，拖动连接形成单键；芳香键、楔线和立体字段仍使用 Ketcher 或属性面板。

当前图片和完整分子快照在 `outputs/editor-seed-check/` 与 `outputs/editor-toolbar-check/`，包括 `c60-methyl-first.png`、`c60-methyl-appended.png`、`toolbar-desktop.png`、`toolbar-narrow.png`。全部功能运行时继续使用本地开源组件，无新增 CDN 或依赖。

## 2026-10-13：局部三维准备、原生 ModelKit 与折叠 C₆₀

本节是先前自动局部 UFF 版本的历史验收；当前版本以页面顶部“独立构象操作、合并键工具与精简界面”为准。

| 验证 | 结果 | 入口及证据 |
|---|---:|---|
| 模型、Ketcher 映射、几何、周期表与新增准备事务 | 84/84 | `outputs/editor-preparation-unit.log`；`tests/editor_preparation_model.test.mjs` 11 项新增边界 |
| Python 几何、文件、SUP 与固定点保护 | 40/40 | service 34、import depiction 3、SUP 3；见 `tests/test_editor_service.py` 等 |
| 实际主页面自动准备 | 10/10 | `scripts/verify_editor_prepare_browser.cjs`；`outputs/editor-prepare-check/browser-report.json` |
| 实际主页面原生菜单、文件与 LAMMPS | 18/18 | `scripts/verify_editor_native_workshop.cjs`；`outputs/editor-native-workshop-check/browser-report.json` |
| 原生 JSmol ModelKit 菜单、鼠标与共同历史 | 19/19 | `frontend/editor/tests/jsmol-native.cjs`；`outputs/editor-jsmol-native-check/report.json` |
| JSmol 宿主回归 / 原始运行库能力 | 18/18、15/15 | `frontend/editor/tests/jsmol-host.cjs`、`jsmol-live.cjs`；对应 `outputs/editor-jsmol-host-check/`、`editor-jsmol-check/` |
| Ketcher 原生 C60 v2 / 集成映射与连接 | 13/13、13/13 | `ketcher-source/probe-native-c60.cjs`、`probe-native-c60-v2.cjs`；`outputs/ketcher-native-c60/`、`ketcher-native-c60-v2/` |

两个主页面脚本均记录实际源码 SHA-256，最终运行 `sourcesChangedDuringRun=[]`、`pageErrors=[]`、`externalRequests=[]`。准备脚本始终 `/api/load=0`，比赛 generation 未变化。原生工坊脚本在所有编辑、文件和几何阶段同样 `/api/load=0`，最后只有两次可见“送入赛道”点击加载。

实际验证内容：

1. 原生二维苯环首次进入三维自动补入 6 个真实 H。返回二维后，每个 H 文字位于其真实原子坐标，最大标签偏移小于 0.00003 px；重复往返不产生额外优化或重复氢。
2. 原生二维添加并连接苯环，成键只替代 1 个编辑器生成的 H。新环 C–C 为 1.39664–1.40795 Å，所有仍存在的旧原子坐标，包括连接锚点，逐值保持。手画 / 导入 H 不被自动删除。
3. TPB 手工二面角 71.23456789° 和全部旧坐标，在新增原生环、自动准备和视图往返后保留。模型专项还验证自动手性失效提示不解冻旧锚点；带电碳及 bracket-H 数量使用统一价态规则。
4. 原生底栏 C₆₀ 模板是真实 SUP，默认化学式下标显示，保持 60C/90 键。明确模板来源和 Pool→UUID 顺序带入球壳，初始半径 3.5187865 Å，1770 个原子对距离保持；首次三维按短边适配完整显示。
5. 原生模板拒绝融合已有原子 / 键。原生单键连接 C60 执行明确的笼双键加成降级，全部碳价态≤4；交易可一次撤销。三维普通拖动移动笼及其直接末端 H，已连接车身不随行。新连接后真实补氢与局部几何准备保持球壳。
6. 原生 ModelKit 菜单的元素替换、H、几何整理、共同 undo/redo、保存及导出在主页面实际点击验证。真原子拖动、Alt 分子旋转和绕键旋转有独立真实鼠标验收。TPB 原生绕键从 -54.4565° 改为约 -33°，原键长和相机保持。选择或双击不能遗留隐藏的键编辑模式。
7. 原生菜单动作在 50 ms 等待下生效；修正仅对工坊 bridge 取消上游 500 ms 菜单延迟。普通移动更新坐标差分后保留 live 模型，无每次拖动重载。
8. 实际下载 `.mcar` 并通过文件选择器重开，逐项精确恢复 C60 稳定 ID、两套坐标、SUP 展开状态、模板来源和 H 生成标记。实际下载 3D V3000 MOL 并重导入，保留 SUP 化学式、完整原子 / 键、连接点和构象，误差在六位小数舍入范围内。没有工坊模板标记的导入 SUP 仍按文件声明成员整体移动，不猜模板身份。
9. 延迟真实 prepare HTTP 结果，同时在原生二维画布继续编辑，旧结果被丢弃。撤销 / 重做、几何预览与取消的已有回归继续通过。
10. 最终连接的 C60–苯加成样例，恢复原生工程精确坐标后通过原生菜单显式整理，再点击送入。HTTP 200 接受 generation 4，真实 LAMMPS 20260930 准备 72 个原子并进入同一 generation 的 `paused`；苯样例 generation 5、12 原子也进入 `paused`。generation 随重跑递增，以 JSON 为准。

此前另一次文件往返构象的 C60 送入被实际引擎拒绝：联合残余力 0.0271224 kcal/mol/Å，要求≤0.0001。页面按相同 generation 显示具体错误，而没有将 HTTP 200 当成成功。该记录提醒几何整理或化学准入通过不保证所有输入构象都能通过联合松弛；并未修改引擎阈值或模拟准备流程。主模拟 agent 随后为准备阶段接入了收敛补救，保持原阈值；最后 18 项联合验收的 C60 样例确已进入 `paused`。该成功来自当前联合运行状态，不能单独归因于编辑器或原生文件精度；不把较早失败记录算成成功。

固定旧端点的新键可保留不合理长度。服务新增实际 10 Å C–C 测试，坐标及原键保持，明确提示固定约束下的收敛不代表键长合理；无 pending、无关区域、正常键长和可动端点四种反例不产生该提示。

当前截图及文件包括 `outputs/editor-prepare-check/` 的氢标签、双苯、TPB、球壳；`outputs/editor-native-workshop-check/` 的 `main-c60-collapsed.png`、`main-c60-sphere.png`、`main-c60-connected.png`、`connected-fullerene-lammps.png`、`benzene-lammps-paused.png`、MOL 与 `.mcar`。每个主页面用例保存完整 canonical 项目和原生视图证据。

本轮仍未完整验收 C60 复杂多点外接、手机触摸和完整剪贴板连续链。图5式空间遮挡球形示意没有替代二维逐原子编辑：默认使用原生 C₆₀ 缩写，展开仍是完整化学图。查询、反应、聚合物、非 SUP S-group、增强 AND/OR 立体组等边界见 README。


## 2026-10-12 历史：原生 Ketcher + JSmol 主工坊

主二维使用本地 Ketcher 3.18.0 的源码构建（含原生 C₆₀ 模板），主三维使用本地 JSmol 16.3.31 ModelKit。本段对应新主页面；下方旧三维工具的 26 项记录不混入本轮统计。

| 验证组 | 结果 | 证据与可复现入口 |
|---|---:|---|
| 成熟工具完整主页面、真实文件和 LAMMPS | **20/20** | `outputs/editor-mature-check/browser-report.json`；`scripts/verify_editor_mature_browser.cjs` |
| JSmol 宿主实际鼠标构造、意图过滤与共同历史 | 18/18 | `outputs/editor-jsmol-host-check/report.json`；`frontend/editor/tests/jsmol-host.cjs` |
| JSmol 原生运行时及显式身份映射 | 15/15 | `outputs/editor-jsmol-check/report.json`；`frontend/editor/tests/jsmol-live.cjs` |
| 原生 C₆₀ 源码构建与真实模板工具 | 13/13 | `outputs/ketcher-native-c60/browser-report.json`；`frontend/editor/ketcher-source/probe-native-c60.cjs` |

完整主页面脚本使用 Windows Chrome 155.0.8059.39，仅访问本任务独立的 **8778** 实例，不启停服务，也不操作 8765 / 8770。最终报告为 `passed=20`、`failed=0`、`pageErrors=[]`、`externalRequests=[]`、`sourcesChangedDuringRun=[]`；每项保存 canonical 分子、修订、两工具编号映射和 JSmol 实际坐标，源码 SHA-256 一并归档。只有最后两个真实“送入赛道”点击调用 `/api/load`。

主页面已实际核对：

1. 原生键工具从空白拖出自由链，一次手势仅一条共同历史；苯模板实际放置得到 6C / 6 键。点击 3D 自动生成真实构象，canonical 和两工具均为 6C + 6H、12 键，原碳 ID 与二维排版保持，新 H 有唯一稳定 ID 和真实连接键。
2. JSmol 的真实鼠标移动改变 `pos3`，相机、拓扑、`pos2` 不变，统一撤销 / 重做恢复精确快照。真实周期表选择元素后点击原子，将 C 替换为 N；所有原有 H、ID、键和两套坐标保留。空白点击添加 O / C、两原子连接、点击已有键改单双键、删除键和原子均经过真实鼠标验证，六次共同撤销恢复整个编辑序列。
3. 原生二维原子拖动只改变 `pos2`；返回三维和连续 2D / 3D 往返不再改变已补氢构象或添加重复 H。四苯基苯实际 JSmol 二面角预览不提交 canonical，显示坐标与预览逐点核对；由 −54.46° 改为 −19.46° 后，应用、两工具往返及共同撤销 / 重做保留构象。
4. 原生 `template-8` 的鼠标放置产生真实 60C / 90 键、30 个双键；每个碳有三个邻接键、键级和为 4，UUID 唯一。主页面中旧赛车片段库和旧自定义片段入口已不存在。原生模板的连出键、MOL 往返和工具历史另在 13 项专项核对。
5. 三维手性 MOL 输入构象保留，并通过两工具切换核对手性字段与坐标。实际下载 3D MOL 后用文件选择器重导入，逐项核对元素、形式电荷、同位素、手性、显式编号映射后的键级与立体字段、3D 标头与三维坐标；本次最大坐标误差为 0 Å。实际 `.mcar` 下载 / 重开精确保留工程 ID、拓扑、两套坐标和两个成熟工具状态。
6. 延迟真实 UFF 返回，在等待期间通过 JSmol 点击选择原子；即使分子修订号未变，旧结果仍因交互版本变化被丢弃。另有明确标记的 **模拟 HTTP 400 补氢失败 UI 用例**：已有三维从 2D 返回 3D，显示具体失败原因，原分子与坐标不变，并可实际保存工程。此项不宣称触发了真实引擎补氢错误。
7. 绘制、几何、文件和切换阶段 `/api/load=0`，比赛 generation 保持 4。断连结构明确送入后 HTTP 200 返回 generation 5，随后同一 generation 进入 `error`，页面显示“赛车必须是一个相互连接的分子；请移除盐、溶剂或未连接的片段”。有效苯分子再次明确送入后 generation 6 经真实 LAMMPS 20260930 准备 12 个原子，进入 `paused` 才显示成功。复跑时 generation 会递增，以 JSON 为准。

本轮实测还发现并修复了 JSmol 的 Java 运行时导出全局 `Error` 后导致界面只显示 “Error” 的问题。应用模块分别保留原生错误构造器，没有修改 J2S 全局运行时；断连失败和补氢失败均严格要求显示具体原因，最终用例没有注入产品修复或模拟比赛成功。

截图和实际下载文件在 `outputs/editor-mature-check/`：`benzene-native-2d.png`、`benzene-jsmol-hydrogens.png`、`jsmol-real-atom-drag.png`、`jsmol-topology-edit.png`、`tpb-jsmol-torsion-preview.png`、`native-c60-template.png`、`chiral-jsmol.png`、`hydrogenation-failure-usable-3d.png`、`disconnected-track-error.png`、`mature-lammps-paused.png`，以及 MOL / `.mcar` 往返样例。

```powershell
node scripts/verify_editor_mature_browser.cjs http://127.0.0.1:8778
node frontend/editor/tests/jsmol-host.cjs
node frontend/editor/tests/jsmol-live.cjs
```

本轮主页面范围为桌面鼠标；没有重新验收手机触摸、完整 Ctrl+C → Ctrl+V 连续操作或原生工具所有化学类型。查询、反应、S-groups、增强 AND / OR 立体组等支持边界见 [README.md](README.md)。旧记录中的 Three.js 手柄、旧片段连接与旧剪贴板通过结果不能作为当前 JSmol 主工坊的验收。

## 2026-10-12 历史：Ketcher 与旧三维主画布、独立 JSmol 实测

本段是主三维尚使用旧 Three.js 构造工具时的记录，JSmol 当时只在独立入口验证。26 项主页面、18 项 Python 以及本段 JSmol 独立探针均为上一版结果；新主工坊结果单独记录，不混合统计。

| 验证组 | 当前结果 | 证据与入口 |
|---|---:|---|
| 统一模型 / Ketcher 适配 / 几何 / 周期表 | 66/66 | `outputs/editor-ketcher-check/unit-summary.json`；adapter 21、model 24、geometry 15、periodic 6 |
| 独立 Python 文件 / 几何 / 准入与导入描画 | 18/18 | `tests/test_editor_service.py` 15 项、`tests/test_editor_import_depiction.py` 3 项 |
| 完整主页面与真实 LAMMPS | 26/26 | `outputs/editor-ketcher-check/browser-report.json`；`scripts/verify_editor_ketcher_browser.cjs` |
| Ketcher 官方发行实际部署 | 7/7 | `outputs/ketcher-deploy/browser-probe.json`；`frontend/editor/research/probe-ketcher-browser.cjs` |
| Ketcher 事件 / 历史 / 生命周期专项 | 16/16 | `outputs/ketcher-deploy/host-probe.json`，同目录 `host-probe.cjs` 可复现 |
| 同位素显示专项 | 已确认天然丰度与显式同位素 | `outputs/ketcher-deploy/labels-probe.json`；省略 isotope 0 消除错误 0C，原坐标 / 拓扑不变 |
| 独立 JSmol ModelKit 真实编辑 | 15 项断言通过 | `outputs/editor-jsmol-check/report.json`；`frontend/editor/tests/jsmol-live.cjs`，含可复现限制断言 |

完整验收只对本任务自己的 **8778** 实例操作。脚本拒绝其他端口，不启停服务。`pageErrors=[]`、`externalRequests=[]`；绘图、几何、导入与保存阶段 `/api/load` 次数为 0，只有两个明确“送入赛道”点击调用加载。后一个成功批次经真实 LAMMPS 20260930 准备 12 个含氢原子并进入 `paused`；断连设计先得到 HTTP 200，再由同一批次的 `error` 给出“赛车必须是一个相互连接的分子”原因。具体 generation 随复跑递增，以归档 JSON 为准。

已实际核对：

1. 原生周期表用真实鼠标添加 Si，再点击同一原子替换 N；稳定 UUID、pos2 / pos3 与撤销重做保持。不依赖元素热键。三维中文周期表的 118 元素同样能真实添加。
2. 原生键工具拖出自由链，一手势一个历史；原生苯模板点击插入真实 6 / 6 分子图，生成三维后用鼠标拖动平移轴和旋转环改变真实坐标。
3. 原生二维选择拖动只改 pos2；四苯基苯二面角预览 / 应用 / 取消、固定侧交换及撤销重做经过 Ketcher / 三维切换，pos3 按快照保留。C60 连接后的 62 个 UUID、92 条键、三维坐标和端口，在原生全选 / 拖动中逐项保持。
4. 片段按显式端口连接，C60 的 60 原子 / 90 键与 60 声明位点真实可核对；饱和位点拒绝且模型不变，明确降低相邻双键后成功连接炔基。
5. 三维手性文件输入构象保持，为二维生成正确楔线；实际二维拖动、3D MOL 下载 / 文件选择器重导入和 `.mcar` 下载 / 重开保留结构、手性和坐标。多记录 SDF 明确拒绝。3D MOL 使用 3D 头，坐标舍入误差在 0.000051 Å 内。
6. 延迟真实 UFF HTTP 响应，同时继续修改模型，旧结果被丢弃。专项还验证原生键盘 change 当场更新交互版本，程序推送期间键盘与指针锁定，避免未提交编辑被异步覆盖。
7. 截图检查发现并修复两个实际显示问题：隐藏 iframe 的零尺寸不能执行 fit，返回二维后再适配可见区域；默认同位素 0 不传 KET，天然碳恢复骨架显示，13C 等明确同位素保持。普通视图切换保留现有 zoom / viewBox，复位不改变坐标。

当前截图 / 下载文件在 `outputs/editor-ketcher-check/`，包括 `workshop-ketcher-2d.png`、`native-periodic-table.png`、`tpb-torsion-preview.png`、`chiral-ketcher.png`、`c60-connected-3d.png`、`race-lammps-paused.png` 与原生 / MOL 样例。JSmol 另有三张实际编辑截图。

JSmol 的 15 项包含原子拖动、35° 片段旋转、增删原子 / 键、键级 / 电荷、显式 ID token、C60 逆序映射、手性 MOL 和坐标安全合并。但实测元素修改可能自动删显式 H / 改键长，原生历史可能丢 ID；因此独立入口与适配器没有被宣称为完整主页面双向同步。见 [详细边界](research/jsmol-live.md)。

剪贴板粘贴已在新隔离浏览器中实际验证：将 Ketcher 自己生成的苯环 MOL 写入剪贴板，真实 Ctrl+V 加鼠标移动 / 放置使原生与工坊从 6 增至 12 原子；新 UUID 唯一、原原子 / pos3 / 键保持，撤销一次恢复。本项计入 26 项。完整 Ctrl+C → Ctrl+V 连续快捷键链尚未稳定验收，不能由适配器复制测试推断该链已通过；具体限制写入 JSON。桌面为本轮主要验收范围，没有重做手机触摸专项。

```powershell
node --disable-warning=MODULE_TYPELESS_PACKAGE_JSON --test tests/editor_ketcher_adapter.test.mjs tests/editor_model.test.mjs tests/editor_ux_model.test.mjs tests/editor_periodic_table.test.mjs
.\.venv\Scripts\python.exe -E -B -m unittest tests.test_editor_service tests.test_editor_import_depiction
node scripts/verify_editor_ketcher_browser.cjs http://127.0.0.1:8778
node frontend/editor/tests/jsmol-live.cjs
```

## 2026-10-11 历史记录

日期：2026-10-11（Asia/Tokyo）。运行环境为项目 CPython 3.11.5、RDKit 2025.9.3、NumPy 1.26.4、本地 Three.js 0.180.0、Windows Chrome、LAMMPS 20260930。

## 本次交互更新

本轮为图标界面、ChemDraw/GaussView 风格操作与完整周期表更新，已重新执行模型、Python 服务、真实浏览器与赛道集成验证。下方首版记录保留为历史，不混入本次统计。

| 验证组 | 结果 | 证据 / 可复现入口 |
|---|---:|---|
| 模型与编辑数据 | 45/45 | `outputs/editor-ux-check/model-report.txt`；`tests/editor_model.test.mjs` 24 项、`tests/editor_ux_model.test.mjs` 15 项、`tests/editor_periodic_table.test.mjs` 6 项 |
| 独立 Python 文件/几何/化学服务 | 15/15 | `outputs/editor-ux-check/backend-report.txt`；`python -E -B -m unittest tests.test_editor_service`，不启动或修改比赛控制器 |
| 独立视图真实浏览器操作 | 21/21 | `outputs/editor-ux-check/views-report.json`；`node frontend/editor/tests/views-browser.cjs`，涵盖框选/套索、相机与三维手柄 |
| 新版完整页面与实际赛道 | 25/25 | `outputs/editor-ux-check/browser-report.json`；`scripts/verify_editor_ux_browser.cjs` |
| 周期表额外隔离 Edge DOM 操作 | 1/1 | `tests/editor_periodic_table.test.mjs` 的可选浏览器测试；`research/reaxys-interaction.md` |

本轮 Chrome 报告 `pageErrors=[]`、`externalRequests=[]`、`loadCalls=2`。两次 `/api/load` 都由明确“送入赛道”触发，普通绘制、几何与文件操作维持原比赛 generation。完整页面验收只在自己的独立 8777 实例运行，未操作 8765 现有服务。真实 LAMMPS 准备后检查同一次 generation 的 `paused`，没有把 HTTP 接受当作准备完成。

已验证的行为：

1. **周期表与替换**：118 元素全部启用；中文、原子序数搜索与键盘选择能真实添加/替换原子。周期表和情境工具条的明确替换保留原子 ID、二维/三维坐标，撤销重做恢复元素。当前模拟不支持 Og 时，编辑仍保留，准入原因清楚显示。
2. **二维自由绘图**：实际指针绘制六条键、闭合苯环；悬停热键改交替双键，一次手势对应一次撤销。普通选择工具拖动只改 pos2，已有 pos3 保留。
3. **真实三维与四苯基苯几何**：实际指针拖动平移轴、旋转环改变真实 pos3，并可一次撤销/重做，拓扑与两套坐标按预期恢复。加载 52 原子示例后当前二面角预填为 −54.46°；改至 −19.46° 的预览没有提交项目，应用后测值一致。交换固定侧、A–D 标签、应用一次历史与撤销重做均验证；取消恢复原显示，2D/3D 往返逐项保留构象。
4. **片段与 C60**：选库先预览不修改图，炔基按位点连接且原有原子坐标保留，撤销/重做一次恢复。自定义片段通过页面勾选唯一连接原子并保存 6 个真实原子和该端口。C60 保持真实 60 原子/90 键和整体操作默认值，编号默认简洁；取消片段插入没有留下编辑。
5. **文件与异步**：实际通过重设计后的菜单下载 `chiral-edited.mol` 与 `chiral-native.mcar`，再通过文件输入导入。三维 MOL 头、结构、手性与坐标往返可核对，本次坐标最大误差为 0；原生工程保留完整图。多记录 SDF 明确拒绝。连续编辑中延迟真实几何 HTTP 响应，确认较旧结果没有覆盖新快照。
6. **赛道与布局**：普通编辑、生成三维、文件与检查未更改比赛 generation。显式提交不连通结构后 HTTP 成功，但本次归档报告的 generation 7 异步进入 error，页面显示“赛车必须是一个相互连接的分子”；随后明确载入苯，generation 8 为 paused，真实 LAMMPS 20260930 准备 12 个含氢原子。桌面及 1024/736/430/320px 布局、周期表检查没有页面横向溢出。

另有三项实际浏览器回归通过：修改已有键后清理绘键旧起点，后续点击原子不产生意外连接；修改几何预览后保存被拦截，清楚提示先应用或取消；周期表 Esc 只关闭弹窗，不清除绘图模式、原子选择或设计。保存、导出、送入赛道均使用同一预览确认检查；本项浏览器用例实际触发并核对的是保存拦截。

本轮报告和截图位于 `outputs/editor-ux-check/`，截图为 `workshop-2d.png`、`workshop-3d.png`、`periodic-table.png`；实际下载文件为 `chiral-edited.mol`、`chiral-native.mcar`。完整页面报告已归档 25 项通过结果，原始首版记录在下方单独保留。

```powershell
node --test tests/editor_model.test.mjs tests/editor_ux_model.test.mjs tests/editor_periodic_table.test.mjs
.\.venv\Scripts\python.exe -E -B -m unittest tests.test_editor_service
node frontend/editor/tests/views-browser.cjs
node scripts/verify_editor_ux_browser.cjs http://127.0.0.1:8777
```

浏览器验收使用自己的独立测试实例，脚本会实际加载赛道。组件/模型测试不访问比赛 API；几何旧结果测试只延迟真实响应，不替换真实 LAMMPS 成功结果。功能限制见 [README.md](README.md) 与 [UX-CHANGES.md](UX-CHANGES.md)。

旧入口 `scripts/verify_editor_browser.cjs` 已兼容转入新版浏览器验收，并实际执行通过；新增运行推荐使用上方 `verify_editor_ux_browser.cjs` 名称。

## 首版交付记录（保留）

### 结果

| 验证组 | 通过 | 证据文件（项目根目录下） |
|---|---:|---|
| 分子图、历史、几何与真实 RDKit 往返 | 24/24 | `outputs/editor-check/model-report.txt` |
| 独立文件/几何/化学服务 | 15/15 | `outputs/editor-check/backend-report.json` |
| 实际浏览器操作和模拟接口 | 16/16 | `outputs/editor-check/browser-report.json` |
| Ketcher 原函数及 JSmol 固定源码技术探针 | 10/10 | `frontend/editor/research/probe-results.json` |

模型用例包含实际四苯基苯与 C60 库，后端用例运行 RDKit；浏览器用例加载本地页面，执行按钮、文件导入导出和真实鼠标拖动。测试没有用假的比赛响应替代 LAMMPS。仅“优化结果过时”用例将真实 HTTP 返回延迟 700 ms，以确定性复现连续编辑。

## 对照验收要求

1. **自由二维绘图与真实三维编辑**：依次点击六个原子、闭合环、改交替双键，生成三维成功。另以片段苯验证 Z 方向平移 0.7 Å、撤销/重做。实际三维鼠标拖动改变 pos3，pos2 和 bonds 逐项不变；实际二维拖动仅改变 pos2。
2. **四苯基苯构象保留**：导入原例的 52 原子/56 键，通过可分离中间键将二面角设为 +73°。切换到二维、整理排版、回到三维后 pos3 逐点完全相同。后端与模型另外交叉核对二面角数值。
3. **真实片段连接**：C60 60 碳/90 键、度数 3/价态 4；饱和位点直接接炔基被拒绝且未改变图。明确将接点相邻双键降为单键后连接成功，得到 62 原子/92 键。模型测试确认默认整体移动。自定义苯环通过页面保存并重新插入真实原子图。
4. **文件往返**：页面导出三维四苯基苯 MOL，检查头部 3D，重新导入图规模和坐标可核对；实际最大坐标分量误差约 0.00004984 Å。页面手性乳酸 MOL 往返保持标签。后端覆盖 R/S、E/Z、同位素、形式电荷、SDF，以及多记录 SDF 的明确拒绝。原生 C60 工程经页面下载再打开，稳定 ID、拓扑、坐标与片段逐项精确相同。
5. **历史与异步**：移动与完整图变换可撤销/重做。实际优化返回被延迟期间修改工程名称，确认旧快照被丢弃。模型用例还覆盖撤销至旧内容时修订号仍单调递增、结构删除后的片段/立体引用恢复。拖动期间发生另一次修订会取消旧拖动。
6. **实际赛道接入**：在我启动的独立服务中，普通编辑/导入/导出/优化/检查均没有改变比赛 generation。未连通双苯环通过 HTTP 接收后异步进入 error，页面显示连接问题而不是成功。随后苯通过显式“送入赛道”，同一 generation 到 paused，真实 LAMMPS 20260930 载入 12 个含氢原子；没有自动开始比赛。

## 额外发现与修复

- 二维虚楔最初可能未传给三维嵌入：已修复为显式解析二维立体意图，回归验证实楔/虚楔分别保持相反且正确的 R/S。
- CW/CCW 依赖邻接序：影响手性中心的拓扑修改现在清除旧指定并显示提示，要求重新指定；避免删除重连同一条键后无提示翻转。
- 部分 E/Z 结构保存为片段会丢立体参考：现在保存前拒绝缺失参考的选择，并在入库前验证模板。
- SVG 隐藏属性必须显式设置：已修复三维视图误叠加二维图，浏览器断言确认三维时 SVG 不可见。
- Windows 同端口双重绑定会导致异常：经主 agent 协调加独占监听；已验证重复启动明确返回 WinError 10048，正常服务仍可访问。

截图：`outputs/editor-check/workshop.png`、`c60-workshop.png`。实际导出文件：`edited-tetraphenylbenzene.mol`、`chiral-roundtrip.mol`、`c60-native.mcar`。这些输出是可复查证据，不是运行依赖。

本记录不表示完整生产验证；GaussView/Ketcher/JSmol 完整 GUI 集成、环内自动构象优化、轴手性、大规模分子性能和跨平台 LAMMPS 尚未验证。完整边界见 README 的“已知限制”。
