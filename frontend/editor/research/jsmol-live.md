# JSmol ModelKit：本地浏览器实测与适配边界

验证日期：2026-10-12。本文件记录初期 **原生 runtime** 的能力与限制。后续主工坊已接入 JSmol，使用意图过滤、显式身份映射和统一历史解决下述副作用；最新集成验收见 [VALIDATION.md](../VALIDATION.md)。

这次使用实际 JSmol 16.3.31 浏览器 runtime，没有用源代码推断代替运行。独立页面能加载四苯基苯、C60、手性乳酸，拖动原子、旋转绕键片段、增删和修改原子/键，并导出 MOL。独立页故意保留未过滤的原生行为用于研究；主工坊由 `jsmol-host.js` 和增强适配器管理实际编辑、双向同步和撤销。

## 打开与操作

现有项目服务运行时，打开：

`http://127.0.0.1:8765/editor/research/jsmol-live.html`

页面独立于工坊和赛道，不调用 `/api/load` 或其它比赛 API。代码和资源全部在本地。无需 Java 插件或 CDN。第一次初始化需要加载 JSmol 脚本。

1. 点击“四苯基苯”，选择“移动原子”并启用，直接拖动一个球。拖动改变分子坐标；选择“相机模式”后才恢复普通视角操作。
2. 保持四苯基苯，旋转轴输入 `0,1`，点击“轴后片段旋转 35°”。这会旋转该键后相连的苯基；页面还提供原生 ModelKit “绕键旋转”模式。
3. “改为氧”“删除原子”“单键”“双键”“删除键”作用于点击的实际原子/键。“添加独立碳原子”在 `{12,12,5}` Å 添加一个不带键的碳；可切换相机模式检查。
4. “读取当前结构”显示真实原子/键数量和对象数据。“导出 JSmol MOL”导出工具内部的当前结构。这个 MOL 是诊断输出，不是工坊原生项目保存格式。

原子索引仅用于这个独立技术验证页中定位操作目标；它不是玩家设计的永久 ID。独立页面没有自动回写工坊，也没有项目撤销历史。

## 实际验收

可重跑命令（项目根目录）：

```powershell
node frontend/editor/tests/jsmol-live.cjs
```

测试使用本机 Chrome、Playwright 和项目现有 `.venv` 中的 RDKit，不增加公共依赖。默认 Playwright 路径是此开发机已有的 Codex runtime；其它机器可设置 `PLAYWRIGHT_MODULE`、`CHROME_PATH` 和 `EDITOR_PYTHON` 到现有安装。脚本启动自己的无头浏览器，以静态请求拦截读取项目文件，禁止外网和 `/api/` 请求，**不启动、不停止、不占用项目服务器**。

结果写入 `outputs/editor-jsmol-check/report.json`；屏幕截图为 `tpb-torsion-live.png`、`c60-live.png` 和 `chiral-live.png`。最近一次 **15 项断言通过**，其中包含对已知限制的复现断言，不表示这些限制已经消失。

| 实测项目 | 结果 |
|---|---|
| runtime 与加载 | 运行时 `_version=1603031`；四苯基苯 52 原子、56 键 |
| 真鼠标拖动 | 鼠标拖动一个 H 原子，位移约 0.842658 Å；相机方向及键图不变 |
| 真实构象编辑 | 四苯基苯绕 0–1 键旋转相连的 11 原子片段；二面角 −54.4565° → −19.4565°，变化恰为 35°；其余原子不变 |
| 工坊历史协作 | 坐标结果只更新 `pos3`；`pos2`、键/立体字段和片段信息保留；通过 `WorkshopModel.applyAsync(baseRevision,...)` 后撤销/重做成立，旧版本再次提交被拒绝 |
| 新增原子 | 在距旧原子仅 0.1 Å 处添加 C，原子 +1，键数不变，不自动补 H |
| 元素与电荷 | C → O、形式电荷 +1 均实际生效，自定义 ID token 保留；原子名称被重置，证明名称不适合作永久身份 |
| 键与删除 | 新建实际键，改单/双/三键，删键、删原子，数量和键级均核对 |
| 中间删除后新增 | 幸存原子的自定义 ID 全部保留；无需按新索引、距离或图匹配恢复身份 |
| JSmol 原生历史限制 | 删除并新增后 `save state` / `restore state` 不能可靠恢复自定义属性；适配器拒绝原生历史命令，需要工坊自己的历史 |
| 氢原子副作用 | 关闭 `addHydrogens` 后，某些元素修改仍删除已存在的显式 H；完整变更提案能看到删除；坐标专用提交拒绝它 |
| 键长副作用 | 修改一个单键相连原子的元素时，JSmol 会把该键长改为 1.5 Å；已记录实际位移 |
| C60 | 60 原子、90 键、每个碳 3 个邻接；将 MOL 原子行完全逆序后，靠显式行映射仍对应正确永久 ID 和原坐标 |
| 手性 MOL 往返 | 手性乳酸整体移动 `{0.3,-0.2,0.4}` Å，MOL 导出标记 3D，RDKit 核对 S 构型保留；MOL 四位小数输出最大误差约 0.0001 Å |
| 本地部署 | 0 外部请求、0 比赛 API 请求、0 浏览器脚本错误；仅初始化时一次无影响的本地目录探测 404 |

## 可复用适配器

`frontend/editor/jsmol-adapter.js` 导出 `JSmolAdapter`，接口是：

```javascript
const adapter = new JSmolAdapter(applet);
adapter.loadProject(project, molText, {
  indexToAtomId, // 必须是导出 MOL 的实际原子行顺序；不能猜测
  revision: model.revision,
});

// 在一次拖动/构象操作结束时读取，不需要每帧重新导入或优化。
const result = adapter.captureGeometry();
// 先由宿主检查手性与当前版本，之后提交一个统一历史步骤。
model.applyAsync(result.baseRevision, result.project, 'JSmol 构象编辑');

// 明确的拓扑动作只能先取得变更提案，不能当成纯坐标提交。
const proposal = adapter.captureEditProposal();
```

稳定 UUID 保存在工坊模型中；JSmol 的 `property_editorid` 保存精确数值 token，外部映射 token → UUID。新增 JSmol 原子的 token 为 0，由适配器分配新 ID；已经删除的 ID 不复用。MOL 导入必须附带行序映射，C60 不使用 SMILES 或对称图匹配。`atomno`、`atomIndex` 和 `atomname` 均不承担永久身份。

`captureGeometry()` 核对原子 ID、元素、电荷、键端点及键级没有变化后，仅合并三维坐标到原模型快照，并带上源编辑版本。它保留手性元数据，但**不会判定新坐标是否已经翻转手性**；正式接入时仍须用现有化学服务验证。主模型负责原生文件、历史、另一视图更新和失效异步结果，JSmol 不另外拥有一套项目历史。

`captureEditProposal()` 返回新增、删除、改变的原子/键、移动的原子和完整快照。它**还没有实现**选择性接纳副作用、二维新原子排版、四面体/E/Z 拓扑变更后的立体信息更新、自定义片段边界同步或主 UI 确认流程。这些工作不能用一次 MOL 导入替代。

## 设置和需要保留的限制

经实际运行验证的关闭方式：

```text
set autobond false;
set modelKitMode true;
modelkit addHydrogens false;
modelkit autobond false;
set picking dragAtom;
```

`getProperty("modelkitInfo")` 中两个开关确为 `false`。它们能阻止新原子自动補氢和距离成键，却不保证所有原生 ModelKit 化学编辑完全没有附带行为。原子元素修改仍可能删除显式 H，单邻接原子的元素修改仍可能调整键长。试图把这些动作一律作为“只换元素”直接同步会破坏玩家构象。

本次 runtime 可靠执行的片段旋转脚本是：

```text
modelkit rotate branch {atomindex=0} {atomindex=1} 35;
```

旧示例中的 `modelkit rotate bond ...` 在此固定版本会报错。选择环内键也不能保证能独立旋转一个无冲突片段，正式界面必须提前检查图拓扑。

原生 JSmol MOL 输出不保存 `property_editorid`，也不保存工坊片段/二维布局。原生 JSmol undo/redo/state restoration 在本轮增删场景中不能充当安全的项目历史；适配器阻止经自身接口发出的这些命令。工具内部菜单或直接调用全局 Jmol API 可以绕过这层保护，正式接入时需收拢这些入口。

以上是最初独立验证阶段的结论。后续已经实现意图过滤宿主并接入主工坊，实际编辑的附带改动不会直接回写共同项目；具体实现、接口和 18 项宿主验收见 [JSmol 主三维宿主](jsmol-host.md)。本独立研究页仍保留原生未过滤行为，便于复现上游限制。

## 版本、来源与许可证

固定上游 runtime：JSmol 16.3.31，官方 2025-07-22 发布的 [jsmol.zip](https://sourceforge.net/projects/jmol/files/Jmol/Version%2016.3/Jmol%2016.3.31/jsmol.zip/download)。压缩包 SHA-256：`4147ba69d31a83ed6c276a1feaf3bd2b6472e82075ec5c96dd03317d0e568655`。

`frontend/editor/vendor/jsmol/` 保留完整 j2s runtime、压缩与未压缩启动代码和上游 README；上游文件共 49,756,122 字节，约 49.8 MB。没有打包完整演示、Java 二进制或 PHP 转发服务。

Jmol/JSmol 为 LGPL-2.1-or-later；启动代码内含 jQuery 3.7.1（MIT）。原许可、版权声明、来源及文件哈希保留在 `vendor/jsmol/NOTICE.md`、`LICENSE-LGPL-2.1.txt`、`LICENSE-jQuery-MIT.txt` 和 `manifest.json`。官方同版本的[完整源码发布包](https://sourceforge.net/projects/jmol/files/Jmol/Version%2016.3/Jmol%2016.3.31/Jmol-16.3.31-full.tar.gz/download)已记录；库保持可替换。2026-10-13 的 [原生菜单接入](jsmol-native.md) 修改了五份菜单相关库文件，原件、补丁、原始与安装哈希均保留。本记录不替整个 Molecular Car 项目选择许可证。

此前 `research/README.md` 的源码探针结论描述较早验证阶段；本文件和真实浏览器 JSON 记录是后续补充。保留独立研究入口是为了复现上述实测副作用，主工坊现已通过宿主管理真实 JSmol 编辑。
