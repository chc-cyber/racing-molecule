# 成熟编辑器最小技术验证（首轮历史）

此文记录 **2026-10-11 首轮**，不是当前运行方案。2026-10-12 已完整本地部署 Ketcher 3.18.0 并接入主二维画布；JSmol 16.3.31 ModelKit 已完成真实浏览器编辑验证。当前方案与入口见 [工坊说明](../README.md)、[实际验收](../VALIDATION.md) 和 [JSmol 实测](jsmol-live.md)。下方保留当时源码探针的范围与决策依据。

验证日期：2026-10-11。结论：**本阶段采用现有 Three.js + RDKit 与项目自有共享分子图，暂不将 Ketcher/JSmol 接入运行页面。** 两者都有真正编辑能力，许可也不是排除原因；主要障碍是当前需求所需的跨工具永久原子身份、立体信息和统一历史记录，不能依靠反复导出/导入 MOL 自动成立。

本记录严格区分实际验证范围：执行了 Ketcher 原始 Pool、KET 原子序列化函数的隔离探针，检查了 JSmol ModelKit 和 MOLWriter 的固定版本源码；**没有完成两个完整浏览器编辑器的端到端集成，也不声称做过完整 JSmol GUI 编辑测试**。发现数据契约障碍后停止扩大依赖安装，转向可直接验证共享图不变量的替代实现。

## 来源、固定版本与复现

| 组件 | 检查版本 | 原始来源 |
|---|---|---|
| Ketcher 源码 | `6fb746d04f86ec368ebcfe020d51367ce354502e` | [EPAM 官方仓库](https://github.com/epam/ketcher/tree/6fb746d04f86ec368ebcfe020d51367ce354502e) |
| Ketcher npm 元数据 | `ketcher-core/react/standalone` 3.18.0 | [官方 npm registry](https://registry.npmjs.org/ketcher-core/3.18.0) |
| Jmol/JSmol 源码 | `9308996bb219424b6cadd3fc8dc2cee04fa58e28`，源码声明 16.5.2 | [维护者官方仓库](https://github.com/BobHanson/Jmol-SwingJS/tree/9308996bb219424b6cadd3fc8dc2cee04fa58e28) |

版本号仅陈述本次读取结果，不表示任意未来版本都有同样行为。19 份源文件与许可证放在 `upstream/`，保持原始文件内容；`upstream-manifest.json` 记录完整 URL、提交号、字节数和 SHA-256，也保存了 npm 固定版本的下载 URL、完整性值与依赖元数据。

从项目根目录执行（Node 24+，无需安装 npm 包）：

```powershell
node --disable-warning=ExperimentalWarning --disable-warning=MODULE_TYPELESS_PACKAGE_JSON frontend/editor/research/probe-upstream.mjs
```

探针重新校验所有源文件哈希，执行 Ketcher 原始 Pool 与抽取的原始 KET 序列化函数，然后针对 JSmol 已检查的源代码契约进行断言，更新 `probe-results.json`。本次 10 项检查全部通过。它是技术决策的可重复证据，**不是替代实际原型验收的 UI 测试**。

## 1. Ketcher 确实适合二维化学绘图，但需要专门 ID 适配

官方 API 提供 `setMolecule()`、`getMolfile()`、`getKet()` 和编辑变更事件；源码具备独立的 AtomAdd、AtomDelete 与化学属性模型。KET 原子转换保留形式电荷、`stereoParity`、`stereoLabel` 及坐标。它是成熟二维化学编辑器，不能因为本次未采用而描述为“缺乏编辑能力”。

实际运行 `Pool` 与 `moleculeToKet()` 得到：

- 添加 3 个原子、删除中间原子、再添加一个：当前池 ID 为 `[0, 2, 3]`，原有 ID 在一次编辑会话中保持。
- 序列化后的 KET 原子记录**没有池 ID 字段**，将其装入新池后 ID 为 `[0, 1, 2]`。
- KET 会保留化学原子映射字段 `mapping = source.aam`，本探针映射 `[101, 103, 104]` 保持不变。因此可构建适配器，但不能把内部池索引或 MOL 行号直接视作项目 UUID。
- 源码序列化也保留 `z`，但这不等于编辑二维布局时能自动保护另外一份玩家三维构象；坐标分离仍须由项目负责。

依据：[Pool](https://github.com/epam/ketcher/blob/6fb746d04f86ec368ebcfe020d51367ce354502e/packages/ketcher-core/src/domain/entities/pool.ts)、[KET 序列化](https://github.com/epam/ketcher/blob/6fb746d04f86ec368ebcfe020d51367ce354502e/packages/ketcher-core/src/domain/serializers/ket/toKet/moleculeToKet.ts)、[KET 导入](https://github.com/epam/ketcher/blob/6fb746d04f86ec368ebcfe020d51367ce354502e/packages/ketcher-core/src/domain/serializers/ket/fromKet/atomBondToStruct.ts)。

成熟集成需要拦截新增、删除、复制、片段插入、导入、undo/redo 等路径，并维护 `项目 UUID ↔ Ketcher 内部池 ID ↔ 临时原子映射号`。不得靠 SMILES/对称图匹配恢复 C60 身份，也不能占用或覆盖用户原本的化学原子映射而不告知。

## 2. JSmol ModelKit 确实能编辑三维，默认操作语义却不能直接照搬

ModelKit 源码存在 `assignAtom`、`assignBond`、实际 `setAtomCoord` 和围绕键轴的 `rotateAboutPointsInternal`。这是原子、键和构象编辑，不只是旋转相机。

同时，普通 ModelKit 分子操作代码会按模式删除附氢/非共价键、调整单邻接原子位置、补氢。源代码还保留可选自动成键路径，并明确说明该旧式自动成键已不是默认、也不推荐。因此不能误称“JSmol 每次移动都会按距离重建键”；准确问题是：必须确定各条编辑路径的配置及输出差异，才能把它用于严格显式拓扑的编辑器。

依据：[ModelKit 实现](https://github.com/BobHanson/Jmol-SwingJS/blob/9308996bb219424b6cadd3fc8dc2cee04fa58e28/src/org/jmol/modelkit/ModelKit.java)。

## 3. 简单 MOL 往返不能保证 ID 与立体标注无损

检查 JSmol `MOLWriter` 的具体行为：

- 按导出集合顺序执行 `atomMap[i] = ++n`，生成新的文件局部编号。
- V3000 原子行的 atom mapping 位写入常量 `0`，不会自动保存项目永久 UUID。
- V2000 `getAtomParity()` 根据当前三维坐标和排序后的内部邻居编号重新计算原子奇偶性。
- 单键 near/far 楔线的输出位于 `is2d` 条件分支；不能承诺三维导出会逐字保留输入楔线标注。
- SDF 另有 `jmol_atom_names` 属性记录，但这与 plain MOL 路径不同，而且原子名可能受 ModelKit 编辑操作影响，不能未经验证就当永久身份。

这些行为可用于正常化学导出，却不满足“透明地把两个编辑器 MOL 来回传递，所有身份、手性、历史都自然不变”的假设。三维坐标也可定义手性；本记录不把重新计算一概视为化学错误，只说明它必须受项目自身的立体信息契约与验证控制。

依据：[MOLWriter 完整源码](https://github.com/BobHanson/Jmol-SwingJS/blob/9308996bb219424b6cadd3fc8dc2cee04fa58e28/src/org/jmol/adapter/writers/MOLWriter.java)。

## 4. 离线部署与许可证

官方 npm 元数据的 unpackedSize：

| Ketcher 3.18.0 包 | 解包字节数 |
|---|---:|
| ketcher-core | 24,803,587 |
| ketcher-react | 29,815,930 |
| ketcher-standalone | 109,835,505 |

这些数字是包的原始解包大小，不是生产 bundle 大小，包之间/内部可能含重复资源，也不含全部传递依赖。Standalone 还声明 `indigo-ketcher` 依赖；它能支持离线运行，但应在专门集成时固定依赖并本地托管 WASM/静态文件。**本次只下载小量源文件，没有下载上述完整包，也没有以体积为由断言无法集成。** JSmol 完整离线 runtime 体积本次未测，不编造数字。

Ketcher Apache-2.0；JSmol 检查的 ModelKit 文件声明 LGPL-2.1-or-later，仓库保留 LGPL 2.1 全文。两者都可以研究和集成，但若分发需遵守各自声明；本次原样保留研究文件的许可证和版权头。没有替整个项目决定许可证。详情见 `../THIRD_PARTY.md`。

## 5. 替代实现的必要约束

选择项目现有 MIT Three.js、BSD RDKit，配合轻量二维 SVG 编辑层与三维显式操作。它的优势是让**一份共享分子图**成为唯一权威数据，避免将工具序列化当同步协议；代价是本轮化学绘图体验和复杂标注覆盖范围不及成熟编辑器，必须在交付中列出真实限制。

实施约束已报告主 agent：稳定 UUID；2D/3D 坐标独立；只有明确拓扑操作才改变键；所有编辑事务进入同一 undo/redo；片段持有真实原子、键、连接位点；异步 RDKit 几何结果必须关联编辑修订号；C60 连接需显式处理饱和价态；进入 `/api/load` 只能由“送入赛道”触发。

将来重新接入 Ketcher/JSmol 时，先在苯环、显式手性样例与 C60 上完成上述 adapter 契约测试，再替换交互层。此次原型的操作验收应以主交付的验证记录为准。
