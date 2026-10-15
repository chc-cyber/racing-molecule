# 分子赛车工坊：依赖、来源与许可证

本文件记录编辑器所用依赖和研究资料，不为项目整体选择许可证。组件在本地托管，不需要公共 CDN。

## 原型运行时组件

| 组件 | 固定版本 | 用途 | 许可证与本地声明 |
|---|---|---|---|
| Ketcher（原生 C₆₀ 扩展） | 3.18.0+molecular-car.c60.3，固定官方 3.18.0 源码 | 主工坊二维绘图、原生环模板及折叠 SUP；内置本地 Indigo WASM | Apache-2.0；`vendor/ketcher-c60/LICENSE`、`NOTICE`、`licenses/` 与 `BUILD-DEPENDENCIES.json`；[官方来源](https://github.com/epam/ketcher/releases/tag/v3.18.0) |
| Ketcher 图标快照 | 3.18.0，固定同一官方提交 | 工坊顶部及三维通用操作的原始 SVG，统一二维/三维外观 | Apache-2.0；原样 SVG、`vendor/ketcher-icons/LICENSE`、`NOTICE`、`manifest.json` 记录来源和 SHA-256；运行内联使用 `currentColor` |
| Three.js（历史编辑器） | 0.180.0 / r180 | 保留旧视图源码与历史测试，主工坊不再加载；赛道仍有自己的 Three.js | MIT；`../vendor/THREE-LICENSE.txt`；[官方源码](https://github.com/mrdoob/three.js/tree/r180) |
| Three.js TransformControls（历史） | r180，配套 Three.js 0.180.0 | 旧视图与历史测试保留，主工坊不再加载 | MIT；`vendor/THREE-LICENSE.txt`；[该版本官方源文件](https://github.com/mrdoob/three.js/blob/r180/examples/jsm/controls/TransformControls.js) |
| JSmol / ModelKit | 16.3.31 | 主工坊三维原生 ModelKit 手势及共同工程桥接；侧栏调用同一原生操作，原生菜单只在独立研究页保留 | LGPL-2.1-or-later；`vendor/jsmol/LICENSE-LGPL-2.1.txt`、jQuery MIT 与 `NOTICE.md`，详见下文 |
| Lucide（历史） | 1.8.0 | 历史界面图标文件保留，当前主工坊不再加载 | ISC；含衍生自 Feather 的部分图标的 MIT 声明，均完整保留在 `vendor/LUCIDE-LICENSE.txt`；[官方仓库](https://github.com/lucide-icons/lucide)、[官网](https://lucide.dev) |
| RDKit | 2025.9.3 | 独立几何准备与化学检查，使用现有 Python 环境 | BSD-3-Clause；已安装分发包保留许可；[官方许可证](https://github.com/rdkit/rdkit/blob/Release_2025_09_3/license.txt) |
| NumPy | 1.26.4 | 现有 RDKit/模拟环境中的数值依赖 | BSD-3-Clause；已安装分发包保留许可；[官方许可证](https://github.com/numpy/numpy/blob/v1.26.4/LICENSE.txt) |

MOL/SDF 与原生项目文件在本地读写。C₆₀ 模板是本项目生成的真实原子与键数据，非视觉贴图。原 `vendor/ketcher/` 官方发行保持原样；主页面使用含原生 C₆₀ 扩展的源码构建，扩展源码、复现步骤和独立来源清单另行保留。旧赛车片段库 UI 已移除。

当前常驻三维图标栏是工坊宿主界面，使用 Ketcher 原始 SVG、相同的白色分组和青绿色激活状态，操作调用 JSmol 原生 picking / ModelKit API。`vendor/ketcher-icons/` 固定 28 份官方图标、1 份已有原创 C₆₀ 扩展图标和 18 份工坊原创三维/专属操作图标；每份来源、修改情况和逐文件 SHA-256 均记录在清单中。官方原始 SVG 不改动，内联渲染只将单色 fill/stroke 归一为 `currentColor`，保留 viewBox 和路径；生成的 `icons.js` 在本地提供字符串，不需要 CDN 或运行时下载图标。Ketcher 的工具栏动作与二维状态并未交给 JSmol 管理，两工具继续使用共同工程。未复制商业工具的侧栏或新增分子编辑引擎。自动三维构建复用现有 RDKit 与 NumPy，不新增依赖，UFF 整理由用户显式触发。新增 Atom、+H、四角适配、UFF 能量曲线与二维排版图标均为原创；绕键旋转原创图标使用带缺口的单环箭头。当前图标清单版本为 `ketcher-3.18.0-workshop-icons.3`；Atom 使用普通字体文字徽标，+H 的 H 使用普通粗体字形，UFF 箭头与曲线分离。`tools/update-icons.cjs` 可重建当前 47 个图标的内联模块和哈希。二维无关工具用原生 hiddenControls 配置隐藏，Ketcher 发行包未因此改动。

Ketcher 来源固定为 tag `v3.18.0`、commit `f4f013e1855b71b0cbd37fc79d55b0cbe51c2371`；原 ZIP SHA-256 `484e7f10a0e74808ae5f43f0e6448db6f4d884bdd2189aa271108474167a3e7e` 与 GitHub 官方 digest 一致。30 份原发行文件约 121 MB，完整保留，无构建步骤。`vendor/ketcher-manifest.json` 记录来源、文件大小和逐文件哈希；`node frontend/editor/research/ketcher-manifest.cjs` 可核对。Ketcher 的 About 入口与版权声明保留，原生功能由同源 iframe 使用。

Lucide 的本地 UMD 文件与许可从 Codex bundled workspace runtime 所带 `lucide` npm 包原样复制。版本由包内 `package.json`（1.8.0）与 UMD 文件头（`lucide v1.8.0 - ISC`）核对。版权声明包括 Copyright (c) 2026 Lucide Icons and Contributors 与 Feather 图标的 Copyright (c) 2013-present Cole Bemis；以保留的完整原文为准。工坊不向 CDN 请求图标。

历史 TransformControls 配合本地 Three.js 使用，主工坊已移除相应 import map。Three.js 原始 MIT 声明 Copyright © 2010-2025 three.js authors 在工坊 `vendor/THREE-LICENSE.txt` 原样保留，赛道已有 Three.js 声明同时保留。使用变换手柄不把原子编号、分子拓扑或撤销状态交给 Three.js 管理。

## 原生 C₆₀ 模板构建

主页面使用 `vendor/ketcher-c60/`，其源代码固定为上述官方提交，扩展原生模板数据、SUP 来源元信息、工具名类型、图标注册及原创 C₆₀ 文字 SVG。TemplateTool 源码加入禁止融合已有原子/键的检查；化学连接通过共同工程的明确加成交易处理。没有 DOM 插入按钮或压缩 bundle 盲改。完整修改补丁为 `ketcher-source/SOURCE-MODIFICATIONS.patch`。

[构建说明](ketcher-source/README.md) 提供下载和复现步骤，`vendor/ketcher-c60-manifest.json` 记录版本、实际依赖版本、来源与 169 个发行文件的 SHA-256。新发行保留 Apache LICENSE / NOTICE、原官方第三方声明，并按实际输出模块收集 144 个依赖包的许可证与版权文件。构建工具的 Less 兼容版本调整只位于独立开发目录，来源和哈希已记录，不修改项目公共依赖配置。

C₆₀ 图来自项目自生成的碳笼数据，构建时复制为 `ketcher-source/c60-source.json`，内部使用 Schlegel 展开二维图，原生 SUP 默认折叠显示 C₆₀；显式来源和顺序映射传递已有球壳三维坐标。保留原始顺序，不靠图匹配识别原子；未连接时完整 60C、90 键、30 双键可核对。运行不需要旧片段库或 npm 安装。`ketcher-source/build/` 是可重建的开发材料，原官方 `vendor/ketcher/` 仍作为原始来源和历史验证保留。

## 交互研究与元素事实数据

ChemDraw、GaussView、Reaxys 和 Marvin JS 仅用于交互研究，没有嵌入其商业组件、图标、商标标识或运行代码。各资料的版权仍归原提供方；[UX-CHANGES.md](UX-CHANGES.md) 与 [research/reaxys-interaction.md](research/reaxys-interaction.md) 记录官方来源以及哪些是工坊自己的改造。

`periodic-table.js` 是本项目自写代码，内嵌 118 个元素的符号、原子序数和中文名等事实数据。符号与常规布局参照 [IUPAC 官方周期表](https://iupac.org/what-we-do/periodic-table-of-elements/)（4 May 2022），中文名参照 [中国化学会官方中文表](https://www.chemsoc.org.cn/library/a2611.html)（2019-01-23）与 [全国科技名词委发布说明](https://cnterm.cn/kxcb/kpwz/201810/t20181029_758620.html)。原表版权属于 IUPAC 与中国化学会；工坊未分发原表 PDF、图形版面或整套原子量数据，仅独立录入所需事实。数据颜色分组用于帮助定位，不表示模拟参数或试跑准入结论。

## JSmol 主三维组件

JSmol 16.3.31 的 `JSmol.min.js`、`j2s/` 和未压缩启动代码在 `vendor/jsmol/` 本地托管，主工坊和独立 ModelKit 研究入口均使用。发行来源为 [Jmol 官方 SourceForge 16.3.31 目录](https://sourceforge.net/projects/jmol/files/Jmol/Version%2016.3/Jmol%2016.3.31/)，归档 SHA-256 `4147ba69d31a83ed6c276a1feaf3bd2b6472e82075ec5c96dd03317d0e568655`。保留 LGPL-2.1、jQuery MIT 文本、原 README 和运行文件内的版权头；声明 / 文件哈希见 `vendor/jsmol/NOTICE.md`、`manifest.json`，官方同版本源码包与替换说明见 [research/jsmol-live.md](research/jsmol-live.md)。JSmol 启动文件未修改；五份库文件带本地原生菜单桥接与即时菜单动作扩展。可读源码、完整原始文件、修改补丁与双份哈希见 [jsmol-source/README.md](jsmol-source/README.md)。原库修改继续采用上游 LGPL-2.1-or-later，适配器为独立模块。

早期 ModelKit 编辑与原生限制见 `outputs/editor-jsmol-check/report.json`。主工坊通过独立宿主与意图过滤控制附带坐标/氢变化及身份映射，原生菜单的几何、历史、保存动作转接共同工程；真实原子/键编辑和坐标手势继续使用 JSmol。合并 Swing 的 packed 文件仅修改与可读源码一致的条件菜单计时表达式，未重新压缩；最新实际集成范围见 `VALIDATION.md`。整包无关 Java 二进制、PHP 代理和示例没有分发，不需要额外 Java / PHP 服务。

## 首轮最小源码验证资料（历史，不由主页面加载）

| 组件 | 研究版本 | 许可 | 来源及保留位置 |
|---|---|---|---|
| Ketcher | 官方提交 `6fb746d04f86ec368ebcfe020d51367ce354502e`；另读取 npm 3.18.0 元数据 | Apache-2.0，Copyright EPAM Systems（年份见各原文件） | [官方仓库](https://github.com/epam/ketcher/tree/6fb746d04f86ec368ebcfe020d51367ce354502e)；`research/upstream/ketcher/LICENSE`；各原文件版权头保持原样 |
| Jmol/JSmol | 官方提交 `9308996bb219424b6cadd3fc8dc2cee04fa58e28`，源码声明 16.5.2 | ModelKit 源文件声明 LGPL-2.1-or-later，Copyright The Jmol Development Team（年份见原文件） | [官方仓库](https://github.com/BobHanson/Jmol-SwingJS/tree/9308996bb219424b6cadd3fc8dc2cee04fa58e28)；`research/upstream/jsmol/LICENSE.txt`；原文件内容保持原样 |

研究目录只包含验证所需的少量未修改源文件、许可证及元数据，不是完整发行包。`research/upstream-manifest.json` 记录每份来源、提交、大小和 SHA-256；文件名中的双下划线对应源路径分隔符，内容未改。上游同名 `package.json` 仅为来源记录，不代表本项目增加 npm 构建步骤。

首轮 Ketcher 函数探针通过 Node 内建 TypeScript 去类型功能执行若干原始函数；首轮 JSmol 只作源码检查。历史依据见 `research/README.md`；2026-10-12 的完整本地 runtime 和真实浏览器验证以上述新记录为准。

项目原有 LAMMPS 与其它第三方说明仍在根目录 `THIRD_PARTY.md`；工坊没有新增或重新分发 LAMMPS 二进制，也没有更改模拟力场。

浏览器整合包的生成阶段将 Ketcher 入口、模块预载及 Indigo Worker 的固定根路径改为相对路径，以支持仓库子目录部署；原生发行目录不改写。对应适配在 `wasm/tools/package_web.py`，主题由本项目独立的 `ketcher-theme.css` 注入。
