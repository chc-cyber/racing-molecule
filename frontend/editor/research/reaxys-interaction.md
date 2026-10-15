# Reaxys / Marvin JS 元素选择研究与本地周期表

核对日期：2026-10-11。只借鉴操作行为；新增模块是本项目独立实现，不加载、复制或分发 Reaxys、Marvin JS 的商业代码、图标或界面资产。

## 官方依据及区分

- Elsevier 的 [Reaxys Reference Guide，May 2025](https://supportcontent.elsevier.com/RightNow%20Next%20Gen/Reaxys/CS-0856_Reaxys%20Reference%20Guide%20-%20PPT_v2%205_25.pdf)，第 13 页明确将 **MarvinJS 或 ChemDrawJS** 用于结构/反应绘制。工坊参考的是这里的结构绘图流程。Reaxys 的分子式检索属于查询输入，不能把它的选元素/计数操作当作结构编辑能力。
- ChemAxon 的 [Marvin JS Periodic Table Dialog](https://docs.chemaxon.com/latest/marvin-js_periodic-table-dialog.html)：从原子工具栏打开周期表，选元素后关闭弹窗，所选元素跟随鼠标，点击画布可连续放置。其 Atom List / NOT List 是检索查询原子，本工坊本轮不引入这类查询结构。
- ChemAxon 的 [Marvin JS Getting Started](https://docs.chemaxon.com/latest/marvin-js_getting-started.html)：常用元素按钮与周期表均可选元素；选中原子后用元素符号/工具栏改变原子类型。这里验证了元素选择与原子替换的两种职责。

工坊新增：符号、简体中文名、原子序数搜索；方向键按表格位置移动；单原子选择时明确的“替换当前选中原子”动作；通用入口选择绘制元素，属性入口可 `open({replaceSelection:true})` 直接替换。中文搜索、按方向键空间导航以及此 API 属于工坊设计，不声称来自 Reaxys 的同一实现。

## 数据来源与归属

- 符号、原子序数、常规 18 列表布局依据 [IUPAC Periodic Table of Elements](https://iupac.org/what-we-do/periodic-table-of-elements/)（页面所列版本 4 May 2022）。
- 中文名参照 [中国化学会 IUPAC 元素周期表（中文版）](https://www.chemsoc.org.cn/library/a2611.html) 的官方下载（2019-01-23 版），并与 [全国科技名词委发布说明](https://cnterm.cn/kxcb/kpwz/201810/t20181029_758620.html) 及 [四个新元素中文定名](https://www.cnterm.cn/news/mtsm/202308/t20230823_760185.html) 核对。原表版权属于中国化学会与 IUPAC；本模块仅独立录入元素符号、序数、中文名等事实数据，不包含原表图形、版面、原子量数据库或 PDF。
- 镧系 57–71、锕系 89–103 独立各 15 格；主表第三列对应范围占位，避免在 La/Ac/Lu/Lr 分属第 3 族的讨论中选定争议形式。元素数据的 `column` 是显示列，不当作 f 区元素的化学族编号。
- 颜色用 `data-category` 提供常见教学分组，只辅助定位，不提供元素反应性、稳定性、模拟参数或价态结论。全部 118 个元素可选择；当前模拟引擎是否接受由独立试跑检查说明。
- 少数超重元素中文字符可能依赖系统字体。按钮同时提供元素符号与原子序数；无障碍名称和搜索别名保留“钅尔”“石田”“气奥”等部件写法，避免罕见字缺字阻塞选取。

## 集成约定

```javascript
import {mountPeriodicTable} from './periodic-table.js';

const picker = mountPeriodicTable({
  dialog: document.querySelector('#element-dialog'),
  getSelection: () => ({count: selectedAtoms.length, element: selectedAtoms[0]?.element}),
  onChoose: (symbol, {replaceSelection}) => {
    if (replaceSelection) replaceSelectedAtomWithHistory(symbol);
    else activateAtomTool(symbol);
  },
});
picker.setCurrent('C');
picker.open(); // 选择绘制元素
picker.open({replaceSelection: true}); // 属性入口：点元素即替换单个所选原子
```

组件只回调选择，不修改分子或原子 ID；宿主负责统一历史、化学提示、所选对象的有效性复核。普通入口里的显式替换按钮使用当前绘制元素，且仅选择单个原子时显示。弹窗取消没有回调，也没有编辑记录。所有元素数据内嵌，不请求 CDN。

CSS 由宿主提供：`.periodic-dialog`、`.periodic-header`、`.periodic-title`、`.periodic-close`、`.periodic-current`、`.periodic-hint`、`.periodic-search-row`、`.periodic-search`、`.periodic-results-status`、`.periodic-preview`、`.periodic-table-scroll`、`.periodic-table`、`.periodic-group`、`.periodic-series`、`.periodic-element`、`.periodic-atomic-number`、`.periodic-symbol`、`.periodic-results`、`.periodic-result`、`.periodic-name`、`.periodic-no-results`、`.periodic-footer`、`.periodic-selection`、`.periodic-replace`。`.periodic-table` 需 `display:grid; grid-template-columns:repeat(18,minmax(...));`，行/列定位由模块写入。

Node 数据/布局/搜索/键盘导航验证：`node --test tests/editor_periodic_table.test.mjs`。页面集成还应验证：打开/取消焦点、搜索并按 Enter、方向键选择、替换单原子且撤销、无选中原子时替换动作不可用。

## 本模块实测

2026-10-11：6 项 Node 测试通过。另用现有本地 Playwright 与隔离的 Edge headless 页面执行第 7 项浏览器测试，通过以下真实 DOM 操作：118 个元素按钮全部启用、搜索“氧”按 Enter 返回 O、原子属性入口选 Og 返回明确替换动作、方向键 C→N 并按 Enter 选择、单原子的显式替换按钮、Esc 无选择回调且恢复打开前焦点、弹窗打开后原选择失效时拒绝替换。测试页面只装载此组件，不访问比赛 API，不启动或停止已有服务。完整工坊的模型替换/撤销由宿主集成验证。

可通过环境变量 `EDITOR_PLAYWRIGHT_MODULE`（本地 Playwright `index.mjs` 绝对路径）与 `EDITOR_BROWSER_EXECUTABLE`（本地浏览器绝对路径）启用同一测试文件的真实浏览器测试，不设置时只运行 6 项无依赖测试。
