# JSmol 原生 ModelKit 菜单与共同工程

日期：2026-10-13。主工坊现在使用 JSmol 16.3.31 的真实 ModelKit
弹出菜单和鼠标操作。右键三维画布，或点击 JSmol 原生左上入口，可以打开
原子、键和其它操作菜单。菜单来自库的资源与动作源代码；没有在画布上
注入仿造工具栏。主页面保留精确键长、键角和二面角面板。

菜单中的移动原子、移动分子、Alt 拖动旋转及绕键旋转，执行原生 JSmol
坐标变换。普通空白拖动转动相机。双击原子选择明确连接的整个分子。
普通点击 C60 后拖动，移动记录中的整个笼；直接连接的末端氢一起移动，
连接在笼外的车身不自动加入。导入的折叠 SUP 基团也按明确成员整体移动，
不靠名称、距离或对称图猜测模板身份。

## 事务边界

原生工具先执行一次用户操作，宿主随后只接纳明确意图：

| 原生入口 | 共同工程如何接纳 |
|---|---|
| 元素、形式电荷、键级、删除 | 只接纳目标字段或明确端点；过滤 ModelKit 附带的删氢和键长调整。化学编辑后重载规范项目。 |
| 添加原子、拖动连接 | 只接纳实际新原子、指定端点与一条明确键；不按距离自动成键，不自动补氢。 |
| 移动、Alt 旋转、绕键旋转 | 核对化学图不变，回收相对原生起点的坐标增量，保留共同工程的精确未动坐标。更新基线后保留原生模型和相机，不每次拖动重载。 |
| 生成三维、补氢、整理几何 | 转接工坊几何服务与版本检查；不执行 JSmol 内部的补氢或最小化脚本。 |
| 撤销、重做、保存项目、导出 MOL | 转接共同工程，所有视图使用同一份完整历史和文件管线。 |

源菜单去掉了“删氢再最小化”、拖动并自动最小化、独立状态保存恢复、
晶体模式、清空模型和退出编辑等可绕过共同工程的入口。产品的全局
JSmol 调试对象仍可在开发者控制台操作；它不是安全沙箱。

JSmol 数字属性 `property_editorid` 对应外部稳定 UUID。MOL 行顺序使用
明确映射；C60 等对称结构也不使用 SMILES 或图匹配恢复身份。原生历史
对自定义属性恢复不可靠，因此菜单撤销转接共同历史。预览只读，旧版本
拖动会取消；异步服务结果由工坊版本检查决定是否接纳。

快速连续拖动保持原生几何手势，避免上游把第二次拖动变成相机操作。
上游 JMenu 的菜单动作默认等待 500 ms；本工坊有 owner bridge 时执行延迟
改为 0。未接入 bridge 的上游研究页保留原时序。视角复位使用
`set zoomLarge false`，按画布短边适配，防止宽屏 C60 上下截断。

## 接入接口

```javascript
new JSmolHost(container, {
  model,
  onMode(mode, options) { /* 更新周边操作说明，不再次设置原生模式 */ },
  onAction(name) { /* 调用共同工程的命令，可返回 Promise */ },
  onSelection(ids, bondId) { /* 同步共同选择和精确几何面板 */ },
  onInteraction() {},
  onError(error) {},
});
```

`onAction` 名称为 `generate`、`hydrogenate`、`relax`、`undo`、`redo`、
`save`、`export`。几何服务不调用 `/api/load`；送入赛道仍须玩家点击。

`onMode` 模式包括 `select`、`measure`、`move`、`molecule`、`rotateBond`、
`atom`、`charge`、`bond`、`delete`。选项含 `native:true`、
`nativePicking`，以及适用的 `element`、`order`、`delta`。
原生增减键级的 `order` 可以是 `p` 或 `m`；父页面只能用有效数值更新
自己的键级下拉框。选择、移动、移动分子、绕键旋转和测量快捷键默认
使用原生手势；精确构造快捷键与按指定轴旋转保留现有共同工程入口。

## 实际验证

```powershell
node frontend/editor/tests/jsmol-native.cjs
node frontend/editor/tests/jsmol-host.cjs
node frontend/editor/tests/jsmol-live.cjs
```

三套脚本各启动自己的本地 Chrome，无头浏览器静态读取项目文件，禁止
外网和比赛 API 请求，不启动、停止或占用现有服务。原生菜单测试以真实
鼠标点选菜单和画布，没有用脚本调用替代菜单动作。机器依赖配置见
[原生运行验证](jsmol-live.md)。

本轮 **19 / 19 项通过**。记录位于 `outputs/editor-jsmol-native-check/report.json`，含当前宿主和
库改动的哈希、动作记录和截图。验证覆盖原生元素添加与替换、拖动连键、
电荷、键级、删除、菜单共同撤销重做、真实拖动与旋转、四苯基苯二面角
变化、C60 与末端氢整体移动、无模板身份的导入 SUP 成员移动、无操作
点击的副作用恢复，以及 50 ms 内生效的菜单动作。没有外部请求或脚本错误。
选择与几何模式显式清除旧的键赋值动作，防止仅关闭 `bondPicking` 显示后，
残留模式仍在双击附近化学键时暗中修改键级；对应宿主回归也已通过。

原始运行测试 15 项、宿主回归 18 项继续保留，分别覆盖原生副作用的复现、
MOL 手性输出、精确坐标、释放时撤销和旧操作版本失效。主工坊真实页面、
文件选择器和 LAMMPS 接入由主页面验收记录覆盖，不能用独立宿主测试替代。

## 源码与依据

五份库文件的修改、上游原件和可逆补丁见
[source extension](../jsmol-source/README.md)。完整 LGPL 和 jQuery MIT
许可保留在 [NOTICE](../vendor/jsmol/NOTICE.md)，项目没有替整体选择许可证。

采用原生鼠标方式的官方依据是
[Jmol Mouse Manual](https://wiki.jmol.org/index.php/Mouse_Manual)；原生
ModelKit 和最小化命令见
[Jmol scripting documentation](https://chemapps.stolaf.edu/jmol/docs/#modelkit)。
当前固定版本的 [16.3.31 发布说明](https://sourceforge.net/projects/jmol/files/Jmol/Version%2016.3/Jmol%2016.3.31/)
仍列有 UNDO 与 FIXED/MODELKIT MINIMIZE 待办，故不把其内部历史与约束优化
当作工坊永久数据管线。原生菜单资源本身使用静态 bundle，通用
[Custom Menus](https://wiki.jmol.org/index.php/Custom_Menus) 的主菜单配置不能
直接替换这个 ModelKit bundle；这里修改真实 native source，保持可审查。

尚存边界：原生自由拖动可能改变立体几何，但适配器保留原拓扑立体标记，
导出与几何服务检查冲突。绕环内键不能保证得到独立刚性片段；精确二面角
面板会根据明确拓扑判断能否分离。首次启动仍需载入本地 Java2Script
运行库，普通三维编辑不运行每帧优化。
