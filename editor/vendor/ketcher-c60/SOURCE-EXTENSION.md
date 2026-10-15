# Ketcher 原生 C60 模板扩展

运行入口：`/editor/vendor/ketcher-c60/index.html`。发行目录已经包含浏览器所需的脚本、样式、Indigo 工作线程和资源；运行时不需要 npm、外部 CDN 或旧的 `fragments.json` 片段库。

这是 EPAM Ketcher 3.18.0 的源码构建。固定来源为 [官方 v3.18.0](https://github.com/epam/ketcher/releases/tag/v3.18.0)，提交 `f4f013e1855b71b0cbd37fc79d55b0cbe51c2371`。原有 `vendor/ketcher/` 官方预构建发行保持不变。源码采用 Apache-2.0；发行中的 `LICENSE`、`NOTICE`、`licenses/` 和 `BUILD-DEPENDENCIES.json` 保留上游及实际打包依赖的版权和许可证。

## 原生扩展

当前扩展版本 `3.18.0+molecular-car.c60.3`。本次仅精简模板放置提示；模板数据及来源映射仍使用版本 2。`prepare-native-c60.cjs` 修改五处上游源码：

- `packages/ketcher-react/src/script/ui/data/templates.js`：在原有八个环模板后追加完整 C60 MOL，使用原生 SUP 缩写基团，默认收缩显示 `C₆₀`；仅该模板附加明确来源标记。
- `packages/ketcher-react/src/components/Icon/utils/iconNameToIcon.ts`：注册 `template-8`。
- `packages/ketcher-react/src/script/ui/action/action.types.ts`：追加原生工具名类型。
- `packages/ketcher-react/src/assets/icons/files/template-8.svg`：原创 C₆₀ 文字图标。
- `packages/ketcher-react/src/script/editor/tool/template.ts`：只对 C60 模板阻止在已有原子或键上的融合，提示先空白放置、再用单键连接。其它模板操作保持原生行为。

底部按钮继续由 Ketcher 的 `TemplatesList` 和 `ActionButton` 渲染，使用原有模板动作及 `TemplateTool`。没有 DOM 注入按钮，没有改写压缩后的 JavaScript。`SOURCE-MODIFICATIONS.patch` 和 `template-8.svg` 可独立审查改动。

`c60-source.json` 保存从项目已有截角二十面体生成器复制的真实化学图，包含 60 个碳、90 条键，其中 30 条双键、60 条单键。每个碳度数为 3、键级和为 4。该快照只用于构建；原旧片段库文件不参与运行。生成器以五边形面做 Schlegel 平面展开，保留原子的输入顺序。收缩时只显示原生 `C₆₀` 名称，展开时显示全部二维原子及键。`c60-template.mol` 和 `extension-data.json` 记录实际模板及数据校验摘要。

同一来源还生成 `../native-c60-data.js`，按明确的模板索引保存真实球壳三维坐标，初始键长均为 1.42 Å，直径约 7.04 Å。此坐标仅是几何初值，后续整理可以形成长短不同的单/双键；不把二维 Schlegel 键长当成物理坐标。

## 操作及边界

单击底部 C₆₀，再点击画布空白插入完整 C60。默认显示原生缩写，右键 `Expand Abbreviation` 展开、`Contract Abbreviation` 收缩；`Remove Abbreviation` 后可逐原子编辑。缩写没有替换或丢弃实际60原子/90键。工坊历史同时保存缩写状态、全部拓扑和构象。

对 C60 已禁止模板融合：点击已有原子或键不会插入，不会把氢偷偷变成笼碳。先在空白放置，再从 `C₆₀` 缩写用原生单键工具拖出连接。工坊把外接单键作为明确加成处理：将连接位点唯一相邻的笼内双键改为单键，另一个笼碳在准备三维时补氢；以同一可撤销事务提交并显示说明。这是原型的通用加成绘制方式，并不承诺合成路径。连接另一端已满价、多个过价键或没有可用于加成的内部双键时拒绝该连接并说明原因。用户仍可保存其他尚未完成的设计。

来源标记 `sgroup.data.molecularCarTemplate={id:'c60',version:2}` 由原生模板源码附加，原生 `fromPaste` 明确映射有序 SUP 成员到新 Pool ID。`ketcher-adapter.js` 再把这些 Pool ID 绑定为稳定 UUID，保存为 `metadata.ketcherGroups` 和 `template.indexToAtomId`。不使用 SMILES、相似图匹配、坐标猜测或用户化学 atom-map 恢复 C60 身份。用户自己建立一个名叫 C60 的缩写不会被当成此模板。原生项目保留来源、缩写状态和三维坐标；普通 MOL 导出保留完整化学图和坐标，但工坊专有缩写/来源信息需要原生项目保存。

## 复现构建

开发环境要求 Node.js 24.14.1 或更高、npm 11、Git，以及首次下载时的网络。所有安装都位于本目录的 `build/`；项目公共依赖无需更改。下列命令从本目录执行：

```powershell
git clone --depth 1 --filter=blob:none --sparse --branch v3.18.0 https://github.com/epam/ketcher.git build/upstream
git -C build/upstream sparse-checkout set packages example
Set-Location build/upstream
npm ci --ignore-scripts --no-audit --no-fund --include-workspace-root --workspace=ketcher-core --workspace=ketcher-react --workspace=ketcher-standalone --workspace=ketcher-macromolecules --workspace=ketcher-react-app
Set-Location ../..
New-Item -ItemType Directory -Force build/tools
Invoke-WebRequest https://registry.npmjs.org/less/-/less-4.2.0.tgz -OutFile build/tools/less-4.2.0.tgz
node install-compatible-less.cjs
node build-native-c60.cjs
```

固定上游 lockfile 的 Less 4.6.4 不导出 `package.json`，与 Vite 8.0.16 的预处理器定位方式不兼容。因此构建目录使用官方 Less 4.2.0，保存原依赖备份并在 manifest 记录下载地址与 SHA-256。没有修改上游 lockfile。`install-compatible-less.cjs` 核对目标位于本目录的 `build/` 后再移动备份。

构建沿用上游 Vite 的源代码解析、SVG 转 React、SDF 导入和 Indigo worker 配置。生成的 `vite-native.config.js` 以固定常量替代仅用于提取常量的 Webpack/Rollup 配置导入，明确构建根目录，并禁用 CSS 压缩，以保留上游旧式 container-query CSS。Indigo 及 Ketcher 核心化学计算不变。`build-native-c60.cjs` 会核对源码提交、预生成官方 KET schema、构建发行、保存依赖声明和写出 `vendor/ketcher-c60-manifest.json`。本机受限环境可能禁止 Vite 原生解析器定位开发依赖，需要允许构建器读取本隔离开发目录；运行页面不需要此权限。

`dependency-notices-plugin.cjs` 以实际输出块的模块清单追踪依赖包及许可证；Indigo 等工作线程的版权还由原官方发行的 `.LICENSE.txt` 完整保留。`build/` 是可重建的开发材料，运行时不需要。

## 实际验证

测试从项目根执行。设置 `PLAYWRIGHT_MODULE`、`CHROME_PATH` 可使用本机已安装的 Playwright 和 Chrome；默认使用本任务的本地运行库。

`node frontend/editor/ketcher-source/probe-superatom.cjs`：最小真实 Ketcher 验证，原生 KET SUP 收缩显示、保持有序60原子/90键及连接点；Unicode 下标名称 C₆₀ 原生 MOL 导出/重导入往返成功。

`node frontend/editor/ketcher-source/probe-native-c60-v2.cjs`：13项真实浏览器集成通过，包含底栏折叠插入、源码来源标记、稳定有序 UUID 和球壳键长、原生右键展开、共用撤销/重做、防融合、原生拖出单键、加成降低指定双键、原生重推同步和一步撤销。记录为 `outputs/ketcher-native-c60-v2/report.json`，截图 `collapsed.png`。页面错误、外部请求及 API 调用均为零。

`node frontend/editor/ketcher-source/probe-native-c60.cjs`：13项纯上游绘图工具检查通过，覆盖原生环栏、60/90拓扑、原生历史和 MOL 重导入。该纯工具页面没有工坊 canonical 化学保护，其连接会复现上游原始的过价状态；合法加成行为由上面的 v2 集成测试验证。

`node --disable-warning=MODULE_TYPELESS_PACKAGE_JSON --test tests/editor_ketcher_adapter.test.mjs`：28项通过，包括 SUP 有序映射、任意同名组不猜测身份、复制产生新 ID、加成价态检查和原构象精确保留。

所有独立浏览器测试通过本地静态文件拦截启动 Chrome，不监听公共端口、不触碰比赛服务。实际主页面切换与几何整理由主页面和后端测试另行验证。
