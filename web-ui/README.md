# Racing Molecule 网页界面

主页、设置、二维／三维工坊与比赛赛道已经合为同一个纯浏览器应用。界面入口由 `dist/app.js` 管理，实际工坊与赛道沿用 `frontend/`，完整可发布目录为 `wasm/dist/`。

## 预览

在项目根目录执行：

```powershell
.\.venv\Scripts\python.exe -E -B wasm/tools/package_web.py
.\.venv\Scripts\python.exe -E -B wasm/tools/serve.py --port 8783
```

打开本机 `http://127.0.0.1:8783/`。也可继续使用 `web-ui/serve.py`，现在它指向同一完整目录。

## 页面与状态

- 首页只保留品牌、加入赛事和设置入口。
- 工坊使用真实 Ketcher / JSmol，共享工程图与撤销历史。文件读写、构象生成、补氢和 UFF 在本机 Worker 内调用 RDKit WASM。
- 保存并进入赛道直接保存当前分子并开始准备，不再额外弹出赛事规则选择；沿用赛道界面当前选择的规则，成功准备后才显示赛道。取消会终止这次准备并阻止延迟结果跳页。
- 页面切换保留工坊和赛道实例；离开赛道自动暂停。刷新页面会终止本次比赛，工坊草稿和设置保留在本机。
- 设置只保存分子模型、画面质量和动画过渡等显示选项；标准／自定义赛事统一在赛道界面选择。已开始或暂停的本场比赛不可切换模式。
- 完整网站需要 HTTP/HTTPS 服务，不能用单 HTML 文件直接启动 WASM。`build-offline.mjs` 仅保留为早期设计原型的历史工具，不用于本版交付。

## 美术

共享色板与控件几何为 `frontend/design-system.css`：暖白、深灰、低饱和朱红；纯色双切角按钮。工坊外壳、内嵌二维绘图工具、比赛控件均套用此主题。元素颜色保持化学含义。

交付 ZIP：`outputs/wasm-preparation/racing-molecule-web-preview.zip`。不包含 Python、SDK、编译缓存或工坊上游源代码；第三方运行资源与许可随包分发。本次只完成本地整合，没有发布到外网。
