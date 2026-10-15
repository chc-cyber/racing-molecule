# 发布与构建

自有代码采用根目录 MIT 许可证。第三方代码、资源及其修改部分保留原许可证；MIT 不覆盖整个发行包。含 LAMMPS 的链接 WASM 模块按 GPL-2.0 分发。详见 THIRD_PARTY.md 及工坊的第三方说明。

## 玩家入口

网站：https://chc-cyber.github.io/racing-molecule/

源码：https://github.com/chc-cyber/racing-molecule

网页的 MD 在玩家设备的 Web Worker 内运行，不需要计算服务器。首次准备需要下载较大的模型与 WASM，请保持网络连接。HTTP 缓存由托管服务和浏览器管理；目前不承诺离线运行。刷新会丢失尚未下载的轨迹，比赛及未保存结果会提示确认。

## Windows 开发构建

需要 Python 3.11、Git、Node.js。首次创建环境：

```powershell
py -3.11 -m venv .venv
.\.venv\Scripts\python.exe -m pip install -r requirements.txt
.\scripts\build-web.ps1
```

构建脚本下载校验固定版本的源码、SDK，在项目内编译并生成静态包；不需要全局安装 LAMMPS。耗时及磁盘需求远高于玩家下载。完整步骤和数值验证见 wasm/README.md。Ketcher 和 JSmol 的可编辑源码、修改及独立构建说明见 frontend/editor/ketcher-source 与 jsmol-source。web-ui/dist 是当前手写界面源码，也纳入版本控制。

只改网页时，已有运行库与势表可复用：

```powershell
.\.venv\Scripts\python.exe wasm/tools/package_web.py
.\.venv\Scripts\python.exe wasm/tools/verify_package.py
.\.venv\Scripts\python.exe wasm/tools/serve.py --port 8783
```

## 发行材料与 Pages

main 保存自有源码、浏览器依赖发行及许可证、测试和构建输入。gh-pages 保存由同一版本生成的网站。Release 配套静态网站 ZIP、源码包及 SHA-256 清单；依赖源码归档包含构建使用的 LAMMPS、RDKit、Boost、better-enums、Catch2、emsdk 和 zlib，不应只提供二进制或把源码下载链接写成将来补交。

发布前运行 package_web.py 和 verify_package.py，保留 source-checksums.json 和来源说明。不要提交 .local、.venv、outputs、用户模型、凭据、个人机器配置或构建缓存。源码构建入口不等于已在全新机器验证，仍需独立安装验收。

## 发布验收

- 浏览器流程：工坊绘制、2D/3D、准备、20 nm 比赛、1×/2×、暂停继续、结果与下载、返回工坊。
- 错误恢复：下载中断重试、存储不足仍可入场、无效模型可读错误、刷新前提示。
- HTTPS：项目子路径中的 Worker/WASM/模型正常加载，冷缓存、重复加载与更新后资源一致。
- 手机真机和 Mac/Safari 由维护者测试；桌面窄屏与自动化 Firefox 不能代替这些设备。

自动化需要 Playwright；设置 PLAYWRIGHT_MODULE 指向自己的安装目录，Firefox 需安装对应浏览器。可设置 CHROME_PATH、EDGE_PATH 覆盖路径。

```powershell
node --test tests/browser-handoff.test.mjs tests/browser-workshop.test.mjs tests/playback.test.mjs tests/race-result.test.mjs wasm/tools/transport.test.mjs
node scripts/verify_web_journey.cjs chrome
node scripts/verify_web_journey.cjs firefox
node scripts/verify_web_recovery.cjs
```

公开部署验收可使用同一脚本：设置 `WEB_TEST_URL=https://chc-cyber.github.io/racing-molecule/` 和 `WEB_TEST_LABEL=public-https`，再运行 `node scripts/verify_web_journey.cjs chrome`。它使用独立的新浏览器环境，不读取玩家账号或浏览器配置；完成真实比赛后重新加载，检查 HTTPS 与资源缓存。
