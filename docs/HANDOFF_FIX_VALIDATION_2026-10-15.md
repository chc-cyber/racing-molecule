# 送入赛道等待问题：发布验证

验证日期：2026-10-15。源码提交 `4fed6a2`，GitHub Pages 提交 `1d6776044251e88c73ff21939c87f90bf9cfe588`。

部署成功：[Pages workflow](https://github.com/chc-cyber/racing-molecule/actions/runs/38739570421)。测试地址：https://chc-cyber.github.io/racing-molecule/

## 原因与修复

公开站点准备赛车时需要下载约 38.3 MB 的共用势能表。旧界面没有下载进度，整个流程使用固定的短超时；慢速下载容易表现为一直停在“准备赛车”。修复前本次测试苯环等待 110.6 秒后成功，未复现结构计算死锁。

现在显示实际下载字节进度，缓存经过长度和 SHA-256 校验的完整资源，并行加载资源及运行时。下载停滞会中止并重试，最终失败会展示原因并恢复编辑。另修复失败后立即重试读取旧界面状态、误判仍在准备的问题。模拟模型和准备计算保持不变。

## 实际验证

`node --test tests/browser-assets.test.mjs tests/browser-handoff.test.mjs wasm/tools/transport.test.mjs`：18 项通过。

`scripts/verify_handoff_progress.cjs`：本地与公开 HTTPS 各 4 项通过。公开站点使用 Chrome 实际浏览器和模拟 Worker：

| 场景 | 结果 |
|---|---|
| 冷缓存苯环送入赛道 | 127.43 秒准备成功，持续显示下载进度 |
| 关闭页面后重新打开，改用乙醇 | 15.05 秒准备成功；拦截所有势能表网络请求，实际请求数为零 |
| 注入资源下载 HTTP 503 | 展示下载失败原因，编辑器解除等待 |
| 解除网络故障后立即重试 | 119.90 秒准备成功 |

本地同一组测试全部通过。公开测试报告 `outputs/handoff-progress-public/report.json` 的 `passed` 为 true。测试脚本保留在仓库内，可重复运行。

首次下载仍受网络速度影响；这些用时不是不同版本的性能对照。缓存复用已验证，不代表全站离线运行。完整打包和资源清单校验通过。
