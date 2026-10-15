# GitHub Pages 线上修复验收 · 2026-10-15

线上入口：https://chc-cyber.github.io/racing-molecule/

- 修复源码：`c4a1fcf`。
- Pages 发布提交：`65c550b6c01243a05bc17bc7f6b2f0faef39c8d0`。
- [Pages 部署任务](https://github.com/chc-cyber/racing-molecule/actions/runs/38704160827)：success。
- 使用全新 Chrome 无头浏览器环境，直接访问公开 HTTPS 工坊；未以本地文件替代线上资源。

| 验收 | 结果 |
|---|---|
| [6]螺烯二维输入 → 三维 | 成功，26 C + 16 H，532 ms |
| 螺烯往返二维/三维、撤销与重做 | 完整工程、坐标、原子 ID 比较通过 |
| Ketcher 原生 C60 按钮放置 → JSmol | 球壳坐标与模板一致，217 ms |
| 中性亚胺氮上实际拖出一根单键 | 原生红色价态提示，隐式 H = 0 |
| 氮价态错误撤销/重做 | 提示随结构正确恢复 |
| 将该氮改为 +1 电荷 | 错误提示消失，可正常生成三维 |

页面错误 0，比赛 API 请求 0。上述计时是本机单次页面操作耗时，不含首次组件下载。自动构建未启用 UFF；本次未重新运行完整比赛。

## 下载文件校验

线上清单与本地构建一致，并对浏览器实际收到的以下文件重新计算 SHA-256：

- `editor/vendor/ketcher-c60/index.html`：`c13fb326f7f79740c309db15bc5a103384e412fa41fe28f89c0d0e9f50a63e89`
- `browser/runtime/molecular-car.wasm`：`cf5c57987dd939c0612cdba0af3adbe40fd813d77490e85ee02a3572c17d4a39`

## 复验

```powershell
$env:WEB_TEST_URL='https://chc-cyber.github.io/racing-molecule'
node scripts/verify_workshop_regressions.cjs
```

记录和实际截图输出到 `outputs/workshop-bugfix-public/`。脚本从公开地址加载 Worker、WASM、Ketcher 和 JSmol，测试数据仅使用测试浏览器内存。
