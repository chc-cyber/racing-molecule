# 赛道后台预加载验证

源码：`0fa1d9e`。Pages：`e93f45cdd4e16784cfabde68b1a1b91cf4148a34`。
部署 [38741814268](https://github.com/chc-cyber/racing-molecule/actions/runs/38741814268) 成功。

页面空闲时启动预加载；直接进入工坊时等待 Ketcher 和三维视图就绪。独立 Worker 逐个低优先级下载、解压、校验势能表，不创建模拟引擎。后台与比赛 Worker 使用同源 Web Locks 和校验缓存，正在下载的表只下载一次。缺少 Web Locks 或缓存不可用时跳过预加载，保留点击准备的正常路径。

## 验证

19 项资源、交接和传输回归测试通过。`scripts/verify_preload.cjs` 在本地与公开 HTTPS 各通过 4 项实际 Chrome 浏览器测试。

公开地址：https://chc-cyber.github.io/racing-molecule/

| 场景 | 线上结果 |
|---|---|
| 后台下载期间通过 Ketcher 模板绘制苯环 | 0.410 秒，原子数量正确；赛道 iframe 尚未创建 |
| 下载期间撤销、重做、再次撤销 | 原子数量按预期恢复 |
| 下载刚开始就送入赛道 | 106.638 秒准备成功；4 份势能表各请求一次 |
| 全新上下文等待预加载完成后送入赛道 | 4.376 秒准备成功；4 份势能表各请求一次，送入时没有重复下载 |

页面异常为零。本地对应准备时间为 2.686 秒与 1.873 秒。公网测试报告保存在 `outputs/preload-public/report.json`，本地为 `outputs/preload/report.json`；测试脚本已提交以供复验。

静态包 2119 个文件的清单、源码来源、资源校验及 ZIP CRC 均通过。ZIP SHA-256：`ad5204689fb76d0381c8de111d9c79bf186951784e0b86e6f033268481e16960`。

预加载将下载等待提前到绘图期间，不减少首次资源总量。上述为本次桌面 Chrome 实测，不是所有网络或设备的速度保证。未新增弹窗或操作按钮，未改变模拟模型、比赛和轨迹。
