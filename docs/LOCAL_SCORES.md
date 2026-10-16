# 本地成绩数据接口

第一版仅保存浏览器本地成绩，不连接全球排行榜服务。入口为赛道“排行榜”和完赛弹窗。

## 记录与排序

- 只记录完成 20 nm 的标准赛事；自定义赛事和超时未完赛不入榜。
- 按模拟时间 `timeFs` 升序；准备、下载和动画用时不计入成绩。
- 同一结构、同一规则及计算版本保留最佳成绩，改名不生成新成绩。
- 使用 IndexedDB，刷新和正常重开后保留。清除网站数据会删除本地成绩；可以导出单条成绩 JSON。
- 第一版结构指纹基于 MOL 原子与键记录，V2000 忽略标题和坐标、保留电荷及立体字段；不是规范图同构算法。重新编号的同一分子可能作为另一条记录，V3000 暂保留坐标参与指纹。

## 版本化接口

`frontend/local-scores.mjs` 提供 `createScore(state, replay, engine)`、`LocalScores.list()`、`LocalScores.save(record)` 和 `scorePayload(record)`。

`racing-molecule-score-v1` 包含：

- `id`、`createdAt`、`name`、`formula`、`verification: "local"`。
- `ruleset`、`engine`、`structureHash`：规则、计算资源及结构指纹。
- `timeFs`、`distanceAngstrom`：最终模拟时间及净位移。
- `competition`：赛事版本、参数、随机种子和请求摆放方式。
- `input`：原始 MOL/SDF、文件名、可选工坊项目 ID 和名称。
- `molecule`、`initialFrame`、`finalFrame`、`runtime`、`startPose`：实际原子和键、准备后的初态、终态及运行元信息。

网页计算版本指纹包括 WASM 下载描述、模型描述及 Fourier 表描述。不同计算资源或赛事参数分榜；摆放姿态属于选手输入，不拆分榜单。原始双精度坐标按 JSON 数值保留。

后续全球服务可消费 `scorePayload()`，但必须在服务端验证版本、结构和参数并复跑。当前没有上传、账户或全球排名 API。

## 命名与摆放

工坊工程增加稳定 `projectId`，网页通过同源父页面同步对应工程与赛道的名称。新建工程使用新的 ID，避免改动另一辆已入场赛车的名称。输入框支持 Enter 确认、失焦确认、Esc 取消。改名不会触发模拟准备或清空轨迹。

标准和自定义赛事都开放现有整体朝向工具（自动、文件朝向、手动旋转、参考原子）。比赛开始和暂停后锁定摆放；重置后可以重新摆放。当前不新增自由平移或高度控制，仍由现有准备流程定位并松弛赛车。

## 验证

- `node scripts/verify_local_scores.cjs`：工坊画苯、双向改名、标准赛事旋转 45°、暂停改名、真实完赛、自动入榜、导出、刷新恢复及再次试跑坐标一致。
- 设置 `WEB_TEST_URL` 可对 GitHub Pages 运行相同流程。
- `tests/local-scores.test.mjs` 覆盖准入、结构指纹、版本分榜、最佳成绩和独立快照。
- `wasm/tools/transport.test.mjs` 覆盖改名不重建运行时、标准姿态和暂停锁定。
