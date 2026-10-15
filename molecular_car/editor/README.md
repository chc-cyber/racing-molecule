# 编辑器独立几何与文件服务

此模块使用项目已有的 RDKit 2025.09.3、NumPy 1.26.4，不新增依赖或独立常驻进程，也不接触 `GameController`、赛道状态和 LAMMPS。HTTP 主服务只需将以下 POST 请求转交 `dispatch(path, data)`，然后将其返回对象输出为 JSON；`EditorError` 继承 `ValueError`，可由现有 400 错误处理捕获。导入及几何计算均不会触发 `/api/load`。

| 路由 | 请求 | 返回 |
|---|---|---|
| `/api/editor/import` | `{filename, content}`，MOL 或单记录 SDF | `{project, warnings, mapping}` |
| `/api/editor/export` | `{project, dimension:2或3, format:"mol"或"sdf"}` | `{filename, content, warnings, mapping}` |
| `/api/editor/geometry` | `{project, operation:"prepare"或"generate"或"hydrogenate"或"relax"或"layout2d", fixedIds:[], revision}` | `{project, revision, warnings, mapping, energy?, construction?, hydrogensAdded?, addedAtomIds?, addedBondIds?}` |
| `/api/editor/validate` | `{project, revision?}` | `{editable, eligible, errors, warnings, revision, atomCount, preparedAtomCount?, maxRaceAtoms, preparationChanges}` |

`revision` 不在服务端递增。请求返回后，前端必须检查版本仍为当前编辑版本，否则丢弃结果。撤销/重做也必须增加编辑版本；请不要仅比较旧项目内容。服务接收并返回完整快照，调用方通过一次 `model.applyAsync` 应用坐标与新增氢图，共用单一撤销栈。旧异步结果不能分批插入原子。若 `hydrogenate` 无需添加氢，返回的项目逐字段不变，模型不应新增空历史。

## 数据约定

`project.format` 为 `molecular-car-workshop-v1`，包含 `name`、`atoms`、`bonds`、`fragments`、`has3d` 和可选 `metadata`。原子序列就是 RDKit 的索引映射；`mapping.indexToAtomId` 从零编号。任何操作均不通过 SMILES 或图匹配猜测原子身份。原生保存应完整保存此对象。普通 MOL/SDF 只保留其格式可表示的信息；重新打开普通文件分配新的编辑器 ID，片段身份应通过原生项目保留。

原子至少包含 `id, element, charge, pos2:[x,y], pos3:[x,y,z]`。可选字段有 `isotope`、`explicitHs`、`noImplicit`、`radicalElectrons`、`chiralTag`、`generatedHydrogen`。新补入的真实 H 设置 `generatedHydrogen:true`，以区别用户手工绘制或文件导入的氢。几何服务不会删除任何已有 H；替代基连接时是否取代自动 H，须由明确的拓扑编辑交易按稳定 ID 处理。允许手性值：`CHI_UNSPECIFIED`、`CHI_TETRAHEDRAL_CW`、`CHI_TETRAHEDRAL_CCW`；CW/CCW 相对于保留的 RDKit 邻接顺序，不是 CIP 的 R/S。

键至少包含 `id, a, b, order`；order 为 1、2、3 或 1.5。可选 `stereo` 值为 `STEREONONE, STEREOANY, STEREOE, STEREOZ, STEREOCIS, STEREOTRANS`。`stereoAtoms` 为空数组或两端邻接参考原子的稳定 ID。`direction` 为 `NONE, BEGINWEDGE, BEGINDASH, ENDUPRIGHT, ENDDOWNRIGHT, EITHERDOUBLE, UNKNOWN`。文件导出后按相同数组索引重新读入检查绝对四面体立体信息和 E/Z，不通过图匹配。

三维 MOL 按 3D 输出。立体标记与坐标冲突时导出失败，并指出原子或键，不自动更改玩家标记。尚未标记的立体化学可能由外部软件从三维坐标推定；服务在这种情况下显式返回警告。二维导出的 E/Z 与布局冲突也会报错，应显式整理二维排版。原生项目允许保存这类待修正状态。

新绘制的楔线/虚楔按二维布局使用 `AssignChiralTypesFromBondDirs` 解析，所有几何操作和导出都执行此检查，包括已经有三维构象的项目。尚无显式手性标签时，几何返回对象会持久化对应标签；楔线与已有标签冲突时拒绝操作。重新整理二维布局时重新生成与同一绝对手性一致的楔线，可以改换楔线所在键及记录端点方向，但不会改变键连接关系。直接镜像二维布局造成手性矛盾不能静默导出；请修正布局或明确更新手性意图。

导入已有三维坐标时不优化、不补氢、不重建。二维文件导入后仅标记 `has3d:false`，进入三维由前端显式调用以下操作。保存、导出、返回二维本身不会暗中触发几何或补氢。

| 操作 | 已有坐标与图 | 新增真实氢 |
|---|---|---|
| `prepare` | 合理初始构建，不执行或检查 UFF；已有三维时仅构建新增原子，全部旧坐标精确保留；没有待处理区域就仅补氢 | 按化学价态补入实际 H 原子和单键 |
| `generate` | ETKDGv3 + UFF 重新生成三维；原有原子/键 ID、二维坐标保留；可固定已有三维原子 | 按化学价态补入实际 H 原子和单键，追加到原数组之后 |
| `hydrogenate` | 需要已有三维；所有原有原子二维/三维坐标逐值保留，不优化、不重建、不删氢 | 仅补缺失 H；没有缺失时项目不变，重复执行幂等 |
| `relax` | 整理三维，原子/键集合保持原样；固定点坐标保留 | 只在临时副本中补氢用于 UFF，不把 H 写入项目 |
| `layout2d` | 只整理二维和一致的楔线表示，三维及原子/键集合保留 | 不补氢 |

推荐切换流程：非空项目从二维进入三维时统一调用 `prepare`。只在响应版本仍匹配当前编辑版本时一次应用并切换视图；失败时保留可编辑草稿并显示原因。用户明确选择“重新生成三维”时仍调用 `generate`，其修改构象的影响应明确提示。

`prepare` 接收 `metadata.geometryPending={atomIds:[新增UUID],changedAtomIds:[拓扑改变的旧UUID],editedAtomIds:[明确修改的旧UUID]}`。所有已有原子的三维坐标，包括旧化学修改中心、用户绘制/文件导入的 H，都按原始数值精确保留。只有真正新增的重原子需要初始构建，之后补齐缺失 H。普通首次结构使用 ETKDGv3；已有构象的新片段单独生成并对齐到旧连接点，避免让整台旧分子参与距离几何计算。不要求 UFF 参数、不创建力场、不执行最小化；二维绘图单位不会直接当作 Å。

原生 C60 的 `metadata.ketcherGroups` 项必须包含 `kind:"c60"` 和按真实模板顺序映射的 60 个稳定 `atomIds`；适配器预先写入真实球壳 `pos3`。服务直接采用此 seed，保持内部几何，用局部连接方向和合理键长放置相邻新结构。二维展开图不会用作球壳 seed；C60–CH3 不再让整颗笼参与 ETKDG。多个固定连接点不能满足时明确失败，不偷偷修改旧原子。服务不依据分子式或对称图匹配重建 template 身份。

成功准备后清除 `metadata.geometryPending` 与 `geometryNotice`，重复切换只补缺失 H。请求输入始终不变；结果包含旧拓扑、两套坐标、SUP 信息和全部新增 H，可一次加入统一撤销记录。前端仍须检查编辑版本及交互版本，旧结果不得覆盖新编辑。

新增氢和键各自分配新的稳定 UUID，并通过 `addedAtomIds`、`addedBondIds` 返回；`hydrogensAdded` 是本次数量。新氢是原子图的实际顶点，能选择、删除、编辑、保存及导出，其键不依赖距离推断。原有原子（含已有氢）、原有键都按原索引保留 ID。新增氢加入其母原子所在片段的 `atomIds`，片段 ID、端口、刚体标记与其他信息保留，后续整体移动会包含这些氢。

`explicitHs` 在 MOL/原生原子属性里有时表示尚未成为图顶点的氢数，例如 `[C@H]`、`[NH3+]`。将这些氢物化为原子后，相应计数归零，避免双计氢；原子形式电荷和同位素信息保留。四面体绝对立体语义与 E/Z 保留；为了让新增第四邻居后的二维表达一致，必要时把楔线改画到新 H 键，原有坐标和键 ID 不变。几何响应仍会拒绝不一致的立体标记，不静默翻转手性。

## 几何与试跑边界

显式 `generate` 使用固定种子的 ETKDGv3，接着执行 UFF；自动 `prepare` 只做合理初始构建和补氢，没有 UFF 能量/收敛字段。显式 `relax` 使用已有构象，支持 `fixedIds`（UFF `AddFixedPoint`），固定点精确保持输入坐标。`hydrogenate` 不运行 UFF；新氢采用 RDKit 的局部几何放置，可能需要用户之后显式整理以消除拥挤。未连接片段仍可独立整理；片段间非键作用不计入此次编辑器整理。显式 UFF 1000 步内未收敛会明确返回警告和 `energy.converged:false`，不是成功收敛提示。拖动过程中不应调用此路由。

`validate` 调用现有 `prepare_molecule` 在副本上执行真实准入与参数检查，然后检查赛道 45 Å 空间跨度。不会启动比赛或替换编辑器构象。准入复用单连通图、补氢后 256 原子、闭壳层、可支持的常规键、元素极化率表、UFF 和 Gasteiger、300 步 UFF 收敛检查。可编辑/保存与可试跑严格分开。当前编辑器序列化上限为 2000 原子、6000 键，交互几何计算暂限 600 个显式原子。

保留玩家比赛构象需要主引擎未来提供显式准备策略，例如 `geometryPolicy:"preserve"`、可选固定原子列表，以及完整的 `editorAtomId`→新补氢后比赛索引映射；否则现有后端仍会补氢、松弛和平移。此模块没有修改准备语义，也没有向现有 `/api/load` 发送未约定字段。

## 片段数据

运行 `python -m molecular_car.editor.fragments` 可重建本地 `frontend/editor/fragments.json`。苯环为 6 原子 6 键，炔基为 2 原子 1 个三键；C60 是数学构造的截角二十面体分子图，60 个真实碳原子、90 条键、30 条 6:6 双键，每个碳度数 3、键级和 4。所有 C60 接点初始满价，必须明确修改邻接键级后再连接。ports.maxOrder=1 表示允许连接的最高键级，不代表当前剩余价态；前端应动态检查。`rigid:true` 要求默认按整个片段操作。

## 限制

- 查询原子/键、增强立体组及非 SUP 的 S-groups 无法无损表达，导入时明确拒绝。SUP 使用真实内部原子及稳定 ID 映射保存 `metadata.ketcherGroups` 的名称、展开状态和连接位点。带 SUP 的 MOL 导出使用 V3000，避免 V2000 SAP 的兼容问题；普通文件重导入分配新 UUID，文件不会保存原生 template 身份。
- 多记录 SDF 明确拒绝；当前没有记录选择器。
- 非四面体/轴手性未实现；原生文件适配器不会静默丢弃未知标记。
- 由于 MOL V2000 坐标写入四位小数，往返允许每分量 0.000051 Å 的舍入误差；原生 JSON 保存完整精度。
- 对价态错误的未完成图允许导入与 MOL 导出，返回警告。无法认证此类图的立体导出，应保留原生项目。
- SDF 属性表不作为原生编辑字段处理；目前只提取分子结构。

## 验证与来源

运行 `.\.venv\Scripts\python.exe -E -B scripts\verify_editor_backend.py`；报告输出到 `outputs/editor-check/backend-report.json`。也可直接运行 `.\.venv\Scripts\python.exe -E -B -m unittest tests.test_editor_service -v`。测试涵盖坐标保留、二面角编辑后二维/三维往返、固定点优化、立体与电荷/同位素往返、冲突拒绝、C60 拓扑及 UFF/准入、多记录拒绝、可编辑与可试跑分离、版本透传和输入不变。还覆盖新绘楔/虚楔在无/已有三维状态下生成相反且正确的 R/S、重新二维排版保留手性，以及镜像二维布局后的冲突拒绝。真实氢测试验证苯环生成 12 个实际原子、每个新 H 只有一条真实单键、片段归属、旧三维/二维坐标精确保留、已有 H 不删除、重复补氢幂等、R/S 与 E/Z 往返，以及错误价态转换失败仍保留草稿。

自动初始构建测试见 `tests/test_editor_initial_construction.py`，用拒绝调用的 UFF/MMFF mock 验证所有自动路径不触及力场，并在 C60–CH3 用例同时禁止 EmbedMolecule。实际无 UFF 参数的 Fe²⁺ 与 He 也可 prepare。覆盖仅新片段＋单 cap 的嵌入、原 60 笼及所有旧氢坐标精确保留、苯/链/新 R/S 官能团、重复往返、模板刚体距离、错误价态及多固定连接拒绝。与原服务、导入和 SUP 检查合计 51 项通过。

主页面当前初始构建回归：独立 8786 服务已启动后，从项目根运行 `node scripts/verify_editor_seed_browser.cjs http://127.0.0.1:8786`。测试强制阻止 `/api/load`，以真实原生模板和鼠标成键验证氢标签、新苯尺度、旧构象保护、TPB 二面角、折叠 C60 球壳往返、C60–CH3 首次构建及已有球壳后添加、共同撤销和旧异步响应丢弃。报告位于 `outputs/editor-seed-check/browser-report.json`，包含当前 source SHA-256 和实际 prepare 响应的 `construction.optimized:false`；历史自动 UFF 报告不能代替当前运行结果。

`frontend/editor/examples/chiral.mol` 是从 `C[C@H](O)C(=O)O` 在本地通过相同固定种子的 ETKDGv3/UFF 生成的三维乳酸测试结构，含一个真实四面体手性中心，用于界面演示与往返核对。

RDKit 保留其安装包中的 BSD 3-Clause 许可证及版权声明；本模块未复制或重新分发 RDKit 源码/二进制。版本与来源沿用项目 `THIRD_PARTY.md`。[RDKit 官方源码及许可证](https://github.com/rdkit/rdkit)、[官方 MOL 立体解释](https://www.rdkit.org/docs/RDKit_Book.html#interpretation-of-the-2d-3d-flag-and-stereochemistry)、[固定点力场接口](https://www.rdkit.org/docs/source/rdkit.ForceField.rdForceField.html)。片段坐标和键图为本模块本地生成，无外部模型文件版权依赖。
