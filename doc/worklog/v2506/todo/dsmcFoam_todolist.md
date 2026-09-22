# dsmcFoam+ Todo List

## Repart + KM remap for replicated-mesh DLB

### 方案来源

- Qiu et al. 2022 的耦合 DSMC/PIC 并行 DLB 方案采用两段式流程：运行中按 cell weight 调用 METIS/ParMETIS 重新划分，然后用 Kuhn-Munkres / Hungarian matching 将新分区重新映射到旧 rank，以降低重平衡后的粒子迁移开销。
- 本项目文献精读报告已记录该思路：`doc/literature/dsmc_parallel_optimization_close_read_report_20260624.md` 中指出，KM remap 可降低 rebalance overhead，尤其在迁移通信较重时更有价值。
- 当前 OFv1706 `hyStrath_dlb` 代码已经具备前半段能力：`src/lagrangian/dsmc/replicatedMesh/dsmcReplicatedMesh.C::reassignByParMetisAdaptiveRepart()` 使用 `ParMETIS_V3_AdaptiveRepart()` 生成新的 `part`，随后 `MPI_Allgatherv()` 得到 `fullPart`，并直接以 `fullPart[cellI]` 更新 `cellOwner_[cellI]`。

### 核心判断

`ParMETIS AdaptiveRepart` 和 `KM remap` 不重复。

- `AdaptiveRepart` 负责决定新的图划分，即每个 cell 属于哪个新 part。
- `KM remap` 负责决定新 part label 到旧 MPI rank label 的一一映射。
- 如果不做 remap，ParMETIS 返回的新 part 编号可能发生全局置换。即使几何分区本质相近，仅因编号对调也会被当前代码统计为大量 `cells changed`，随后触发不必要的 `migrateParticlesByCellOwner()`。
- 因此完整流程应为：`old cellOwner_` -> `ParMETIS AdaptiveRepart fullPart` -> `KM/greedy part-to-rank remap` -> `remapped fullPart` -> 更新 `cellOwner_` -> 必要时迁移粒子。

### 收益预期

- 降低 DLB 后的 `nChanged`，减少 owner label 置换造成的伪迁移。
- 降低 `autoRebalanceMigrationWallTime_` 和 `totalCellsChanged_`，尤其对 repeated DLB、rank 数较多、粒子迁移较贵的 case 更明显。
- 提高 DLB 稳定性：相近的连续两次 repartition 会倾向保持同一 rank ownership，减少分区编号抖动。
- 对 ParMETIS 本体时间影响很小。匹配规模是 `nProcs x nProcs`，而不是 `nCells x nCells`；在当前 MPI8、MPI16、MPI64 这类规模下，KM 或 greedy matching 的成本通常远小于一次真实粒子迁移。
- 如果当前 ParMETIS 已经高度稳定，收益可能主要体现在诊断指标上；如果日志出现类似 `Phase C ParMETIS AdaptiveRepart: 0.6-0.9 cells changed ratio`，则该方案可能直接减少大量不必要迁移。

### 可行性分析

- 代码接入点清晰：应放在 `MPI_Allgatherv(part -> fullPart)` 之后、`Update cellOwner_` 之前。
- 所需数据已经齐全：旧 `cellOwner_`、新 `fullPart`、`nProcs_`、全局 cell 数、可选 `globalCellParticles`。
- 代价矩阵可按 overlap 或迁移代价构造：
  - 基础版：`cost[newPart][oldRank] = - count(cells where fullPart == newPart && old cellOwner_ == oldRank)`，最大化旧新重叠。
  - 粒子权重版：用 `globalCellParticles[cellI]` 替代 cell 数，优先减少实际 parcel 迁移。
  - 混合版：`cell overlap + particle overlap`，先做基础版验证，再考虑加权。
- 可先实现 greedy overlap matching，复杂度低、代码短、足够验证收益；若 greedy 效果不足，再实现完整 Hungarian/Kuhn-Munkres。
- 风险主要在映射方向：需要明确 `remap[newPart] = oldRank`，然后用 `newOwner = remap[fullPart[cellI]]` 更新 `cellOwner_`。方向写反会破坏 owner 分布。
- 该方案不改变 ParMETIS 分区质量，只改变 part label 到 rank label 的绑定；因此正确性风险低，但必须验证每个 remap 后 owner 仍在 `[0, nProcs_)` 且每个 rank 只被匹配一次。

### 实施方案

1. 增加运行时开关：

   - `replicatedMeshDLBRemap true/false`，默认建议先设为 `false`，用于 A/B 测试。
   - `replicatedMeshDLBRemapMode greedyOverlap|hungarianOverlap|greedyParticle|hungarianParticle`，第一阶段只实现 `greedyOverlap` 即可。
2. 在 `reassignByParMetisAdaptiveRepart()` 中保留旧 owner：

   - 在更新 `cellOwner_` 前复制 `oldOwner = cellOwner_`。
   - 用 `fullPart` 与 `oldOwner` 构造 `nProcs_ x nProcs_` overlap 矩阵。
3. 实现 remap：

   - 基础 greedy：每轮选择当前最大 overlap 的 `(newPart, oldRank)`，锁定这一对，直到所有 part 匹配完成。
   - 对未出现 cell 的 part 或 tie case，按 rank id 做确定性 fallback，保证结果可复现。
   - 输出 `remap[newPart] = oldRank`。
4. 应用 remap：

   - 在 changed-cell 统计前，将 `fullPart[cellI]` 转为 `remappedPart[cellI] = remap[fullPart[cellI]]`。
   - 用 `remappedPart` 与旧 `cellOwner_` 比较，更新 `cellOwner_` 并计算 `nChanged`。
5. 增加诊断输出：

   - remap 前 changed cells。
   - remap 后 changed cells。
   - remap saved cells = before - after。
   - overlap 总收益或 estimated saved particle count。
   - 每次 DLB 打印 `Phase C remap: mode=..., changedBefore=..., changedAfter=...`。
6. 验证流程：

   - 用同一个 case 分别跑 `replicatedMeshDLBRemap false` 和 `true`。
   - 对比 `Phase C ParMETIS AdaptiveRepart` 的 changed ratio、`autoRebalanceMigrationWallTime_`、`Phase C auto DLB wall max`、总 wall time。
   - 检查粒子数守恒、`Total Iterations`、`End main`、采样输出是否正常。
   - 先在 MPI8 做 smoke test，再扩展到 MPI16/MPI64 或混合 MPI+OMP。

### 阶段性验收标准

- 功能验收：remap 开启后程序稳定跑完，owner 范围合法，无粒子丢失。
- 性能验收：至少在一个 DLB 触发 case 中，`changedAfter < changedBefore`，并且 migration wall time 不上升。
- 诊断验收：日志能清楚区分 ParMETIS 原始 changed cells、KM/greedy remap 后 changed cells、实际 migration 时间。
- 决策验收：若 greedy overlap 已能显著降低 changed cells，暂不实现 Hungarian；若 greedy 在 rank 数较大或 tie 较多时效果不稳定，再补完整 KM。

## Owner-gated inflow and boundary parcel generation

### 问题来源

- 当前 `dsmcCloud::evolve()` 在 move 之前调用 `controllers_.controlBeforeMove()` 和 `boundaries_.controlBeforeMove()`，随后才执行 replicated-mesh 的初始分发或周期性 `migrateParticlesByCellOwner()`。
- `dsmcCloud::addNewParcel()` 和 `addNewStuckParcel()` 已经有 replicated-mesh owner 过滤：若 `cellOwner()[cellI] != myRank()`，直接 `return`，因此当前不是“所有 rank 真正生成 parcel 再迁移”的模式。
- 但多个 inlet/general boundary 模型仍在调用 `addNewParcel()` 之前执行完整 face/species 循环、随机粒子数采样、位置采样、速度采样、转动/振动/电子能级采样。非 owner rank 最终不会插入 parcel，但前面的 CPU 和 RNG 开销已经发生。
- 典型位置包括：
  - `src/lagrangian/dsmc/boundaries/derived/generalBoundaries/dsmcFreeStreamInflowPatch/dsmcFreeStreamInflowPatch.C`
  - `src/lagrangian/dsmc/boundaries/derived/generalBoundaries/dsmcFreeStreamInflowFieldPatch/dsmcFreeStreamInflowFieldPatch.C`
  - `src/lagrangian/dsmc/boundaries/derived/generalBoundaries/dsmcChapmanEnskogFreeStreamInflowPatch/dsmcChapmanEnskogFreeStreamInflowPatch.C`
  - `src/lagrangian/dsmc/boundaries/derived/generalBoundaries/dsmcMassFlowRateInlet/dsmcMassFlowRateInlet.C`
  - `src/lagrangian/dsmc/boundaries/derived/generalBoundaries/dsmcLiouFangPressureInlet/dsmcLiouFangPressureInlet.C`
  - `src/lagrangian/dsmc/boundaries/derived/generalBoundaries/dsmcWangPressureInlet/dsmcWangPressureInlet.C`

### 核心判断

这项优化不是减少已生成 parcel 的 MPI 迁移，而是减少 replicated-mesh 模式下非 owner rank 的无效 boundary/inflow 生成计算。

- 已有 `addNewParcel()` owner gate 保护正确性，避免非 owner rank 插入真实 parcel。
- 剩余低效发生在 owner gate 之前：每个 rank 都可能为不属于自己的入口 face 消耗随机数、构造 tet 面三角形、采样速度和能级。
- 如果将 owner gate 前移到 boundary face 循环早期，可跳过大部分无效计算。
- 需要谨慎处理 `accumulatedParcelsToInsert_` 这类 per-face 累计余数，否则 DLB 改变 `cellOwner_` 后，新 owner 可能继承不到合理的 fractional remainder，影响入口通量统计连续性。

### 收益预期

- 对入口面较多、inflow 粒子生成密集、MPI rank 数较高的 case，可降低 `boundaries_.controlBeforeMove()` 相关开销。
- 可减少非 owner rank 的 RNG 消耗和复杂采样开销，改善 move 前阶段的 rank-wall 不均衡。
- 不增加通信，不改变 `migrateParticlesByCellOwner()` 路径。
- 对当前 `addNewParcel()` 已经过滤非 owner 插入的正确性逻辑影响较小；主要风险集中在累计余数和随机序列变化。

### 可行性分析

- 接入点清晰：在每个 general boundary 的 face 循环中，拿到 `cellI = cells_[f]` 后即可判断该 face 是否由本 rank owner 负责。
- 需要在 `dsmcCloud` 或 `dsmcGeneralBoundary` 层提供轻量 helper，例如：
  - `cloud_.replicatedMeshActive()`
  - `cloud_.replicatedMesh().cellOwner()[cellI]`
  - `cloud_.replicatedMesh().myRank()`
- 不建议只依赖 `addNewParcel()` 的 late return，因为它无法避免前置采样开销。
- 不建议第一版改变通量累计公式或跨 rank 同步 accumulator；先做保守 owner-only 生成路径，并用粒子数、质量通量和采样场验证。

### 实施方案

1. 增加一个小型 helper：

   - 可放在 `dsmcCloud` 中，例如 `bool ownsCellForInjection(label cellI) const`。
   - 非 replicated-mesh 模式返回 `true`。
   - replicated-mesh 模式返回 `cellOwner()[cellI] == myRank()`。
2. 在主要入口模型中前移 owner gate：

   - 在 face loop 中读取 `cellI = cells_[f]` 后立即判断。
   - 对非 owner face，跳过 tet 分解、随机位置、速度、能级采样和 `addNewParcel()`。
3. 保守处理 `accumulatedParcelsToInsert_`：

   - 第一版建议仍允许所有 rank 更新 deterministic accumulator 增量，避免 DLB 后 owner 切换导致入口通量余数突变。
   - 将随机取整、Poisson 采样、`faceAccumulator -= nInserted` 放到 owner-only 分支内，避免非 owner rank 消耗随机数。
   - 若发现 owner 切换后的 accumulator 状态不一致，再考虑在 DLB 后按 owner 广播/重置 boundary accumulator。
4. 增加诊断：

   - 每个 boundary 或全局累计 `injectionFacesSkippedNonOwner`。
   - `injectionCandidatesSampled`、`injectionParcelsInserted`。
   - 可选统计 owner gate 前后 `boundaries_.controlBeforeMove()` wall time。
5. 验证流程：

   - 先选单一 inflow case 做 MPI8 no-write smoke test。
   - 对比优化前后总 parcel 数、每步新增 parcel 数、入口通量、`rhoN`、`dsmcNMean`、壁面热流。
   - 对比 `move pre-control/partition` 或 boundary 相关 profile 时间。
   - 再扩展到包含 pressure inlet/outlet calculated molar fraction 的复杂 case。

### 阶段性验收标准

- 功能验收：replicated-mesh 模式下非 owner face 不执行昂贵随机采样，真实插入 parcel 仍只发生在 owner rank。
- 正确性验收：总粒子数、入口质量流率和关键采样场与基线在随机误差范围内一致。
- 性能验收：入口生成阶段 wall time 或 move 前控制阶段时间下降；migration 时间不应因该改动上升。
- 诊断验收：日志能报告跳过的非 owner face 数和实际插入粒子数，用于判断收益上限。

## Tier-2: 构造期按 cellOwner 过滤初始 cloud（读取时过滤）

### 方案来源

- 2026-09-14 800w-rcj 崩溃调查的两个直接产物：
  - `detail_mix/thread_offset_int32_overflow_fix_20260914.md`：构造期全量
    cloud 导致线程偏移 int32 溢出（已修，但 1-rank 大 N 配置仍依赖该修复）；
  - `sacct -j 4616685` 实测构造期 **MaxRSS 350.56 GB/rank**（75.7M 初始
    粒子全量物化 + ASCII 读入瞬态缓冲），远超理论估算 ~80 GB——构造期
    峰值内存是大 case 的第一短板。
- 架构讨论结论（2026-09-14）：replicated mesh 启动设计为"每 rank 全量
  构造、首步 evolve 才删除非 owned"
  （`distributeInitialParticles()`，dsmcReplicatedMesh.C:3367，纯删除零
  通信）。Tier-1（构造器内提前过滤）不省内存；Tier-2（读取时过滤）才
  是内存收益的来源。

### 核心判断

构造期各 rank 只需要最终 owned 的 N/R 个粒子（以及瞬态读取缓冲），
全量物化 N 个 parcel 对象纯属浪费。读取时过滤把构造期峰值从

```
现状  ≈ mesh 22G + 750B×N(全量 parcel) + ASCII 瞬态(~350G 的主要成分)
Tier2 ≈ mesh 22G + 750B×N/R(owned) + 分块流式瞬态(每块常数)
```

对 256 GB 单节点、R=4：N_max 从 ~73M（800w-rcj 恰好卡死在这条线上）
提升到 ~192M（×2.6）。对理想形态（1 节点 1 rank）无收益——它不解决
mesh 复制地板和 750B/粒子稳态成本，那些属于 P1 粒子数据布局。

前置事实（已核实）：

1. `dsmcParcel::readFields`（dsmcParcelIO.C:133）bulk 模式：positions
   流式建 parcel，U/typeId/ERot/vibLevel 等是整场 IOField 读入后按索引
   赋值。若在 append 时就丢弃非 owned parcel，`c.size() < N` 会触发
   `checkFieldIOobject` 尺寸失配 FATAL——**这是 Tier-2 的主要工程量**。
2. `replicatedMesh_->initialize()`（dsmcReplicatedMesh.C:505）依赖闭包
   干净：仅 MPI + controlDict + mesh（METIS 按 cell 数分解），不依赖
   粒子/occupancy/boundaries，可前移到 readFields 之前。
3. `noTimeCounter::initialConfiguration()` 为空函数；
   `buildCollisionSelectionRemainderFromScratch()` 只用 rndGen+nCells
   ——构造器中位于 buildCellOccupancy 之后的消费者对"过滤后建
   occupancy"无感。
4. processor checkpoint restart 路径天然兼容（cloud 本来就是 owned
   子集，过滤器为 no-op）。

### 实施方案

1. **dsmcCloud 构造器时序重排**（fresh start 分支）：

   ```
   readFields（改造为 owned-filtered 读取，见 2）
   → replicatedMesh_->initialize()        // 前移：先有 cellOwner_
   → distributeInitialParticles()         // 变 no-op 兜底（保留）
   → buildCellOccupancyFromScratch()      // 只 bin N/R
   → 其余初始化不变
   ```

   首步 evolve 的 `migrationCalls()==0` 分支保留作兜底（对旧格式初始
   目录/无过滤路径仍有效）。
2. **owned-filtered readFields**（核心改动，建议分两步走）：

   - 第一步（保守）：bulk 数组照读（瞬态全量），构造 parcel 对象时跳过
     `cellOwner_[cellI] != myRank_` 的索引。`checkFieldIOobject` 尺寸
     检查改为"每个字段 size == N 全量"（读的是全量数组），赋值循环里
     只物化 owned——消除 750B×N 的 parcel 对象层，保留 IOField 瞬态。
     第二步（完全体）：positions 改分块流式（每块 ~1M 行），块内先解析
     cellI（需要 cell 查找或预存的 cell 映射），只物化 owned，其余字段
     按"owned 索引列表"稀疏抽取——消除全部全量瞬态。
   - 注意 cellOwner_ 前置可用性：`initialize()` 前移后
     `computeCellOwnerScotch()` 产出与现状完全一致（同 METIS 输入），
     RNG 流不变（过滤不消耗随机数）。
3. **诊断输出**：

   ```
   Replicated mesh: rank R read N particles, materialized M owned
       (skipped S non-owned during read), read transient peak = X MB
   ```

   配合 slurm 采样器（`/proc/<pid>/status` VmHWM，30s 间隔）记录逐 rank
   峰值，形成 before/after 内存曲线。
4. **验证流程**：

   - 单测口径：`zb-cylinder-react-validate/mpi4omp2` 300 步，对比 final
     particles/Total energy/碰撞计数与历史基线（1,957,830 /
     1.9622e-3，`sigmaTcRMax_value_ctor_zero_collision_fix_20260914.md`
     的表）；
   - `kept/deleted` 日志改为 read-time skipped 统计，全局守恒断言：
     Σ rank materialized == N；
   - checkpoint restart 冒烟（过滤器 no-op 路径）；
   - 超算 800w-rcj 复跑：`sacct -j <id> --units=G -o MaxRSS` 对比
     350.56 GB 基线，目标是 < 120 GB（R=12）；
   - DLB 触发 case（palphd 3.3.1）跑通首步 + 首次 rebalance。

### 收益预期（256 GB 单节点，800w 级 case）

| 部署                    | 现状构造峰值/rank |  Tier-2 后 |   N_max 提升 |
| ----------------------- | ----------------: | ---------: | -----------: |
| R=4                     | ~350G（实测口径） | < 60G 量级 | 73M → ~192M |
| R=12（800w-rcj 现配置） |              350G | < 35G 量级 |     余量充足 |
| R=1（理想形态）         |              不变 |       不变 |       无收益 |

（第一步改造后 IOField 瞬态仍在，预计 ~150G 量级；完全体才是上表数值。）

### 边界与不做的事

- 不动 mesh 复制（22G/rank 地板）与稳态 750B/粒子——属 P1 粒子数据
  布局（SoA/线格式瘦身）领地；
- 不改 DLB/迁移语义：过滤只发生在 fresh-start 读取，时间循环内一切
  照旧；
- int64 偏移修复独立保留（1-rank 大 N 仍需要）；
- `decomposePar`/标准 decomposed 路径不受影响（无 replicated mesh）。

### 阶段性验收标准

- 功能验收：所有验证 case `Total Iterations`/`End main`/stuck=0/碰撞
  非零（含 `collisions` 门）。
- 正确性验收：末态 particles/Total energy 与未过滤基线在 RNG 噪声带内
  一致（参考 6/25 消融散布 0.04%/0.035%）。
- 内存验收：超算 800w-rcj MaxRSS 从 350G 降至目标值（第一步 < 200G，
  完全体 < 120G @ R=12）。
- 性能验收：构造 wall time 不劣化（bin N/R 应更快）；时间循环各阶段
  计时与基线一致。
- 诊断验收：read-time skipped 统计 + VmHWM 采样曲线归档。

### 实施状态更新（2026-09-16）

Tier-2 已实现并投入实战，**并行为主路径**（保守串行版降为小文件回退）：
OMP 分块解析 + 条目自带 cellI 精确预过滤（免 findCell，93% 非owned
零构造成本）+ mmap 共享页缓存 + 进度打印。实现与坑的完整记录见
`todo/tier2_parallel_filtered_read_20260916.md`。

实战（230wcell-1bparticle，116M parcels，m24o16）：串行版 3.5h+ 无输出
→ 并行版 ~10 分钟完成解析；每 rank RSS 87GB（预测）→ **实测 15.6GB**；
24 ranks 总 RSS 205GB（无过滤路径 ≈ 2.2TB，物理不可行）。
**过滤开关（replicatedMeshFilterInitialRead）在该规模为必需项。**

剩余：bulk 字段读取串行 I/O（17GB/rank，下一瓶颈）；时间循环运行观察中。

### 优化 3 设计文档已出（2026-09-17）

processor 分布式初始 cloud + OMP 填充的完整设计：
`todo/distributed_initial_cloud_omp_fill_design_20260917.md`。
要点：两级并行（MPI ÷24 × OMP ÷16）、两形态（3a 分布式写 / 3b in-solver
填充，推荐后者）、per-cell RNG 确定性（与 rank/线程无关）、sigmaTcRMax
内算、守恒按构造精确。目标：230w/1.16 亿初始化 6.5h → <5 min。
工作量 ~1-2 周（P1-P4），风险与验证门齐备。决策：1 立即做（独立收益），
3 排期做（正解），2 被 3 作废可跳过。
