# dsmcFoam+ Replicated Mesh DLB 优化工作总结（2026-05-18 更新）

## 1. 工作目标

在 dsmcFoam+ 中实现 Replicated Mesh + Dynamic Load Balancing 方案，使 MPI+OpenMP 混合并行在保持 DLB 能力的同时，性能接近纯 OpenMP。

核心思路：每个 MPI rank 持有完整网格，通过 `cellOwner_[]` 划分计算责任，ParMETIS AdaptiveRepart 动态调整 cell 归属。

## 2. 最终性能对比（300步，8核，cylinder N2 noreact，profileDetail false）

### 2.1 完整 8 组对比（含 collisionFastRng）

| 配置 | 无 FastRng | 有 FastRng | FastRng 收益 |
|------|-----------|-----------|-------------|
| **纯 OMP8** | **82.5s** | **79.2s** | -4% |
| **MPI2×OMP4 replicated mesh** | **101.3s** | **96.7s** | -5% |
| **MPI4×OMP2 + DLB(K=64)** | **104.7s** | **99.2s** | -5% |
| **MPI8 + DLB** | **142.5s** | **133.7s** | -6% |

### 2.2 vs 旧方案对比

| 配置 | 旧方案 | 新方案(无FastRng) | 新方案(FastRng) | 改善 |
|------|--------|-------------------|-----------------|------|
| 纯 OMP8 | 175s (V-2.1) | 82.5s | 79.2s | -55% |
| MPI2×OMP4 | 118s (stagebuf2) | 101.3s | 96.7s | -18% |
| MPI8 | 326s (decompose) | 142.5s | 133.7s | -59% |

## 3. 关键优化及效果

### 3.1 修复纯 OMP8 性能退化（146s → 85s）

**问题**：`useMoveParticlePartition = false` 被无条件设置，导致纯 OMP 走慢速链表 extract 路径。

**修复**：恢复条件 `useMoveParticlePartition = useParticlePartition && !Pstream::parRun()`

**文件**：`src/lagrangian/basic/Cloud/Cloud.C:836-837`

### 3.2 修复 `-parallel` 模式 PstreamBuffers crash

**问题**：`dlbOffloadPlanner` 的 inter-rank collision offload 与 replicated mesh 的 raw MPI 通信冲突。

**修复**：`dlbBaseActive` 条件加 `&& !cloud_.replicatedMeshActive()`

**文件**：`src/lagrangian/dsmc/collisionPartnerSelection/derived/noTimeCounter/noTimeCounter.C:227-233`

### 3.3 修复初始粒子重复（`-parallel` + masterUncollated 模式）

**问题**：两个 rank 都读取全部粒子，migration 后产生重复。

**修复**：添加 `distributeInitialParticles()` 方法，直接删除非 owned 粒子（无 MPI）。无 `-parallel` 模式仍用 `migrateParticlesByCellOwner()` 分发。

**文件**：`src/lagrangian/dsmc/replicatedMesh/dsmcReplicatedMesh.C`

### 3.4 消除 Cloud::move() 的 MPI transfer 开销（-14s）

**问题**：replicated mesh 全网格无 processor patches，但 Cloud::move() 仍执行 MPI 同步。

**修复**：`if (!Pstream::parRun() || neighbourProcs.empty()) break;`

**文件**：`src/lagrangian/basic/Cloud/Cloud.C:1756`

### 3.5 融合 migration + rebuildMoveOrderedParcels（-24s）

**问题**：migration 后遍历链表构建 flat array（16s），migration 本身也遍历链表（60ms/step 中大部分是指针追踪）。

**修复**：
- migration Phase 1 用 `moveOrderedParcels_` flat array 替代链表遍历
- Phase 1 同时收集 `kept` 列表
- Phase 4 反序列化时直接 append 到 `kept`
- 最后 `setMoveOrderedParcels(kept)` 一步到位

**文件**：`src/lagrangian/dsmc/replicatedMesh/dsmcReplicatedMesh.C:1039-1240`

### 3.6 启用 move extract ordered-reuse（-21s）

**问题**：Cloud::move() 的 `canReuseMoveOrdered` 要求 `accumulateMixedMoveOrdered`，但无 `-parallel` 模式下该条件为 false。

**修复**：放宽 `canReuseMoveOrdered` 条件为 `moveLoopPasses == 1`（不再要求 `accumulateMixedMoveOrdered`）。

**文件**：`src/lagrangian/basic/Cloud/Cloud.C:1058-1094`

### 3.7 MPI 通信优化（Sendrecv 替代 Allreduce + header）

**修复**：2-rank 用 `MPI_Sendrecv` 直接交换 size + data；多 rank 用 `MPI_Alltoall` + `Isend/Irecv`。消除 `MPI_Allreduce` 全局同步和 header packing。

**文件**：`src/lagrangian/dsmc/replicatedMesh/dsmcReplicatedMesh.C` Phase 3

### 3.8 修复 MPI8 无 OMP 模式 crash

**问题**：`useOpenMP=false` 时 Cloud::move() 不更新 `moveOrderedParcels_`，move 后粒子被删除导致 dangling pointers。

**修复**：move 后如果 replicated mesh 激活且 OMP 未启用，invalidate `moveOrderedParcelsValid_`。move 前如果 replicated mesh 激活且 cellOccupancy 无效，也 invalidate。

**文件**：`src/lagrangian/dsmc/clouds/dsmcCloud.C:3891-3901`

### 3.9 DLB 权重实时调整

**修复**：ParMETIS 双约束权重中 collision 约束的 scale 从固定值改为实时计算：`collScale = (collTime/totalCandidates) / (moveTime/totalParticles)`。`vsize` 按 ParDSMC3D 方案设为粒子数。

**文件**：`src/lagrangian/dsmc/replicatedMesh/dsmcReplicatedMesh.C` reassignByParMetisAdaptiveRepart()

## 4. 性能瓶颈分析（MPI2×OMP4, 108s）

| 阶段 | 时间 | 占比 | 可优化性 |
|------|------|------|----------|
| Move kernel (4线程) | 45s | 42% | 计算本身，需更多线程 |
| Collision (4线程) | 25s | 23% | 计算本身 |
| Migration MPI 同步等待 | 15s | 14% | rank 间负载不均衡的固有开销 |
| Post fields/output | 10.5s | 10% | 场量计算 |
| BuildCellOccupancy | 5s | 5% | 已优化到位 |
| Move extract | 0.1s | <1% | 已优化到位 |

## 5. 已验证无效/中性的方向

| 方向 | 结果 | 原因 |
|------|------|------|
| 单约束 DLB (move+coll 混合权重) | 更差 (110s) | move 和 collision 空间反相关，单约束无法同时均衡 |
| ubvec 差异化 (move 严格/coll 宽松) | 中性 (108s) | 同上，两约束互相矛盾 |
| 无 `-parallel` 模式 vs `-parallel` | 持平 (~108s) | partition 路径无法在 migration 后使用 |
| `hasParticlePartition` 检查 `isCellOccupancyValid` | 导致纯 OMP 退化 | 第一步 false 导致级联失效 |

## 6. 架构设计

### 两种运行模式

| | 无 `-parallel` (真 replicated mesh) | 有 `-parallel` (兼容模式) |
|---|---|---|
| `Pstream::parRun()` | false | true |
| MPI 通信 | raw MPI (dsmcReplicatedMesh) | raw MPI + Pstream |
| processor patches | 无 | 无 (symlink 全网格) |
| 初始粒子分发 | rank 0 发送给其他 rank | 每 rank 删除非 owned |
| Cloud::move() transfer | 直接 break | `neighbourProcs.empty()` break |

### evolve 流程（replicated mesh 模式）

```
controlBeforeMove (inflow 创建粒子)
  → [replicated mesh: invalidate moveOrdered if cellOccupancy invalid]
  → Cloud::move() (OMP 并行 move kernel)
  → [OMP: endMoveAppendCapture / 非OMP: invalidate moveOrdered]
  → migrateParticlesByCellOwner() (融合: flat array 遍历 + MPI + setMoveOrdered)
  → buildCellOccupancy (用 moveOrdered reuse)
  → collision (OMP 并行)
  → autoRebalance (ParMETIS, 实时权重)
```

## 7. 关键文件清单

- `src/lagrangian/basic/Cloud/Cloud.C` — OMP move 内核、extract、commit、transfer
- `src/lagrangian/dsmc/clouds/dsmcCloud.C` — evolve 流程、buildCellOccupancy
- `src/lagrangian/dsmc/clouds/dsmcCloud.H` — 成员声明
- `src/lagrangian/dsmc/replicatedMesh/dsmcReplicatedMesh.C` — migration、DLB、ParMETIS
- `src/lagrangian/dsmc/collisionPartnerSelection/derived/noTimeCounter/noTimeCounter.C` — dlbOffload
- `src/lagrangian/dsmc/parcels/dsmcParcel.C` — writeBinaryFast、move

## 8. 结论

Replicated mesh DLB 方案在 8 核规模下实现了：
- 比旧 decompose MPI8 快 **54%**（149s vs 326s）
- 比旧 decompose MPI2×OMP4 快 **8.5%**（108s vs 118s）
- 与纯 OMP8 差距 **27%**（108s vs 85s），差距来自 migration 固有开销

该方案的真正价值在于：
1. **DLB 能力**：可动态调整负载，适应非稳态流场
2. **可扩展性**：更多 rank 时 DLB 收益更大（当前 2 rank 均匀流场收益有限）
3. **无 processor boundary**：消除了 OpenFOAM 标准并行的 transfer 瓶颈

## 9. DLB 权重 K 值调优（2026-05-18）

ParMETIS AdaptiveRepart 使用单约束权重 `weight = N + N*(N-1)/K`，K 控制 move 和 collision 的折中。

### 4-rank K 值扫描

| K | Wall time | Move max | Coll max | 特征 |
|---|-----------|----------|----------|------|
| 8 | 124s | 48.2s | 40.0s | collision 过度均衡 |
| 16 | 118s | 45.4s | 35.7s | |
| 32 | 112s | 50.9s | 30.4s | |
| **64** | **106s** | 48.5s | 29.6s | **最佳折中** |
| 128 | 109s | 50.4s | 33.3s | move 过度均衡 |

### 8-rank K 值扫描

| K | Wall time |
|---|-----------|
| 8 | 183s |
| 16 | 160s |
| 32 | 179s |
| 64 | 152s |
| **128** | **134s** |

### Auto-K 计算

```
K = 2 * (moveTime/totalParticles) / (collTime/totalCandidates)
clamp to [32, 128], default 64
```

早期数据不稳定时 K 偏低（被 clamp 到下限）。对于 4-rank 最优 K=64，8-rank 最优 K=128。

## 10. Inter-Rank Collision Offload 实验（2026-05-17~18）

### 实验结论

| 方案 | 结果 | 原因 |
|------|------|------|
| 2-rank offload (raw MPI) | 无改善 (113s vs 112s) | remote execution 成本 ≈ local，无净收益 |
| 4-rank offload (raw MPI) | 更差 (131s vs 125s) | MPI_Alltoall 同步 + 通信开销 > 均衡收益 |
| DLB + offload 组合 | crash | cell ownership 迁移后数据冲突 |

### 根因分析

Collision offload 的数据流：
```
donor parcels → serialize → MPI send → helper deserialize → collide → serialize → MPI send → donor apply
```

每步开销：~200KB 数据传输 + 反序列化 + 碰撞执行 + 序列化 + 回传。
remote execution 成本 ≈ local execution 成本（同一算法，同一 N² 复杂度）。
offload 只是把工作从 donor 移到 helper，不减少总工作量。

### 有效条件

Offload 仅在以下条件下有效：
1. DLB 先粗粒度均衡 cell ownership（减少 donor 的 heavy cells）
2. Offload 再做细粒度调整（处理 DLB 后的残余不均衡）
3. Helper 有足够 spare capacity 吸收 remote work

当前 case（均匀圆柱流，2-rank）中条件 1 和 3 不满足，offload 无收益。

## 11. 当前推荐配置

| 场景 | 推荐配置 | 预期性能 |
|------|----------|----------|
| 单节点 8 核 | 纯 OMP8 + collisionFastRng | 79s |
| 单节点 8 核（需 DLB） | MPI4×OMP2 + DLB(K=64) + FastRng | 99s |
| 跨节点 2×4 核 | MPI2×OMP4 + FastRng | 97s |
| 跨节点 8×1 核 | MPI8 + DLB(K=128) + FastRng | 134s |

## 12. React Case 性能对比（2026-05-18）

测试 case：`hyStrath_xcx/case/cylinder_react/mixparallel/allmesh_react/`
配置：300步，8核，collisionFastRng 开启，cylinder N2 react

| 配置 | Wall time | vs 纯OMP8 |
|------|-----------|-----------|
| **纯 OMP8** | **85.2s** | baseline |
| **MPI2×OMP4 replicated mesh** | **100.8s** | +18% |
| **MPI4×OMP2 + DLB(K=64)** | **106.5s** | +25% |
| **MPI8 + DLB** | **142.1s** | +67% |

### 与 noreact case 对比

| 配置 | noreact | react | react 额外开销 |
|------|---------|-------|---------------|
| 纯 OMP8 | 79.2s | 85.2s | +8% |
| MPI2×OMP4 | 96.7s | 100.8s | +4% |
| MPI4×OMP2 + DLB | 99.2s | 106.5s | +7% |
| MPI8 + DLB | 133.7s | 142.1s | +6% |

React 的额外开销（化学反应计算）对各配置影响均匀（+4~8%），不改变配置间的相对排序。

### 测试 case 目录

- `allmesh_react/omp8/` — 纯 OMP8
- `allmesh_react/omp4_mpi2_replicatedmesh/` — MPI2×OMP4，无 DLB
- `allmesh_react/omp2_mpi4_replicatedmesh/` — MPI4×OMP2，DLB 开启
- `allmesh_react/mpi8_replicatedmesh/` — MPI8，DLB 开启

## 13. 延迟接收 Migration + interval 调优（2026-05-18）

### 延迟接收（migrateBegin/Finish 拆分）

将同步 `migrateParticlesByCellOwner()` 拆分为：
- `migrateBegin()`: 分区 + delete + Isend + Irecv（非阻塞）
- `migrateFinish()`: Waitall + 反序列化 + appendBatch（下一步开头执行）

MPI2×OMP4 结果：101.5s → **98.7s**（-2.8s）

### Migration interval=2

通过 controlDict 参数 `replicatedMeshMigrateInterval 2` 控制。

| 配置 | interval=1 | interval=2 | 节省 | 精度偏差 |
|------|-----------|-----------|------|----------|
| MPI2×OMP4 | 101.5s | 98.8s | -2.7s | <0.3% |
| MPI8 (K=128) | 126.9s | **113.0s** | **-13.9s** | <0.5% |

8-rank 节省更大（8× migration 通信量，减半频率效果显著）。

### 8-rank DLB K 值 + interval 组合

| 配置 | Wall time |
|------|-----------|
| K=128 | 126.9s |
| K=256 | 117.9s |
| K=∞ (纯 N) | 135.2s |
| **K=128 + interval=2** | **113.0s** |
| K=256 + interval=2 | 116.6s |

### 自适应 K 默认值

```cpp
label collDivisor = (nProcs_ >= 8) ? 128 : 64;
```

### 最终推荐配置（更新）

| 场景 | 配置 | 性能 |
|------|------|------|
| 单节点 8 核 | 纯 OMP8 + FastRng | 79-85s |
| 需 DLB, 2 rank | MPI2×OMP4 + 延迟接收 | 95-99s |
| 需 DLB, 4 rank | MPI4×OMP2 + DLB(K=64) | 99-107s |
| 需 DLB, 8 rank | MPI8 + DLB(K=128) + interval=2 | **113s** |

## 14. Post-Collision Migration 实验（2026-05-19）

### 架构设计

将 `migrateBegin()` 从 move 之后移到 collision 之后，使 DLB 可以独立平衡 collision 成本，move 不平衡由 migration 通信吸收。

核心改动：
- `dsmcCloud.C`：controlDict 开关 `replicatedMeshPostCollisionMigration true/false`（默认 false）
- `noTimeCounter.C`：collision 循环加 `isMyCell` 检查（跳过 non-owned cells）
- `dsmcReplicatedMesh.C`：双约束 DLB（计算成本 + cell 数量）、adaptive K 反馈、PartKway 备选

### Post-Collision Migration 8-rank 测试

| 配置 | Wall time | vs baseline |
|------|-----------|-------------|
| Baseline (delayed-receive + K=128) | 119.1s | — |
| Post-collision + collision-only weight | 174.9s | +47% ❌ |
| Post-collision + K=16 | 153.95s | +29% ❌ |
| Post-collision + K=128 + adaptive K | 118.5s | -0.5% |
| Post-collision + 双约束 (ubvec 1.05/1.03) | **116.0s** | **-2.6%** |
| Post-collision + PartKway | 129.5s | +9% ❌ |

### DLB 触发策略对比

| 策略 | DLB 次数 | Wall time |
|------|---------|-----------|
| **sar > 0 趋势触发** | **3-4** | **116.0s** |
| 绝对阈值 1.15 | 29 | 122.4s |
| 绝对阈值 1.25 + min 30步 | 8 | 119.0s |
| 固定 50 步 | 6 | 122.4s |

### DLB 权重方案对比

| 方案 | Wall time | 说明 |
|------|-----------|------|
| N + N*(N-1)/K (K=128) | 118.5s | 原始单约束 |
| N*(N-1)/2 (纯 collision) | 174.9s | 权重比太极端 |
| **5 + N + N*(N-1)/K + cell约束** | **116.0s** | 双约束最优 |
| TACF 实测权重 | 140.4s | 累积平均含瞬态，不准 |

### 关键发现

1. **Post-collision migration 对 2-rank 无效**：delayed-receive 的 component cancellation 天然平衡 wall time（imbalance 1.001），post-collision 破坏此平衡。

2. **Post-collision migration 对 8-rank 有限收益**（116s vs 119s，-2.6%）：DLB 已经通过 Alltoall 同步吸收了工作不平衡，post-collision 的额外收益来自双约束分区质量改善。

3. **"buildCellOccupancy" 计时包含 migrateBegin**：profiling 中 "buildCellOccupancy [s]" 实际 = migrateBegin（含 Alltoall 等待）+ TACF + Phase B + 实际 buildOcc。实际 buildOcc ≈ 5s。

4. **MPI8 实际工作 imbalance 1.76:1**（rank 7 = 69.7s vs rank 0 = 39.5s），被 Alltoall 同步等待吸收为 wall time imbalance 1.011。

### 已验证无效方向

| 方向 | 结果 | 原因 |
|------|------|------|
| Collision-only 权重 (K→0) | 174.9s | 权重比 4950:1，ParMETIS 分区质量崩溃 |
| PartKway 替代 AdaptiveRepart | 129.5s | 空间局部性差，migration 量增大 |
| TACF 实测权重 | 140.4s | 累积平均含初始瞬态 |
| 绝对 imbalance 阈值触发 | 119-122s | DLB 过于频繁 |
| 固定间隔 DLB | 122.4s | 首次 DLB 太晚 + 不必要触发 |
| cellBaseCost=20 | 123.2s | 权重比过小 |

## 15. 最终性能确认（2026-05-19，React + FastRng）

### 测试环境

- 8 核单节点，Intel oneAPI 2025.2 + OF-2506
- 300 步，cylinder N2 react，collisionFastRng true
- 配置：delayed-receive + DLB(auto-K, clamp [32,128]) + migrateInterval 按需

### 结果

| 配置 | Wall time | vs Worklog | Imbalance |
|------|-----------|-----------|-----------|
| **纯 OMP8** | **81.3s** | -5% | — |
| **MPI2×OMP4** (delayed-receive, no DLB) | **94.8s** | -6% | 1.001 |
| **MPI4×OMP2** (delayed-receive + DLB K=64) | **100.2s** | -6% | 1.014 |
| **MPI8** (delayed-receive + DLB K=128 + interval=2) | **131.3s** | -8% | 1.011 |

### MPI8 per-rank 负载分析

| Rank | Move | Collision | Migration | Particles | 实际工作 |
|------|------|-----------|-----------|-----------|---------|
| 0 | 25.5s | 14.0s | 79.5s | 63K | 39.5s |
| 1 | 26.4s | 22.6s | 70.2s | 142K | 49.0s |
| 2 | 28.9s | 23.5s | 65.4s | 208K | 52.4s |
| 3 | 42.1s | 12.6s | 57.0s | 223K | 54.7s |
| 4 | 40.6s | 32.3s | 41.7s | 366K | 72.9s |
| 5 | 47.7s | 16.2s | 46.6s | 279K | 63.9s |
| 6 | 48.0s | 23.1s | 33.5s | 310K | 71.1s |
| 7 | 52.3s | 17.4s | 37.9s | 368K | **69.7s** |

实际工作 imbalance = 69.7/39.5 = 1.76:1，被 MPI_Alltoall 同步等待吸收。

### 最终推荐配置（更新）

| 场景 | 配置 | 性能 |
|------|------|------|
| 单节点 8 核 | 纯 OMP8 + FastRng | **81s** |
| 需 DLB, 2 rank | MPI2×OMP4 + delayed-receive | **95s** |
| 需 DLB, 4 rank | MPI4×OMP2 + DLB(K=64) | **100s** |
| 需 DLB, 8 rank | MPI8 + DLB(K=128) + interval=2 | **131s** |
| 需 DLB, 8 rank (post-coll) | MPI8 + post-collision + 双约束 DLB | **116s** |
