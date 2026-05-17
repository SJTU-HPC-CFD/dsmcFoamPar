# dsmcFoam+ Replicated Mesh DLB 优化工作总结（2026-05-16）

## 1. 工作目标

在 dsmcFoam+ 中实现 Replicated Mesh + Dynamic Load Balancing 方案，使 MPI+OpenMP 混合并行在保持 DLB 能力的同时，性能接近纯 OpenMP。

核心思路：每个 MPI rank 持有完整网格，通过 `cellOwner_[]` 划分计算责任，ParMETIS AdaptiveRepart 动态调整 cell 归属。

## 2. 最终性能对比（300步，8核，cylinder N2 noreact）

| 配置 | Wall time | vs 纯OMP8 | 备注 |
|------|-----------|-----------|------|
| 旧 decompose MPI8 | 326s | +284% | 标准 OpenFOAM 分区 |
| 旧 decompose MPI2×OMP4 | 118s | +39% | stagebuf2 优化后 |
| **纯 OMP8** | **85s** | baseline | dlbOffload 开启 |
| **MPI2×OMP4 replicated mesh** | **108s** | +27% | 新方案，有 DLB |
| **MPI4×OMP2 replicated mesh** | **118s** | +39% | rank 越多 migration 越大 |
| **MPI8 replicated mesh (无OMP)** | **149s** | +75% | migration 开销主导 |

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
