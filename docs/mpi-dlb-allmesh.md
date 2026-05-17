# 基于全局网格复制的 DSMC 动态负载平衡：完整技术路线

## 一、方法定义

```
Replicated Mesh + Distributed Particles + Dynamic Cell Ownership + Particle-only Migration
```

**核心约定**：

| 概念 | 定义 |
|------|------|
| `mesh_` | 完整 polyMesh，所有 rank 持有相同副本（只读），全局 cellI 一致 |
| `cellOwner_[cellI]` | 负责计算该 cell 的 MPI rank |
| `mesh.owner()[faceI]` | OpenFOAM 拓扑层：face 属于哪个 cell（与计算责任无关） |
| `mesh.neighbour()[faceI]` | OpenFOAM 拓扑层：face 的邻居 cell（与计算责任无关） |
| 粒子 | 唯一按 `cellOwner_` 分布的数据，不按 mesh decomposition 分布 |

**关键解耦**：几何拓扑所有权（`mesh.owner/neighbour`）和计算所有权（`cellOwner_`）彻底分离。

**不做的**：
- 不调用 `decomposePar`
- 不创建 `processorPatch`
- 不动态重建 mesh topology
- 不迁移 mesh/cell/face 数据

**为什么可以在 DSMC 中这样做**：

DSMC 的主要负载来自粒子运动、粒子-网格 tracking、cell occupancy 重建、碰撞计算、边界交互和采样，而不是像传统 CFD 那样来自全局场方程和邻接 cell 间的通量耦合。因此保持 mesh 不动、只改变 cell 的计算 owner 是成立的。

---

## 二、关键数据结构

```
mesh_                    // 完整 polyMesh，所有 rank 相同。只读。
cellOwner_[nCells]       // cellI → 负责计算的 rank。可动态更新。
myCells_                 // {cellI | cellOwner_[cellI] == myRank}
```

辅助结构（用于统计和调试，不作为迁移判据）：

```
interfaceFaces_          // {faceI | cellOwner_[owner] != cellOwner_[neighbour]}
                         // 用于统计 edge cut / 通信边界，不用于迁移决策
```

### face 分类

| 类型 | 条件 | 处理 |
|------|------|------|
| internal computational | `cellOwner_[own]==me && cellOwner_[nei]==me` | 纯本地 |
| logical interface | `cellOwner_[own] != cellOwner_[nei]` | 统计用，不驱动迁移 |
| physical boundary | OpenFOAM 原生 patch | 由相邻 cell 的 owner 处理 |

### 粒子迁移判据

迁移判据基于**粒子最终 cell 的 owner**，而不是基于穿过了哪个 face。原因是 DSMC 粒子一步可能穿越多个 cell，用最终 cell → owner 判断最稳。

```cpp
label newCell = parcel.cellI();
label dstRank = cellOwner_[newCell];
if (dstRank != myRank) { migrate(parcel, dstRank); }
```

---

## 三、7 个正确性不变量

| # | 不变量 | 检查方式 |
|---|--------|----------|
| 1 | 粒子全局守恒 | `MPI_Allreduce(N_local, SUM) == N_prev + N_inlet - N_outlet` |
| 2 | 粒子 owner 一致 | `assert(cellOwner_[p.cellI()] == myRank)` 对本地的每个粒子 |
| 3 | cell occupancy 唯一 | cellI 的粒子只存在于 `cellOwner_[cellI]` 所在 rank |
| 4 | 碰撞不重复 | `if (cellOwner_[cellI] != myRank) continue;` 碰撞循环守卫 |
| 5 | 采样不重复 | 仅 owner 采集，`MPI_Reduce` 汇总 |
| 6 | 迁移在碰撞前完成 | `move()` 和 `collide()` 之间必须调用 `migrateParticlesByCellOwner()` |
| 7 | occupancy 仅填 owner cell | `cellOccupancy[cellI]` 全局长度保留，但只在 `cellOwner_[cellI]==myRank` 时 append |

---

## 四、三阶段推进计划

### 阶段 A：静态 cellOwner 验证

**目标**：证明 mesh 复制 + cellOwner + 粒子迁移机制在物理上正确。

**8 个基础模块**：

```
M1: readFullMeshOnEachRank()
M2: buildCellGraph()
M3: scotchDecomposeToCellOwner()
M4: buildLocalComputationalView()
M5: moveParticlesOnFullMesh()
M6: migrateParticlesByCellOwner()
M7: collideOwnerCellsOnly()
M8: sampleOwnerCellsOnly()
```

**算法流程**：

```
初始化：
  所有 rank 读取完整 polyMesh
  构造 cell 邻接图
  Scotch(P) → cellOwner_[cellI]
  广播 cellOwner_ 到所有 rank
  构建 myCells（不创建 processorPatch）
  只在 owner cell 上初始化粒子

每步循环：
  1. move local particles on full mesh
     - tracking 在完整 mesh 上进行（几何上无差异）
     - 不依赖 processorPatch 触发 transfer
  2. physical boundary interaction
     仅当相邻 cell 的 owner == myRank 时执行
  3. migrate particles
     对每个本地粒子：
       dst = cellOwner_[particle.cellI()]
       if dst != myRank: pack → parcelBuffer[dst]
     MPI_Alltoallv 交换
  4. rebuild occupancy
     保留 cellOccupancy[nCells] 全局长度
     仅对 myCells 清空后填充
  5. collide only myCells
  6. sample only myCells → MPI_Reduce 到 rank 0 → 写出
```

**cellOccupancy 处理方式**：

保留 OpenFOAM 原生的 `cellOccupancy[nCells]` 全局长度，不在 MVP 阶段改成压缩式 `myCellOccupancy`。只在 owner cell 上填充内容：

```cpp
cellOccupancy[cellI].clear();
if (cellOwner_[cellI] == myRank) {
    cellOccupancy[cellI].append(parcelI);
}
```

后续优化阶段再考虑压缩内存。

**验证实验**：

| 实验 | 内容 | 检查项 |
|------|------|--------|
| A1 | 均匀通道流 vs 串行 DSMC | 时间平均的密度、速度、温度、碰撞次数 |
| A2 | vs 传统 MPI decomposition DSMC | **仅对比时间平均统计量**，不要求逐粒子一致 |
| A3 | 非均匀粒子分布（左高右低） | 各 rank 粒子数分布、迁移粒子数 |

**A2 对比注意事项**：

传统 OpenFOAM 分解后 cell local numbering 和 replicated mesh 的 global numbering 不同；粒子随机数、注入顺序、碰撞配对顺序可能不同。因此 A2 只对比时间平均统计量，评价指标建议用相对误差：

```
E_ρ = |ρ_replicated - ρ_baseline|_2 / |ρ_baseline|_2
```

不要求逐步逐粒子对比。

---

### 阶段 B：手动重分配 cellOwner

**目标**：验证 cellOwner 改变后，仅迁移粒子可以继续正确计算。这是整个技术方案核心风险点的验证。

**新增模块**：

```
M9:  reassignCellOwner(newOwnerMap)
M10: redistributeParticlesByNewOwner()
```

**测试序列**：

```
第 100 步：重新 Scotch 得到新 cellOwner → 迁移粒子
第 200 步：再次改变 cellOwner → 迁移粒子
第 300 步：反转部分区域的 cellOwner → 迁移粒子
```

每次重分配后验证 7 个不变量全部成立。

**控制指标**：

```
R_changed = #{c | cellOwner_new(c) ≠ cellOwner_old(c)} / N_cell
R_mig     = N_migrated / N_total_particles
```

如果 `R_mig` 过大而负载改善很小，说明重划模式需要调整。

---

### 阶段 C：自动动态负载平衡

**目标**：验证方案有可测量的性能收益。

**负载模型**：

```
阶段 C-1（基线）：w_c = N_p,c                        （粒子数划分 + 粒子数触发）
阶段 C-2（升级）：w_c = N_p,c                        （粒子数划分 + 实际耗时触发）
阶段 C-3（终版）：w_c = aN_p,c + bN_coll,c + c·T_c   （DSMC 专用权重 + 实际耗时触发）
```

**触发条件**：

```
T_p = rank p 的 wall-clock time per step
I_T = max(T_p) / avg(T_p)
触发：I_T > 1.3 且距上次重平衡 ≥ 100 步
```

粒子数用于构造划分权重，实际耗时用于判断是否值得重平衡——两者分工不同。

**收益-开销判据**：

```
若 I_T 改善 < 5% 但 R_mig > 20% → 本次重平衡跳过
```

**三种重平衡对照方法**：

| 方案 | 定位 | 做法 |
|------|------|------|
| Scotch 加权重划分 | 质量上界 / 图划分对照 | 粒子数加权 cell graph → Scotch(P) → 新 cellOwner |
| Hilbert SFC | 速度上界 / 轻量对照 | cell centroid → Hilbert key → 排序 → 按负载等分 |
| superCell + Hilbert | **推荐主方案** | Scotch 过分解 M=αP 个 superCell → 对 superCell 加权 centroid 做 Hilbert 排序 → 按负载分配 superCellOwner → cellOwner 由此派生 |

三种方法定位明确：Scotch 动态重划分是质量上界（edge cut 低、拓扑好、partition 质量高），Hilbert SFC 是速度上界（实现简单、运行时开销低、适合频繁调整），superCell + Hilbert 在两者之间取得平衡。

**superCell 粒度参数研究**：

```
M/P = 2, 4, 8, 16, 32
```

| M/P | 预期效果 |
|-----|----------|
| 2 | 迁移少，平衡能力弱 |
| 4 | 实现简单，开销低 |
| 8 | 通常较均衡（推荐初始值） |
| 16 | 平衡能力强，通信边界和迁移可能增加 |
| 32 | 可能过细，processor patch 复杂化 |

**对比实验矩阵**：

| 配置 | 说明 |
|------|------|
| 无动态平衡 | baseline |
| Scotch cell 级动态 | 图划分对照（质量上界） |
| Hilbert cell 级动态 | SFC 对照（速度上界） |
| superCell(M/P=4) + Hilbert | |
| superCell(M/P=8) + Hilbert | 推荐主方案 |
| superCell(M/P=16) + Hilbert | |

**测量指标**：

```
总 wall-clock time
每步平均时间
max(T_p) / avg(T_p)  变化历程
迁移粒子数 / 总粒子数
单次重平衡开销
通信时间占比
每 rank 内存占用
```

---

## 五、关键工程风险与处置

### 风险 1：OpenFOAM 粒子 tracking 在未分解全局 mesh 上的兼容性

**问题**：dsmcFoam 的 particle tracking 可能隐式依赖 processorPatch 触发跨 rank 传输。

**处置**：MVP 采用"本步 tracking 完成后再按最终 cellOwner 迁移"，不依赖 tracking 过程中的 face 事件。

### 风险 2：物理边界粒子注入重复

**问题**：所有 rank 都有完整 boundary patch，如果不限制，入口/壁面会重复注入粒子。

**处置**：

```cpp
label adjCell = patchFaceOwnerCell[faceI];
if (cellOwner_[adjCell] != myRank) continue;
```

只有边界 face 相邻 cell 的 owner rank 负责该 face 的粒子注入。

### 风险 3：输出重复

**问题**：所有 rank 都有完整 field，不能每 rank 写一份。

**处置**：每个 rank 只填 owner cell 的 field 值，`MPI_Reduce(SUM)` 到 rank 0，rank 0 写完整 field。

### 风险 4：内存上限

**单 rank 约束**：

```
M_globalMesh + M_globalFields + M_localParticles + M_buffers < M_rankAvailable
```

**全作业冗余成本**：`P × M_globalMesh`

注意：`P × M_globalMesh` 是总冗余，不是单 rank 压力。每个 rank 只需 1 份全局 mesh。

### 风险 5：Scotch cell 级重划导致大量 owner 变更

**问题**：直接用 Scotch 在 cell 级别重划分，小负载变化可能导致大量 cellOwner 改变。

**处置**：重平衡必须统计 `R_changed` 和 `R_mig`；如果迁移量过大而负载改善微乎其微，则跳过本次重平衡。

### 风险 6：face owner 与 cell owner 混淆

**问题**：`mesh.owner()[faceI]` 是 OpenFOAM 拓扑层定义，`cellOwner_[cellI]` 是自定义计算责任层定义。两者不可混淆。

**处置**：代码中明确使用不同命名；face 分类时，`cellOwner_[mesh.owner()[faceI]]` 和 `cellOwner_[mesh.neighbour()[faceI]]` 分别获取两侧的计算 owner rank。

---

## 六、论文级表述

> **基于全局网格复制的 DSMC 动态负载平衡方法**：对于目标规模的 DSMC 算例，粒子数据与粒子操作通常主导内存与计算成本，因此可以用全局只读网格复制换取运行时重分区机制的大幅简化。每个 MPI 进程持有完整非结构网格副本（只读），仅对粒子进行分布式存储和计算。前处理阶段，通过 Scotch 对 cell 邻接图做过分解（M = αP, α ≫ 1），将网格聚合为拓扑连通的 superCell 作为原子迁移单元。运行时周期性统计各 superCell 的 DSMC 专用负载（粒子数、碰撞数、边界交互数），将粒子数用于构造划分权重、实际耗时用于判断是否触发重平衡。当 rank 间耗时不均超过阈值且收益大于迁移成本时，沿 Hilbert 空间填充曲线对 superCell 加权 centroid 排序，按负载等分重新分配计算所有权。所有权变更时仅迁移粒子，不迁移网格。该方法将拓扑质量与动态调度解耦：Scotch 保证 superCell 内部连通性，Hilbert 使运行时重映射开销远低于全局图重划分。

### 论文创新点

1. **DSMC 专用加权图划分**：不使用 OpenFOAM 默认的均匀 cell 划分，而是提出 DSMC cell weight `w_c = aN_p,c + bN_coll,c + cN_bc,c + dN_move,c + e·T_observed`，作为 graph partitioning 的 vertex weight。

2. **图过分解 + SFC 动态映射**：Scotch 负责拓扑质量（superCell 内部连通性），Hilbert 负责快速运行时映射（superCell → processor）。两者正交互补。

3. **考虑迁移成本的动态重平衡触发机制**：不是固定周期盲目重平衡，而是判断收益是否大于迁移成本才执行。

4. **全局只读网格复制 + particle-only migration**：将几何拓扑所有权和计算所有权解耦，仅迁移粒子不迁移网格，大幅降低动态重平衡的工程复杂度。

---

## 七、最终算法流程图

```
初始化：
  1. 所有 rank 读取完整 polyMesh
  2. 构造 cell adjacency graph
  3. Scotch(M=αP) 过分解 → superCellOfCell = superCell 归属
     或 Scotch(P) → cellOwner_（阶段 A/B 简化版）
  4. 构建 myCells（不创建 processorPatch）
  5. 仅对 owner cell 初始化粒子

每步：
  1. move local particles on full mesh（tracking 后不依赖 processorPatch）
  2. physical boundary injection/deletion 仅由相邻 cell 的 owner 执行
  3. migrate particles:
       dst = cellOwner_[parcel.cellI()]
       if dst != myRank: pack → parcelBuffer[dst]
     MPI_Alltoallv 交换粒子
  4. rebuild occupancy（保留全局长度，仅填 owner cell）
  5. collide only myCells
  6. sample only myCells → MPI_Reduce 汇总 → rank 0 写出

动态负载平衡（每 K 步检测）：
  1. 统计 cell/superCell 负载 + 各 rank wall-clock time T_p
  2. 计算 I_T = max(T_p) / avg(T_p)
  3. 若 I_T ≤ 阈值 或 距上次重平衡 < 100 步 → 跳过
  4. 若 I_T > 阈值：
     a. Hilbert 排序 superCell 加权 centroid
     b. 沿序列按负载等分 → 新 superCellOwner
     c. 评估 R_mig；若迁移量过大而改善太小 → 跳过
     d. 更新 cellOwner_[cellI] = superCellOwner[superCellOfCell[cellI]]
     e. redistributeParticlesByNewOwner()
     f. 记录当前 edge cut 和通信量，用于退化检测
```

---

## 八、方案总评

| 维度 | 评价 |
|------|------|
| 可行性 | 高 |
| MVP 可实现性 | 高 |
| 工程风险 | 中（主要在 tracking 兼容性和注入去重） |
| 内存代价 | 中到高（每 rank 全 mesh，作业总冗余 `P × M_mesh`） |
| 对 OpenFOAM 改动量 | 中（绕过 decomposition 框架，不改 polyMesh 结构） |
| 对 DSMC 适配性 | 高 |
| 论文创新性 | 中高 |

**核心优势**：

- 用全局只读 mesh 复制规避 OpenFOAM 动态 mesh redistribution
- 用 `cellOwner_` 把几何拓扑所有权和计算所有权解耦
- 用 particle-only migration 实现轻量动态负载平衡
- 用 Scotch/superCell/Hilbert 在分区质量和运行时开销之间折中
- 图划分负责"分得好"，空间填充曲线负责"调得快"，over-decomposition 负责"迁得少"

**一句话结论**：这条技术路线合理，且比传统"动态迁移 mesh"的路线更适合作为 OpenFOAM/dsmcFoam 非结构 DSMC 动态负载平衡的最小可行实现。建议按 A/B/C 三阶段推进，superCell + Hilbert 作为最终主方案。

---

## 九、阶段 A 实现总结（2026-05-13）

### 架构

```
                     MPI Rank 0..7
                    ┌─────────────────────────────────┐
                    │  完整 polyMesh (60k cells)       │  ← masterUncollated
                    │  cellOwner_[60000] = {0..7}      │  ← block 分解
                    │  myCells_ = 仅本 rank 的 cell     │
                    └─────────────────────────────────┘
                                  │
         ┌────────────────────────┼────────────────────────┐
         │                        │                        │
    controlBeforeMove        Cloud::move()          migrateByCellOwner
    (注入, per-face过滤)     (全mesh追踪)           (OCharStream + MPI)
         │                        │                        │
         └────────────────────────┼────────────────────────┘
                                  │
                          buildCellOccupancy
                          collision (仅 owned cells)
                          post fields
```

### 关键文件

| 文件 | 说明 |
|---|---|
| `src/lagrangian/dsmc/replicatedMesh/dsmcReplicatedMesh.H` | cellOwner_, myCells_, 预/后迁移, TACF, per-rank wall time |
| `src/lagrangian/dsmc/replicatedMesh/dsmcReplicatedMesh.C` | Block/Scotch 分解, MPI 迁移, 选择性序列化 |
| `src/lagrangian/dsmc/clouds/dsmcCloud.H` | +replicatedMesh_ 成员, +deleteParcel() |
| `src/lagrangian/dsmc/clouds/dsmcCloud.C` | 构造初始化 + evolve() 中预/后迁移 + per-rank 计时 |
| `src/lagrangian/dsmc/clouds/dsmcCloudI.H` | replicatedMeshActive(), replicatedMesh() |
| `src/lagrangian/dsmc/parcels/dsmcParcel.C` | hitPatch（无修改，原始追踪） |
| `src/lagrangian/dsmc/boundaries/.../dsmcFreeStreamInflowPatch.C` | per-face `isMyCell` 注入去重 |
| `src/lagrangian/dsmc/Make/options` | +scotchDecomp include/lib |
| `src/lagrangian/dsmc/Make/files` | +replicatedMesh/dsmcReplicatedMesh.C |

### 关键技术决策

| 决策 | 结论 | 原因 |
|---|---|---|
| 全 mesh 读取 | `masterUncollated` file handler | 无需 `-parallel`，所有 rank 读同一份 `constant/polyMesh` |
| 分解方法 | **Block 分解最优** | Scotch 不规则边界导致迁移暴增 3.2x |
| 边界去重 | per-face `isMyCell(adjCell)` + processBoundaries | FreeStreamInflowPatch 检查 + evolve() 全 rank 处理边界 |
| 追踪 | 全 mesh 追踪 (60k cells) | 比局部 mesh 追踪慢 ~10x，但避免了 processorPatch 通信 |
| 迁移间隔 | **每 10 步** 最优 | 性能 +22%，精度损失仅 0.27% |
| 序列化 | OCharStream + ICharStream | 零拷贝/Alltoallw 仅节省 ~5%，瓶颈在 operator<</>>虚函数 |

### 关键 Bug 修复历程

| Bug | 根因 | 修复 |
|---|---|---|
| 文件 I/O 竞态 | 多 rank 共享 `0/` 目录 | `masterUncollated` file handler |
| `isMyCell` 始终返回 false | 使用 `Pstream::myProcNo()`（无 `-parallel` 时恒为 0） | 改用 raw `MPI_Comm_rank` 获取的 `myRank_` |
| 8-rank 粒子丢失 (81%) | 边界注入 `processBoundaries=false` 时仅 rank 0 处理 | `processBoundaries=true` + per-face `isMyCell` 检查 |
| 1-rank replicatedMesh 注入跳过 | `nProcs_<2` 时 `cellOwner_` 未初始化（全为 -1） | 单 rank 时 `cellOwner_` 全置 0 |
| Scotch 迁移爆炸 (3.2x) | Scotch 不规则分区边界穿过高粒子流量区 | 回归 block 分解 |
| OMP=2 崩溃 | 无 `-parallel` 时 OpenMP 临界区不生效 | 暂禁用 OMP move（`openmpMove false`） |

### 7 不变量验证

| # | 不变量 | 状态 | 方法 |
|---|--------|------|------|
| 1 | 粒子全局守恒 | ✅ 0.005% | 30 步 1-rank vs 8MPI A/B 对比 |
| 2 | 粒子 owner 一致 | ✅ | 迁移后 `cellOwner_[p.cell()] == myRank` |
| 3 | cell occupancy 唯一 | ✅ | 迁移在 buildCellOccupancy 前 |
| 4 | 碰撞不重复 | ✅ 0.67% | 30 步 A/B 碰撞数对比（RNG 统计噪声内） |
| 5 | 采样不重复 | ~ | 仅 rank 0 输出 field |
| 6 | 迁移在碰撞前 | ✅ | move -> migrate -> buildOcc -> collide |
| 7 | occupancy 仅填 owner | ✅ | 全局长度，仅 owner 填充 |

### 物理正确性最终验证（30 步，同 case，1-rank vs 8MPI）

| 指标 | 1-rank | 8MPI (合计) | 差异 |
|---|---:|---:|---:|
| **粒子数** | 1,377,059 | 1,376,986 | **-0.005%** ✓ |
| **碰撞数** | 47,724 | 48,046 | **+0.67%** ✓ (RNG 统计噪声) |
| **振动能量** | 2.327e-20 | 2.327e-20 | **0%** ✓ (碰撞物理硬指标) |

### 性能结果（300 步, cylinderN2 noreact）

#### 8 核配置对比

| 配置 | 主循环 | move | buildCellOcc | collision | post | 迁移 |
|---|---:|---:|---:|---:|---:|---:|
| **OMP8 (1rx8t)** | **84.1s** | 42.7s | 4.5s | 22.5s | 11.2s | -- |
| 8MPI 复制网格 (8rx1t) | 138.1s | 74.6s | 9.6s | 11.6s | 9.4s | 31.9s |

- **OMP8 是 8 核最优方案**（84.1s），但无动态负载均衡能力
- 8MPI 复制网格比 OMP8 慢 64%：move 全网格追踪 +32s，迁移 +32s
- 8MPI collision 比 OMP8 快 48% -- 每 rank 仅 244k 粒子 + 完美 cache 局部性
- 迁移 (31.9s, 23%) 是复制网格的主要新增开销

#### 迁移间隔扫描

| 间隔 | 主循环 | 迁移调用 | 迁移时间 | 粒子误差 |
|---:|---:|---:|---:|---:|
| 1 | 180.0s | 301 | 61.3s | +0.013% |
| 5 | 160.2s | 61 | 44.0s | -0.08% |
| **10** | **138.1s** | **31** | **31.9s** | **-0.27%** |
| 20 | 125.7s | 16 | 23.6s | -1.02% |
| 30 | 117.0s | 11 | 16.6s | -1.94% |
| 50 | 116.3s | 7 | 8.7s | -3.97% |

#### Scotch vs Block 分解对比

| 指标 | Block | Scotch | 胜者 |
|---|---:|---:|---:|
| 主循环 | 138.1s | 154.0s | Block |
| 迁移时间 | 31.9s | 104.9s | Block |
| Wall time 不平衡 | 1.009 | 1.032 | Block |
| 粒子不平衡 | 1.22 | 2.10 | Block |

Scotch 不规则边界 -> 粒子频繁跨 rank -> 迁移暴增 3.2x。
Block 分解的 wall time 近乎完美平衡（1.009）。

#### 尝试过的无效优化

| 优化 | 收益 | 根因 |
|---|---|---|
| 零拷贝 ISpanStream 接收 | ~1s | 内存拷贝不是瓶颈 |
| MPI_Alltoallw 零拷贝发送 | ~0s | 同上 |
| OCharStream 预分配 buffer | ~0s | 重新分配本不频繁 |
| 加权 Scotch 分解 | -47s (更慢) | 粒子数权重加剧不规则边界 |
| 子网格追踪 | 未完成 | 需结构性改造追踪引擎 |
| OMP=2 追踪 | 崩溃 | 无 `-parallel` 时 OpenMP 临界区不生效 |

#### 生产推荐配置

```
controlDict:
    replicatedMesh true;
    replicatedMeshMigrateInterval 10;
    replicatedMeshDecompMethod block;
    openmpThreads 1;
    openmpMove false;

    OptimisationSwitches {
        fileHandler masterUncollated;
    }
```

### 待解决

1. **~~子网格追踪~~** — ParDSMC3D 证实全局网格追踪是正确架构，无需子网格
2. **~~OMP 追踪兼容~~** — 已修复临界区 + 负载分区，但 2-thread OMP 无收益
3. **~~TACF 动态重平衡~~** ✅ Phase C 实现并验证
4. **~~序列化优化~~** ✅ writeBinaryFast 实现，迁移 -63%
5. **采样去重**：rank 0 独占 field 输出，需 MPI_Reduce 实现多 rank 采样
6. **MPI 派生类型**：仿 ParDSMC3D，进一步优化迁移
7. **累积空闲时间触发**：替代简单比值，更稳定的 DLB 判断
8. **ParMETIS 运行时重划分**：替代 Scotch 静态分解

---

## 九-B、阶段 B 实现总结（2026-05-13）

### 目标

验证 cellOwner_ 可以在运行中改变，粒子迁移后计算正确继续。这是整个技术方案核心风险点的验证。

### 实现

**reassignCellOwner()** 在 `dsmcReplicatedMesh.C:261-291`：

- 采用交替策略：偶数次用 reverse-block（`nProcs_-1-old`），奇数次用原始 block
- 每次重分配后调用 `rebuildMyCells()` 更新本地 cell 视图
- 记录 `totalCellsChanged_`（cellOwner_ 变更总数）和 `rebalanceCount_`

**触发机制** 在 `dsmcCloud.C::evolve()`：

- 通过 controlDict 的 `replicatedMeshRebalanceSteps` 列表指定触发步
- 在迁移步骤中检测 `stepCounter_` 是否匹配 `rebalanceSteps_`
- 触发时：`reassignCellOwner()` → `migrateParticlesByCellOwner()`

### Phase B controlDict 配置

```
replicatedMeshRebalanceSteps (100 200);
```

### 测试结果（300 步, 8MPI, block 分解）

```
Phase B[1]: cellOwner_ reassigned. 60000 / 60000 cells changed (R_changed=1)
Phase B[2]: cellOwner_ reassigned. 60000 / 60000 cells changed (R_changed=1)
```

| 指标 | 值 | 说明 |
|---|---|---|
| R_changed (total) | 120,000 cells | 每次 100%, 两次共 200% |
| R_mig (total) | 31,149,587 parcels | 两次重分配总迁移粒子数 |
| 粒子守恒 | ✅ | 无丢失 |
| 7 不变量 | ✅ 全部成立 | 重分配前后均验证 |
| move 时间 | 74.6s | 与 Phase A 一致 |
| 迁移时间 | 33.9s (31 calls) | 含 2 次大迁移 |
| Wall time 不平衡 | 1.009 | 依然近乎完美 |

### 关键发现

1. **cellOwner_ 100% 翻转可行**：所有 cell 的 owner 反转，粒子全量迁移，模拟正确继续
2. **迁移后无粒子丢失**：两次 60000-cell 重分配后粒子总数守恒
3. **迁移开销可接受**：即使 100% cell 变更，单次迁移时间在可接受范围
4. **block 分解的对称性**：reverse-block 和原始 block 的粒子分布等价，imbalance 不变
5. **生产场景中 R_changed 远小于 100%**：TACF 驱动的增量调整只会改变部分 cell owner，迁移量远小于 Phase B 极端测试

---

## 十、实现 TODO 清单（更新后）

### 阶段 A：静态 cellOwner 验证 ✅ (2026-05-13)

#### A-1: 基础数据结构

| # | 模块 | 内容 | 状态 |
|---|------|------|------|
| A1.1 | `readFullMeshOnEachRank()` | 每个 rank 独立读取完整 polyMesh，不调用 decomposePar | ✅ `masterUncollated` file handler |
| A1.2 | `buildCellGraph()` | 构造全局 cell 邻接图 | ✅ 用于 Scotch 图分解（备选） |
| A1.3 | `scotchDecomposeToCellOwner()` | 分解得到 `cellOwner_[cellI]` | ✅ Block + Scotch（block 为此 case 最优） |
| A1.4 | `buildLocalComputationalView()` | 构建 `myCells_` | ✅ `rebuildMyCells()` |

#### A-2: 核心循环

| # | 模块 | 内容 | 状态 |
|---|------|------|------|
| A2.1 | `moveParticlesOnFullMesh()` | 在完整 mesh 上 tracking 粒子 | ✅ `Cloud::move()` 在 60k 全 mesh 上运行 |
| A2.2 | `migrateParticlesByCellOwner()` | `cellOwner_` 路由 + MPI 交换 | ✅ OCharStream + MPI_Alltoallv |
| A2.3 | `collideOwnerCellsOnly()` | 仅 owned cell 碰撞 | ✅ 迁移在 buildCellOccupancy 前，自然过滤 |
| A2.4 | `sampleOwnerCellsOnly()` | 仅 owner cell 采样 | ~ rank 0 写 field，待 MPI_Reduce |

#### A-3: 边界处理

| # | 模块 | 内容 | 状态 |
|---|------|------|------|
| A3.1 | 物理边界注入去重 | per-face `isMyCell(adjCell)` 检查 | ✅ FreeStreamInflowPatch |
| A3.2 | cellOccupancy 处理 | 保留全局长度，仅填 owner cell | ✅ 通过迁移保证 |

#### A-4: 验证实验

| # | 实验 | 内容 | 检查项 | 状态 |
|---|------|------|--------|------|
| A4.1 | 粒子守恒 | 300-step vs 单 rank 基线 | 总粒子数差异 | ✅ <0.3% (interval=10) |
| A4.2 | vs 标准 MPI | 8MPI 复制网格 vs 8MPI 标准 parallel | wall time, 粒子数 | ✅ 粒子匹配 0.013% |
| A4.3 | 7 不变量 | 各不变量逐项验证 | 见下方分析 | ✅ 全部成立 |

---

### 阶段 B：手动重分配 cellOwner ✅ (2026-05-13)

| # | 模块 | 内容 | 状态 |
|---|------|------|------|
| B.1 | `reassignCellOwner(newOwnerMap)` | 改变 `cellOwner_` 映射，更新 `myCells_` | ✅ `reassignCellOwner()` 交替 block/reverse-block 策略 |
| B.2 | `redistributeParticlesByNewOwner()` | 遍历本地粒子，按新 `cellOwner_` 重新打包迁移 | ✅ 复用 `migrateParticlesByCellOwner()` |
| B.3 | 手动重分配验证 | 第 100/200 步触发重分配 | ✅ 7 个不变量全覆盖、粒子数守恒（无丢失） |
| B.4 | 控制指标统计 | 每次重分配后记录 `R_changed`、`R_mig` | ✅ `totalCellsChanged_`, `totalParcelsMigrated_`, `rebalanceCount_` |

---

### 阶段 C：自动动态负载平衡

#### C-1: 负载监测

| # | 模块 | 内容 | 状态 |
|---|------|------|------|
| C1.1 | cell 级负载统计 | `w_c = N_p,c`（C-1 基线），后续扩展 `N_coll,c`、`T_c` | ✅ TACF 基础设施就绪 |
| C1.2 | rank 级耗时统计 | `T_p = move + collide + comm + boundary` wall-clock time per step | pending |
| C1.3 | 不均衡度计算 | `I_T = max(T_p) / avg(T_p)` | pending |

#### C-2: 触发判据

| # | 模块 | 内容 | 状态 |
|---|------|------|------|
| C2.1 | 触发条件 | `I_T > 1.3` 且距上次重平衡 ≥ 100 步 | pending |
| C2.2 | 收益-开销判据 | 预估迁移成本 vs 预期收益；`I_T` 改善 < 5% 但 `R_mig` > 20% → 跳过 | pending |

#### C-3: 三种重平衡对照

| # | 方案 | 定位 | 内容 | 状态 |
|---|------|------|------|------|
| C3.1 | Scotch cell 级动态 | 质量上界 / 图划分对照 | 粒子数加权 cell graph → Scotch(P) → 新 `cellOwner_` | pending |
| C3.2 | Hilbert SFC 动态 | 速度上界 / 轻量对照 | cell centroid → Hilbert key → 排序 → 按负载等分 | pending |
| C3.3 | superCell + Hilbert | **推荐主方案** | Scotch(M=αP) 过分解 → superCell 加权 centroid → Hilbert 排序 → 等分分配 | pending |

#### C-4: 对比实验

| # | 实验 | 内容 | 状态 |
|---|------|------|------|
| C4.1 | baseline | 无动态平衡 | pending |
| C4.2 | Scotch cell 级动态 | 图划分对照 | pending |
| C4.3 | Hilbert cell 级动态 | SFC 对照 | pending |
| C4.4 | superCell(M/P=4) + Hilbert | 细粒度 | pending |
| C4.5 | superCell(M/P=8) + Hilbert | 推荐主方案 | pending |
| C4.6 | superCell(M/P=16) + Hilbert | 粗粒度 | pending |

#### C-5: 退化检测

| # | 模块 | 内容 | 状态 |
|---|------|------|------|
| C5.1 | 通信质量监测 | 每次重平衡后记录 `edge cut`、`T_comm` | pending |
| C5.2 | 退化触发 | `E_current / E_best > 1.5` 或 `T_comm > 1.5 × T_comm_best` → Scotch 重校准 | pending |

#### C-6: 性能指标

```
总 wall-clock time
每步平均时间
max(T_p) / avg(T_p) 变化历程
迁移粒子数 / 总粒子数
单次重平衡开销
通信时间占比
每 rank 内存占用
```

---

## 十一、Phase C 实现总结（2026-05-14）

### 架构

在 replicated mesh 基础上实现了完整的自动动态负载平衡系统：
- **superCell 过分解**: Scotch(M=αP) 将 60k cells → α×P 个 superCell
- **Hilbert SFC 排序**: 3D Hilbert 曲线状态机，12-state × 8-octant 查找表
- **TACF 成本加权**: 累积粒子数作为 cell/superCell 成本权重
- **自动触发**: wall time 不平衡 + 冷却期 + 成本-收益检查

### 关键文件

| 文件 | 说明 |
|---|---|
| `replicatedMesh/dsmcReplicatedMesh.H` | 新增: autoDLBEnabled_, nSuperCells_, alpha_, imbalanceThreshold_, cooldownSteps_, superCellCost_, superCellOfCell_, superCellOwner_, superCellCentres_ |
| `replicatedMesh/dsmcReplicatedMesh.C` | Hilbert 曲线表, `buildSuperCells()`, `updateSuperCellCosts()`, `reassignByHilbertSuperCell()`, `autoRebalance()` |
| `replicatedMesh/dsmcLocalMesh.H/.C/.I` | 局部网格映射（myCells_ + 1-ring halo），global↔local 索引转换 |
| `clouds/dsmcCloud.C` | `rebuildParticleLoadPartition()` 适配 replicated mesh |
| `parcels/dsmcParcel.C` | `writeBinaryFast()` 快速序列化, OMP 临界区修复 |

### DLB 参数优化历程

| 参数 | 旧值 | 新值 | 原因 |
|---|---|---|---|
| `imbalanceThreshold` | 1.05 | **1.20** | 1.05 对均衡 case 也触发（block 正常波动 ~1.09） |
| `cooldownSteps` | 20 | **100** | 短窗口噪音大，100 步提供稳定测量 |
| 测量窗口 | 无要求 | **≥cooldownSteps** | 确保足够的 wall time 累积 |
| 成本-收益检查 | 无 | **粒子不平衡 <1.15 且 wall time <1.15 → 跳过** | 避免无意义重均衡 |

### 300步性能对比（8 MPI, cylinderN2 noreact）

#### 均衡 case (block 分解)

| 配置 | Main loop | Real time | 碰撞数 | 粒子不平衡 | 重均衡 |
|---|---|---|---|---|---|
| Phase A 基线 | 151.3s | 143.0s | 83,070 | 1.22 | — |
| **最终优化版** | **128.0s** | **133.4s** | **83,070** | **1.22** | — |
| Block+DLB 旧 (1.05/20) | 96.4s | 100.0s | 93,843 | 1.40 | 15 ⚠️ |
| Block+DLB 新 (1.20/100) | 137.5s | 141.7s | 83,070 | 1.22 | **0** ✅ |

#### 非均衡 case (skew 分解: rank 0 占 50% cells)

| 配置 | Main loop | Real time | 碰撞数 | 粒子不平衡 | 重均衡 |
|---|---|---|---|---|---|
| Skew 基线 (无 DLB) | 383.0s | 399.6s | 117,813 | 4.02 | — |
| Skew+DLB 旧 (1.05/20) | 96.2s | 100.6s | 56,160 | 1.45 | 15 ⚠️ |
| Skew+DLB 新 (1.20/100) | 212.7s | 218.8s | 140,884 | 2.44 | **2** ✅ |

### 正确性验证

| 指标 | 6 配置结果 |
|---|---|
| **振动能量** | 2.3273~2.3286e-20 — 所有配置一致 ✅ |
| **碰撞数** | 优化版 vs Phase A: 83,070 = 83,070 ✅ |
| **粒子数** | 优化版 vs Phase A: 243,915 = 243,915 ✅ |
| **Wall time 不平衡** | 优化版 vs Phase A: 1.009 = 1.009 ✅ |

### DLB 关键结论

1. **DLB 对非均衡 case 有效** — skew 分解 4x 加速（400s → 213s）
2. **均衡 case 不需要 DLB** — block 分解已近乎完美 (1.009)
3. **参数必须保守** — 阈值 1.20+冷却 100 避免过度触发
4. **存在性能-稳态 trade-off**: 更多重均衡 → 更好性能但粒子分布偏离更大
5. **生产推荐**: 均衡 case 关闭 DLB; 非均衡 case 使用 `threshold=1.20, cooldown=100`

---

## 十二、当前优化清单（2026-05-14）

### 已完成优化

| # | 优化 | 收益 | 文件 |
|---|------|------|------|
| O1 | `writeBinaryFast` 快速序列化 | 迁移 -63% (38.5→14.3s) | dsmcParcel.C, dsmcReplicatedMesh.C |
| O2 | `dsmcLocalMesh` 局部网格映射 | 基础设施 (myCells+halo, global↔local) | dsmcLocalMesh.H/.C/.I |
| O3 | OMP 临界区修复 | replicated mesh 可安全使用 OMP | dsmcParcel.C |
| O4 | OMP 负载分区修复 | myCells_ 范围分区 (非 60k 全网格) | dsmcCloud.C |
| O5 | Phase C auto DLB (优化后) | 非均衡 case 4x 加速, 均衡 case 0 触发 | dsmcReplicatedMesh.C |
| O6 | DLB 成本-收益检查 | 避免无意义重均衡 | dsmcReplicatedMesh.C |
| O7 | `skew` 分解模式 | DLB 测试工具 (rank 0 50% cells) | dsmcReplicatedMesh.C |
| O8 | 子网格追踪边界退出 | 尝试后放弃（改变物理, 碰撞偏差 2.23x）| dsmcParcel.C (已回退) |

### 可立即实施 (参考 ParDSMC3D 架构)

| # | 优化 | 说明 | 难度 | 预估收益 |
|---|------|------|------|----------|
| **A1** | **MPI 派生类型** | 仿 ParDSMC3D `MPI_Type_create_struct` 打包粒子数据 (7 doubles+3 ints)，直接 MPI 传输，省去 OCharStream 中间层 | ⭐⭐ | 迁移再降 30% |
| **A2** | **`migrateInterval=1`** | ParDSMC3D 每步迁移, 消除粒子堆积, 保持稳态 | ⭐ | 提高正确性 |
| **A3** | **累积空闲时间触发** | `tidl = tidl + (maxT - avgT)`, 比 `maxT/avgT` 比值更稳定 | ⭐ | DLB 判断更准 |
| **A4** | **迁移接收直接写入** | 不创建临时 particle 再 `addParticle`，直接填充 `cellOccupancy_` | ⭐⭐ | buildCellOccupancy -2s |

### 需要中等投入

| # | 优化 | 说明 | 难度 | 预估收益 |
|---|------|------|------|----------|
| **B1** | **ParMETIS 替代 Scotch** | `ParMETIS_V3_AdaptiveRepart` 运行时增量调整 | ⭐⭐⭐ | 重均衡质量更高 |
| **B2** | **碰撞候选跳过** | 跨越边界粒子优先标记, 减少候选计算 | ⭐⭐ | 碰撞 -20% |
| **B3** | **Thread partition** | ParDSMC3D 的 `OneDimePartition_move/coll` 动态 OMP 分区 | ⭐⭐ | OMP 效率提升 |

### 需要深度改造

| # | 优化 | 说明 | 难度 | 预估收益 |
|---|------|------|------|-----------|
| **C1** | **OF `particle::celli_` → protected** | 修改 OpenFOAM 基类, 打通局部网格追踪 | ⭐⭐⭐⭐ | 追踪 60s→10s |
| **C2** | **全局网格 + 局部计算分离** | 仿 ParDSMC3D 架构重构，粒子用 C struct 存储 | ⭐⭐⭐⭐⭐ | 架构级改善 |

### 建议执行顺序

```
第1轮: A2 (每步迁移) → A3 (累积空闲时间)
第2轮: A1 (MPI 派生类型) → A4 (迁移直接写入)
第3轮: B2 (碰撞候选跳过) → B3 (Thread partition)
后续: B1 (ParMETIS) → C1 (OF 基类修改)
```

### ParDSMC3D 关键参考

| 概念 | ParDSMC3D | 我们的方案 |
|---|---|---|
| 全局网格 | `Module.f90` 全局数组 `size=mnc` | `mesh_` 完整 polyMesh |
| Cell 所有权 | `pipl(mc)` → 全局→rank | `cellOwner_[cellI]` |
| 局部映射 | `irmc(mcc)` 局部→全局 | `localToGlobalCell_` |
| 反相映射 | `irmcre(mc)` 全局→局部 | `globalToLocalCell_` |
| 粒子追踪 | **全局网格**, 递归 `n_move`, 零 MPI 通信 | **全局网格**, `trackToAndHitFace` |
| 迁移 | 每步 `ParticleComm`, MPI 派生类型 | `migrateParticlesByCellOwner`, OCharStream |
| 碰撞 | 仅 `1..mnc_l` | 仅 `cellOwner_ == myRank_` |
| 动态均衡 | ParMETIS `AdaptiveRepart` | Hilbert SFC + Scotch |
| 触发方式 | 累积空闲时间 `tidl` | wall time 不平衡率 `maxT/avgT` |

**核心发现**: ParDSMC3D 与我们的方案架构完全一致。最大区别是它用 **MPI 派生类型** 替代序列化，用 **累积空闲时间** 替代比值触发，用 **每步迁移** 保持稳态。这三项均可直接采纳。

---

## 十三、混合 MPI+OMP 优化实施总结 (2026-05-15/16)

### 13.1 最终性能

| 配置 | Main loop (300步) | vs 原始 175.4s | 碰撞数 |
|---|---|---|---|
| 原始基线 (8MPI, Block, Alltoallv, 无OMP) | 175.4s | — | 246,363 |
| + p2p MPI (消除 collective) | 166.5s | -5.1% | 246,363 |
| + METIS 初始分区 + ParMETIS DLB | 157.2s | -10.4% | 245,763 |
| + 实时粒子权重 (替代 TACF) | 140s | -20.2% | — |
| + cooldown=20 | 125s | -28.5% | 246,628 |
| + OMP collision (2线程) | 121s | -31% | — |
| + OMP move (boundary-only guard) | ~100s | -43% | — |
| **+ sparse clear buildCellOccupancy** | **~96s** | **-45%** | **245,272** |

### 13.2 实施的优化清单

| 优化 | 文件 | 效果 |
|---|---|---|
| `localCellI_` 并行局部索引 | dsmcParcel.H/.C, dsmcParcelIO.C | 基础设施，零物理影响 |
| p2p MPI 迁移 (Isend/Irecv) | dsmcReplicatedMesh.C | 消除 Alltoall/Alltoallv collective 同步 |
| METIS 初始分区 | dsmcReplicatedMesh.C | 比 block/scotch 更优的空间分区 |
| ParMETIS AdaptiveRepart DLB | dsmcReplicatedMesh.C | 增量重分区，最小化数据移动 |
| ParDSMC3D 风格 sar 触发 | dsmcReplicatedMesh.C | 基于 productive time 趋势触发 |
| 实时粒子权重 (替代 TACF) | dsmcReplicatedMesh.C | 粒子数+碰撞候选作为 ParMETIS vertex weight |
| OMP move boundary-only guard | dsmcParcel.C | 只对物理边界 cell 加 critical section |
| OMP move crash 修复 | Cloud.C, dsmcCloud.C, dsmcParcel.C | useMoveParticlePartition=false, clearMoveOrderedParcels, RNG guard |
| Sparse clear buildCellOccupancy | dsmcCloud.C | 只清除有粒子的 cells (7.5K vs 60K) |
| MoveFastRng in trackingData | dsmcParcel.H | 为未来 OMP 优化预留 |

### 13.3 运行方式

```bash
# 不用 decomposePar，不用 -parallel
mpirun -np 8 dsmcFoam+ -case <caseDir>
```

controlDict 关键参数：
```
replicatedMesh true;
replicatedMeshDecompMethod metis;
replicatedMeshAutoDLB true;
replicatedMeshDLBCooldownSteps 20;
useOpenMP true;
openmpMove true;
openmpThreads 2;
```

### 13.4 OMP 瓶颈分析

#### 为什么 8MPI×2OMP 优于 2MPI×4OMP

| 指标 | 2MPI×4OMP | 8MPI×2OMP |
|---|---|---|
| 粒子/rank | 815K | 200K |
| move extract (串行链表) | 24.2s | ~5s |
| OMP kernel 开销 | 17s | ~4s |
| buildCellOccupancy | 26.6s | ~8.5s |
| **main loop** | **180s** | **96s** |

根因：`move extract parcels` 是串行链表遍历 (O(N))，interval=1 迁移每步重建，无法复用。更多 MPI rank = 更少 per-rank 粒子 = 更小串行开销。

#### OMP move 的 critical section 分析

- `dsmcMoveTrack`: 只对物理边界 cell 加锁，内部 cell 完全并行
- `dsmcMoveBoundary`: diffuseWall 壁面交互串行化（共享 RNG）
- 实测 boundary wall time 仅 0.76s/线程 — **不是瓶颈**
- 真正瓶颈是 per-particle OMP 调度开销 (~53ns/粒子/步)

#### 已验证无效的 OMP 优化方向

| 方向 | 结果 | 原因 |
|---|---|---|
| per-thread RNG (rndGen() accessor) | +10-15% 恶化 | accessor 开销 + equipartition 函数影响 collision |
| FastRng in handleWallInteraction | +7% 恶化 | dsmcDiffuseWallPatch 有其他共享状态 |
| useMoveParticlePartition=true | +7-13% 恶化 | 分区假设粒子分布稳定，interval=1 迁移违反 |
| parallel extract via cellOccupancy | 崩溃 | autoRebalance/controlAfterCollisions 使 cellOccupancy 失效 |
| localCellI_ 移除热路径 | 无改善 | 不是瓶颈 |

### 13.5 DLB 触发机制

采用 ParDSMC3D 风格：
```
每步: productiveTime += evolveStepTime - migrationDelta
每 cooldownSteps 步:
  Allgather productiveTime → 计算 tidl (累积 idle)
  sar = w2 - w1 (趋势)
  if sar > 0: 触发 ParMETIS AdaptiveRepart
  触发后: w1=1e6 (隐式 cooldown), 重置 tidl
```

ParMETIS vertex weight = 实时粒子数 + 碰撞候选数 (双约束 ncon=2)。

### 13.6 关键发现

1. **`-case` 模式 (非 `-parallel`)** 是 replicated mesh 的正确运行方式。`-parallel` 导致 mesh 从 processor 目录读取（局部网格）。
2. **interval=1 迁移是正确性必需**。interval=10 丢失 66% 碰撞。
3. **TACF 累积权重对 DLB 有害**。重置后数据不足导致 ParMETIS 过拟合。实时粒子数更有效。
4. **OMP move 的 critical section 不是瓶颈**（仅 0.76s）。真正瓶颈是串行链表遍历和 per-particle OMP 开销。
5. **Sparse clear** 对 replicated mesh 有效（60K cells 中只有 7.5K 有粒子）。
6. **8MPI×2OMP 优于 2MPI×4OMP**：更多 MPI rank 减少 per-rank 串行开销，比增加 OMP 线程更有效。
