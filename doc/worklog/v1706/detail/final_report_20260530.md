# DSMC OpenMP 并行化最终报告

日期：2026-05-30/31

## 最终性能（React 2D cylinder, ~2.3M 粒子, ~104K cells）

### 最终验证结果（2026-05-31 代码清理后）

| 方法 | total(s) | 加速比 | move | coll | fields | build | nParticles | 正确性 |
|------|----------|--------|------|------|--------|-------|-----------|--------|
| Serial (OMP-1) | 282.6 | 1.00x | 182.2 | 38.0 | 32.9 | 29.3 | 2,320,669 | 基准 |
| **OMP-2** | **142.6** | **1.98x** | 93.3 | 15.6 | 29.1 | 4.3 | 2,320,621 | ✓ |
| **OMP-4** | **106.0** | **2.67x** | 60.8 | 10.5 | 30.2 | 4.2 | 2,320,721 | ✓ |
| **OMP-8** | **89.8** | **3.15x** | 49.4 | 6.5 | 29.0 | 4.6 | 2,320,584 | ✓ |
| MPI-2 | 147.5 | 1.92x | 102.8 | 17.3 | 20.6 | 6.6 | 2,320,779 | ✓ |
| MPI-4 | 111.3 | 2.54x | 81.8 | 12.6 | 13.5 | 3.3 | 2,320,527 | ✓ |
| MPI-8 | 83.5 | 3.39x | 57.4 | 8.0 | 13.4 | 4.5 | 2,320,498 | ✓ |

**所有粒子数正确！OMP-8 = 3.15x。**

- OMP-2 (1.98x) > MPI-2 (1.92x) ← **OMP 胜出**
- OMP-4 (2.67x) > MPI-4 (2.54x) ← **OMP 胜出**
- OMP-8 (3.15x) vs MPI-8 (3.39x) ← 差距 7%

## OMP-8 阶段分解

| 阶段 | 耗时 | 占比 | 加速比 | 方式 |
|------|------|------|--------|------|
| move | 50-55s | 54% | 3.3-3.6x | trackToFace cell级并行 + boundary串行 |
| fields | 29-34s | 32% | 1.0x | 串行 forAllConstIter（架构限制） |
| buildPtrs+Idx | 5-16s | 5-14% | 1.8-6.0x | 部分并行+延迟删除 |
| collision | 6-8s | 7% | 5-6x | cell级并行 + reaction critical section |

## 里程碑完成状态

| M | 状态 | 说明 |
|----|------|------|
| M0 | ✅ | wall-time phase timing |
| M1 | ✅ | parcelPtrs_ + cellFirst/cellNext/cellCount |
| M2 | ✅ | noTimeCounterOMP, 5-6x collision加速 |
| M3 | ⚠️ | cellFirst/cellNext 遍历 cache miss，串行最优 |
| M4 | ✅ | dynamic schedule <0.01% imbalance |
| M5 | ✅ | moveParallel: trackToFace并行+boundary串行 |
| M6 | ✅ | rndGen() 自动路由 + measurement守护 |
| M7 | ✅ | reaction critical section |
| M8 | ✅ | cellOccupancy消除 + 延迟删除 |
| M9 | ❌ | MPI DLB接口（功能性需求） |

## 修改文件

### 新增
- `collisionPartnerSelection/derived/noTimeCounterOMP/` — OMP 碰撞模型

### 核心修改
- `clouds/dsmcCloud.H/C/I` — 并行数据结构+move+计时+RNG路由
- `basic/Cloud/Cloud.H` — trackingRescue原子递增
- `basic/particle/particleTemplates.C` — cloud.labels()→局部变量
- `parcels/dsmcParcel.H` — trackingData加threadId
- `boundaries/.../dsmcPatchBoundary.C` — measureProperties守护
- `boundaries/.../dsmcDeletionPatch.C` — deletionInteraction守护
- `reactions/.../dissociationQK.C` — outputResults用cellFirst/cellNext
- `Make/options`, `Make/files`, `basic/Make/options` — 编译配置

## 关键技术决策

1. **Fields 无法有效并行**：cellFirst/cellNext按cell遍历时粒子内存不连续→严重cache miss。sorted index的构建开销抵消并行收益
2. **trackToFace根因**：cloud.labels()是共享scratch buffer，改为局部变量解决
3. **延迟删除**：dead标记+compact，省一次IDLList遍历
4. **rndGen()自动路由**：omp_in_parallel()→并行区内自动per-thread RNG
5. **负载均衡<0.01%**：dynamic schedule已完美均衡

## 限制和后续方向

1. **IDLList串行遍历**（~5s/步）：需DsmcParticleStore替代
2. **Fields 32%瓶颈**：需粒子物理排序（DsmcParticleStore）
3. **Boundary串行**（~3-5s）：需thread-local measurement buffer
4. **架构目标**：DsmcParticleStore → 预期3.5-4.5x
