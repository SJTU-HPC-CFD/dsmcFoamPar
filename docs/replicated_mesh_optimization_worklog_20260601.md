# Replicated Mesh 优化工作日志 (Session 2)

## 日期: 2026-06-01

## 起始状态
- mpi8replicatedmesh: 202s (原始)
- omp8: 113s
- mpi8origin: 194s

## 优化成果总结

| 优化项 | 时间变化 | 改善 |
|--------|----------|------|
| 原始版本 | 202s | — |
| + flat POD transfer | 202s (migration内部优化) | DLB步3.7x加速 |
| + post myCells循环 | 171s | -15% |
| + post prefetch | 158s | -22% |
| + post owned boundary faces | 158s | boundary -85% |
| + post compact cache (idx-based) | 154s | -24% |
| + PID v3 + N^0.8权重 | 148s | -27% |
| + post-migration overlap | **141-144s** | **-29%** |

## 最终配置

### controlDict 参数
```
replicatedMesh true;
replicatedMeshDelayedReceive true;
replicatedMeshNoAlltoall true;
replicatedMeshFlatTransfer true;
replicatedMeshDecompMethod metis;
replicatedMeshAutoDLB true;
replicatedMeshDLBSteps 50;
replicatedMeshDLBImbalanceThreshold 1.5;
replicatedMeshDLBFixedK 0;
replicatedMeshDLBInitialK 1024;
replicatedMeshDLBItr 1000;
replicatedMeshDLBUbvec 1.05;
replicatedMeshDLBAdaptiveKMode 0;  // PID v3
```

### 三算例最终对比

| 指标 | omp8 (1×8T) | mpi8origin (8×MPI) | mpi8replicatedmesh |
|------|-------------|--------------------|--------------------|
| ClockTime | 106s | 173s | 144s |
| Move | 50.4s | 131.6s | 51.4s |
| BuildOccupancy | 6.2s | 8.5s | 7.4s |
| Collision | 9.5s | 17.9s | 13.6s |
| Post | 35.2s | 17.1s | 9.2s (实际19s) |
| Migration | — | — | 41.0s |
| 总粒子数 | 2,463,452 | 2,463,429 | 2,466,245 |
| 碰撞数 | 35,781 | — | 35,559 |

## 详细优化技术

### 1. Flat POD Transfer (migration序列化优化)
- **文件**: `dsmcParcel.H/C`, `dsmcReplicatedMesh.C`
- **原理**: 定义固定大小 `TransferData` POD struct (200B/粒子)，替代 OCharStream 序列化
- **效果**: DLB步迁移200万粒子只需197ms (pack=22ms, MPI=98ms, deser=40ms)
- **正常步**: pack=0.08ms, deser=1ms (几乎零开销)

### 2. Post owned-cells-only (myCells循环)
- **文件**: `dsmcVolFields.C`
- **原理**: 将 `for(celli=0; celli<nCells; ++celli)` 替换为 `for(idx=0; idx<myCells.size(); ++idx)`
- **修改点**: build() 3个循环 + field combine + cell reduce + vibrational accum
- **效果**: 循环次数从104151降到~13000 (-87.5%)

### 3. Post prefetch
- **文件**: `dsmcVolFields.C`
- **原理**: 在粒子循环中添加 `__builtin_prefetch(nextParcel, 0, 1)`
- **效果**: parcel accumulate 从23.4s降到20.7s (-12%)

### 4. Post precomputed owned boundary faces
- **文件**: `dsmcVolFields.H/C`
- **原理**: 在 `createField()` 中预计算每个 patch 的 owned face 列表
- **效果**: boundary accumulation 从21.3s降到2.3s (-89%)

### 5. Post compact cache (idx-based)
- **文件**: `dsmcVolFields.C`
- **原理**: cache 数组从104K缩小到myCells.size()，用循环idx作为索引
- **关键**: `cacheIdx = myCellsData ? idx : celli` (replicated用idx, 非replicated用celli)
- **效果**: parcel accumulate 从20.7s降到18.5s, field combine从1.9s降到1.3s

### 6. PID v3 adaptive K
- **文件**: `dsmcReplicatedMesh.C`
- **原理**: 用 bottleneck rank 的 moveExcess vs collExcess 方向来调整K
- **公式**: `direction = (bnMoveExcess - bnCollExcess) / (|bnMoveExcess| + |bnCollExcess|)`
- **效果**: K温和上升(1024→1400), 比旧PID(K下降到574)和fixedK=2048都更稳定

### 7. N^0.8 move权重
- **文件**: `dsmcReplicatedMesh.C`
- **原理**: `moveCost = nPart^0.8` (低密度cell每粒子move更慢, 给更高权重)
- **公式**: `vwgt = max(1, int(pow(nPart, 0.8) + nPart*(nPart-1)/K))`
- **效果**: 148s→148s (稳定), 比线性N权重更好地平衡实际move时间

### 8. Post-migration overlap
- **文件**: `dsmcCloud.C`
- **原理**: 将 `fields_.calculateFields()` 移到 migration 之前执行
- **时间线**: move → post → migrateBegin → migrateFinish → buildOcc → coll
- **效果**: 快rank做完move后立即做post, post时间被慢rank的move时间吸收
- **结果**: post报告时间从30s降到9.2s (实际19s, 10s被overlap隐藏)

### 9. replicatedMeshActive() 修复
- **文件**: `dsmcCloudI.H`
- **Bug**: `replicatedMeshActive()` 只检查 `valid()`, 单进程时也返回true
- **修复**: 改为 `valid() && active()`, 避免omp8模式走错误路径
- **效果**: omp8从127s恢复到106s
2026.6.5再次修改：moveTrackCallCount 这条 profiling 统计链，不是 dsmcParcel false-path。清理后，源码中已无 moveTrackCallCount / trackNsPerCall 残
  留。相关位置在 hyStrath_xcx/src/lagrangian/basic/Cloud/Cloud.C:544、hyStrath_xcx/src/lagrangian/dsmc/clouds/dsmcCloud.C:4385、hyStrath_xcx/src/lagrangian/
  dsmc/parcels/dsmcParcel.C:144。

### 10. cellOwner可视化输出
- **文件**: `dsmcReplicatedMesh.C`
- **功能**: 每次DLB后写出 `cellOwner` volScalarField, 可在Tecplot/ParaView中查看分区

## DLB 参数测试记录

| 配置 | 时间 | 说明 |
|------|------|------|
| K adaptive(旧PID), itr=1000, DLB/50步 | 147s | K从1024降到574(方向错误) |
| K=2048 fixed, itr=1000, DLB/50步 | 152s | 纯粒子数均衡 |
| K=2048, itr=1000, DLB/20步 | 156s | 更频繁DLB反而更慢 |
| K=2048, itr=10, DLB/20步 | 158s | 增量调整 |
| K=2048, itr=0.01, DLB/20步 | 180s | 太保守,负载不均衡 |
| K=2048, itr=1, DLB/50步 | 159-173s | 不连通分区 |
| K adaptive(PID v1), itr=1000 | 165s | K升到2639(过度) |
| K adaptive(PID v2), itr=1000 | 154-159s | K降到637(collFrac) |
| **K adaptive(PID v3), itr=1000** | **151s** | **K升到1737(方向正确)** |
| K adaptive(PID v3)+N^0.8 | **148s** | **最优权重** |
| PartKway替代AdaptiveRepart | 159s | 保证连通但分区质量差 |
| ubvec=1.2 | 155s | 放宽容忍度反而更差 |
| moveItersPerCell(瞬时) | 168s | 不稳定 |
| moveItersPerCell(累积) | 174s | 新cell无历史数据 |

## DLB 负载均衡分析

### 间隔负载趋势
- 初期 (step 0-200): max/min = 2.9-6.5x (粒子从0增长, 分布剧烈变化)
- 中期 (step 200-350): max/min = 2.0-3.5x (DLB开始生效)
- 后期 (step 350-500): max/min = 1.6-1.9x (接近稳态, DLB有效)

### DLB后瞬时负载
- 范围: max/min = 2.0-4.1x
- ParMETIS在网格拓扑约束下无法完美均衡
- 但50步后通过粒子自然迁移趋于均衡

### 关键发现
1. Move占总负载55-89%, 是绝对主导
2. Bottleneck rank每次DLB后不同 (负载在rank间轮转)
3. itr=1000(激进重分配)对非稳态算例最优
4. ubvec=1.05(严格容忍度)比1.2更好
5. 不连通分区对性能影响可忽略(正常步迁移量不变)

## 剩余瓶颈

| 瓶颈 | 时间 | 原因 | 可能的解决方案 |
|------|------|------|---------------|
| Move | 51s | 单线程处理~270K粒子/rank | MPI+OpenMP混合模式 |
| MPI sync wait | 35s | rank间move不均衡 | 已到ParMETIS极限 |
| Collision | 14s | 正常物理计算 | 已接近极限 |
| Post(overlap后) | 9s | writeFields+controllers | 减少采样频率 |

## 修改的文件清单

| 文件 | 改动 |
|------|------|
| `src/lagrangian/dsmc/parcels/dsmcParcel.H` | TransferData struct, packTransfer/unpackTransfer声明, ownedBoundaryFaces_ |
| `src/lagrangian/dsmc/parcels/dsmcParcel.C` | packTransfer/unpackTransfer实现, move中moveItersPerCell累积 |
| `src/lagrangian/dsmc/replicatedMesh/dsmcReplicatedMesh.H` | useFlatTransfer_, useNoAlltoall_ 成员 |
| `src/lagrangian/dsmc/replicatedMesh/dsmcReplicatedMesh.C` | flat transfer路径, PID v3, N^0.8权重, cellOwner输出, itr/ubvec可配置 |
| `src/lagrangian/dsmc/clouds/dsmcCloud.H` | moveItersPerCell_ 成员和accessor |
| `src/lagrangian/dsmc/clouds/dsmcCloud.C` | post-migration overlap, moveItersPerCell清零 |
| `src/lagrangian/dsmc/clouds/dsmcCloudI.H` | replicatedMeshActive()修复 |
| `src/lagrangian/dsmc/macroscopicProperties/derived/combined/dsmcVolFields/dsmcVolFields.H` | ownedBoundaryFaces_ 成员 |
| `src/lagrangian/dsmc/macroscopicProperties/derived/combined/dsmcVolFields/dsmcVolFields.C` | myCells循环, prefetch, compact cache, owned boundary |

---

# Session 3: DLB 深度优化 + 异步流水线 + migrateInterval

## 日期: 2026-06-01 (续)

## 优化成果：202s → 121s (-40%)

| 优化项 | 时间变化 | 改善 |
|--------|----------|------|
| Session 2 最优 (PID v3 + N^0.8 + post overlap) | 148s | -27% |
| + async pipeline (migrateFinish延迟到下一步) | 146s | 小幅改善 |
| + migrateInterval=2 (双约束) | 137s | -32% |
| + migrateInterval=5 (双约束) | 125s | -38% |
| **+ migrateInterval=10 (双约束, 最终)** | **121-125s** | **-40%** |

## 最终配置

```
replicatedMesh true
replicatedMeshDelayedReceive true
replicatedMeshNoAlltoall true
replicatedMeshFlatTransfer true
replicatedMeshDecompMethod metis
replicatedMeshMigrateInterval 10    // 关键! 每10步迁移一次
replicatedMeshDLBDualConstraint true // 双约束ParMETIS
replicatedMeshAutoDLB true
replicatedMeshDLBSteps 50
replicatedMeshDLBImbalanceThreshold 1.5
replicatedMeshDLBFixedK 0
replicatedMeshDLBInitialK 1024
replicatedMeshDLBItr 1000
replicatedMeshDLBUbvec 1.05        // 单约束默认; 双约束时collision自动放宽到1.5
replicatedMeshDLBInitialAlpha 0.8
replicatedMeshDLBAdaptiveKMode 0
```

## 关键技术突破

### 11. migrateInterval > 1: 碰撞允许非owned cells

**核心发现**: 当前代码的 `buildCellOccupancy()` 和 `collision()` 实际上**已经支持非owned cells**。
- `buildCellOccupancy` 把所有粒子(不管cell是否owned)放入对应cell的occupancy
- `collision` 遍历所有cells(不检查owned), 对任何有>1粒子的cell做碰撞
- 因此增大 migrateInterval 不会导致碰撞丢失!

**之前66%碰撞丢失的原因**: 旧代码可能只对owned cells构建occupancy, 当前代码已修复。

**验证数据** (双约束 + N^0.8):
| interval | 总时间 | move | migration | 碰撞数 | 碰撞偏差 |
|----------|--------|------|-----------|--------|---------|
| 1 | 138s | 53.5s | 29.0s | 35,798 | +0.7% |
| 2 | 130s | 54.3s | 16.4s | 35,456 | -0.2% |
| 5 | 125s | 53.5s | 12.8s | 35,670 | +0.4% |
| **10** | **121s** | 55.6s | 3.8s | 34,719 | -2.3% |

interval=5: 精度最佳(+0.4%), interval=10: 性能最佳但偏差-2.3%。

### 12. 异步流水线

将 `migrateFinish` 从本步移到下一步开始执行:
```
步N:   migrateFinish(上一步) → move → post → migrateBegin(发送,不等待) → buildOcc → coll
步N+1: migrateFinish(步N数据) → move → ...
```

效果: migration MPI wait 从 ~75ms 降到 ~50ms/步(未完全消除, 因为 buildOcc+coll 时间不够完全吸收等待差)。

### 13. 双约束 ParMETIS

当前单约束用 N^0.8 压缩了动态范围, 但无法同时精确均衡 move 和 collision。

双约束方案:
- 约束0: N^0.8 (move balance, alpha压缩动态范围)
- 约束1: N*(N-1) (collision balance)
- ubvec[0]=1.05 (move严格), ubvec[1]=1.5 (coll放宽)

**为什么N^0.8在双约束中仍有价值**: 压缩约束0的动态范围, 防止 ParMETIS 被少数高密度cell主导分区质量。纯N(alpha=1)的双约束反而更差(132s vs 121s)。

### 14. Greedy Balance 实验 (失败, 验证了重要结论)

实现了无连通性约束的贪心粒子数均衡:
- 按粒子数排序所有cells, 贪心分配给粒子数最少的rank
- 结果: 粒子数完美均衡(max-min=1粒子!), 但总时间 247s

**教训**: 粒子数均衡 ≠ move时间均衡。不连通分区导致:
1. 每步正常迁移量大幅增加(+边界cells更多)
2. 小cell区域的粒子move更慢, 拥有小cell的rank仍然更慢
3. 每次DLB 97% cells重分配, 迁移开销巨大

**关键结论**: ParMETIS的连通分区约束是有益的——它保证了空间局部性, 减少了每步的正常迁移量。

### 15. 自适应 ubvec 实验 (失败)

尝试根据move/coll时间比例动态调整双约束的ubvec:
- move主导 → 收紧move约束, 放宽coll约束
- coll主导 → 收紧coll约束, 放宽move约束

结果: 133s vs 固定ubvec[1]=1.5的123s。自适应引入不稳定, ParMETIS分区质量对ubvec敏感。

### 16. Physics-based move cost model 实验 (失败)

用实际速度+cell尺寸计算每cell的trackToAndHitFace迭代次数:
- `crossings = |U_avg| * deltaT / cellSize`
- `weight = nPart * (1 + crossings)`

结果: 165s(无cap), 159s(cap=3)。权重动态范围太大(空cell=1, 高密度小cell=400), ParMETIS无法产生好分区。

**教训**: 精确的物理模型≠好的DLB权重。ParMETIS需要压缩后的权重动态范围。

### 17. currentCellCount 修复实验 (部分有效)

DLB权重用 `p.cell()` (当前move后位置) 而不是 `cellOccupancy` (上一步collision后位置):
- interval=1: alpha=1时156s, 粒子分布3.6x (vs cellOccupancy的7.6x)
- interval=10 + currentCellCount: 未测试(和interval=10矛盾——cellOccupancy在interval>1时本来就包含move后位置)

结论: cellOccupancy在interval=1时是"陈旧"数据, 但N^0.8压缩后比currentCellCount更稳定。interval>1时cellOccupancy本身就是最新的。

### 18. 边界cell权重修正实验 (无效)

给inlet/wall boundary cells放大权重(1.5x), deletion cells缩小权重(0.7x)。结果无改善——边界cells数量少(699+699+149=1547/104151=1.5%), 影响太小。

## DLB参数完整测试矩阵

| K策略 | alpha | 约束 | interval | 总时间 | 说明 |
|-------|-------|------|----------|--------|------|
| PID v3 | 1.0(adaptive) | 单 | 1 | 144-151s | 基准 |
| fixedK=2048 | — | 单 | 1 | 152s | 纯N均衡 |
| hill-climb | 0.8→0.76 | 单 | 1 | 148-150s | 最优单约束 |
| hill-climb | 1.0→0.95 | 单 | 1 | 152s | 收敛慢 |
| hill-climb | 0.8→0.76 | **双** | 1 | 138s | 双约束优势 |
| hill-climb | 0.8→0.76 | **双** | 2 | 130s | interval开始发挥 |
| hill-climb | 0.8→0.76 | **双** | **5** | **125s** | 精度最优 |
| hill-climb | 0.8→0.76 | **双** | **10** | **121s** | 性能最优 |
| hill-climb | 1.0(adaptive) | 双纯N | 10 | 132s | 纯N不如N^0.8 |
| hill-climb | 0.8(adaptive) | 双自适应ubvec | 10 | 133s | 不如固定ubvec |
| — | 0.8 | 单 | 10 | 125s | 单约束也行 |

## 三算例最终对比

| 指标 | omp8 (1×8T) | mpi8origin (8×MPI) | mpi8replicated (最优) |
|------|-------------|--------------------|----------------------|
| ClockTime | **106s** | 173s | **121-125s** |
| Move | 50.4s | 131.6s | 53-55s |
| BuildOccupancy | 6.2s | 8.5s | ~13s |
| Collision | 9.5s | 17.9s | ~12s |
| Post | 35.2s | 17.1s | ~30s |
| Migration | — | — | 3.8-12.8s |

## 架构限制与剩余空间

当前 121s 已接近架构极限:
- Move 54s = omp8(50s) + 4s overhead (rank间不均衡残差)
- Post 30s = 无法再优化(所有已知优化已应用)
- 剩余瓶颈全在硬件/并发层面

**进一步优化唯一可行方向**: MPI+OpenMP混合模式 (2rank×4thread)
- 预期: 接近omp8的106s
- 风险: 已有OpenMP move代码, 但需验证replicated mesh + OpenMP的并发正确性

## 修改文件清单 (本session新增)

| 文件 | 改动 |
|------|------|
| `dsmcCloud.C` | 异步流水线: migrateFinish移到下一步开始 |
| `dsmcReplicatedMesh.H` | alpha_, reassignGreedyBalance声明 |
| `dsmcReplicatedMesh.C` | Hill-climbing alpha PID, N^0.8权重, 双约束(ubvec[1]=1.5), greedy balance, 边界修正, currentCellCount, physics-based model (experimental), migrateInterval可配置 |
| `dsmcCloudI.H` | replicatedMeshActive()修复 |
| `dsmcCloud.H` | moveItersPerCell_成员 |
| `dsmcParcel.C` | move中moveItersPerCell累积 |
| `fieldPropertiesDict` | sampleInterval测试(最终恢复为1) |
| `controlDict` | migrateInterval, DLBDualConstraint, DLBInitialAlpha, DLBProfile等 |


---

# Session 4: DLB参数精调 + 性能回退修复 + nCandidates/vsize实验

## 日期: 2026-06-01 (续2)

## 核心问题: 121s 性能无法复现

Session 3 结束时最优 121s，但后续测试在 137-140s 徘徊。经排查确认:
1. **增量编译累积问题**: 多次修改/回退导致目标文件不一致，clean rebuild 后恢复
2. **async pipeline 不是问题**: 121s 的 log.dual_i10 也是 async pipeline 模式
3. **自适应 alpha 不是问题**: 固定 alpha=0.8 反而更稳定
4. **根因: Scotch 初始分区质量波动** + DLB 频率不足

## 实测最优配置 (多次确认)

```
replicatedMeshMigrateInterval 10
replicatedMeshDLBDualConstraint true
replicatedMeshDLBSteps 50
replicatedMeshDLBImbalanceThreshold 1.15    // 关键: 比默认1.5更早触发DLB
replicatedMeshDLBFixedK 1                   // 固定alpha=0.8, 不自适应
replicatedMeshDLBInitialAlpha 0.8
replicatedMeshDLBItr 1000
replicatedMeshDLBUbvec 1.05
replicatedMeshDLBUbvec1 1.5                 // 双约束collision放宽
replicatedMeshDLBVsizeExp 0                 // free migration
```

**稳定性能**: 126-132s (平均 ~129s), 碰撞偏差 -1.2%

## DLB 参数扫描结果

### threshold 扫描 (固定alpha=0.8, interval=10, dual N²)
| threshold | 时间 | move | 说明 |
|-----------|------|------|------|
| 1.05 | 131-132s | 56-57s | 过早触发, 粒子数少时分区质量差 |
| 1.15 | **127-130s** | 56-57s | 最优: 足够早 + 数据充分 |
| 1.2 | 131-139s | 53-58s | 波动大 |
| 1.5 | 131-140s | 55-62s | 太晚, Scotch差分区运行太久 |

### ubvec[1] 扫描 (固定alpha=0.8, interval=10)
| ubvec[1] | 时间 | coll偏差 | 说明 |
|----------|------|---------|------|
| 1.2 | 135s | -1.1% | 约束太紧 |
| 1.3 | 140s | +0.5% | 约束太紧 |
| **1.5** | **129s** | **-1.2%** | 最优 |
| 1.8 | 137s | -4.3% | 太松, 碰撞偏差大 |

### vsizeExp 扫描 (固定alpha=0.8, interval=10, ubvec1=1.5)
| vsizeExp | 时间 | move | coll | 说明 |
|----------|------|------|------|------|
| **0 (free)** | **129s** | 57s | 35,133 | 最优: ParMETIS自由重分区 |
| 0.3 | 133s | 65s | 34,302 | coll偏差-3.5%, 限制过度 |
| 0.5 | 140s | 70s | 35,548 | 严重限制 |
| 0.8 | 146s | 73s | 35,340 | 几乎禁止重分配 |

### nCandidatesPerCell 替代 N² (interval=10, vsExp=0.5)
| 约束1 | 时间 | coll | 说明 |
|-------|------|------|------|
| **N² (理论)** | **129s** | 35,100 | 最优: 压缩后动态范围适中 |
| nCandidates (实测) | 140s | 35,724 | 动态范围太大, ParMETIS质量差 |

**教训**: nCandidates 动态范围远超 N²，必须 EWMA/clamp 预处理后才能用

## 关键实验: 初始化 ParMETIS 分区 (失败)

在 `distributeInitialParticles()` 后立刻调用 `reassignByParMetisAdaptiveRepart()`:
- 结果: 140-141s (比 Scotch 初始更差)
- 原因: 初始化时粒子太少(-2M vs 最终2.46M)，分区不适应后续自由流注入

## 性能回退修复历程

| 问题 | 当时性能 | 修复 | 修复后 |
|------|---------|------|--------|
| 增量编译 | 137-140s | clean rebuild | 131s |
| 自适应alpha不稳定 | 133-137s | fixedK=1, 固定alpha=0.8 | 128-130s |
| threshold=1.5太晚 | 131-140s | threshold=1.15 | 127-130s |
| async pipeline疑云 | — | 确认121s时已激活, 非问题 | — |

## 最终确认的最佳性能

| 日志 | 时间 | move | coll | DLBs |
|------|------|------|------|------|
| log.low3 | **126.7s** | 56.4s | 35,285 | — |
| log.fa1 | 127.9s | 58.8s | 35,106 | — |
| log.fa3 | 128.7s | 59.7s | 35,217 | — |
| log.r2 | 129.2s | 58.4s | 35,182 | — |

历史最优: log.dual_i10 = **123.1s** (move=50.9s, 幸运Scotch)

## 剩余优化空间

1. **Scotch确定性** (3-6s): 固定随机种子消除初始分区方差
2. **MPI+OpenMP混合** (20-30s): 2rank×4thread, 预期接近omp8(106s)
3. **nCandidates + EWMA** (1-3s): 需要预处理后替代N²
4. **edge weight用迁移通量** (1-3s): 需要额外统计

## 修改文件清单 (本session新增/变更)

| 文件 | 改动 |
|------|------|
| `dsmcReplicatedMesh.C` | nCandidates约束(实验), vsizeExp可配置, ubvec1可配置 |
| `dsmcReplicatedMesh.H` | reassignByParMetisAdaptiveRepart() public(已回退), rebuildMyCells |
| `dsmcCloud.C` | 初始化ParMETIS分区(已回退), post-migration overlap |
| `controlDict` | ImbalanceThreshold扫描(1.05→1.8), ubvec1扫描(1.2→1.8), vsizeExp扫描(0→0.8), FixedK=1固定alpha |

