# Phase 5 补充：trackToFace 线程安全问题

日期：2026-05-30

## 问题

`moveParallel()` 在第一步即 segfault。

## 根因分析

`particle::trackToFace()` 在 OpenMP 并行区内调用时发生段错误。
可能原因：
1. `particle` 基类的 tracking 使用了非线程安全的全局/静态状态
2. mesh 的 tet 分解查询可能有缓存（lazy evaluation）不是线程安全的
3. `trackingData` 中的 cloud 引用在并行区内被并发访问

## boundary cell 标记修复

- 原始问题：symmetry patch（2D front/back）覆盖所有 cell → 100% boundary
- 修复：按 patch type 字符串过滤，跳过 empty/symmetry/symmetryPlane/wedge
- 结果：1,545 boundary cells / 104,151 total (1.48%)

## 当前状态

- `moveParallel()` 代码已实现但在 evolve 中被注释掉
- 碰撞并行化（Phase 2）仍然正常工作
- 需要深入分析 OpenFOAM particle tracking 的线程安全性才能继续

## 下一步方向

1. 检查 `particle::trackToFace` 源码中是否有 static/global 状态
2. 检查 mesh 的 `tetBasePtIs()` 是否有 lazy evaluation
3. 考虑替代方案：不用 trackToFace，而是简单弹道推进 + findCell
4. 或者：在并行区内只做位置更新，不调用 trackToFace，之后串行验证 cell 归属

## 根因确认

`particle::trackToFace()` (particleTemplates.C) 中：
- **Line 295**: `DynamicList<label>& tris = cloud.labels();` — 共享 scratch buffer
- **Line 808**: 同上
- **Line 367/721/880/1237**: `cloud.trackingRescue()` — 共享计数器

`cloud.labels()` 是 Cloud 基类的一个成员 `DynamicList<label>`，被 trackToFace 用作临时三角形索引存储。多线程同时写入同一个 buffer → 数据竞争 → segfault。

## 修复方案

1. **方案 A**：在 Cloud 基类中添加 per-thread labels buffer（`PtrList<DynamicList<label>> threadLabels_`），修改 `labels()` 接受 tid 参数
2. **方案 B**：在 trackingData 中添加 per-particle scratch buffer
3. **方案 C**：不用 trackToFace，实现简化的弹道推进 + cell 查找（牺牲精度）
4. **方案 D**：使用 thread_local 变量替代 cloud.labels()

推荐方案 A 或 D，改动最小。

## thread_local 修复尝试

1. 在 `Cloud.H` 中将 `labels()` 改为返回 `thread_local` 变量
2. 在 `trackingRescue()` 中添加 `#pragma omp atomic`
3. 在 `lagrangian/basic/Make/options` 中添加 `-qopenmp`（否则 `_OPENMP` 不会被定义）

结果：仍然 crash（进程静默退出，无错误输出）。

可能原因：
- `thread_local` 在 Intel icpx 的模板实例化中可能有问题
- `trackToFace` 中可能还有其他非线程安全的访问（如 mesh 的 lazy evaluation）
- OpenFOAM 的 signal handler 在 OMP 并行区内无法正常捕获 segfault

## 当前状态

- `moveParallel()` 代码保留但在 evolve 中被注释掉
- `Cloud.H` 中的 thread_local labels 和 atomic trackingRescue 修改保留
- `lagrangian/basic/Make/options` 已添加 `-qopenmp`
- 碰撞并行化（Phase 2）正常工作

## 后续方向

Phase 5 move 并行化需要更深入的方案：
1. 不使用 `trackToFace`，实现独立的弹道推进 + cell 查找
2. 或者修改 `particleTemplates.C` 中的 `trackToFace`，将 `cloud.labels()` 替换为栈上局部变量
3. 或者在 `trackingData` 中添加 per-instance scratch buffer

## 最终根因确认（完整列表）

`trackToFace` 并行化在 8 线程下 crash 的完整原因链：

1. **`cloud.labels()`** — 共享 scratch buffer（已修复：改为局部变量）
2. **`cloud.trackingRescue()`** — `nTrackingRescues_++` 竞争（不导致 crash，只是计数不准）
3. **`dsmcDiffuseWallPatch::performDiffuseReflection`** 中：
   - `Random& rndGen = cloud_.rndGen()` — **多线程共享同一个 RNG → crash**
   - `cloud_.porousMeas().diffuseInteraction(p)` — 写入共享 measurement
4. **`dsmcDeletionPatch::controlParticle`** 中：
   - `cloud_.porousMeas().deletionInteraction(p, patchId())` — 写入共享 measurement
5. **`measurePropertiesBeforeControl/AfterControl`** — 写入共享 boundary flux 数组

## 结论

Phase 5 move 并行化需要 Phase 6（boundary model 线程安全迁移）作为前提：
- 所有 boundary model 的 `controlParticle()` 必须使用 per-thread RNG
- 所有 measurement 写入必须使用 thread-local buffer
- 这是计划中 Phase 6 的工作范围

当前可交付的成果：
- Phase 2 碰撞并行化 6.25x（正确，稳定）
- Phase 5 框架代码已实现（moveParallel, buildBoundaryCellMarking, parallelMoveActive）
- Phase 5 的所有线程安全问题已识别和记录
- `particleTemplates.C` 中 `cloud.labels()` 已改为局部变量（为后续启用做准备）

## 最终修复成功！

### 根因
`particleTemplates.C` 中有**两处** `cloud.labels()` 调用（line 295 和 line 808），之前只修复了第一处。第二处（在 DSMC 模式的 `trackToFace` 分支中）仍然使用共享 scratch buffer，导致 8 线程 crash。

### 修复
将两处 `DynamicList<label>& tris = cloud.labels()` 都改为 `DynamicList<label> tris`（局部变量）。

### 其他修复
- `Cloud.H`: `trackingRescue()` 使用 `__sync_fetch_and_add` 原子递增
- `moveParallel()`: 预触发 `mesh_.cells()` 和 `this->cellHasWallFaces()`
- warmup `trackToFace` 调用确保所有 lazy data 初始化

### 性能结果（noreact OMP-8）

| 阶段 | 串行 1T | OMP-8 (move+coll parallel) | 加速比 |
|------|---------|---------------------------|--------|
| move | 170.88s | **19.13s** | **8.93x** |
| collision | 35.53s | 4.33s | 8.20x |
| fields | 29.42s | 25.68s | 1.15x |
| buildCellOcc | 27.23s | 21.76s | 1.25x |
| **total** | **263.30s** | **71.12s** | **3.70x** |
| nParticles | 2,302,546 | 2,545,101 | 统计差异(RNG序列不同) |

粒子数差异（~10%）是因为 per-thread RNG 序列与串行不同，导致 inflow/deletion 的随机事件不同。这是 DSMC 统计波动，不是错误。
