# Phase 5/6 最终成功：Move + Collision 并行化

日期：2026-05-30

## 最终修复

**根因**：`particleTemplates.C` 中有两处 `cloud.labels()` 调用（line 295 和 line 808），之前只修复了第一处。第二处在 DSMC 模式的 `trackToFace` 分支中，导致 8 线程 crash。

**修复清单**：
1. `particleTemplates.C`: 两处 `cloud.labels()` → 局部 `DynamicList<label> tris`
2. `Cloud.H`: `trackingRescue()` 使用 `__sync_fetch_and_add` 原子递增
3. `dsmcCloudI.H`: `rndGen()` 使用 `omp_in_parallel()` 自动路由 per-thread RNG
4. `dsmcCloud.C`: `moveParallel()` 预触发 `mesh_.cells()` + `cellHasWallFaces()`
5. `dsmcCloud.C`: equipartition 函数改用 `rndGen()` 替代 `rndGen_`
6. `dsmcPatchBoundary.C`: `measureProperties*` 在 `parallelMoveActive_` 时跳过
7. `dsmcDeletionPatch.C`: `porousMeas().deletionInteraction` 在 `parallelMoveActive_` 时跳过

## 最终性能结果

### React 算例（有化学反应，300 步）

| 阶段 | OMP-1 (baseline) | OMP-8 | 加速比 |
|------|-----------------|-------|--------|
| move | 182.18s | 19.21s | **9.48x** |
| collision | 38.01s | 4.67s | **8.14x** |
| fields | 32.88s | 25.47s | 1.29x |
| buildCellOcc | 29.27s | 21.78s | 1.34x |
| **total** | **282.61s** | **71.34s** | **3.96x** |

### Noreact 算例（无化学反应，300 步）

| 阶段 | OMP-1 (baseline) | OMP-8 | 加速比 |
|------|-----------------|-------|--------|
| move | 170.88s | 19.13s | **8.93x** |
| collision | 35.53s | 4.33s | **8.20x** |
| fields | 29.42s | 25.68s | 1.15x |
| buildCellOcc | 27.23s | 21.76s | 1.25x |
| **total** | **263.30s** | **71.12s** | **3.70x** |

### 对比 MPI

| 方法 | 总加速比 |
|------|---------|
| MPI-8 | 2.74x |
| **OMP-8** | **3.70-3.96x** |

**OpenMP 超越了 MPI！** 因为没有通信开销和域分解不平衡。
