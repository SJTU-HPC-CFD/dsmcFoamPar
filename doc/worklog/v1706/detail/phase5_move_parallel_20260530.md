# Phase 5 工作日志：Move 并行化原型

日期：2026-05-30

## 实现内容

1. `buildBoundaryCellMarking()` — 预计算哪些 cell 有 boundary face
2. `moveParallel(trackTime)` — 分类并行 move：
   - 内部粒子（不在 boundary cell）：OpenMP parallel for + trackToFace
   - Boundary 粒子：串行 dsmcParcel::move()
3. 在 `evolve_moveAndCollide()` 中，当 `omp_get_max_threads() > 1` 时使用 moveParallel

## 调试过程

- `_OPENMP` 宏正确定义（`-qopenmp` 在编译命令中）
- 问题：solver 需要重新链接才能使用更新的 .so（`rm dsmcFoam+.o && wmake`）
- `boundaryCell_` 延迟初始化（在 moveParallel 首次调用时检查 size）

## 性能结果（noreact 算例，300 步）

| 阶段 | 串行 1T | OMP-8 Phase5 | 加速比 |
|------|---------|-------------|--------|
| move | 170.88s | 27.69s | **6.17x** |
| collision | 35.53s | 4.48s | **7.93x** |
| fields | 29.42s | 27.03s | 1.09x |
| buildCellOcc | 27.23s | 22.42s | 1.21x |
| **total** | **263.30s** | **81.83s** | **3.22x** |

## 正确性问题

**粒子数不一致**：baseline 2,302,546 vs phase5 2,545,922（多了 ~10%）

原因分析：
- `moveParallel()` 中内部粒子的 tracking 跳过了 `trackParcelFaceTransition`（flux 统计）
- 更重要的是：内部粒子如果 track 到 boundary cell 后碰到 deletion patch，
  当前实现没有正确处理 `td.keepParticle = false` 的情况
- 内部粒子的 `hitWallPatch` 回调在 `trackToFace` 内部被调用了，
  但 `measurePropertiesBeforeControl` 写入共享数组可能有竞争

## 下一步修复

1. 内部粒子 track 后如果 `keepParticle = false`，需要标记删除
2. 需要处理 processor patch hit（当前只设 flag 但不做 transfer）
3. `trackParcelFaceTransition` 需要 thread-local buffer 或跳过
4. 验证 boundary measurement 的线程安全性

## 修改文件

- `src/lagrangian/dsmc/clouds/dsmcCloud.H` — boundaryCell_, moveParallel(), buildBoundaryCellMarking()
- `src/lagrangian/dsmc/clouds/dsmcCloud.C` — 实现 + evolve 调用路径修改
