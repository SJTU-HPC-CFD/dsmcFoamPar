# Phase 0 工作日志：基线计时与 Profiling

日期：2026-05-30

## 完成内容

1. 在 `Make/options` 中添加 `-qopenmp` 编译选项（Intel icpx 编译器）
2. 在 `dsmcCloud.H` 中添加 `#include <omp.h>`（`_OPENMP` 守护）和计时成员变量
3. 在 `dsmcCloud.C` 中实现 `evolve_moveAndCollide()` 和 `evolve_fields()` 的 wall-time 阶段计时
4. 在 `evolve()` 中添加周期性计时输出（`timingOutputInterval` 控制，从 `controlDict` 读取）
5. 编译通过，物理结果正确

## Profiling 结果

### 串行 Baseline（OMP_NUM_THREADS=8，无并行 kernel）

算例：`pal-phd3.3.1-2dcylinder/react/omp8`，300 步

| 阶段 | 时间(s) | 占比 |
|------|---------|------|
| move | 158.54 | 69.4% |
| collisions | 27.88 | 12.2% |
| fields | 29.85 | 13.1% |
| buildCellOccupancy | 11.87 | 5.2% |
| other | 0.23 | 0.1% |
| **total** | **228.37** | 100% |

per-step average: 0.761 s, nParticles: 2,320,691, nCells: 104,151

### MPI 基准

| 指标 | MPI-2 | MPI-4 | MPI-8 |
|------|-------|-------|-------|
| ExecutionTime (s) | 157 | 117 | 89 |
| per-step (s) | 0.492 | 0.371 | 0.278 |
| move 占比 | 69.7% | 73.5% | 68.8% |
| collision 占比 | 11.7% | 11.3% | 9.6% |
| fields 占比 | 14.0% | 12.1% | 16.1% |
| buildCellOcc 占比 | 4.5% | 3.0% | 5.4% |
| nParticles (total) | 2,320,779 | 2,320,527 | 2,320,498 |

MPI 加速比（相对串行 228s）：MPI-2=1.55x, MPI-4=2.05x, MPI-8=2.74x

## 关键结论

1. **Move 占 69-74%**，是绝对瓶颈，远超原计划假设的 50%
2. Collision 仅占 10-12%，fields 占 12-16%
3. MVP 路线（只并行 collision + fields）的 Amdahl 上限：`1/0.694 = 1.44x`
4. **Phase 5 move 并行化必须尽早启动**，否则无法获得有意义的总加速
5. MPI-8 加速比 2.74x 说明 move 阶段本身是可并行的（MPI 域分解有效）
6. 测试算例有 reactions（dissociationQK），Phase 2 collision OMP 会整体回退到串行

## 对计划的影响

- 原计划执行顺序 `M0→M1→M2→M3→M4→M5` 需要调整
- 推荐：`M0→M1→M2(最小)→M5(立即启动)` 与 `M3/M4` 并行
- Phase 2 collision 对本算例收益有限（~12% × 加速比），但仍有框架价值
- Phase 5 move 是获得 >2x 加速的唯一路径

## 修改文件

- `src/lagrangian/dsmc/Make/options` — 添加 `-qopenmp`
- `src/lagrangian/dsmc/clouds/dsmcCloud.H` — 添加 omp.h、计时成员
- `src/lagrangian/dsmc/clouds/dsmcCloud.C` — 实现阶段计时和输出
- `run/.../react/omp8/system/controlDict` — 添加 `timingOutputInterval 300`
- `run/.../react/mpi2,mpi4,mpi8/` — 创建 MPI 基准算例

## 日志文件位置

- `run/.../react/omp8/log.dsmcFoam+` — 串行 baseline
- `run/.../react/mpi2/log.dsmcFoam+` — MPI-2 基准
- `run/.../react/mpi4/log.dsmcFoam+` — MPI-4 基准
- `run/.../react/mpi8/log.dsmcFoam+` — MPI-8 基准
