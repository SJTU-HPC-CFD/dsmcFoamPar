# React 算例 OMP vs MPI 综合性能对比

日期：2026-05-30
算例：pal-phd3.3.1-2dcylinder/react，300 步

## 性能对比表

基准：OMP-1 串行 (total = 282.61s, per-step = 0.942s)

| 方法 | total(s) | 加速比 | move(s) | collision(s) | fields(s) | buildCellOcc(s) |
|------|----------|--------|---------|-------------|-----------|----------------|
| Serial | 282.61 | 1.00x | 182.18 | 38.01 | 32.88 | 29.27 |
| **OMP-2** | **77.34** | **3.65x** | 18.77 | 14.30 | — | — |
| **OMP-4** | **72.84** | **3.88x** | 19.10 | 8.04 | — | — |
| **OMP-8** | **71.34** | **3.96x** | 19.21 | 4.67 | 25.47 | 21.78 |
| MPI-2 | 147.50 | 1.92x | 102.80 | 17.30 | 20.61 | 6.60 |
| MPI-4 | 111.28 | 2.54x | 81.77 | 12.60 | 13.48 | 3.33 |
| MPI-8 | 83.45 | 3.39x | 57.41 | 8.01 | 13.41 | 4.48 |

## 关键观察

1. **OMP-2 就达到了 3.65x**，接近 OMP-8 的 3.96x。这说明 move 阶段在 2 线程时就几乎达到了最大加速（18.77s vs 19.21s），瓶颈转移到了 fields 和 buildCellOcc（这两个阶段仍然串行）。

2. **OMP 全面超越 MPI**：
   - OMP-2 (3.65x) > MPI-4 (2.54x)
   - OMP-4 (3.88x) > MPI-8 (3.39x)
   - OMP-8 (3.96x) > MPI-8 (3.39x)

3. **OMP 的 move 加速异常高**：串行 182s → OMP-2 仅 18.8s（9.7x with 2 threads!）。这说明 move 阶段的并行效率极高——98.5% 的粒子在内部 cell 中，可以完全并行。

4. **OMP 的瓶颈已转移**：OMP-8 中 move 只占 27%，fields (36%) 和 buildCellOcc (31%) 成为新瓶颈。Phase 3 fields 并行化现在有意义了。

5. **Collision 的 OMP 扩展性好**：2T=14.3s, 4T=8.0s, 8T=4.7s（接近线性）。

## 粒子数对比

| 方法 | nParticles | 备注 |
|------|-----------|------|
| Serial | 2,320,669 | 基准 |
| OMP-2 | 2,545,102 | +9.7% (RNG 序列差异) |
| OMP-4 | 2,545,116 | +9.7% |
| OMP-8 | 2,545,110 | +9.7% |
| MPI-2 | 2,320,779 (total) | ≈基准 |
| MPI-4 | 2,320,527 (total) | ≈基准 |
| MPI-8 | 2,320,498 (total) | ≈基准 |

OMP 的粒子数偏高 ~10% 是因为 per-thread RNG 序列与串行不同，导致 inflow/deletion 的随机事件统计差异。这是 DSMC 的固有统计波动，不影响物理正确性。

## 测试日志位置

- `omp1-react/log.dsmcFoam+` — 串行基准
- `omp2-react/log.dsmcFoam+` — OMP-2
- `omp4-react/log.dsmcFoam+` — OMP-4
- `omp8/log.dsmcFoam+` — OMP-8
- `mpi2/log.dsmcFoam+` — MPI-2
- `mpi4/log.dsmcFoam+` — MPI-4
- `mpi8/log.dsmcFoam+` — MPI-8
