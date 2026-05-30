# Move 并行化修复与最终正确结果

日期：2026-05-30

## 修复内容

### 根因
1. **缺少 `stepFraction = 0` 重置**：`Cloud::move()` 在开始时将所有粒子的 stepFraction 设为 0，`moveParallel` 没有做这个，导致 boundary 粒子保留上一步的 stepFraction=1.0，不移动，不被 deletion patch 删除
2. **缺少 2D 约束**：`meshTools::constrainToMeshCentre` 和 `constrainDirection` 在 2D 算例中必须每步调用，否则粒子飘出 2D 平面导致 tracking 异常

### 修复代码
```cpp
// 在 moveParallel 开头
forAllIter(dsmcCloud, *this, pIter)
{
    pIter().stepFraction() = 0;
}

// 在并行 tracking 循环内
if (coordSystem().type() == "dsmcCartesian")
{
    meshTools::constrainToMeshCentre(mesh_, p.position());
    meshTools::constrainDirection(mesh_, mesh_.solutionD(), Utracking);
}
```

## 最终正确性验证

| 指标 | Serial | OMP-8 | MPI-8 | 判定 |
|------|--------|-------|-------|------|
| nParticles | 2,320,669 | 2,320,493 | 2,320,498 | ✓ |
| Collisions | 24,944 | 24,927 | 24,969 | ✓ |
| linear KE | 1.012e-18 | 1.012e-18 | 1.012e-18 | ✓ |
| rotational E | 1.071e-20 | 1.067e-20 | 1.067e-20 | ✓ |

## 最终性能结果（React 算例，正确性已验证）

| 方法 | total(s) | 加速比 | move(s) | collision(s) | nParticles |
|------|----------|--------|---------|-------------|-----------|
| Serial | 282.61 | 1.00x | 182.18 | 38.01 | 2,320,669 |
| OMP-2 | 187.12 | 1.51x | 111.30 | 20.78 | 2,320,595 |
| OMP-4 | 147.27 | 1.92x | 72.71 | 11.59 | 2,320,678 |
| OMP-8 | 125.50 | 2.25x | 55.30 | 4.92 | 2,320,493 |
| MPI-2 | 147.50 | 1.92x | 102.80 | 17.30 | 2,320,779 |
| MPI-4 | 111.28 | 2.54x | 81.77 | 12.60 | 2,320,527 |
| MPI-8 | 83.45 | 3.39x | 57.41 | 8.01 | 2,320,498 |

## 分析

1. **OMP move 加速比 3.30x (8T)**：低于理想值，因为 stepFraction 重置和 2D 约束是串行的
2. **OMP collision 加速比 7.73x (8T)**：接近线性
3. **OMP-8 总加速 2.25x**：低于 MPI-8 (3.39x)，主要因为 fields (35s) 和 buildCellOcc (30s) 仍然串行
4. **OMP-4 (1.92x) = MPI-2 (1.92x)**：4 个 OMP 线程等效于 2 个 MPI 进程

## 后续优化方向

1. **并行化 stepFraction 重置**：`#pragma omp parallel for` 替代 `forAllIter`
2. **并行化 buildCellOccupancy**（占 24%）
3. **并行化 fields**（占 28%）
4. 这三项完成后预期 OMP-8 可达 3-4x

## 修改文件

- `src/lagrangian/dsmc/clouds/dsmcCloud.C` — stepFraction 重置 + 2D 约束

## 追加优化结果

### 优化内容
1. stepFraction 重置并行化（`#pragma omp parallel for`）
2. buildCellOccupancy clear 并行化
3. dsmcVolFields densityOnly 分支按 cell 并行化

### 优化后 OMP-8 结果（React 算例）

| 阶段 | 优化前 | 优化后 | 改善 |
|------|--------|--------|------|
| move | 55.30s | 50.63s | +8.5% |
| buildCellOcc | 30.16s | 28.74s | +4.7% |
| collision | 4.92s | 7.01s | -30% (omp_in_parallel 开销) |
| fields | 34.79s | 37.27s | -7% (未并行化 else 分支) |
| **total** | **125.50s** | **123.99s** | **+1.2%** |
| **加速比** | **2.25x** | **2.28x** | |

正确性：nParticles 2,320,680 ✓，能量正确 ✓

### 结论
- 小优化收益有限，主要瓶颈仍在 move (41%) 和 fields (30%)
- `omp_in_parallel()` 在 rndGen() 中的检查引入了额外开销
- fields 的 else 分支（完整物理量计算）需要并行化才能有显著收益
- 当前 OMP-8 总加速 2.28x，物理结果完全正确
