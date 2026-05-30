# Move 优化：延迟删除 + parcelPtrs_ compact

日期：2026-05-30

## 优化内容

在 `moveParallel` 中：
1. 不立即调用 `deleteParticle()`，而是标记 dead 粒子（`stepFraction = -1`）
2. Boundary 粒子串行 move 后也用同样的 dead 标记
3. 最后一次性遍历 `parcelPtrs_`：删除 dead 粒子 + compact 数组
4. `evolve_moveAndCollide` 中 move 后跳过 `buildParcelPtrs()`（已经 compact 了），只重建 `buildCellIndex()`

## 效果

省掉了 move 后的 `buildParcelPtrs()` IDLList 遍历（~7s/300步）。

| 阶段 | 之前 | 优化后 | 改善 |
|------|------|--------|------|
| buildCellOcc | 12.3s | **4.4s** | **+64%** |
| **total** | **88.2s** | **85.8s** | **+3%** |
| **加速比** | **3.20x** | **3.29x** | |

nParticles: 2,320,577 ✓ 正确
