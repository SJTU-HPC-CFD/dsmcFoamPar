# Phase 5 分析：Move 并行化设计

日期：2026-05-30

## 代码分析

### Cloud::move() 结构
- while 循环处理 processor transfer（多次迭代直到无 transfer）
- 每次迭代遍历所有粒子调用 `dsmcParcel::move(td, trackTime)`
- 粒子 move 后检查 `switchProcessor`，准备 MPI transfer

### dsmcParcel::move() 结构
- 设置 `stepFraction`（新粒子随机化）
- while 循环：`trackToFace()` → 检查 boundary → 更新 stepFraction
- `trackToFace(position + dt*U, td, true)` 第三个参数 `true` 表示 DSMC 模式
  - 内部调用 `hitWallPatch()` → boundary model `controlParticle()`
  - 内部调用 `hitProcessorPatch()` → 设置 `switchProcessor = true`
- `tracker().trackParcelFaceTransition()` — 写入共享 flux 统计

### 线程安全分析

| 操作 | 线程安全？ | 原因 |
|------|-----------|------|
| trackToFace 几何追踪 | ✓ | 只读 mesh |
| 修改粒子 position/cell/tet | ✓ | per-particle |
| hitWallPatch → controlParticle | △ | 取决于 boundary model |
| hitProcessorPatch | ✓ | 只设 flag |
| trackParcelFaceTransition | ✗ | 写共享 flux 数组 |
| rndGen_.sample01 (新粒子) | ✗ | 共享 RNG |
| cyclic boundary model | ✗ | 可能有共享状态 |

### 关键发现

1. `trackToFace()` 不能被拆分为"先追踪再处理 boundary"——boundary hit 发生在追踪过程中
2. 但 `hitWallPatch` 对于 specular/diffuse wall 只修改当前粒子的速度——实际上是线程安全的
3. `trackParcelFaceTransition` 是唯一确定的共享写入
4. processor transfer 需要 MPI 同步，必须串行

### 可行方案

**方案 A：直接并行 dsmcParcel::move()（最大收益，中等风险）**
- 在 `parcelPtrs_` 上 `#pragma omp parallel for`
- 用 thread-local buffer 替代 `trackParcelFaceTransition` 的直接写入
- 用 `cloud_.rng(tid)` 替代 `rndGen_` 用于新粒子 stepFraction
- Processor patch 粒子标记后串行处理 MPI transfer
- 需要验证所有 boundary model 的 `controlParticle()` 是否线程安全

**方案 B：预筛选 + 弹道推进（保守，低风险）**
- 预计算哪些粒子在当前步不会碰到任何 face（`position + dt*U` 仍在当前 cell 内）
- 这些粒子直接更新 position，不调用 trackToFace
- 剩余粒子走串行 Cloud::move()
- 预期覆盖率：取决于 cell 大小 vs 粒子速度*dt

**方案 C：绕过 Cloud::move()，自建并行 move 框架（最大改动）**
- 不使用 Cloud::move() 的 IDLList 遍历
- 在 parcelPtrs_ 上并行调用 tracking
- 自行处理 processor transfer（收集 → MPI → 分发）
- 风险最高但收益最大

### 推荐

先实现方案 A 的原型：
1. 在 `evolve_moveAndCollide()` 中，用并行版本替代 `Cloud<dsmcParcel>::move(td, dt)`
2. 需要修改 `dsmcParcel::move()` 接受 thread-local tracker
3. 需要验证 `dsmcDiffuseWallPatch::controlParticle()` 线程安全性

## 下一步

- 先完成 Phase 3 fields 并行化（低风险确定收益）
- 然后回来实现 Phase 5 方案 A
