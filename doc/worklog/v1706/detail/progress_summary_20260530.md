# 工作进展总结 2026-05-30

## 已完成

### Phase 0：基线计时
- OpenMP 编译支持（`-qopenmp`）
- Wall-time 阶段计时（move/collision/fields/buildCellOcc/coordSystem/controllers）
- 并行诊断输出（nParticles, nCells）
- Baseline profiling: **move=69%, collision=12%, fields=13%, buildCellOcc=5%**
- MPI 基准：MPI-2=1.55x, MPI-4=2.05x, MPI-8=2.74x

### Phase 1：parcelPtrs_ 桥接索引
- `DynamicList<dsmcParcel*> parcelPtrs_` — 随机访问指针数组
- `cellFirst_/cellNext_/cellCount_` — SPARTA 风格 cell 链表索引
- `PtrList<Random> threadRng_` — per-thread RNG
- `rng(tid)` — 统一 RNG 入口
- 验证通过

### Phase 2：noTimeCounterOMP 碰撞并行
- Runtime selection 注册成功
- Reaction fallback 正确工作（测试算例有 reactions）
- **无 reaction 算例碰撞加速比：4.85x（8 线程）**
- Zero-copy backend（直接在 dsmcParcel* 上操作）

### Phase 5：Move 分析
- `trackToFace()` 内部调用 `hitWallPatch()` 回调
- `measurePropertiesBeforeControl()` 写入共享 boundary flux 数组 — 非线程安全
- 方案确定：分类并行（内部粒子并行 + boundary 粒子串行）
- 需要预计算 boundary cell 标记

## 当前瓶颈

Move 占 69%，是获得有意义总加速的唯一路径。
即使 collision + fields 完美并行，总加速上限仅 1.44x。

## 下一步

1. Phase 5 实现：预计算 boundary cell → 内部粒子并行 move → boundary 粒子串行
2. Phase 3 fields 并行化（低风险，13% 收益）
3. Phase 7 reaction delayed event（使 react 算例也能用 collision OMP）

## 文件变更清单

- `src/lagrangian/dsmc/Make/options` — 添加 `-qopenmp`
- `src/lagrangian/dsmc/Make/files` — 添加 noTimeCounterOMP
- `src/lagrangian/dsmc/clouds/dsmcCloud.H` — omp.h, 计时成员, parcelPtrs_, cellIndex, threadRng_, rng()
- `src/lagrangian/dsmc/clouds/dsmcCloud.C` — 计时实现, buildParcelPtrs, buildCellIndex, validateCellIndex, rng, threadRng 初始化
- 新建 `collisionPartnerSelection/derived/noTimeCounterOMP/noTimeCounterOMP.H`
- 新建 `collisionPartnerSelection/derived/noTimeCounterOMP/noTimeCounterOMP.C`
- `run/.../react/omp8/constant/dsmcProperties` — 切换到 noTimeCounterOMP
- `run/.../react/omp8/system/controlDict` — 添加 timingOutputInterval
- 新建 `run/.../react/mpi2/`, `mpi4/` — MPI 基准算例
- 新建 `run/.../react/omp8-noreact/` — 无 reaction 验证算例
