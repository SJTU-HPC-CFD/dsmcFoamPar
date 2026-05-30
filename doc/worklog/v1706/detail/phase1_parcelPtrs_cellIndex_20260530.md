# Phase 1 工作日志：parcelPtrs_ 桥接索引 + Cell Index + Thread RNG

日期：2026-05-30

## 完成内容

1. 在 `dsmcCloud.H` 中添加：
   - `DynamicList<dsmcParcel*> parcelPtrs_` — 随机访问指针数组
   - `labelList cellFirst_`, `cellNext_`, `cellCount_` — SPARTA 风格 cell 链表索引
   - `PtrList<Random> threadRng_` — per-thread RNG
   - 公共方法：`buildParcelPtrs()`, `buildCellIndex()`, `validateCellIndex()`, `rng(tid)`
   - 访问器：`cellFirst()`, `cellNext()`, `cellCount()`, `parcelPtrs()`

2. 在 `dsmcCloud.C` 中实现：
   - `buildParcelPtrs()` — 从 IDLList 构建指针数组
   - `buildCellIndex()` — 反向遍历构建 cell 链表（O(n)）
   - `validateCellIndex()` — 对比 cellOccupancy 验证正确性
   - `rng(tid)` — 统一 RNG 入口
   - 构造函数中初始化 threadRng_（8 线程）
   - 在 move 后的 buildCellOccupancy 同步点重建 parcelPtrs + cellIndex
   - autoMap（load balance）后也重建

3. 编译修复：
   - `Random::integer()` 不存在 → 改用 `Random::sample01<scalar>()`
   - `Random::scalar01()` 是 private → 改用 `sample01<scalar>()`

## 测试结果

OMP_NUM_THREADS=8，300 步：
- 粒子数：2,320,550（baseline: 2,320,691，统计等价）
- 碰撞数：25,103（baseline: 25,214，统计等价）
- 线性动能：1.012e-18（一致）
- threadRng 初始化：8 thread RNGs

计时对比（300 步 total phase time）：
- Phase 0 baseline: 228.37s (buildCellOcc 5.2%)
- Phase 1: 249.42s (buildCellOcc 11.2%)
- 新增开销：buildParcelPtrs + buildCellIndex 约 +6% (~16s/300步)

## 备注

buildParcelPtrs + buildCellIndex 的额外开销（~53ms/步）在预期范围内。
后续 Phase 2 碰撞并行化后，这个开销会被碰撞加速抵消。
如果需要进一步优化，可以考虑只在碰撞前构建一次（而非每次 buildCellOccupancy 都重建）。

## 修改文件

- `src/lagrangian/dsmc/clouds/dsmcCloud.H`
- `src/lagrangian/dsmc/clouds/dsmcCloud.C`
