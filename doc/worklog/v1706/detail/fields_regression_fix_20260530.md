# Fields 性能回退修复

日期：2026-05-30

## 问题
OMP-2 的 fields 时间 67-104s，远超串行 33s。

## 根因
`dsmcVolFields::calculateField()` 的 densityOnly 和 else 分支被改为 `cellFirst/cellNext` 遍历 + `#pragma omp parallel for`。
`cellFirst/cellNext` 遍历粒子时，粒子在 `parcelPtrs_` 中的物理内存位置与 cell 归属无关，导致严重的 **cache miss**（每次访问粒子数据都是随机内存跳跃）。
原始的 `forAllConstIter(dsmcCloud, cloud_, iter)` 虽然也是链表遍历，但 OpenFOAM 的 parcel 分配器有更好的内存局部性。

## 修复
恢复 `dsmcVolFields.C` 的粒子遍历为原始的 `forAllConstIter`（串行但 cache-friendly）。
移除所有 `#pragma omp parallel for` from dsmcVolFields.C。

## 结论
Fields 的粒子累加循环**不适合用 cellFirst/cellNext 并行化**——cache miss 的代价远超并行收益。
正确的并行化方案需要先按 cell 排序粒子（`sortedParticleIds_`），但这需要额外的排序开销。
对于当前算例（fields 占 ~30%），保持串行是最优选择。

## 修复后结果

| 方法 | total(s) | 加速比 | fields(s) |
|------|----------|--------|-----------|
| Serial | 282.6 | 1.00x | 32.9 |
| OMP-2 | 149.5 | 1.89x | 29.2 ✓ |
| OMP-4 | 108.5 | 2.60x | 31.3 ✓ |
| OMP-8 | 95.8 | 2.95x | 34.6 ✓ |
