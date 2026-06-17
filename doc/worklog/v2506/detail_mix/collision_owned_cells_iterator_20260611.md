# collision owned-cells iterator 实现与测试 — OFv1706 hyStrath_dlb

日期：2026-06-11

## 1. 问题

MPI replicated mesh 模式下，collision 循环遍历全部 nCells（OMP 路径）或全部
`occupancyCollisionCells_`（串行路径），导致每 rank 迭代量 = nCells，总迭代量
= nRanks × nCells。OMP8 单进程对所有 nCells 做 8 线程平分，每线程 ~nCells/8；
MPI 模式下每线程迭代量从 nCells/8 膨胀到 nCells/2 ~ nCells，collision 从 3.52s
膨胀到 16.87~32.80s。

根因：collision 循环没有用 `myCells_`（replicated mesh 提供的 owned-cell 列表）
裁剪迭代空间。

## 2. 修改

### 2.1 `dsmcCloud.H` — 新增成员和访问器

- 新增 `labelList occupancyOwnedCollisionCells_`（owned + 有 >1 粒子的 cell 列表）
- 新增 `occupancyOwnedCollisionCells()` inline 访问器

### 2.2 `dsmcCloud.C` — `buildCellOccupancy()` 中构建列表

在 OMP 和串行两条路径的 `buildCellOccupancy()` 末尾，分别插入：

```cpp
if (replicatedMeshActive())
{
    DynamicList<label> ownedCC(occupancyCollisionCells_.size());
    forAll(occupancyCollisionCells_, i)
    {
        const label cellI = occupancyCollisionCells_[i];
        if (replicatedMesh_->isMyCell(cellI))
            ownedCC.append(cellI);
    }
    occupancyOwnedCollisionCells_.transfer(ownedCC);
}
else
{
    occupancyOwnedCollisionCells_ = occupancyCollisionCells_;
}
```

构造函数初始化列表同步追加 `occupancyOwnedCollisionCells_()`。

### 2.3 `noTimeCounter.C` — 碰撞循环改用 owned 列表

**OMP 路径**（static/guided/dynamic 三种 schedule）：

```cpp
// 旧：for (label cellI = 0; cellI < nCells; ++cellI)
// 新：
const labelList& ownedCC = cloud_.occupancyOwnedCollisionCells();
const label nOwnedCC = ownedCC.size();
for (label idx = 0; idx < nOwnedCC; ++idx)
{
    const label cellI = ownedCC[idx];
    localCollisions += processCell(cellI, threadI);
    localCandidates += nCandidatesPerCell[cellI];
}
```

**串行路径**：

- `nCandidatesPerCell` 清零从 `forAll(nCandidatesPerCell, cellI)` 改为
  只清零 `ownedCollCells` 中的 cell
- 碰撞循环从 `occupancyCollisionCells()` 改为 `occupancyOwnedCollisionCells()`

## 3. 构建

```bash
source doc/scripts/env.sh
cd src/lagrangian/dsmc && wmake -j libso     # libdsmcFoam+.so
cd applications/.../dsmcFoam+ && wmake -j     # dsmcFoam+ binary
```

构建通过（2026-06-11 16:49）。

## 4. 测试

### 4.1 Smoke test — OMP8 10 步

- exit 0, End main, 0 stuck particles
- 粒子数 2,065,844，总能量 1.1506 — 与历史 OMP8 一致
- collision phase = 0.019s / 10 步

### 4.2 快速性能对比 — ourmesh 500 步 no-write

| 模式 | real (s) | move | buildOcc | collision | post | vs 旧基准 |
|---|---|---|---|---|---|---|
| OMP8 | 62.31 | 47.75 | 6.58 | **3.49** | 0.69 | -0.9% |
| MPI2xOMP4 | 69.89 | 52.17 | 7.58 | **16.04** | 1.06 | **-4.9%** |
| MPI8 | 97.59 | 68.05 | 10.22 | **24.82** | 4.87 | **-24.3%** |

旧基准取自 `dsmcFoam_plus_mpi_omp_mixed_technical_report_20260610.md` 的
no-write 正式对比（ourmesh 3 次均值）。

### 4.3 Correctness gates

| 模式 | final particles | collisions | stuck | total energy |
|---|---|---|---|---|
| OMP8 | 2,463,761 | 35,601 | 0 | 1.243103838 |
| MPI2xOMP4 | 2,463,507 | 35,531 | 0 | 1.243108864 |
| MPI8 | 2,463,680 | 35,262 | 0 | 1.242212945 |

全部通过，粒子数和能量在 DSMC 随机波动范围内一致。

## 5. 分析

### OMP8 — 持平（-0.9%，噪声级）

OMP8 不使用 replicated mesh，`occupancyOwnedCollisionCells_` = `occupancyCollisionCells_`，
循环量不变。3.52s → 3.49s 在正常波动范围。

### MPI2xOMP4 — 温和改善（-4.9%）

OMP 碰撞路径从 `for (cellI=0; cellI<nCells)` 改为只遍历 owned collision cells。
每 rank 循环量从 104,151 → ~52,000。但 `processCell` 对空 cell 返回极快
（仅做 occupancy count 检查和 early return），所以实际 CPU 节省有限。

### MPI8 — 显著改善（-24.3%）

串行碰撞路径双重受益：
1. 循环量从 `occupancyCollisionCells_`（~13,000 个 owned cells with >1 parcel）
   改为 `occupancyOwnedCollisionCells_`（显式过滤，更安全）
2. `nCandidatesPerCell` 清零从 `forAll(104151)` 改为只清零 owned collision
   cells（~13,000），500 步累计节省可观

## 6. 文件清单

| 文件 | 修改内容 |
|---|---|
| `src/lagrangian/dsmc/clouds/dsmcCloud.H` | +`occupancyOwnedCollisionCells_` 成员 + 访问器 |
| `src/lagrangian/dsmc/clouds/dsmcCloud.C` | `buildCellOccupancy()` 两处构建 owned 列表；构造函数初始化列表 |
| `src/lagrangian/dsmc/collisionPartnerSelection/derived/noTimeCounter/noTimeCounter.C` | OMP 路径 3 处 schedule 循环 + 串行路径改用 owned 列表 |

## 7. 后续建议

1. **补做 MPI4xOMP2 测试** — 预期 collision 改善幅度介于 MPI2xOMP4 和 MPI8 之间
2. **zb-cylinder-react 验证** — 反应流 case 的 collision/reaction 路径也需要验证
3. **多次重复** — 单次 run 有波动，建议每模式 3 次重复取均值
4. **其他 collision 模型** — `noTimeCounterSubCycled`、`simplifiedBernoulliTrials`
   等模型仍有 `forAll(cellOccupancy, cellI)` 循环，可按相同方式修改
