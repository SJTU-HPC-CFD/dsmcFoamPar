# DSMC 完整 OpenMP 并行化与负载平衡：实施计划（修订版）

日期：2026-05-30
基于：dsmc_sparta_plan_20260528.md（碰撞并行化计划）
参考：refcode/sparta-24Sep2025/（SPARTA Kokkos 并行化模式）

## 0. 对原计划的代码验证与修正

在详细审查 `work/v1706-sparta` 实际代码后，发现原计划的以下假设需要修正：

### 0.1 代码结构事实

| 原计划假设 | 实际代码 | 影响 |
|-----------|---------|------|
| cellOccupancy 是某种 map | `DynamicList<DynamicList<dsmcParcel*>>` 按 cellI 索引 | 替换模式不变，但类型签名需匹配 |
| BinaryCollisionModel::collide 接受 `dsmcParcel&` | 签名为 `collide(dsmcParcel& pP, dsmcParcel& pQ, const label cellI, scalar cR=-1)` | 重载需匹配此签名 |
| 反应模型独立调用 | 反应在 noTimeCounter::collide() **内部**调用：`cloud_.reactions().reactions()[rMId]->reaction(parcelP, parcelQ)` | 反应必须在碰撞循环内处理，不能分离 |
| move 阶段可简化为弹道推进 | 使用 `trackToFace()` + tet 分解追踪，涉及 patch 回调（壁面模型、cyclic、processor） | move 并行化比预想复杂得多 |
| controllers 只做统计 | temperatureController 等**直接修改粒子速度** `p->U()` | controller 并行化需要处理写冲突 |
| field 计算是简单 reduction | 通过 `cellOccupancy` 迭代，每个 field 独立调用 `calculateField()` | 可并行化但需适配新索引 |
| coordSystem().evolve() 无副作用 | axisymmetric 模式调用 `addNewParcel()`/`deleteParticle()` 修改粒子拓扑 | 必须在 move 并行化之外串行处理 |

### 0.2 关键代码路径（实际）

```
dsmcFoam+.C main loop:
  while (runTime.loop())
      dsmc.evolve()
          evolve_moveAndCollide()
              controllers_.controlBeforeMove()      ← 直接修改粒子 U_
              boundaries_.controlBeforeMove()
              [removeElectrons + buildCellOccupancy]
              Cloud<dsmcParcel>::move(td, dt)        ← trackToFace, patch 回调
              buildCellOccupancy()
              [addElectrons + buildCellOccupancy]
              coordSystem().evolve()                 ← axisym: clone/delete particles
              controllers_.controlBeforeCollisions()
              boundaries_.controlBeforeCollisions()
              collisions()
                  collisionPartnerSelectionModel_->collide()  ← noTimeCounter
                      for each cell in cellOccupancy:
                          build 8 sub-cells
                          NTC pair selection
                          reactions().reactions()[rMId]->reaction(pP, pQ)
                          binaryCollision().collide(pP, pQ, cellI)
              [buildCellOccupancy if reactions]
              controllers_.controlAfterCollisions()
              boundaries_.controlAfterCollisions()
          evolve_fields()
              fields_.calculateFields()             ← per-field calculateField()
              controllers_.calculateProps()
              boundaries_.calculateProps()
      dsmc.loadBalanceCheck()                       ← MPI 级重分区
```

### 0.3 dsmcParcel 实际字段

```cpp
class dsmcParcel : public particle  // particle 提供 position_, cellI_, tetFaceI_, tetPtI_
{
    vector U_;                    // 速度
    scalar RWF_;                  // 径向权重因子
    scalar ERot_;                 // 转动能
    label ELevel_;                // 电子能级
    label typeId_;                // 组分 ID
    label newParcel_;             // 新粒子标记（patch ID 或 -1）
    TrackedParcel tracked_;       // 追踪信息
    label classification_;        // 激波分类
    StuckParcel* stuck_;          // 粘附信息（指针，可为 nullptr）
    labelList vibLevel_;          // 振动量子数（动态大小！）
};
```

**关键发现**：`vibLevel_` 是 `labelList`（动态大小），不是固定数组。原计划的 `MAXVIBMODE=4` 假设需要运行时验证。

## 1. 扩展目标

### 1.1 目标加速比

| 阶段 | 典型占比 | 目标加速比(8T) | 难度 | 可行性评估 |
|------|---------|---------------|------|-----------|
| Collide | 40-70% | 5-7x | 中 | 高：per-cell 独立 |
| Move (Cloud::move) | 20-40% | 2-4x | 极高 | 中：trackToFace 涉及 patch 回调 |
| buildCellOccupancy | 5-10% | 3-5x | 低 | 高：纯索引构建 |
| Field calculation | 5-15% | 5-7x | 低 | 高：per-cell reduction |
| Controllers | 2-5% | 2-3x | 中-高 | 中：直接修改粒子 |

整体目标：8 线程 **4-5x 总加速比**（保守估计，考虑 move 阶段的限制）。

### 1.2 SPARTA 参考的适用性

| SPARTA 模式 | 适用于本项目？ | 原因 |
|------------|--------------|------|
| Kokkos parallel_for over cells (collide) | ✓ 直接适用 | per-cell 碰撞独立 |
| Kokkos parallel_reduce over particles (move) | △ 部分适用 | OF 追踪比笛卡尔复杂，但粒子间仍独立 |
| parallel_for + scan (sort) | ✓ 可简化适用 | 串行 O(n) 已足够快，大规模时可并行 |
| RCB load balance | ✗ 不直接适用 | 这是 MPI 级；我们需要线程级 |
| Kokkos ScatterView (field reduction) | ✓ 等价于 thread-local buffer | 避免 atomic |

## 2. 架构设计

### 2.1 核心数据结构

```
┌──────────────────────────────────────────────────────────────┐
│  dsmcCloud (扩展)                                            │
│                                                              │
│  ┌────────────────────────────────────────────────────────┐ │
│  │  Cloud<dsmcParcel> IDLList (保留)                       │ │
│  │  用于：OF 追踪 / IO / MPI 迁移 / patch 回调            │ │
│  └────────────────────────────────────────────────────────┘ │
│                                                              │
│  ┌────────────────────────────────────────────────────────┐ │
│  │  DynamicList<dsmcParcel*> parcelPtrs_                   │ │
│  │  [0][1][2]...[nLocal-1]  — 指针数组，O(1) 随机访问     │ │
│  └────────────────────────────────────────────────────────┘ │
│                                                              │
│  ┌──────────────┐  ┌──────────────┐  ┌──────────────┐      │
│  │ cellFirst_[] │  │ cellNext_[]  │  │ cellCount_[] │      │
│  │ nCells       │  │ nLocal       │  │ nCells       │      │
│  └──────────────┘  └──────────────┘  └──────────────┘      │
│                                                              │
│  ┌────────────────────────────────────────────────────────┐ │
│  │  PtrList<Random> threadRng_ (per-thread RNG)            │ │
│  └────────────────────────────────────────────────────────┘ │
│                                                              │
│  ┌────────────────────────────────────────────────────────┐ │
│  │  ThreadLoadBalancer threadBalance_ (线程负载平衡器)      │ │
│  └────────────────────────────────────────────────────────┘ │
└──────────────────────────────────────────────────────────────┘
```

**设计决策**：不引入独立的 `DsmcParticleData` 结构体用于 move 阶段。

原因：
1. Move 阶段使用 `trackToFace()` 需要完整的 `particle` 基类状态（position_, cellI_, tetFaceI_, tetPtI_, facei_, stepFraction_）
2. 拷贝这些字段到平坦结构体再拷回的开销可能抵消并行收益
3. 直接在 `dsmcParcel*` 指针数组上并行操作更高效

**碰撞阶段**使用 `DsmcCollisionData` 结构体（碰撞不需要追踪状态，平坦数组对 cache 更友好）。

### 2.2 碰撞粒子结构体（仅碰撞阶段使用）

```cpp
// DsmcCollisionData.H
enum { MAXVIBMODE = 5 };  // 覆盖 N2(1), O2(1), CO2(3), CH4(4), SF6(5)

struct DsmcCollisionData
{
    vector U;               // 24B  速度（直接用 OF vector 类型）
    scalar ERot;            //  8B  转动能
    scalar RWF;             //  8B  径向权重因子
    label  ELevel;          //  8B  电子能级
    label  typeId;          //  8B  组分索引
    label  cellI;           //  8B  所在单元
    label  origIdx;         //  8B  在 parcelPtrs_ 中的索引
    label  vibLevel[MAXVIBMODE]; // 40B 振动量子数
    label  nVibModes;       //  8B  实际振动模态数
    label  flag;            //  8B  NORMAL/DELETED/NEW
    point  position;        // 24B  位置（sub-cell 分配需要）
};
// ~160B per particle
```

### 2.3 线程级负载平衡策略

**策略 A：Dynamic Scheduling（默认）**
- `#pragma omp for schedule(dynamic, chunk)`
- chunk = max(1, nCells / (nThreads * 16))

**策略 B：加权贪心装箱（高不均匀度时启用）**
- 按 cellCount_[cellI]² 估算碰撞工作量
- 贪心分配 cells 到线程

**自适应切换**：`max(cellCount) / avg(cellCount) > 4` 时启用策略 B。

## 3. 实施阶段

### Phase 0：基础设施（1-2 周）

#### 目标

添加 OpenMP 编译支持，定义碰撞粒子 POD 结构体和线程 RNG，不改变任何现有行为。

#### 具体步骤

1. **修改编译选项** — `src/lagrangian/dsmc/Make/options` 添加 `-fopenmp`
2. **定义 DsmcCollisionData** — 新建 `collisionPartnerSelection/basic/DsmcCollisionData.H`（见 2.2 节）
3. **添加线程 RNG** — `dsmcCloud.H` 新增 `PtrList<Random> threadRng_`，构造时按 `omp_get_max_threads()` 初始化
4. **添加编译宏** — `#ifdef DSMC_OPENMP` 控制 OpenMP 代码路径，便于串行回退

#### 修改文件

- `src/lagrangian/dsmc/Make/options`
- 新建 `src/lagrangian/dsmc/collisionPartnerSelection/basic/DsmcCollisionData.H`
- `src/lagrangian/dsmc/clouds/dsmcCloud.H`（声明 threadRng_）
- `src/lagrangian/dsmc/clouds/dsmcCloud.C`（初始化 threadRng_）

### Phase 1：平坦数组 + Cell 索引（2-3 周）

#### 目标

在碰撞前从 IDLList 构建 parcelPtrs_ 指针数组和 SPARTA 风格 cellFirst_/cellNext_ 索引，碰撞后回写修改。

#### 具体步骤

1. **在 dsmcCloud 中添加成员**

```cpp
DynamicList<dsmcParcel*> parcelPtrs_;   // 索引 → parcel 指针
labelList cellFirst_;                    // 每个 cell 的第一个粒子索引（-1 表示空）
labelList cellNext_;                     // 同 cell 下一个粒子索引（-1 表示末尾）
labelList cellCount_;                    // 每个 cell 的粒子数
```

2. **实现 buildParcelPtrs()**

```cpp
void dsmcCloud::buildParcelPtrs()
{
    parcelPtrs_.clear();
    forAllIter(typename Cloud<dsmcParcel>, *this, iter)
    {
        parcelPtrs_.append(&iter());
    }
}
```

3. **实现 buildCellIndex()**（SPARTA particle::sort() 的等价物）

```cpp
void dsmcCloud::buildCellIndex()
{
    const label nCells = mesh_.nCells();
    const label nPart = parcelPtrs_.size();

    cellFirst_.setSize(nCells, -1);
    cellNext_.setSize(nPart, -1);
    cellCount_.setSize(nCells, 0);

    // 反向遍历构建链表（O(n)，结果为正序）
    for (label i = nPart - 1; i >= 0; i--)
    {
        const label ic = parcelPtrs_[i]->cell();
        cellNext_[i] = cellFirst_[ic];
        cellFirst_[ic] = i;
        cellCount_[ic]++;
    }
}
```

4. **验证**：buildCellIndex() 后 assert 每个 cell 的粒子数与 cellOccupancy_ 一致

#### 修改文件

- `src/lagrangian/dsmc/clouds/dsmcCloud.H`
- `src/lagrangian/dsmc/clouds/dsmcCloud.C`

### Phase 2：OpenMP 碰撞循环 + 线程负载平衡（3-4 周）

#### 目标

实现新的碰撞选择模型 `noTimeCounterOMP`，在平坦数组上用 OpenMP 并行执行碰撞，含线程级负载平衡。

#### 2.1 碰撞数据构建与回写

```cpp
void dsmcCloud::buildCollisionData()
{
    const label nPart = parcelPtrs_.size();
    collisionData_.setSize(nPart);

    forAll(parcelPtrs_, i)
    {
        dsmcParcel& p = *parcelPtrs_[i];
        DsmcCollisionData& cd = collisionData_[i];
        cd.U = p.U();
        cd.ERot = p.ERot();
        cd.RWF = p.RWF();
        cd.ELevel = p.ELevel();
        cd.typeId = p.typeId();
        cd.cellI = p.cell();
        cd.origIdx = i;
        cd.position = p.position();
        cd.nVibModes = min(p.vibLevel().size(), label(MAXVIBMODE));
        for (label m = 0; m < cd.nVibModes; m++)
            cd.vibLevel[m] = p.vibLevel()[m];
        cd.flag = 0; // NORMAL
    }
}

void dsmcCloud::reconcileCollisionData()
{
    forAll(collisionData_, i)
    {
        DsmcCollisionData& cd = collisionData_[i];
        if (cd.flag == PARTICLE_DELETED) { /* 延迟删除 */ continue; }
        dsmcParcel& p = *parcelPtrs_[cd.origIdx];
        p.U() = cd.U;
        p.ERot() = cd.ERot;
        p.ELevel() = cd.ELevel;
        for (label m = 0; m < cd.nVibModes; m++)
            p.vibLevel()[m] = cd.vibLevel[m];
    }
    // 处理 DELETED 和 NEW 粒子（串行）
}
```

#### 2.2 并行碰撞核心循环

```cpp
void noTimeCounterOMP::collide()
{
    cloud_.buildParcelPtrs();
    cloud_.buildCellIndex();
    cloud_.buildCollisionData();

    const label nCells = mesh_.nCells();
    label totalCollisions = 0;

    // 线程局部缓冲
    const int nThreads = omp_get_max_threads();
    List<DynamicList<DsmcCollisionData>> threadNewParticles(nThreads);
    List<DynamicList<label>> threadDelList(nThreads);

    #pragma omp parallel reduction(+:totalCollisions)
    {
        const int tid = omp_get_thread_num();
        Random& rng = cloud_.threadRng(tid);
        DynamicList<label> subCells[8];

        #pragma omp for schedule(dynamic, 64)
        for (label cellI = 0; cellI < nCells; cellI++)
        {
            if (cloud_.cellCount()[cellI] < 2) continue;

            // 构建 sub-cell 分配（同 noTimeCounter 的 8 octant 方式）
            // NTC 碰撞对选择
            // 调用碰撞模型（DsmcCollisionData 重载版本）
            // 反应产生新粒子 → threadNewParticles[tid]
            // 反应删除粒子 → threadDelList[tid]
            totalCollisions += cellCollisions;
        }
    }

    // 串行合并：处理新粒子和删除
    cloud_.reconcileCollisionData();
    mergeNewParticles(threadNewParticles);
    processDeleteList(threadDelList);
}
```

#### 2.3 线程安全分析

| 数据 | 访问模式 | 安全性 |
|------|---------|--------|
| collisionData_[plist[i]] | 每 cell 只被一个线程处理 | ✓ 安全 |
| sigmaTcRMax_[cellI] | per-cell 写入，一 cell 一线程 | ✓ 安全 |
| collisionSelectionRemainder_[cellI] | 同上 | ✓ 安全 |
| threadRng_[tid] | per-thread | ✓ 安全 |
| cloud_.constProps(typeId) | 只读 | ✓ 安全 |
| 新粒子创建 | thread-local buffer | ✓ 安全 |
| 粒子删除 | thread-local 标记 | ✓ 安全 |

#### 2.4 线程级负载平衡

```cpp
void noTimeCounterOMP::buildWorkAssignment()
{
    const label nCells = mesh_.nCells();
    const int nThreads = omp_get_max_threads();

    // 估算每 cell 工作量：碰撞对数 ~ n*(n-1)/2
    scalarList cellWork(nCells);
    forAll(cellWork, c)
    {
        label n = cloud_.cellCount()[c];
        cellWork[c] = scalar(n) * scalar(max(n-1, label(0)));
    }

    // 贪心装箱分配
    scalarList threadLoad(nThreads, 0.0);
    threadWorkList_.setSize(nThreads);

    // 按工作量降序排列 cells，逐个分配给负载最小的线程
    // ...（标准贪心装箱算法）
}
```

自适应选择：
```cpp
scalar maxN = max(cloud_.cellCount());
scalar avgN = scalar(parcelPtrs_.size()) / scalar(nCells);
if (maxN / max(avgN, SMALL) > 4.0)
    useGreedyBinPacking_ = true;  // 高不均匀度
else
    useGreedyBinPacking_ = false; // dynamic scheduling 足够
```

#### 2.5 碰撞模型适配

- `BinaryCollisionModel` 新增重载：`collide(DsmcCollisionData& pP, DsmcCollisionData& pQ, label cellI, scalar cR=-1)`
- `LarsenBorgnakkeVariableHardSphere` 实现该重载
- 碰撞后能量分配函数接受 thread RNG 参数
- **不修改现有接口**，仅添加重载

#### 2.6 反应模型适配

- 反应在碰撞循环内调用（与 noTimeCounter 相同）
- 新增接口：`reaction(DsmcCollisionData& pP, DsmcCollisionData& pQ, DynamicList<DsmcCollisionData>& newParticles)`
- 新粒子写入 thread-local buffer，碰撞结束后串行合并到 cloud
- 删除标记写入 `cd.flag = PARTICLE_DELETED`

#### 修改文件

- 新建 `collisionPartnerSelection/derived/noTimeCounterOMP/noTimeCounterOMP.H`
- 新建 `collisionPartnerSelection/derived/noTimeCounterOMP/noTimeCounterOMP.C`
- `collisions/derived/LarsenBorgnakkeVariableHardSphere/` — 添加 DsmcCollisionData 重载
- `reactions/derived/dissociationQK/` — 延迟创建接口
- `reactions/derived/ionisationQK/` — 延迟创建接口
- `reactions/derived/exchangeQK/` — 延迟创建接口
- `clouds/dsmcCloud.H` — threadRng() 访问器、collisionData_ 成员
- `Make/files` — 添加新源文件

### Phase 3：Move 阶段并行化（4-6 周）

#### 目标

对粒子 move 阶段实现 OpenMP 并行化。这是难度最高的部分。

#### 3.1 可行性分析

`Cloud<dsmcParcel>::move()` 遍历 IDLList，对每个粒子调用 `dsmcParcel::move(td, dt)`。
`dsmcParcel::move()` 内部调用 `trackToFace()` 进行 tet 分解追踪。

| 操作 | 线程安全？ | 原因 |
|------|-----------|------|
| trackToFace() 读取 mesh topology | ✓ | mesh 在 move 期间不变 |
| 修改粒子自身 position/cellI/tetFaceI | ✓ | per-particle 独立 |
| 粒子到达 processor patch → switchProcessor | ✗ | 需要加入迁移列表 |
| 粒子命中壁面 → patch model 回调 | ✗ | 壁面模型可能有共享状态 |
| 粒子命中 cyclic patch | ✗ | cyclic model 可能有状态 |
| 粒子被删除（出域） | ✗ | 修改 IDLList |
| faceTracker 统计通量 | ✗ | 共享计数器 |

#### 3.2 方案：分类并行

**核心思路**：大部分粒子（80-95%）在内部 cell 间自由飞行，不触及任何 boundary patch。只有少数粒子需要 patch 交互。

```cpp
void dsmcCloud::moveParticlesParallel(scalar dt)
{
    buildParcelPtrs();
    const label nPart = parcelPtrs_.size();

    // Step 1: 标记哪些粒子可能触及 boundary
    //   方法：检查粒子所在 cell 是否有 boundary face
    //   预计算 boundaryCell_[cellI] 标记（每步不变，可缓存）
    boolList nearBoundary(nPart, false);
    forAll(parcelPtrs_, i)
    {
        if (boundaryCell_[parcelPtrs_[i]->cell()])
            nearBoundary[i] = true;
    }

    // Step 2: 并行移动内部粒子（无 patch 交互）
    const int nThreads = omp_get_max_threads();
    List<DynamicList<label>> threadMigrateList(nThreads);

    #pragma omp parallel
    {
        const int tid = omp_get_thread_num();

        #pragma omp for schedule(dynamic, 256)
        for (label i = 0; i < nPart; i++)
        {
            if (nearBoundary[i]) continue;  // 跳过，留给串行处理

            dsmcParcel& p = *parcelPtrs_[i];
            // 简化 move：弹道推进 + cell 穿越（仅内部 face）
            moveInternalParticle(p, dt, tid, threadMigrateList[tid]);
        }
    }

    // Step 3: 串行移动 boundary 附近粒子（使用完整 OF 追踪）
    forAll(parcelPtrs_, i)
    {
        if (!nearBoundary[i]) continue;
        dsmcParcel& p = *parcelPtrs_[i];
        p.move(td, dt);  // 完整追踪，含 patch 回调
    }
}
```

#### 3.3 内部粒子追踪的简化实现

```cpp
void dsmcCloud::moveInternalParticle(
    dsmcParcel& p, scalar dt, int tid,
    DynamicList<label>& migrateList
)
{
    // 弹道推进
    const point xNew = p.position() + p.U() * dt;

    // 快速检查：是否仍在当前 cell
    if (mesh_.pointInCell(xNew, p.cell()))
    {
        p.position() = xNew;
        return;
    }

    // 穿越 cell：使用邻居遍历找到新 cell
    // （比 findCell 快，利用 mesh topology）
    label newCell = findCellByNeighbourWalk(p.cell(), xNew);

    if (newCell >= 0 && !boundaryCell_[newCell])
    {
        p.position() = xNew;
        p.cell() = newCell;
        // 更新 tetFaceI_, tetPtI_（需要）
        mesh_.findTetFacePt(newCell, xNew, p.tetFaceI(), p.tetPtI());
    }
    else
    {
        // 到达 boundary cell 或找不到 → 标记为需要精确追踪
        p.flag_ = NEEDS_FULL_TRACKING;
    }
}
```

#### 3.4 替代方案：直接并行化 Cloud::move()

如果简化追踪精度不够，可以直接在 `parcelPtrs_` 上并行调用完整的 `trackToFace()`：

```cpp
#pragma omp parallel
{
    const int tid = omp_get_thread_num();
    DynamicList<label>& myMigrates = threadMigrateList[tid];
    DynamicList<label>& myDeletes = threadDeleteList[tid];

    #pragma omp for schedule(dynamic, 128)
    for (label i = 0; i < nPart; i++)
    {
        dsmcParcel& p = *parcelPtrs_[i];
        // trackToFace 本身是线程安全的（只读 mesh）
        // 但 hitPatch 回调不是 → 需要延迟处理
        bool hitBoundary = moveWithDeferredPatches(p, dt, tid);
        if (hitBoundary) myMigrates.append(i);
    }
}

// 串行处理所有 hit boundary 的粒子
forAll(threadMigrateList, tid) { ... }
```

**推荐**：先实现 3.2 的分类方案（简单、安全），性能不足时再尝试 3.4。

#### 3.5 coordSystem().evolve() 处理

axisymmetric 的 `evolve()` 调用 `addNewParcel()`/`deleteParticle()` 修改粒子拓扑，**必须串行执行**。
它在 move 之后、collide 之前调用，不影响 move 和 collide 的并行化。

#### 修改文件

- `src/lagrangian/dsmc/clouds/dsmcCloud.H` — moveParticlesParallel 声明
- `src/lagrangian/dsmc/clouds/dsmcCloud.C` — 实现 + evolve_moveAndCollide 调用路径
- 可能需要修改 `src/lagrangian/basic/Cloud/Cloud.C`（或绕过，直接在 dsmcCloud 层实现）

### Phase 4：Field 计算与 Controller 并行化（2-3 周）

#### 目标

并行化 `evolve_fields()` 中的场计算和 controller 的粒子循环。

#### 4.1 Field 计算并行化

`dsmcFieldProperties::calculateFields()` 调用每个 field 的 `calculateField()`。
每个 field 独立遍历 cellOccupancy（或新的 cellFirst_/cellNext_）做 per-cell reduction。

```cpp
// 典型 field 计算模式（如 dsmcDensityFluctuationsZone）
void dsmcSomeField::calculateFieldParallel()
{
    const labelList& cells = zone_.cells();
    const label nZoneCells = cells.size();

    #pragma omp parallel for schedule(static)
    for (label c = 0; c < nZoneCells; c++)
    {
        const label cellI = cells[c];
        scalar rhoN = 0;
        vector momentum = vector::zero;

        for (label ip = cellFirst_[cellI]; ip >= 0; ip = cellNext_[ip])
        {
            const dsmcParcel& p = *parcelPtrs_[ip];
            rhoN += p.RWF();
            momentum += p.RWF() * p.U();
        }

        rhoNField_[c] = rhoN / cellVol[cellI];
        momentumField_[c] = momentum / cellVol[cellI];
    }
}
```

**线程安全**：每个 cell 的 reduction 独立写入 per-cell 数组，无冲突。

#### 4.2 Controller 并行化

Controllers 在 `controlBeforeMove()` / `controlBeforeCollisions()` 中直接修改粒子速度。

**关键观察**：每个 controller 操作特定 zone 的粒子，不同 controller 的 zone 通常不重叠。

**方案 A：per-controller 内部并行化**

```cpp
void temperatureController::controlParcelsBeforeMove()
{
    const labelList& cells = controlZone();

    // Step 1: 并行计算统计量（reduction）
    #pragma omp parallel for schedule(static) reduction(+:sumMass,sumMom)
    for (label c = 0; c < cells.size(); c++)
    {
        // 累加质量、动量...
    }

    // Step 2: 计算修正因子 chi_（串行，很快）
    computeChi();

    // Step 3: 并行应用修正
    #pragma omp parallel for schedule(static)
    for (label c = 0; c < cells.size(); c++)
    {
        const label cellI = cells[c];
        for (label ip = cellFirst_[cellI]; ip >= 0; ip = cellNext_[ip])
        {
            dsmcParcel& p = *parcelPtrs_[ip];
            if (findIndex(typeIds_, p.typeId()) != -1)
            {
                p.U() -= UMean_[c];
                p.U() *= chi_[c];
                p.U() += UMean_[c];
            }
        }
    }
}
```

**方案 B：多 controller 间并行（仅当 zone 不重叠时）**

```cpp
// 如果能证明 zone 不重叠，可以并行执行多个 controller
// 但这需要运行时检查，复杂度高，收益有限
// 推荐：先实现方案 A
```

#### 4.3 buildCellOccupancy 并行化

当前 `buildCellOccupancy()` 被频繁调用（move 后、electron 操作后、reaction 后）。
用 `buildCellIndex()` 替代时，可以并行化统计步骤：

```cpp
void dsmcCloud::buildCellIndexParallel()
{
    const label nPart = parcelPtrs_.size();
    const label nCells = mesh_.nCells();

    cellFirst_.setSize(nCells, -1);
    cellNext_.setSize(nPart, -1);
    cellCount_.setSize(nCells, 0);

    // 对于百万粒子级别，串行 O(n) 反向遍历已足够快（~2ms）
    // 仅在 nPart > 5M 时考虑并行版本
    for (label i = nPart - 1; i >= 0; i--)
    {
        const label ic = parcelPtrs_[i]->cell();
        cellNext_[i] = cellFirst_[ic];
        cellFirst_[ic] = i;
        cellCount_[ic]++;
    }
}
```

#### 修改文件

- `macroscopicProperties/derived/` 下各 field 类 — 添加并行版本
- `controllers/derived/temperature/` 等 — 内部循环并行化
- `clouds/dsmcCloud.C` — buildCellIndexParallel

### Phase 5：消除双重存储 + 性能优化（4-6 周，可选）

#### 目标

将所有 cellOccupancy_ 消费者迁移到 cellFirst_/cellNext_，减少内存开销，统一数据路径。

#### 5.1 迁移清单

需要替换 cellOccupancy 访问的模块：

- `collisionPartnerSelection/derived/` — 所有碰撞选择模型（~6 个类）
- `reactions/derived/` — 所有反应模型（~5 个类）
- `coordinateSystem/derived/dsmcAxisymmetric/` — axisymmetricWeighting
- `macroscopicProperties/derived/` — 所有 field 计算类（~10 个类）
- `controllers/derived/` — 所有 controller（~10 个类）
- `boundaries/derived/` — 所有 boundary 模型

替换模式统一：
```cpp
// 旧：
const DynamicList<dsmcParcel*>& cellParcels = cellOccupancy[cellI];
forAll(cellParcels, i) { dsmcParcel& p = *cellParcels[i]; ... }

// 新：
for (label ip = cellFirst_[cellI]; ip >= 0; ip = cellNext_[ip])
{
    dsmcParcel& p = *parcelPtrs_[ip]; ...
}
```

#### 5.2 内存节省

- 移除 `DynamicList<DynamicList<dsmcParcel*>> cellOccupancy_`
  - 当前开销：nCells 个 DynamicList 对象 + nParticles 个指针 ≈ nCells*48B + nPart*8B
- 新结构开销：cellFirst_(nCells*8B) + cellNext_(nPart*8B) + cellCount_(nCells*8B)
- 净节省：约 nCells*32B（DynamicList 管理开销）

#### 5.3 性能优化（可选）

- **NUMA first-touch**：按线程亲和性初始化 parcelPtrs_ 内存页
- **Cache line padding**：避免 sigmaTcRMax_ 等 per-cell 数组的 false sharing
- **粒子排序**：按 cellI 排序 parcelPtrs_ 提升碰撞阶段的 cache 命中率

## 4. 验证策略

| 阶段 | 验证方法 | 通过标准 |
|------|---------|---------|
| Phase 1 | buildCellIndex() 后对比 cellOccupancy_ | 每 cell 粒子数完全一致 |
| Phase 2 | OMP_NUM_THREADS=1 运行 noTimeCounterOMP | 碰撞计数与 noTimeCounter 一致（相同 RNG 种子） |
| Phase 2 | 多线程 heatBath-5species | 稳态温度/密度在统计误差内（<1%） |
| Phase 3 | 并行 move vs 串行 move | 粒子最终位置一致（tolerance 1e-12） |
| Phase 3 | 壁面反射计数 | 并行/串行一致 |
| Phase 4 | 并行 field vs 串行 field | 数值 bitwise 一致 |
| 全局 | hypersonicCorner 1/2/4/8 线程 | 加速比测量 + 结果统计等价 |
| 全局 | 长时间运行（1000+ 步） | 无 crash、无 drift、无内存泄漏 |

## 5. 风险与缓解

| 风险 | 严重度 | 缓解措施 |
|------|--------|----------|
| 反应模型在并行区内调用 addNewParcel | 高 | thread-local buffer + 串行合并 |
| vibLevel_ 动态大小超过 MAXVIBMODE=5 | 中 | 运行时检查 + FatalError 提示 |
| trackToFace 的 patch 回调非线程安全 | 高 | 分类方案：内部粒子并行，boundary 粒子串行 |
| controller 直接修改粒子 U_ | 中 | per-controller 内部并行（zone 内 cell 独立） |
| axisymmetric evolve() 修改粒子拓扑 | 中 | 保持串行，在 move 并行之后执行 |
| 双重存储增加内存 ~160B/粒子 | 低 | Phase 5 消除；或仅碰撞阶段临时分配 |
| RNG 序列改变导致结果不可复现 | 低 | 接受统计等价；保留串行模式用于调试 |
| mesh_.pointInCell / findCell 性能差 | 中 | 使用邻居遍历替代全局搜索 |
| OpenMP fork/join 开销 | 低 | persistent thread pool（OMP 默认行为） |

## 6. 工作量与优先级

| 阶段 | 内容 | 工作量 | 风险 | 覆盖计算占比 |
|------|------|--------|------|-------------|
| Phase 0 | 基础设施 | 1-2 周 | 低 | — |
| Phase 1 | 平坦数组 + Cell 索引 | 2-3 周 | 低 | — |
| Phase 2 | 碰撞并行 + 负载平衡 | 3-4 周 | 中 | 40-70% |
| Phase 3 | Move 并行 | 4-6 周 | 高 | 20-40% |
| Phase 4 | Field/Controller 并行 | 2-3 周 | 低-中 | 5-15% |
| Phase 5 | 消除双重存储 + 优化 | 4-6 周 | 中 | 内存/cache 优化 |

**推荐实施顺序**：0 → 1 → 2 → 4 → 3 → 5

理由：
- Phase 4 难度低、收益确定，在攻克 Move 之前先积累经验
- Phase 3 风险最高，放在有经验后再做
- Phase 5 是优化，非功能性必需

**MVP（最小可行产品）**：Phase 0-2 = 6-9 周，覆盖最大热点
**推荐目标**：Phase 0-2 + 4 = 8-12 周，覆盖 50-85%
**完整目标**：Phase 0-4 = 12-18 周，覆盖 ~90%

## 7. 与 SPARTA 的关键差异

| 方面 | SPARTA | 本方案 |
|------|--------|--------|
| 网格类型 | 笛卡尔/AMR | 非结构化多面体 |
| 粒子追踪 | 直接索引 O(1) | tet 分解 trackToFace |
| 并行模型 | Kokkos (CUDA/OMP/HIP) | 纯 OpenMP |
| 数据结构 | 纯平坦数组 | 混合（IDLList + 指针数组 + 碰撞平坦数组） |
| 负载平衡 | MPI 级 RCB | MPI 级重分区 + 线程级贪心装箱 |
| 壁面交互 | 解析几何求交 | OF patch 回调模型 |
| 反应处理 | 碰撞循环内 | 碰撞循环内（相同） |

## 8. 参考

- SPARTA 源码：`refcode/sparta-24Sep2025/src/`
  - `particle.cpp`（sort）、`update.cpp`（move）、`collide.cpp`（collide）
  - `KOKKOS/particle_kokkos.cpp`、`KOKKOS/update_kokkos.cpp`（并行模式）
  - `fix_balance.cpp`（负载平衡）
- 原始碰撞并行化计划：`doc/worklog/v1706/dsmc_sparta_plan_20260528.md`
- 数据结构选型：`doc/worklog/v1706/dsmc_load_balancing_data_structure_20260528.md`
