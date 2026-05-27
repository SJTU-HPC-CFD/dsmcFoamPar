# DSMC OpenMP 碰撞并行化：实施计划

日期：2026-05-28

## 1. 背景与目标

hyStrath_dlb 的 DSMC 碰撞阶段占总计算时间 40-70%，且粒子分布不均导致严重的负载不平衡。当前粒子存储采用侵入式双向链表（IDLList），无法 OpenMP 并行化。目标是引入 SPARTA 风格的平坦数组 + OpenMP dynamic scheduling，在最小改动下实现碰撞阶段的线程级负载平衡。

### 核心矛盾

OpenFOAM 的粒子追踪（Cloud::move）深度耦合于 IDLList + tet 分解追踪架构，完全替换风险极高。

### 解决方案

采用**混合方案**：保留 Cloud<dsmcParcel> 用于追踪/IO/MPI，新增平坦碰撞数组用于碰撞阶段的 OpenMP 并行。

## 2. 方案概述

采用分阶段实施，每个阶段可独立测试验证：

- Phase 0：基础设施（OpenMP 编译支持、碰撞粒子结构体定义）
- Phase 1：构建平坦碰撞数组 + SPARTA 风格 cell 索引
- Phase 2：实现 OpenMP 碰撞循环（核心目标）
- Phase 3（可选）：消除双重存储，统一数据结构

Phase 0-2 为必须完成的工作，预计 6-9 周。Phase 3 为后续优化。

## 3. Phase 0：基础设施（1-2 周）

### 目标

添加 OpenMP 编译支持，定义碰撞粒子 POD 结构体，不改变任何现有行为。

### 具体步骤

1. **修改编译选项**
   - `src/lagrangian/dsmc/Make/options`：添加 `-fopenmp` 到 EXE_INC 和 LIB_LIBS

2. **定义碰撞粒子结构体** — 新建 `src/lagrangian/dsmc/collisionPartnerSelection/basic/DsmcCollisionData.H`

```cpp
enum { MAXVIBMODE = 4 };

struct DsmcCollisionData
{
    double x[3];        // 24B  位置（反应创建新粒子需要）
    double v[3];        // 24B  速度
    double erot;        //  8B  转动能
    double evib;        //  8B  振动能（标量，兼容 SPARTA）
    double RWF;         //  8B  径向权重因子
    label  ELevel;      //  8B  电子能级
    label  typeId;      //  8B  组分索引
    label  cellI;       //  8B  所在单元
    label  origIdx;     //  8B  在 parcelPtrs_ 中的索引（用于回写）
    label  vibLevel[MAXVIBMODE]; // 32B 振动量子数（固定大小）
    label  nVibModes;   //  8B  实际振动模态数
    label  flag;        //  8B  状态标志（正常/待删除/新创建）
};
// 约 152B per particle
```

3. **添加线程 RNG** — 在 dsmcCloud 中添加 `PtrList<Random> threadRng_`，构造时按线程数初始化

### 修改文件

- `src/lagrangian/dsmc/Make/options`
- 新建 `src/lagrangian/dsmc/collisionPartnerSelection/basic/DsmcCollisionData.H`
- `src/lagrangian/dsmc/clouds/dsmcCloud.H`（声明 threadRng_）
- `src/lagrangian/dsmc/clouds/dsmcCloud.C`（初始化 threadRng_）

## 4. Phase 1：平坦碰撞数组 + Cell 索引（2-3 周）

### 目标

在碰撞前从 IDLList 构建平坦数组和 SPARTA 风格 first[]/next[] 索引，碰撞后回写修改。与现有 cellOccupancy_ 并行运行用于验证。

### 具体步骤

1. **在 dsmcCloud 中添加成员**

```cpp
label nCollisionParticles_;
label maxCollisionParticles_;
DsmcCollisionData* collisionData_;       // 平坦碰撞数组
DynamicList<dsmcParcel*> parcelPtrs_;    // 索引 → IDLList 指针映射
labelList cellFirst_;                     // 每个 cell 的第一个粒子索引
labelList cellNext_;                      // 同 cell 下一个粒子索引
```

2. **实现 buildCollisionData()**
   - 遍历 IDLList 一次（forAllIter）
   - 将每个粒子的碰撞相关字段拷贝到 collisionData_[i]
   - 存储 parcelPtrs_[i] = &parcel（用于回写）
   - 设置 collisionData_[i].origIdx = i
   - 处理 vibLevel：拷贝 min(parcel.vibLevel().size(), MAXVIBMODE) 个元素

3. **实现 buildCellIndex()**（SPARTA sort() 的等价物）
   - 初始化 cellFirst_[all] = -1
   - 反向遍历 collisionData_：

```cpp
for (label i = nCollisionParticles_-1; i >= 0; i--) {
    label ic = collisionData_[i].cellI;
    cellNext_[i] = cellFirst_[ic];
    cellFirst_[ic] = i;
}
```

4. **实现 reconcileCollisionData()**
   - 遍历 collisionData_，将修改后的 v, erot, evib, ELevel, vibLevel 写回 parcelPtrs_[origIdx]
   - 处理标记为"待删除"的粒子：调用 deleteParticle()
   - 处理新创建的粒子（来自反应）：调用 addNewParcel()

5. **验证**：在 evolve_moveAndCollide() 中，buildCollisionData() 后对比 cellFirst_/cellNext_ 与 cellOccupancy_ 的每 cell 粒子数，assert 一致。

### 修改文件

- `src/lagrangian/dsmc/clouds/dsmcCloud.H`
- `src/lagrangian/dsmc/clouds/dsmcCloud.C`

## 5. Phase 2：OpenMP 碰撞循环（3-4 周）

### 目标

实现新的碰撞选择模型 `noTimeCounterOMP`，在平坦数组上用 OpenMP 并行执行碰撞。

### 具体步骤

1. **新建 noTimeCounterOMP 类**（继承 collisionPartnerSelection）
   - 运行时通过 dsmcProperties 字典选择，不影响现有 noTimeCounter
   - 核心 collide() 方法结构：

```cpp
void noTimeCounterOMP::collide()
{
    cloud_.buildCollisionData();
    cloud_.buildCellIndex();

    const label nCells = mesh_.nCells();
    label totalCollisions = 0;

    #pragma omp parallel reduction(+:totalCollisions)
    {
        const int tid = omp_get_thread_num();
        Random& rng = cloud_.threadRng(tid);
        DynamicList<label> localDelList;
        DynamicList<DsmcCollisionData> localNewParticles;
        List<DynamicList<label>> subCells(8);

        #pragma omp for schedule(dynamic, 64)
        for (label cellI = 0; cellI < nCells; cellI++)
        {
            label ip = cloud_.cellFirst()[cellI];
            if (ip < 0) continue;

            // 构建 plist（同 SPARTA collisions_one）
            // NTC 碰撞逻辑...
            // 反应产生新粒子 → localNewParticles
            // 反应删除粒子 → localDelList
        }
    }

    // 串行合并：处理新粒子和删除
    cloud_.reconcileCollisionData();
}
```

2. **线程安全分析**
   - `collisionData_[plist[i]]`：每个 cell 只被一个线程处理，cell 内粒子不会被其他线程访问 → 安全
   - `sigmaTcRMax_[cellI]`：per-cell 写入，一个 cell 一个线程 → 安全
   - `collisionSelectionRemainder_[cellI]`：同上 → 安全
   - RNG：per-thread → 安全
   - 新粒子创建：thread-local buffer → 安全
   - `cloud_.constProps(typeId)`：只读 → 安全

3. **碰撞模型适配**
   - `BinaryCollisionModel::collide()` 和 `sigmaTcR()` 当前接受 `dsmcParcel&`
   - 新增重载版本接受 `DsmcCollisionData&`，或创建轻量适配器
   - 主要碰撞模型：`LarsenBorgnakkeVariableHardSphere`
   - 碰撞后能量分配函数（`postCollisionRotationalEnergy`, `postCollisionVibrationalEnergyLevel`）需要接受 thread RNG

4. **反应模型适配**
   - 反应模型（dissociationQK, ionisationQK, exchangeQK 等）在碰撞循环内被调用
   - 当前调用 `cloud_.addNewParcel()` 创建新粒子 — 非线程安全
   - 解决方案：反应模型新增接口，将新粒子数据写入 thread-local `DynamicList<DsmcCollisionData>`
   - 碰撞循环结束后串行合并

5. **运行时选择**
   - 在 dsmcProperties 中：`collisionPartnerSelectionModel noTimeCounterOMP;`
   - 保留 noTimeCounter 作为串行回退

### 修改文件

- 新建 `src/lagrangian/dsmc/collisionPartnerSelection/derived/noTimeCounterOMP/noTimeCounterOMP.H`
- 新建 `src/lagrangian/dsmc/collisionPartnerSelection/derived/noTimeCounterOMP/noTimeCounterOMP.C`
- `src/lagrangian/dsmc/collisions/derived/LarsenBorgnakkeVariableHardSphere/` — 添加 DsmcCollisionData 重载
- `src/lagrangian/dsmc/reactions/derived/dissociationQK/dissociationQK.C` — 延迟创建接口
- `src/lagrangian/dsmc/reactions/derived/ionisationQK/ionisationQK.C` — 延迟创建接口
- `src/lagrangian/dsmc/clouds/dsmcCloud.H` — threadRng() 访问器
- `src/lagrangian/dsmc/Make/files` — 添加新源文件

## 6. Phase 3（可选）：消除双重存储（4-6 周）

### 目标

将 cellOccupancy_ 的所有消费者迁移到 cellFirst_/cellNext_，减少内存开销。

### 步骤概要

- 替换所有 controller 中的 cellOccupancy 访问（~10 个类）
- 替换所有 reaction 中的 cellOccupancy 访问（~5 个类）
- 替换 coordSystem 中的 cellOccupancy 访问
- 替换 macroscopicProperties 中的 cellOccupancy 访问
- 最终移除 cellOccupancy_ 成员和 buildCellOccupancy()

此阶段改动面广但机械性强，每个文件的修改模式相同：

```cpp
// 旧：
const DynamicList<dsmcParcel*>& cellParcels = cellOccupancy[cellI];
forAll(cellParcels, i) { dsmcParcel& p = *cellParcels[i]; ... }

// 新：
for (label ip = cellFirst_[cellI]; ip >= 0; ip = cellNext_[ip]) {
    DsmcCollisionData& p = collisionData_[ip]; ...
}
```

## 7. 验证策略

1. **Phase 1 验证**：buildCollisionData() 后 assert 每个 cell 的粒子数与 cellOccupancy_ 一致
2. **Phase 2 单线程验证**：OMP_NUM_THREADS=1 运行 noTimeCounterOMP，对比 noTimeCounter 的碰撞计数（相同 RNG 种子应完全一致）
3. **Phase 2 多线程统计验证**：多线程运行 heatBath-5species 算例，验证稳态温度/密度与串行结果在统计误差内一致
4. **性能测试**：在 hypersonicCorner 或 orion107kmNR 算例上测量 1/2/4/8 线程的碰撞阶段加速比

## 8. 风险与缓解

| 风险 | 缓解措施 |
|------|----------|
| 反应模型在并行区内调用 addNewParcel | 延迟创建：thread-local buffer + 串行合并 |
| vibLevel 超过 MAXVIBMODE=4 | 编译期 static_assert + 运行时检查 |
| 双重存储增加内存（~152B/粒子额外） | Phase 3 消除；或接受（碰撞阶段临时分配） |
| 随机数序列改变导致结果不可复现 | 接受统计等价；保留串行模式用于调试 |
| BinaryCollisionModel 接口改动影响其他碰撞模型 | 仅添加重载，不修改现有接口 |

## 9. 工作量估算

| 阶段 | 工作量 | 风险等级 |
|------|--------|----------|
| Phase 0 | 1-2 周 | 低 |
| Phase 1 | 2-3 周 | 低-中 |
| Phase 2 | 3-4 周 | 中 |
| Phase 3 | 4-6 周 | 中-高 |
| **总计（Phase 0-2）** | **6-9 周** | |

## 10. 参考

- SPARTA 源码：`refcode/sparta-24Sep2025/src/particle.h`, `collide.cpp`
- 数据结构选型分析：`doc/worklog/v1706/dsmc_load_balancing_data_structure_20260528.md`

