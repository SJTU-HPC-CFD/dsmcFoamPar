# DSMC OpenMP 全流程并行与负载均衡：详细工作计划

日期：2026-05-30

状态：详细实施计划草案

相关文档：

- `doc/worklog/v1706/dsmc_sparta_plan_20260528.md`
- `doc/worklog/v1706/dsmc_load_balancing_data_structure_20260528.md`
- `doc/worklog/v1706/dsmc_full_openmp_plan_20260530.md`

参考代码：

- `work/v1706-sparta/src/lagrangian/dsmc/`
- `refcode/sparta-24Sep2025/src/particle.h`
- `refcode/sparta-24Sep2025/src/particle.cpp`
- `refcode/sparta-24Sep2025/src/collide.cpp`
- `refcode/sparta-24Sep2025/src/update.cpp`
- `refcode/sparta-24Sep2025/src/KOKKOS/update_kokkos.cpp`
- `refcode/sparta-24Sep2025/src/KOKKOS/collide_vss_kokkos.cpp`
- `refcode/sparta-24Sep2025/src/balance_grid.cpp`
- `refcode/sparta-24Sep2025/src/fix_balance.cpp`

## 1. 背景

### 1.1 原计划定位

`dsmc_sparta_plan_20260528.md` 的目标是 DSMC 碰撞阶段 OpenMP 并行化。该计划的核心判断是正确的：

1. 当前 `dsmcCloud` 继承 `Cloud<dsmcParcel>`，粒子主体是 OpenFOAM Lagrangian 链表结构，不适合 OpenMP 随机访问和分块调度。
2. 当前 `cellOccupancy_` 是 `DynamicList<DynamicList<dsmcParcel*>>`，每步从链表重建，适合串行模型但不适合作为完整并行主数据结构。
3. SPARTA 使用连续粒子数组 + cell 索引，碰撞阶段天然可以按 cell 并行。

但该计划只覆盖 collision partner selection 和 binary collision，尚未覆盖 DSMC 时间步中同样重要的 move、boundary、coordSystem、field、controller、reaction 创建/删除和统计输出。若目标扩大为完整 OpenMP 并行和负载均衡，原计划应作为第一阶段原型，而不是最终架构。

### 1.2 当前 `work/v1706-sparta` 主链路

当前 DSMC 时间步主流程位于：

- `src/lagrangian/dsmc/clouds/dsmcCloud.C`
- `dsmcCloud::evolve()`
- `dsmcCloud::evolve_moveAndCollide()`
- `dsmcCloud::evolve_fields()`

核心顺序如下：

```text
controllers_.controlBeforeMove()
boundaries_.controlBeforeMove()
removeElectrons() + buildCellOccupancy()            // 如果有电子
Cloud<dsmcParcel>::move(td, deltaTValue())          // 串行链表追踪
buildCellOccupancy()
addElectrons() + buildCellOccupancy()               // 如果有电子
coordSystem().evolve()
controllers_.controlBeforeCollisions()
boundaries_.controlBeforeCollisions()
collisions()
buildCellOccupancy()                                // 如果有 reactions
controllers_.controlAfterCollisions()
boundaries_.controlAfterCollisions()
reactions_.outputData()
fields_.calculateFields()
fields_.writeFields()
controllers_.calculateProps()
boundaries_.calculateProps()
boundaryMeas_.outputResults()
```

这说明完整 OpenMP 并行不能只替换 `collisionPartnerSelection`。真正的目标必须覆盖：

1. 粒子存储与 cell 索引。
2. 粒子运动与边界交互。
3. 碰撞与反应。
4. 坐标系权重与 clone/delete。
5. controllers、boundaries、fields、measurements。
6. 每阶段的线程级负载均衡。

### 1.3 实际代码结构

当前 DSMC 库的主要模块位于 `src/lagrangian/dsmc/`：

| 模块 | 关键文件 | 当前并行化相关问题 |
|------|----------|--------------------|
| Cloud 主体 | `clouds/dsmcCloud.H`, `clouds/dsmcCloud.C` | 持有 `cellOccupancy_`、`rndGen_`、sub-models、主时间步流程 |
| Parcel | `parcels/dsmcParcel.H`, `parcels/dsmcParcel.C` | `move()` 深度耦合 OpenFOAM tracking、boundary、tracker |
| Collision selection | `collisionPartnerSelection/basic`, `collisionPartnerSelection/derived/noTimeCounter` | 当前遍历 `cellOccupancy()`，使用 `cloud_.randomLabel()` |
| Binary collision | `collisions/basic/BinaryCollisionModel`, `collisions/derived/*` | 接口接受 `dsmcParcel&`，内部使用 `cloud_.rndGen()` |
| Reactions | `reactions/basic`, `reactions/derived/*QK` | 接口接受 `dsmcParcel&`，可能调用 `addNewParcel()` 和修改粒子类型 |
| Coordinate system | `coordinateSystem/derived/axisymmetric`, `coordinateSystem/derived/spherical` | 遍历 `cellOccupancy()`，可能 clone/delete |
| Boundaries | `boundaries/basic`, `boundaries/derived/*` | patch hit 直接调用 `controlParticle()`，模型内使用随机数、测量、add/delete |
| Controllers | `controllers/basic`, `controllers/derived/*` | 多数直接访问 `cellOccupancy()`，部分 add/delete 粒子 |
| Fields | `macroscopicProperties/basic/dsmcFieldProperties`, `dsmcStandardFields` | 当前遍历 `dsmcCloud` 链表并写 volFields |
| Measurements | `faceTracker`, `boundaryMeasurements`, `porousMeasurements`, `cellMeasurements` | 存在共享累加状态 |
| Build | `Make/files`, `Make/options` | 需要新增 OpenMP 编译选项和新源文件 |

### 1.4 吸收 Claude 版计划的有用修正

Claude 版计划最有价值的地方，是把完整重构拆成了更贴近 `work/v1706-sparta` 现状的两层路线：短期先桥接现有 `Cloud<dsmcParcel>`，长期再逐步把连续粒子结构变成时间步主结构。这个修正应纳入本计划：

1. 第一轮不应直接把 move 阶段 POD 化。`dsmcParcel::move()` 深度依赖 OpenFOAM `particle` 基类、tet tracking、`trackToFace()`、patch 回调、`trackingData` 和 `stepFraction()`，直接改成独立平坦结构风险过高。
2. 短期先在 `dsmcCloud` 中增加 `parcelPtrs_`、`cellFirst_`、`cellNext_`、`cellCount_` 和 per-thread RNG，以指针数组方式获得 O(1) 随机访问和 SPARTA 风格 cell 链表索引。
3. 碰撞阶段可以使用临时 `DsmcCollisionData`，只复制碰撞需要的字段；真正的 IO、MPI、tracking、boundary 回调仍以 `Cloud<dsmcParcel>` 为准。
4. `cellOccupancy_` 的消费者应渐进迁移，优先 collision 和 fields；不能假设所有 controller、boundary、coordSystem 都能一次切到新 store。
5. move 的“内部粒子并行、边界粒子串行”思路可作为原型方向，但不能把简化的 `pointInCell + neighbourWalk` 当成正式替代 `trackToFace()` 的可靠实现。
6. 推荐实施顺序调整为：Phase 0 基线和开关，Phase 1 `parcelPtrs_` 桥接索引，Phase 2 碰撞最小并行，Phase 3 fields/tally，Phase 4 线程级负载均衡，之后再进入 move 原型。

因此，本计划后续将同时保留两条线：MVP 桥接线负责尽快落地和验证收益，长期 Store 线负责最终完整 OpenMP 架构。

## 2. 目标

### 2.1 总目标

实现 `work/v1706-sparta` DSMC 求解器的 rank 内完整 OpenMP 并行化，并建立可扩展的负载均衡机制。

完整目标包括：

1. 建立 SPARTA 风格连续粒子主数据结构。
2. 将碰撞阶段迁移到 OpenMP cell 级并行。
3. 将 cell index 构建迁移到平坦数组路径。
4. 将 move 阶段逐步迁移到平坦数组路径。
5. 将常用 boundary/controller/coordSystem/field 迁移到线程安全路径。
6. 对创建、删除、迁移、测量、统计使用 thread-local buffer 和同步点合并。
7. 建立 per-cell work history，用于 OpenMP 动态调度和未来 MPI 级负载均衡。

### 2.2 性能目标

目标按 8 OpenMP threads 估算，但必须区分 MVP 路线和完整路线。MVP 路线指 M0-M4：move、复杂 boundary/controller、未迁移 reaction 仍可能串行；完整路线指 M0-M9：move、boundary/controller、reaction delayed event 和 legacy cleanup 完成。

| 阶段 | 当前问题 | 目标加速 |
|------|----------|----------|
| cell index / occupancy | 链表遍历 + pointer lists | 3-6x |
| collision，无 reaction | cell 串行循环 | 4-7x |
| collision，reaction 未迁移 | 整体回退旧 `noTimeCounter` | 接近 1x |
| move | 链表串行 tracking | 2-5x |
| fields / tally | 链表遍历 + shared field 写入 | 4-7x |
| controller / boundary | 模型分散、共享状态多 | 1.5-4x |
| MVP 总步进，move 串行 | Amdahl 受限 | 1.5-3x |
| 当前 profiling 算例 MVP，总步进 | `move≈50%`、`collision≈30%` | 1.3-1.7x |
| 完整总步进 | 多阶段串行瓶颈基本消除 | 3-6x |

性能目标是统计意义上的工程目标，不承诺逐步 bitwise reproducibility。

当前 profiling 目标算例已知大致为 `move≈50%`、`collision≈30%`。若进一步假设 `fields=10%`、`buildIndex=5%`、`other=5%`，即使 collide 和 fields 都达到 6x、buildIndex 达到 4x，MVP 总时间仍为：

```text
T_new = 0.50 + 0.30/6 + 0.10/6 + 0.05/4 + 0.05 = 0.6292
Speedup = 1.59x
```

若只并行 collision 而 fields/index 仍未充分优化，则上限更接近：

```text
T_new = 0.50 + 0.30/6 + 0.20 = 0.75
Speedup = 1.33x
```

因此，本项目当前 profiling 算例下，MVP 的总加速验收不能使用完整路线的 3-6x 目标，也不应把 Phase 5 视为后置可选项。Phase 5 move prototype 应在 M2 最小 collision 原型后立即启动，并与 M3/M4 并行推进。

### 2.3 正确性目标

1. `OMP_NUM_THREADS=1` 时，新路径应与旧路径尽可能接近，至少保持统计一致。
2. `OMP_NUM_THREADS>1` 时，密度、温度、速度、热流、反应率等宏观量应落在 DSMC 统计误差内。
3. 粒子数、cell count、species count、质量、动量、能量、boundary event 数应有明确检查工具。
4. 原串行路径必须可回退，便于定位问题。

### 2.4 非目标

第一轮不承诺：

1. GPU/Kokkos/Cabana 后端。
2. 在线 MPI cell 迁移。
3. 所有 WIP boundary/controller/field 模型立即线程安全。
4. 多线程下随机数序列与串行完全一致。
5. 对单个超大 cell 做 cell 内碰撞并行。

## 3. 可行性分析

### 3.1 数据结构可行性

SPARTA 的 `Particle::OnePart` 是连续数组，保存 `x/v/erot/evib/dtremain/weight/ispecies/icell/flag`。`Particle::sort()` 通过 `cinfo.first/count` 和 `particle.next` 建立 cell 索引。该模型适合 `work/v1706-sparta`：

1. 当前 DSMC 主要算法仍是 cell-local。
2. `noTimeCounter` 已按 cell 遍历 `cellOccupancy()`，迁移到 `cellFirst/cellNext` 自然。
3. `sigmaTcRMax_` 和 `collisionSelectionRemainder_` 是 per-cell 数据，适合 cell 并行。
4. fields、controllers、coordSystem 多数也按 cell 遍历，后续可统一迁移。

结论：连续数组 + cell 索引可行，且是完整 OpenMP 的必要前提。

### 3.2 碰撞阶段可行性

当前 `noTimeCounter.C` 的碰撞逻辑：

1. 遍历 `cloud_.cellOccupancy()`。
2. 每个 cell 构建 subcell lists。
3. 基于 `nC`、`sigmaTcRMax[cellI]`、`collisionSelectionRemainder[cellI]` 计算 candidate 数。
4. 随机选粒子对。
5. 调用 `binaryCollision().sigmaTcR()` 和 `binaryCollision().collide()`。
6. 调用 reaction model。

并行化可行原因：

1. 外层 cell 独立。
2. 同一粒子只属于一个 cell。
3. per-cell 状态写入无交叉。
4. 新增/删除粒子可延迟到同步点。

主要改造点：

1. `BinaryCollisionModel` 需要增加 `DsmcParticle` 或 `DsmcParticleView` 接口。
2. reaction model 需要增加延迟创建/删除接口。
3. RNG 需要线程安全。
4. 统计输出需要 reduction。

结论：碰撞阶段是最适合第一阶段实施的模块。

### 3.3 Move 阶段可行性

当前 `dsmcParcel::move()` 依赖：

1. `particle::trackToFace()` 和 tet 分解 tracking。
2. `dsmcParcel::trackingData`。
3. `boundaries().patchBoundaryModels()[...]->controlParticle()`
4. `boundaries().cyclicBoundaryModels()[...]->controlMol()`
5. `tracker().trackParcelFaceTransition()`
6. `rndGen()` 初始化新 parcel 的 `stepFraction()`。

并行化难度高，但不是不可行。关键是不能在并行 tracking 中直接写共享对象，而应拆成：

1. 粒子几何推进 kernel。
2. boundary hit event 记录。
3. event 合并与 boundary model 应用。
4. 粒子 flag 更新与 compact。

短期可行策略：

1. 第一轮保留 `Cloud<dsmcParcel>::move()` 串行，只并行 collision/index/fields。
2. 第二轮基于 `parcelPtrs_` 和 OpenFOAM tracking 状态抽取无复杂 boundary 的受限 move 原型。
3. 第三轮迁移常用 wall boundary。
4. 最后处理 cyclic、sticking、catalytic、field patch 等复杂模型。

结论：move 可行但风险最高，必须分阶段，不应阻塞 collision 原型。

### 3.4 Fields 与 tally 可行性

当前 `dsmcStandardFields::calculateFields()` 遍历 `dsmcCloud` 链表，并写 `rhoN/rhoM/linearKE/...` 等 volFields。并行化方式明确：

1. 遍历连续粒子数组。
2. 使用 per-thread cell field buffers 或按 cell 独占写。
3. 最后 reduce 到 OpenFOAM fields。

结论：fields 阶段中低风险，适合作为 collision 后的第二批并行目标。

### 3.5 Boundary/controller 可行性

当前 boundary 和 controller 模型数量多，且接口都是面向 `dsmcParcel&`：

- `dsmcPatchBoundary::controlParticle(dsmcParcel&, trackingData&)`
- `dsmcStateController::controlParcelsBeforeMove()`
- `dsmcStateController::controlParcelsBeforeCollisions()`
- `dsmcStateController::controlParcelsAfterCollisions()`

可行策略是分级迁移：

| 等级 | 含义 | 处理方式 |
|------|------|----------|
| Level 0 | 不改，串行回退 | WIP、复杂全局状态模型 |
| Level 1 | 只读/只改当前粒子 | 直接 `DsmcParticleView` 并行 |
| Level 2 | 会新增/删除粒子 | thread-local event buffer |
| Level 3 | 会写测量/输出/全局统计 | per-thread tally + reduce |

结论：不能一次性全部并行，但常用模型可逐步迁移。

## 4. 总体架构设计

### 4.1 核心原则

1. 新并行路径以连续数组为主，`Cloud<dsmcParcel>` 作为 IO/MPI/legacy 镜像。
2. 不在 OpenMP 并行区直接调用 `addParticle()`、`deleteParticle()`、`addNewParcel()`、`deleteParticle()`。
3. 不在 OpenMP 并行区直接写共享 measurement、volField、controller 全局状态。
4. 所有跨线程共享写入通过 thread-local buffer、reduction、atomic 或同步点完成。
5. 每个阶段都保留串行回退开关。

### 4.2 新增核心数据结构

本计划采用“双层数据结构”路线。第一层是 MVP 桥接结构，直接挂在 `dsmcCloud` 上，优先支撑 collision、cell index、fields 和负载统计；第二层才是长期 `DsmcParticleStore`，用于逐步替换 legacy 链表消费者。

第一轮建议先在 `dsmcCloud` 中增加：

```cpp
DynamicList<dsmcParcel*> parcelPtrs_;
labelList cellFirst_;
labelList cellNext_;
labelList cellCount_;
PtrList<Random> threadRng_;
```

职责划分：

1. `Cloud<dsmcParcel>` 仍是真实粒子容器，负责 IO、MPI、OpenFOAM tracking、patch 回调和 legacy 模型。
2. `parcelPtrs_` 是每步可重建的随机访问桥接数组，避免在碰撞/统计中反复遍历 IDLList。
3. `cellFirst_`、`cellNext_`、`cellCount_` 是 SPARTA 风格 cell 链表索引，用于替代高成本的 `cellOccupancy_` 遍历路径。
4. `threadRng_` 为 collision、field sampling、future event path 提供线程本地随机数入口。
5. 碰撞阶段用临时 `DsmcCollisionData` 做 cache-friendly 计算，结束后只回写碰撞改变的字段。

这个 MVP 层不要求立刻移除 `cellOccupancy_`，而是允许在每个阶段结束后调用 `buildCellOccupancy()` 供未迁移模型继续使用。这样可以把风险限制在单个模块内。

长期 Store 层仍然建议保留，作为第二阶段和第三阶段重构目标。

建议新增目录：

```text
src/lagrangian/dsmc/particleStore/
```

建议新增文件：

```text
particleStore/DsmcParticle.H
particleStore/DsmcParticleStore.H
particleStore/DsmcParticleStore.C
particleStore/DsmcCellIndex.H
particleStore/DsmcEventBuffers.H
particleStore/DsmcThreadContext.H
particleStore/DsmcLoadBalancer.H
particleStore/DsmcLoadBalancer.C
```

长期版 `DsmcParticle`：

```cpp
enum DsmcParticleFlag
{
    DSMC_ACTIVE = 0,
    DSMC_DELETE = 1,
    DSMC_NEW = 2,
    DSMC_MIGRATE = 3,
    DSMC_EXIT = 4,
    DSMC_STUCK = 5
};

enum { MAXVIBMODE = 5 };

struct DsmcParticle
{
    scalar x[3];
    scalar U[3];
    scalar ERot;
    scalar RWF;
    scalar stepFraction;

    label ELevel;
    label typeId;
    label cellI;
    label tetFaceI;
    label tetPtI;
    label newParcel;
    label classification;
    label flag;

    label vibLevel[MAXVIBMODE];
    label nVibModes;

    label origId;
};
```

说明：

1. 这是完整时间步所需字段，不只是碰撞字段；第一轮不强制把 move 改到这个结构上。
2. `origId` 用于和 legacy parcel 映射，早期可以是 `parcelPtrs_` index。
3. stuck particle 的 wallTemperature/wallVectors 建议先用 side arrays，避免主结构膨胀。
4. 当前 `dsmcParcel::vibLevel_` 是动态 `labelList`，`MAXVIBMODE` 只是长期 POD 结构的初始固定上限；若任一 species 的 vib modes 超出上限，第一版应报错并回退串行路径，后续再实现 overflow side storage。

### 4.3 `DsmcParticleStore`

职责：

1. 管理 `DsmcParticle* particles_`、`nLocal_`、`capacity_`。
2. 从 `Cloud<dsmcParcel>` 构建 store。
3. 将 store 回写到 `Cloud<dsmcParcel>`。
4. 管理 `parcelPtrs_` 兼容映射。
5. 提供 append、markDelete、compact、reserve。
6. 管理 cell index。

第一版可挂在 `dsmcCloud` 内：

```cpp
autoPtr<DsmcParticleStore> particleStore_;
Switch useOpenMPDsmc_;
Switch syncStoreEveryStep_;
```

### 4.4 Cell index

保留两套视图：

1. SPARTA 风格：

```cpp
labelList cellFirst_;
labelList cellNext_;
labelList cellCount_;
```

2. 排序 range 风格：

```cpp
labelList cellStart_;
labelList sortedParticleIds_;
```

第一阶段实现 `cellFirst/cellNext/cellCount` 即可。后续 fields/move 优化时再引入 `cellStart/sortedParticleIds`。

### 4.5 Thread context

每个 OpenMP thread 持有：

```cpp
struct DsmcThreadContext
{
    Random rng;
    DynamicList<label> deleteList;
    DynamicList<DsmcParticle> newParticles;
    DynamicList<DsmcBoundaryEvent> boundaryEvents;
    DynamicList<DsmcTallyEvent> tallyEvents;
    label collisionAttempts;
    label acceptedCollisions;
    label reactions;
};
```

初期使用 per-thread `Random`。后续引入 counter/hash RNG，以减少调度依赖。

### 4.6 Event buffers

需要明确所有并行阶段的延迟提交事件：

| 事件 | 来源 | 提交动作 |
|------|------|----------|
| delete particle | reaction、boundary、controller、coordSystem | mark delete + compact |
| create particle | reaction、inlet、clone | append to store or addNewParcel |
| boundary hit | move | 调用对应 boundary model 或 tally |
| face transition | move | reduce 到 `dsmcFaceTracker` |
| field contribution | fields | reduce 到 volFields |
| porous interaction | add/delete/boundary | reduce 到 porousMeasurements |

## 5. 详细实施计划

### Phase 0：基线、编译与计时基础设施

目标：建立安全开关、计时和回归基线，不改变默认行为。

预计工作量：1 周

修改文件：

- `src/lagrangian/dsmc/Make/options`
- `src/lagrangian/dsmc/clouds/dsmcCloud.H`
- `src/lagrangian/dsmc/clouds/dsmcCloud.C`
- `doc/scripts/build-dsmcFoam.sh`
- `doc/scripts/run-*.sh`

任务：

1. 增加可选 OpenMP 编译选项。
2. 增加 `DSMC_OPENMP` 或 dictionary switch，例如：

```text
openMPDsmc
{
    enabled false;
    nThreads 8;
    parallelCollision false;
    parallelFields false;
    parallelMove false;
}
```

3. 在 `evolve_moveAndCollide()` 和 `evolve_fields()` 增加阶段 wall-time 计时。OpenMP 性能分析必须以 wall time 为主，CPU time 只能作为辅助输出。
4. 输出每阶段耗时和并行诊断量：

```text
DSMC timings:
  controllersBeforeMove
  boundariesBeforeMove
  move
  buildParcelPtrs
  buildCellIndex
  buildCollisionData
  collisionKernel
  reconcileCollisionData
  buildCellOccupancy
  coordSystem
  collisions
  fieldsKernel
  fieldsReduce
  totalStep

DSMC parallel diagnostics:
  nParticles
  nCells
  nNonEmptyCells
  maxCellCount / avgNonEmptyCellCount
  collisionCandidates
  acceptedCollisions
  fallbackCells
  fallbackParticles
  collisionDataBytesCopied
  fieldBufferBytes
  fieldMode
  threadWorkMax / threadWorkAvg
```

5. 建立 baseline cases：

- heatBath-5species
- 2D cylinder / wall-heavy case
- hypersonicCorner 或 orion107kmNR

6. 在 Phase 0 期间冻结 4 个实现决策，避免 Phase 1/2 开始后接口反复变化：

| 决策 | 推荐方案 | 理由 |
|------|----------|------|
| vibLevel 超限 fallback 粒度 | 优先只将含超限粒子的 cell 回退到串行碰撞，其余 cell 继续并行；若旧 `noTimeCounter` 尚未拆出 per-cell helper，MVP 可先整段 collision 回退 | 避免一个稀有 species 让整步丧失 OpenMP 收益，同时保留第一版可实现性 |
| reaction fallback 触发点 | `noTimeCounterOMP::collide()` 入口检查 `cloud_.reactions().nReactions() > 0`，若 reaction 未迁移则整体回退旧 `noTimeCounter` | 当前 reaction 是 pair 内部调用，在并行区逐对判断会复杂且容易破坏一致性 |
| RNG 入口 | 在 `dsmcCloud` 增加 `Random& rng(const label tid = -1)`，串行阶段返回 `rndGen_`，并行阶段返回 `threadRng_[tid]` | 避免各模型直接混用 `cloud_.rndGen()` 和 thread RNG |
| fields 内存策略 | `nCells * nThreads * nAccumFields * sizeof(scalar) < 100MB` 时使用 full per-thread fields，否则使用 cell-chunk exclusive 或 sparse events | 明确内存上限，避免大算例因 thread-local field 爆内存 |
| collision backend | Phase 2 同时保留 zero-copy parcel backend 和 packed collision backend 的设计入口 | 通过 profiling 决定拷贝成本和 cache locality 的 trade-off |

7. Phase 0 profiling 后执行优先级复核。当前目标 profiling 算例已知 `moveTime/totalStepTime≈50%`、`collisionTime/totalStepTime≈30%`，因此默认按 move-heavy 路线规划：

- Phase 5 move prototype 在 M2 最小 collision 原型后立即启动，与 Phase 3/4 并行推进。
- Phase 3 fields/tally 仍应实施，但不应阻塞 Phase 5，因为 fields 即使高加速也无法突破 move 串行瓶颈。
- 若 `moveTime / totalStepTime < 15%`，继续按 M0-M4 优先 collision、fields、load balance。
- 若 reaction case 占目标算例主体，Phase 7 reaction delayed event 不能长期后置，否则 Phase 2 collision OMP 对这些算例无收益。
- 若 `buildCollisionData + reconcileCollisionData` 超过串行 collision kernel 时间的 20-30%，Phase 2 优先 zero-copy backend 或并行 copy/reconcile。

验收：

1. 默认配置行为不变。
2. 无 OpenMP 时可编译。
3. 开启 OpenMP 编译但禁用并行路径时结果不变。
4. 每个基准算例有阶段耗时输出。
5. 上述 4 个接口/策略在设计文档和代码入口中都有明确落点。
6. 输出 Amdahl 预测：基于 baseline 阶段占比估算 M0-M4 和完整路线的理论加速上限。

### Phase 1：`parcelPtrs_` 桥接索引、cell index 与 Store shell

目标：先引入轻量桥接结构和 SPARTA 风格 cell 索引，不改变物理算法；同时建立长期 `DsmcParticleStore` 的最小 shell，但不让第一轮依赖完整 Store。

预计工作量：2-3 周

新增文件：

- `particleStore/DsmcParticle.H`
- `particleStore/DsmcParticleStore.H`
- `particleStore/DsmcParticleStore.C`
- `particleStore/DsmcCellIndex.H`

修改文件：

- `clouds/dsmcCloud.H`
- `clouds/dsmcCloud.C`
- `Make/files`

任务：

1. 在 `dsmcCloud` 中增加 MVP 桥接成员：

```cpp
DynamicList<dsmcParcel*> parcelPtrs_;
labelList cellFirst_;
labelList cellNext_;
labelList cellCount_;
PtrList<Random> threadRng_;
```

2. 实现 `buildParcelPtrs()`，从 OpenFOAM 链表生成随机访问指针数组：

```cpp
void dsmcCloud::buildParcelPtrs()
{
    parcelPtrs_.clear();
    parcelPtrs_.reserve(this->size());

    forAllIter(dsmcCloud, *this, iter)
    {
        parcelPtrs_.append(&iter());
    }
}
```

3. 实现 `buildCellIndex()`，以 `parcelPtrs_` 为输入构建 SPARTA 风格 cell 链表：

```cpp
void dsmcCloud::buildCellIndex()
{
    const label nCells = mesh_.nCells();
    const label nPart = parcelPtrs_.size();

    cellFirst_.setSize(nCells);
    cellCount_.setSize(nCells);
    cellNext_.setSize(nPart);

    cellFirst_ = -1;
    cellCount_ = 0;
    cellNext_ = -1;

    for (label i = nPart - 1; i >= 0; --i)
    {
        const label cellI = parcelPtrs_[i]->cell();
        cellNext_[i] = cellFirst_[cellI];
        cellFirst_[cellI] = i;
        ++cellCount_[cellI];
    }
}
```

注意：不要只依赖 `setSize(n, value)` 完成初始化，因为 OpenFOAM List 在 size 未变化时可能不会按预期重置旧值。这里显式执行 `cellFirst_ = -1`、`cellCount_ = 0`、`cellNext_ = -1`，避免上一步残留索引污染新一轮 cell index。

4. 增加 debug 校验工具：

```cpp
void validateParcelPtrs() const;
void validateCellIndexAgainstOccupancy() const;
void dumpCellIndexStats(Ostream& os) const;
```

5. 增加统一 RNG 入口，后续新并行代码禁止直接访问 `cloud_.rndGen()`：

```cpp
Random& dsmcCloud::rng(const label tid = -1)
{
    if (tid < 0)
    {
        return rndGen_;
    }

    return threadRng_[tid];
}
```

说明：

1. 串行阶段和 legacy 模型继续走 `rng(-1)` 或现有 `rndGen()` wrapper。
2. OpenMP parallel region 内必须显式传入 `tid`，从 `threadRng_[tid]` 取随机数。
3. 第一批迁移 collision 相关模型；boundary、controller、reaction 后续按 capability 逐步迁移。

6. `DsmcParticleStore` 第一轮只做 shell：

```cpp
class DsmcParticleStore
{
public:
    void reserve(const label n);
    void clear();
    label size() const;
};
```

完整 `buildFromCloud()`、`syncToCloud()`、`compact()` 后移到 Phase 8 或 move/reaction 需要时再补齐。这样 Phase 1 的成功不依赖一次性复制全部 parcel 状态。

7. 在这些同步点重建桥接索引：

- 初始进入 OpenMP collision/field 前。
- `Cloud<dsmcParcel>::move()` 后。
- `coordSystem().evolve()` 后。
- reaction、controller、boundary 发生 add/delete 后。
- legacy `buildCellOccupancy()` 后的 debug 校验阶段。

验收：

1. `parcelPtrs_.size()` 等于 cloud 粒子数。
2. 每个 `parcelPtrs_[i]` 非空，且 `parcelPtrs_[i]->cell()` 合法。
3. 每 cell 的 `cellCount_` 与 `cellOccupancy_[cell].size()` 一致。
4. `cellFirst_/cellNext_` 遍历得到的总粒子数等于 `parcelPtrs_.size()`，无环、无越界。
5. `rng(-1)` 与旧 `rndGen()` 行为一致，OpenMP 路径能通过 `rng(tid)` 使用线程本地随机数。
6. 不启用后续并行算法时，物理结果不变。

### Phase 2：OpenMP 碰撞原型

目标：实现 `noTimeCounterOMP`，第一版基于 `parcelPtrs_` + `cellFirst_/cellNext_` 做 cell 级并行，并同时保留两种 collision backend 的设计入口：zero-copy parcel backend 和 packed collision backend。不要在未测量拷贝成本前只押注单一 packed `DsmcCollisionData` 路线。

预计工作量：3-5 周

新增文件：

- `collisionPartnerSelection/basic/DsmcCollisionData.H`
- `collisionPartnerSelection/derived/noTimeCounterOMP/noTimeCounterOMP.H`
- `collisionPartnerSelection/derived/noTimeCounterOMP/noTimeCounterOMP.C`

修改文件：

- `collisionPartnerSelection/basic/collisionPartnerSelection.H`
- `collisions/basic/BinaryCollisionModel/BinaryCollisionModel.H`
- `collisions/derived/VariableHardSphere/*`
- `collisions/derived/VariableSoftSphere/*`
- `collisions/derived/LarsenBorgnakkeVariableHardSphere/*`
- `collisions/derived/LarsenBorgnakkeVariableSoftSphere/*`
- `reactions/basic/dsmcReaction/dsmcReaction.H`
- `reactions/derived/*QK/*`
- `clouds/dsmcCloud.H`
- `clouds/dsmcCloud.C`
- `Make/files`

任务：

1. 新增两种 collision backend：

| Backend | 数据路径 | 优点 | 风险 | 推荐使用场景 |
|---------|----------|------|------|--------------|
| zero-copy parcel backend | 直接通过 `parcelPtrs_` 操作 `dsmcParcel&` | 无 `buildCollisionData()`/`reconcileCollisionData()` 拷贝，最快验证并行正确性 | cell 内粒子内存不连续，cache miss 较多，接口与 legacy parcel 耦合更深 | 中小规模、碰撞计算量不高、copy/reconcile 占比高的算例 |
| packed collision backend | 拷贝到 `DsmcCollisionData` 或后续 sorted collision buffer | cache locality 更好，接口干净，利于长期 Store 化 | gather/scatter 成本和 origIdx 映射复杂，若串行 copy 会成为瓶颈 | 大规模、碰撞计算量高、cell locality 收益明显的算例 |

第一版推荐先实现 zero-copy backend，建立 OpenMP cell loop、RNG、统计、fallback 和验证框架；随后实现 packed backend，并用 Phase 0/2 profiling 决定默认路径。

2. packed backend 使用 collision-only 平坦数据结构：

```cpp
enum { MAXVIBMODE_COLLISION = 5 };

struct DsmcCollisionData
{
    point position;
    vector U;
    scalar ERot;
    scalar RWF;
    label ELevel;
    label typeId;
    label cellI;
    label origIdx;
    label vibLevel[MAXVIBMODE_COLLISION];
    label nVibModes;
    label flag;
};
```

说明：

1. `MAXVIBMODE_COLLISION` 只是第一版固定上限，必须运行时检查 `dsmcParcel::vibLevel().size()`。
2. 若某 species 超过上限，优先只将含超限粒子的 cell 标记为 serial fallback，其余 cell 继续并行；若旧 `noTimeCounter` 尚未拆出 per-cell helper，MVP 可以先整段 collision 回退。
3. `position` 用于 subcell 分配；`cellI`、`origIdx` 用于回写与 debug。
4. 第一版只回写碰撞会改变的字段：`U`、`ERot`、`ELevel`、`vibLevel`、可能变化的 `typeId`。position、cell、tet、stepFraction 仍由 move/Cloud 路径维护。

3. packed backend 的 build/reconcile 不能保持串行瓶颈：

```text
parallel buildCollisionData or buildSortedCollisionData
parallel collisionKernel
parallel reconcileCollisionData for dirty particles
```

性能要求：

1. `buildCollisionData()` 必须支持 OpenMP parallel for。
2. `reconcileCollisionData()` 只回写 dirty particles 或碰撞可能改变的字段，避免全字段写回。
3. 若采用 packed backend，优先评估按 cell 顺序写入：

```text
cellStart_[cellI] = writePos
collisionData_[writePos++] = extractFrom(*parcelPtrs_[ip])
origIdx_[writePos] = ip
```

这样 collision loop 可以按 `[cellStart_[cellI], cellStart_[cellI+1])` 连续访问；代价是回写时需要通过 `origIdx_` 间接寻址。

4. 新增 runtime selection model：`collisionPartnerSelectionModel noTimeCounterOMP;`
5. `noTimeCounterOMP::collide()` MVP 流程：

```text
if cloud.reactions().nReactions() > 0 and reactions not OMP-ready:
  fallback to legacy noTimeCounter
cloud.buildParcelPtrs()
cloud.buildCellIndex()
if backend == packed:
  parallel cloud.buildCollisionData() or buildSortedCollisionData()
mark cells requiring serial fallback
parallel over cells
  skip serial-fallback cells
  build local plist
  calculate selected pairs
  choose candidates
  sigmaTcR
  binary collision
merge thread contexts
run serial fallback cells if any
if backend == packed:
  parallel cloud.reconcileCollisionData()
cloud.buildCellOccupancy()         // 兼容 legacy consumers
```

6. 给 `BinaryCollisionModel` 增加并行接口，不破坏原接口。注意现有真实签名为 `collide(dsmcParcel& pP, dsmcParcel& pQ, const label cellI, scalar cR=-1)`，新接口必须显式区分。

zero-copy backend 接口：

```cpp
virtual scalar sigmaTcR
(
    const dsmcParcel& pP,
    const dsmcParcel& pQ
) const;

virtual void collide
(
    dsmcParcel& pP,
    dsmcParcel& pQ,
    const label cellI,
    DsmcThreadContext& ctx
);
```

packed backend 接口：

```cpp
virtual scalar sigmaTcR
(
    const DsmcCollisionData& pP,
    const DsmcCollisionData& pQ
) const;

virtual void collide
(
    DsmcCollisionData& pP,
    DsmcCollisionData& pQ,
    const label cellI,
    DsmcThreadContext& ctx
);
```

7. 第一批优先迁移：

- `VariableHardSphere`
- `VariableSoftSphere`
- `LarsenBorgnakkeVariableHardSphere`
- `LarsenBorgnakkeVariableSoftSphere`

8. reaction 第一版策略：

- 无 reaction 时先完成碰撞并行。
- `noTimeCounterOMP::collide()` 入口检查 `cloud_.reactions().nReactions() > 0`。
- 若 reaction 未迁移且 `nReactions() > 0`，整体回退旧 `noTimeCounter`，不要在 OpenMP pair loop 内逐对调用 `returnModelId()` 后再决定。
- 原因是 reaction 当前在 `noTimeCounter::collide()` 内部直接调用 `cloud_.reactions().reactions()[rMId]->reaction(parcelP, parcelQ)`，会涉及随机数、add/delete、共享计数和粒子拓扑修改。
- 第二步实现 QK reaction 的 delayed event 接口。

9. `sigmaTcRMax_[cellI]` 和 `collisionSelectionRemainder_[cellI]` 按 cell 更新，无需锁。
10. collision/reaction 计数使用 OpenMP reduction 或 thread context 汇总。
11. 每个 OpenMP thread 通过 `cloud.rng(tid)` 使用独立 RNG。第一版接受调度导致的随机序列变化，后续再引入 deterministic hash RNG。
12. Phase 2 即引入最低限度负载统计：`cellWork = nC*(nC - 1)`、`collisionAttempts`、`threadWorkMax/threadWorkAvg`。默认调度可用 `schedule(dynamic, chunk)`，但 chunk 和 work imbalance 必须输出。
13. 使用 thread-local scratch buffers 构建 subcell lists，避免在 cell loop 内频繁分配 `DynamicList` 造成 allocator 竞争。

验收：

1. `OMP_NUM_THREADS=1` 下 noTimeCounterOMP 与 noTimeCounter 尽可能接近，至少统计等价。
2. `OMP_NUM_THREADS=2/4/8` 下 heatBath-5species 宏观量统计等价。
3. 无 reaction case 先通过。
4. reaction case 在未迁移 reaction 时于入口整体回退，已迁移 reaction 才允许进入并行 pair loop。
5. vibLevel 超限 cell 有明确 serial fallback 或整段 fallback 日志，不允许静默截断。
6. zero-copy backend 和 packed backend 至少实现一个可运行路径；若只实现 packed backend，build/reconcile 必须并行化。
7. 输出 collision 性能分解：`buildParcelPtrs`、`buildCellIndex`、`buildCollisionData`、`collisionKernel`、`reconcileCollisionData`、`buildCellOccupancy`。
8. 若 `buildCollisionData + reconcileCollisionData` 超过串行 collision kernel 时间的 20-30%，必须优先切换 zero-copy backend 或优化 packed copy/reconcile。
9. 碰撞 kernel 阶段有可测加速，collision 总阶段加速单独报告。

### Phase 3：并行 fields / tally

目标：将宏观场计算从链表遍历迁移到 `parcelPtrs_`/cell index 或后续 store，并行 reduce。

预计工作量：2-3 周

修改文件：

- `macroscopicProperties/basic/dsmcFieldProperties/*`
- `macroscopicProperties/basic/dsmcStandardFields/*`
- `macroscopicProperties/derived/combined/dsmcVolFields/*`
- `boundaryMeasurements/*`
- `faceTracker/*`
- `cellMeasurements/*`

任务：

1. 给 `dsmcFieldProperties` 增加可选 `parcelPtrs_`/store 计算路径。
2. `dsmcStandardFields::calculateFields()` 改造为：

```text
zero thread-local fields
parallel over non-empty cells or particles
  accumulate rhoN/rhoM/linearKE/rotE/vibE/momentum
reduce thread-local fields into volFields
correctBoundaryConditions()
```

3. fields 阶段通常是 memory-bandwidth bound，调度策略要同时考虑写冲突和读局部性：

- 优先按 non-empty cell 遍历，减少空 cell 开销。
- 如果使用 `cellFirst_/cellNext_`，同 cell 粒子仍可能内存不连续；若 fields 成为瓶颈，应提前引入 `cellStart_/sortedParticleIds_` 或 sorted collision buffer 复用。
- 对 full per-thread fields，reduce 成本是 `O(nThreads * nCells * nAccumFields)`；稀疏网格中应优先 sparse events 或 non-empty cell chunk。

4. 对大 mesh 避免 `nThreads * nCells * nFields` 内存过大：

- 估算 `fieldBufferBytes = nCells * nThreads * nAccumFields * sizeof(scalar)`。
- 当 `fieldBufferBytes < 100MB` 时使用 full per-thread field，最后 reduce 到 volFields。
- 当 `fieldBufferBytes >= 100MB` 时使用 cell-chunk exclusive 写，或使用 sparse thread-local events。
- `100MB` 作为第一版默认阈值，后续可放入 dictionary，例如 `openMPDsmc.fieldBufferLimit 100MB`。

5. face/boundary/porous measurements 使用 event buffer。

验收：

1. fields 与串行统计等价。
2. `fields_.calculateFields()` 阶段有明确加速。
3. 内存策略日志输出 `fieldBufferBytes`、`nAccumFields`、`nNonEmptyCells/nCells` 和所选模式。
4. 内存占用在大算例中可控。

### Phase 4：线程级负载均衡

目标：为 collision、fields、boundary events、move 原型提供统一负载估计和调度。

预计工作量：2 周

新增文件：

- `particleStore/DsmcLoadBalancer.H`
- `particleStore/DsmcLoadBalancer.C`

任务：

1. 定义 per-cell work：

```cpp
struct DsmcCellWork
{
    scalar nParticles;
    scalar collisionAttempts;
    scalar boundaryEvents;
    scalar reactions;
    scalar moveCost;
    scalar collisionCost;
    scalar fieldCost;
    scalar totalCost;
};
```

2. 每步更新指数滑动平均：

```cpp
workAvg[cellI] = alpha*workNow[cellI] + (1-alpha)*workAvg[cellI];
```

3. collision 的默认工作量估计应从 Phase 2 就开始记录：

```cpp
collisionWork[cellI] = cellCount[cellI]*(cellCount[cellI] - 1);
```

NTC candidate 数近似随 `nC^2` 增长，因此调度不能只看 cell 数。`collisionAttempts` 和实测 cell 耗时应逐步替代估算值。

4. 提供三种调度：

| 策略 | 用途 |
|------|------|
| OpenMP dynamic | 默认低风险路径 |
| heavy-first task list | 粒子分布强不均 |
| static by contiguous cell range | fields/cache-friendly |

5. heavy-first task list 的重建频率：

- 不建议每步都完整排序所有 cells，除非 profiling 证明开销可忽略。
- 默认每 `10-50` 步重建一次 task list。
- 中间步使用上一轮 task list，并用 `workAvg` 做指数滑动更新。
- 若 `threadWorkMax/threadWorkAvg` 连续超过阈值，例如 `1.5`，提前重建。

6. 输出线程负载统计：

- 每线程 cell 数
- 每线程粒子数
- 每线程 collision attempts
- 每线程 boundary events
- 每线程耗时

验收：

1. 非均匀算例中线程 idle 时间下降。
2. heavy-first 比 simple dynamic 在热点算例中更稳定。
3. 输出 `threadWorkMax/threadWorkAvg`、task-list rebuild 周期和调度策略。
4. 输出 per-cell work，可用于后续 MPI DLB。

### Phase 5：Move 原型

目标：建立受限的并行 move 原型，先覆盖可以证明安全的简单场景；正式路径仍以 OpenFOAM `dsmcParcel::move()`/`trackToFace()` 语义为基准。

预计工作量：4-8 周

修改文件：

- `parcels/dsmcParcel.C`
- `parcels/dsmcParcel.H`
- `clouds/dsmcCloud.C`
- `boundaries/basic/*`
- `boundaries/derived/patchBoundaries/dsmcSpecularWallPatch/*`
- `boundaries/derived/patchBoundaries/dsmcDiffuseWallPatch/*`
- `boundaries/derived/patchBoundaries/dsmcAbsorbingWallPatch/*`

任务：

1. 第一版不把简化 `pointInCell + neighbourWalk` 作为正式替代 `trackToFace()` 的实现。若做这个方向，只能作为实验分支，必须单独验证 `stepFraction`、tet 状态、face transition、reduced-D 约束和 boundary hit 行为。
2. 优先抽取 `dsmcParcel::move()` 中可复用的几何推进和状态更新边界，明确哪些字段必须继续由 `particle` 基类维护。
3. 定义 `DsmcTrackingState`，只承载并行路径确实需要且可验证的字段，不试图一次替代完整 `dsmcParcel::trackingData`。
4. 使用 Claude 版建议的桥接思路：通过 `parcelPtrs_` 识别可安全并行推进的候选粒子，只有能保证不会触发 patch/cyclic/processor/outlet 回调的粒子才进入并行 move。
5. 并行遍历安全候选粒子：

```text
parallel over particles
  compute track time
  track with validated OpenFOAM-compatible state update
  update x/cell/tet/stepFraction
  if boundary hit: append boundary event
  if processor/outlet: set flag
```

6. 若任一粒子进入未支持路径，应回退到完整串行 move，而不是在并行区留下半更新状态。
7. 第一版支持：

- Cartesian
- no processor migration within rank test
- specular/diffuse/absorbing wall
- no sticking/catalytic/cyclic complex behavior

8. boundary event 合并后再调用或模拟 boundary model。对于会写 measurement、add/delete 或依赖 shared RNG 的 boundary，必须先串行回退。
9. 输出 move prototype 覆盖率和收益诊断：

```text
parallelMoveParticles / totalParticles
parallelMoveTime / totalMoveTime
serialMoveFallbackParticles / totalParticles
serialMoveFallbackTime / totalMoveTime
fallback reasons:
  boundary-risk
  cyclic/processor/outlet
  reduced-D
  unsupported patch model
  tracking-state mismatch
```

10. 若 `parallelMoveParticles / totalParticles < 70%`，需要复核安全粒子分类策略。可选方向包括放宽内部粒子判定、扩大支持的简单 wall patch、或把更多 boundary hit 转为 event buffer，而不是继续只优化已并行的 move kernel。

验收：

1. 无 wall case 轨迹与串行一致或统计等价，并通过 cell/tet/stepFraction 校验。
2. 简单 wall case boundary hit 数一致。
3. move 阶段有可测加速。
4. 不支持模型自动回退串行 move。
5. 对 reduced-D、cyclic、processor patch、sticking/catalytic 等复杂路径有显式禁用或 fallback 日志。
6. 输出 `parallelMoveParticles/totalParticles`、`parallelMoveTime/totalMoveTime`、fallback 粒子比例和 fallback 原因分类。
7. 若可并行粒子比例低于 70%，必须给出下一步提高覆盖率的策略，而不是只报告 move kernel 加速。

### Phase 6：Boundary、controller、coordSystem 迁移

目标：逐步消除 move 前后主要串行热点。

预计工作量：4-8 周

优先迁移：

1. `dsmcCartesian`：基本无额外粒子操作。
2. `dsmcAxisymmetric` 和 `dsmcSpherical`：RWF 更新、clone/delete 延迟提交。
3. 常用 wall patch：

- `dsmcSpecularWallPatch`
- `dsmcDiffuseWallPatch`
- `dsmcAbsorbingWallPatch`
- `dsmcCLLWallPatch`

4. 常用 inflow/outlet：

- `dsmcFreeStreamInflowPatch`
- `dsmcMassFlowRateInlet`
- pressure inlet/outlet 系列

5. 常用 temperature/gravity/forcing controllers。

任务：

1. 给模型增加 capability 标记：

```cpp
enum DsmcOMPCompatibility
{
    SerialOnly,
    ParticleLocal,
    ThreadLocalEvents,
    ParallelReduce
};
```

2. 基类增加默认 `SerialOnly`，逐个模型声明能力。
3. 支持模型走 store 并行路径，不支持模型自动串行回退。
4. add/delete 统一走 event buffer。
5. measurement 写入统一走 reduce。

验收：

1. 常用边界组合可以全并行。
2. 不支持模型不会错误并行执行。
3. coordSystem clone/delete 后 cell index 正确。

### Phase 7：Reaction delayed event 完整化

目标：让主要 QK reactions 在线程安全路径工作。

预计工作量：3-5 周

修改文件：

- `reactions/basic/dsmcReaction/*`
- `reactions/basic/dsmcReactions/*`
- `reactions/derived/dissociationQK/*`
- `reactions/derived/ionisationQK/*`
- `reactions/derived/associativeIonisationQK/*`
- `reactions/derived/exchangeQK/*`
- `reactions/derived/chargeExchangeQK/*`
- `reactions/derived/mixed/*`

任务：

1. 新增 `DsmcReactionContext`：

```cpp
struct DsmcReactionContext
{
    Random& rng;
    DynamicList<DsmcParticle>& newParticles;
    DynamicList<label>& deleteList;
    label cellI;
};
```

2. reaction 不直接调用 `cloud_.addNewParcel()`。
3. reaction 不直接删除 parcel。
4. reaction rates 使用 per-thread counters。
5. mixed reaction 模型后迁移，先确保单一 reaction 模型。

验收：

1. dissociation/ionisation/exchange 单模型统计等价。
2. reaction 创建/删除后 store compact 正确。
3. `reactions_.outputData()` 输出正确。

### Phase 8：Cloud 镜像最小化与 legacy 消费者清理

目标：减少每步 store 与链表之间的双向拷贝。

预计工作量：3-6 周

任务：

1. 统计所有 `cellOccupancy()` 使用点并分类：

| 类别 | 处理方式 |
|------|----------|
| 只读按 cell 遍历 | 替换为 store cell iterator |
| 修改粒子局部字段 | 替换为 `DsmcParticleView` |
| add/delete | event buffer |
| output/IO | 同步到 Cloud |
| WIP/复杂模型 | 保留 legacy 回退 |

2. 提供 adapter：

```cpp
forAllDsmcParticlesInCell(store, cellI, p)
```

3. 将已迁移模型从 `cellOccupancy()` 切到 store iterator。
4. 只在 IO、MPI、legacy model 需要时同步 `Cloud<dsmcParcel>`。

验收：

1. 常用路径每步不再强制完整 syncToCloud。
2. `buildCellOccupancy()` 调用次数减少。
3. 总步进性能进一步提升。

### Phase 9：MPI 级负载均衡接口

目标：为未来跨 rank 动态负载均衡准备数据和接口。

预计工作量：2-4 周

任务：

1. 输出 per-cell work field：

- particle count
- collision attempts
- boundary events
- measured time

2. 与现有 `dynamicLoadBalancing/dsmcDynamicLoadBalancing.C` 对接。
3. 先实现离线建议：生成 decomposition weights。
4. 后续评估在线 cell+particle migration。

验收：

1. 可输出 SPARTA 风格 particle/time weights。
2. 可用于现有重分区流程。

## 6. 具体文件改动清单

### 6.1 第一批必须修改

| 文件 | 改动 |
|------|------|
| `Make/options` | OpenMP 编译/链接选项 |
| `Make/files` | 新增 particleStore 和 noTimeCounterOMP 源文件 |
| `clouds/dsmcCloud.H` | 增加 `parcelPtrs_`、cell index、Store shell、thread contexts、开关、访问器 |
| `clouds/dsmcCloud.C` | 初始化桥接索引、阶段计时、并行路径入口 |
| `collisionPartnerSelection/basic/collisionPartnerSelection.H` | 保持接口，新增并行模型不破坏旧模型 |
| `collisionPartnerSelection/derived/noTimeCounterOMP/*` | 新 OMP 碰撞模型 |
| `collisions/basic/BinaryCollisionModel.H` | 新增 `DsmcCollisionData` 重载或 adapter |
| `collisions/derived/LarsenBorgnakke*` | 第一批并行碰撞实现 |

### 6.2 第二批修改

| 文件 | 改动 |
|------|------|
| `macroscopicProperties/basic/dsmcStandardFields/*` | `parcelPtrs_`/store fields 路径 |
| `macroscopicProperties/basic/dsmcFieldProperties/*` | 调度并行 field |
| `faceTracker/*` | event buffer reduce |
| `boundaryMeasurements/*` | thread-local measurement |
| `cellMeasurements/*` | cell reduce |

### 6.3 第三批修改

| 文件 | 改动 |
|------|------|
| `parcels/dsmcParcel.*` | 抽取 tracking 核心或 adapter |
| `boundaries/basic/*` | OMP compatibility 标记 |
| `boundaries/derived/patchBoundaries/*` | 常用 wall 模型 particle-view/store 接口 |
| `boundaries/derived/generalBoundaries/*` | inlet/outlet event buffer |
| `controllers/basic/*` | OMP compatibility 标记 |
| `controllers/derived/*` | 常用 controller particle-view/store 接口 |
| `coordinateSystem/derived/*` | RWF/clone/delete 并行路径 |

### 6.4 第四批修改

| 文件 | 改动 |
|------|------|
| `reactions/basic/*` | delayed event 接口 |
| `reactions/derived/*QK/*` | 线程安全 reaction |
| `dynamicLoadBalancing/*` | per-cell work 对接 |

## 7. 难点评估

### 7.1 `Cloud::move` 与 tet tracking

难度：高

原因：

1. `dsmcParcel::move()` 使用 OpenFOAM particle tracking，状态分散在 parcel、trackingData、mesh、boundary 中。
2. wall/cyclic/processor patch 处理会调用模型并写共享对象。
3. 新 parcel 的 `stepFraction()` 使用 `rndGen()`。

缓解：

1. 不把 move 作为第一阶段。
2. 先保留串行 move，建立 collision/index/field 并行收益。
3. 再抽取最小 tracking kernel。
4. boundary event 延迟提交。

### 7.2 Reaction 创建/删除粒子

难度：高

原因：

1. reaction 当前直接作用于 `dsmcParcel&`。
2. 可能改变 species、internal energy、创建产品、删除 reactant。
3. reaction rates 和输出存在共享计数。

缓解：

1. 无 reaction case 先通过。
2. reaction 新接口使用 context。
3. new/delete 延迟合并。
4. reaction counters thread-local。

### 7.3 Boundary/controller 数量多

难度：中高

原因：

1. 模型多，行为差异大。
2. 部分模型 WIP 或使用全局状态。
3. 部分模型 add/delete parcel 或写测量。

缓解：

1. capability 分级。
2. 默认 `SerialOnly`。
3. 常用模型优先迁移。
4. 不支持模型显式日志提示并回退。

### 7.4 RNG 可复现性

难度：中

问题：

1. per-thread RNG 会随调度顺序改变随机序列。
2. DSMC 只要求统计一致，但 debug 需要稳定复现。

缓解：

1. 初期 per-thread RNG，接受统计等价。
2. 增加 deterministic RNG mode：

```text
seed = hash(globalSeed, timeStep, cellI, particleI, eventI, streamId)
```

3. 串行路径保留用于对比。

### 7.5 `cellOccupancy_` 消费者迁移

难度：中高

问题：

1. `cellOccupancy()` 在 collision、reaction、boundary、controller、fields、coordSystem 中广泛使用。
2. 一次性替换风险高。
3. 若 Phase 2 后仍每步强制 `buildCellOccupancy()` 供 legacy consumer 使用，它会成为 collision 加速后的兼容层开销。

缓解：

1. `parcelPtrs_`/store 与 `cellOccupancy_` 双轨运行。
2. 每迁移一个模块就增加验证。
3. 提供 cell iterator adapter，降低修改量。
4. Phase 2 验收后立即统计剩余 `cellOccupancy()` consumer；若 Phase 3 fields 已迁移，应评估取消 collision 后的无条件 `buildCellOccupancy()`。

### 7.6 内存占用

难度：中

问题：

1. 早期 store + Cloud 双存储会增加内存。
2. per-thread field buffers 也可能占用大。
3. packed collision backend 的 `DsmcCollisionData` 会引入额外 gather/scatter 数据量；若 build/reconcile 串行，会成为新的 Amdahl 瓶颈。

缓解：

1. store 初期仅在 OpenMP 模式启用。
2. fields 根据 mesh 大小选择 full buffer 或 chunk reduce。
3. packed collision backend 的 build/reconcile 必须支持并行化，并记录 `collisionDataBytesCopied`。
4. 对中小规模或 copy 占比高的算例优先使用 zero-copy backend。
5. 后期减少 syncToCloud 频率。
6. side arrays 按需分配。

### 7.7 单 cell 热点

难度：高但可后置

问题：

1. cell 级并行无法拆分单个超大 cell。
2. cell 内并行会产生粒子写冲突。

缓解：

1. 首选网格/subcell 改善。
2. heavy-first 调度减少普通不均衡。
3. cell 内 candidate graph coloring 作为研究项，不纳入第一轮。

## 8. 验证计划

### 8.1 单元级验证

| 验证项 | 方法 |
|--------|------|
| `parcelPtrs_` build | 粒子数、species count、cell count 对比 |
| cell index | 与 `cellOccupancy_` 每 cell 数量对比 |
| collision data reconcile | 碰撞字段 round-trip 对比 |
| compact | delete/new 后粒子数和 cell index 对比 |
| RNG | thread count 变化下统计稳定性 |

### 8.2 物理回归

| 算例 | 目的 |
|------|------|
| heatBath-5species | 碰撞/反应统计 |
| 2D cylinder | wall boundary、move、load imbalance |
| hypersonicCorner | 非均匀高负载 |
| orion107kmNR | 大规模性能 |

### 8.3 阶段验收

每个 phase 必须满足：

1. 可编译。
2. 默认串行路径不变。
3. 新路径可通过 `OMP_NUM_THREADS=1`。
4. 新路径可通过至少一个多线程算例。
5. 有阶段耗时数据。
6. 若结果不是 bitwise 一致，必须给出统计对比。

补充说明：

1. 碰撞计数在使用 per-thread RNG 和 dynamic schedule 后不要求逐步完全一致，应比较宏观统计量、接受率范围和长期均值。
2. fields/tally 因 reduction 顺序不同，浮点结果不要求 bitwise 一致，但必须在可解释阈值内。
3. `OMP_NUM_THREADS=1` 应作为 debug 近似等价路径，多线程结果以统计等价为标准。
4. 对 add/delete、reaction、boundary event，必须额外检查粒子数、species count、质量/动量/能量统计和 event count。

### 8.4 性能测试矩阵

线程数：

```text
1, 2, 4, 8, 16
```

指标：

1. 总 wall time。
2. move time。
3. cell index time。
4. build/copy time：`buildParcelPtrs`、`buildCellIndex`、`buildCollisionData`。
5. collision kernel time。
6. collision total time：包含 copy/reconcile/compat overhead。
7. reconcile/writeback time。
8. field kernel time。
9. field reduce time。
10. sync/compact/buildCellOccupancy time。
11. 每线程 work imbalance：`threadWorkMax/threadWorkAvg`。
12. fallbackCells、fallbackParticles、reaction fallback 次数。
13. move 覆盖率：`parallelMoveParticles/totalParticles`、`parallelMoveTime/totalMoveTime`、`serialMoveFallbackParticles/totalParticles`。
14. move fallback 原因分类。
15. 内存占用、`fieldBufferBytes`、`collisionDataBytesCopied`。
16. 基于 baseline 阶段占比的 Amdahl speedup prediction。

## 9. 时间与里程碑

| 里程碑 | 内容 | 预计周期 | 风险 |
|--------|------|----------|------|
| M0 | 基线计时与 OpenMP 开关 | 1 周 | 低 |
| M1 | `parcelPtrs_` 桥接 + cell index + Store shell | 2-3 周 | 中 |
| M2 | noTimeCounterOMP 无 reaction | 3-5 周 | 中 |
| M3 | fields/tally 并行 | 2-3 周 | 中 |
| M4 | 线程级负载均衡 | 2 周 | 中 |
| M5 | move 原型 | 4-8 周 | 高 |
| M6 | boundary/controller/coordSystem 迁移 | 4-8 周 | 高 |
| M7 | reaction delayed event | 3-5 周 | 高 |
| M8 | legacy cleanup / sync 最小化 | 3-6 周 | 中高 |
| M9 | MPI DLB 接口 | 2-4 周 | 中 |

建议执行顺序：

```text
M0 -> M1 -> M2 -> M3 -> M4 -> M5 -> M6 -> M7 -> M8 -> M9
```

当前 profiling 算例为 move-heavy，推荐执行顺序调整为：

```text
M0 -> M1 -> M2 -> M5
              \-> M3 -> M4
M5 -> M6 -> M7 -> M8 -> M9
```

解释：

- M2 仍有价值，因为 collision 占比约 30%，且它会建立 `parcelPtrs_`、thread RNG、cell work 和 fallback 框架。
- M5 不应等到 M3/M4 全部完成后再开始，否则总步进会长期被 50% 串行 move 封顶。
- M3/M4 与 M5 并行推进：fields 和 weighted scheduler 提供增量收益，move prototype 解决主要瓶颈。

Phase 0 profiling 后仍允许微调：

- 若 `moveTime/totalStepTime > 25%`，执行 `M0 -> M1 -> M2` 后应将 M5 提前，至少与 M3/M4 并行推进。
- 若目标算例主要启用 reaction，且 reaction 未迁移导致 M2 整体回退，则 M7 应提前到 M3/M4 之前或并行推进。
- 若 `buildCollisionData + reconcileCollisionData` 成为 M2 主瓶颈，应优先做 zero-copy backend 或 packed backend 的并行 copy/reconcile，不应继续优化 collision kernel。

如果项目目标更偏短期性能，可以先执行：

```text
M0 -> M1 -> M2 -> M3 -> M4
```

这条路线可在不触碰高风险 move 的情况下获得碰撞和 field 加速，但对当前 move-heavy profiling 算例，总步进目标应按 1.3-1.7x 评估；若需要超过约 2x，必须推进 Phase 5。

## 10. 决策点

### 10.1 Store 是主结构还是影子结构

推荐决策：短期影子，长期主结构。

原因：

1. 短期先用 `parcelPtrs_` 影子索引能降低风险。
2. 长期若每阶段都 sync，性能会被拷贝吞掉。
3. 完整 OpenMP 必须让 store 成为时间步内主结构。

### 10.2 RNG 策略

推荐决策：

1. Phase 1 增加 `dsmcCloud::rng(label tid = -1)` 统一入口。
2. Phase 2 用 `cloud.rng(tid)` 访问 per-thread RNG，不在 OpenMP 区直接调用 `cloud_.rndGen()`。
3. 串行阶段和 legacy 模型继续使用 `rng(-1)` 或现有 `rndGen()` wrapper。
4. Phase 4 或 Phase 5 增加 deterministic hash RNG。

### 10.3 Reaction 是否进入第一版 OMP collision

推荐决策：

1. 第一版 noTimeCounterOMP 支持无 reaction。
2. `noTimeCounterOMP::collide()` 入口检查 `cloud_.reactions().nReactions() > 0`。
3. 有 reaction 且 reaction 未迁移时整体回退旧 `noTimeCounter`，不要在并行 pair loop 内逐对 fallback。
4. delayed reaction event 在 Phase 7 完整化。

### 10.4 Move 是否立即并行

推荐决策：对当前 profiling 算例，M2 后立即启动 Phase 5 move prototype，并与 Phase 3/4 并行推进。

理由：

1. move 风险最高。
2. collision/index/field 可先贡献可观加速。
3. 但当前目标算例 `move≈50%`，若 move 串行，MVP 总加速上限约 1.3-1.7x。
4. 因此不能长期只优化 collision kernel；M5 至少需要启动受限 prototype，先覆盖内部粒子或简单 boundary 场景。
5. Phase 5 仍不能绕过 `trackToFace()` 语义，必须保留 unsupported path fallback。

### 10.5 Collision backend 选择

推荐决策：Phase 2 同时保留 zero-copy 和 packed 两个后端入口，用 profiling 选择默认路径。

理由：

1. zero-copy backend 避免 `DsmcCollisionData` 的 gather/scatter 成本，适合先验证 OpenMP cell loop 和中小规模算例。
2. packed backend 可能改善 cache locality，但 build/reconcile 必须并行化，最好支持按 cell 排序存储。
3. 若 copy/reconcile 开销超过串行 collision kernel 时间的 20-30%，packed backend 不应作为默认路径。
4. 后续 `cellStart_/sortedParticleIds_` 可同时服务 collision 和 fields，不应无限后置。

## 11. 推荐近期任务

接下来建议按以下顺序推进：

1. 在本计划基础上补充 `cellOccupancy()` 使用点清单，按模块分类。
2. 实施 Phase 0：OpenMP 开关、阶段 wall-time 计时、并行诊断指标和 Amdahl 预测。
3. 实施 Phase 1 MVP：`parcelPtrs_`、显式 reset 的 `buildCellIndex()`、per-thread RNG 和校验工具。
4. 实施 Phase 2 最小版本：优先 zero-copy `parcelPtrs_ + noTimeCounterOMP`，或 packed backend 但必须并行 build/reconcile；无 reaction，`OMP_NUM_THREADS=1/2/4/8` 测试。
5. Phase 2 同时输出 collision 总阶段和 collision kernel 阶段的加速，判断拷贝/兼容层是否成为瓶颈。
6. M2 最小 collision 原型通过后，立即启动 Phase 5 move prototype 设计和风险隔离，不等待 M3/M4 全部完成。
7. 同步推进 Phase 3 fields/tally，并根据 `nNonEmptyCells/nCells` 和 `fieldBufferBytes` 选择 full buffer、chunk 或 sparse events。
8. Phase 4 weighted scheduler 应尽早接入 Phase 2 collision，至少用 `nC*(nC-1)` 和 `collisionAttempts` 作为 work 输入。
9. move prototype 第一目标不是全覆盖，而是量化可安全并行的粒子比例、fallback 比例和 move kernel 潜在加速。

## 12. 成功标准

第一轮成功标准：

1. 串行默认路径完全保留。
2. `parcelPtrs_` 桥接索引与 cell index 稳定通过校验。
3. `noTimeCounterOMP` 在无 reaction case 上统计正确。
4. collision kernel 阶段 8 线程获得 4x 以上加速；collision 总阶段必须单独报告 copy/index/compat overhead 后的实际加速。
5. 对当前 move-heavy profiling 算例，MVP 总步进获得 1.3-1.7x 即可视为 collision/fields 路线阶段性有效；若要超过约 2x，必须有 Phase 5 move 收益。
6. Phase 5 move prototype 至少输出可并行粒子比例、fallback 原因分类、move kernel time 和串行 fallback time。
7. 若 move 仍占总时间超过 40%，后续优先级继续放在 move/boundary，而不是继续微调 collision kernel。

完整成功标准：

1. 常用 DSMC 算例在 move/collision/fields/controller/boundary 阶段均可 OpenMP 并行。
2. 8 线程总步进达到 3-6x 加速。
3. 非均匀算例中线程负载不均显著降低。
4. 常用 reaction 和 boundary 模型支持 delayed event 并统计正确。
5. per-cell work 可输出，并能服务未来 MPI 动态负载均衡。
