# processor 分布式初始 cloud + OMP 填充（优化 3）— 完整设计 2026-09-17

承接：`dsmcFoam_todolist.md` 初始化优化三项之 3（正解）；
前置实现：Tier-2 过滤读取（`tier2_parallel_filtered_read_20260916.md`）。

## 1. 问题定义（实测数据）

初始化 6.5h 的成本结构（230wcell-1bparticle，job 4641127 实测）：

| 阶段 | 耗时 | 性质 |
|---|---:|---|
| dsmcInitialise+ 串行生成 + 写 17GB | ~6-10 min | ExecutionTime 外，但随规模线性恶化（1B → ~50 min + 150GB） |
| positions 并行解析（Tier-2 ✓） | ~10-15 min | 已优化 |
| **bulk 字段全量读 + ASCII token 解析** | **~2.5-3.5h** | **单核；24 ranks × 17GB = 408GB Lustre 重复读** |
| 字段重映射 / occupancy / 首步 | ~0.5h | |

结构性根因：**全局初始 cloud 文件被 24 个 rank 完整重复读取**（replicated
mesh 架构的固有代价），且 ASCII 格式的 Foam token 解析极慢（CPU≈wall，
实测无 I/O 等待——dsmcInitialise+ 刚写完、page cache 全命中）。

规模外推：1B 粒子 → 字段 ~150GB/rank × 24 = **3.6TB Lustre** + 串行生成
~50 min——**不可运行**。

## 2. 设计目标

1. 初始化（生成→可计算）总耗时：6.5h → **< 5 min**（230w/1.16 亿粒子级）
2. 消除全量文件往返：**零全量读、零（或可选）全量写**
3. 确定性：同 mesh + 同 dsmcInitialiseDict + 同种子 → **逐 cell 初始状态
   与 rank 数/线程数无关**（per-cell RNG 流）
4. 守恒：Σ rank 填充数 == dsmcInitialiseDict 密度积分期望值（按构造精确）
5. 重启语义保持：首写点之后 checkpoint 路径不变；0/ 物化可选
6. 复用已有设施：ctor-2 填充路径、processorN 写格式（checkpoint）、
   Tier-2 的 cellOwner 前移

## 3. 总体架构：两级并行 + 零文件往返

```
MPI 层：24 ranks，每 rank 只填充 cellOwner==myRank 的 cells（÷24）
  └─ OMP 层：rank 内 16 线程按 owned cells 并行采样（÷16）
```

cell 划分完备且不相交（METIS 分区）→ **Σ rank 填充 = 全域填充（按构造
精确守恒，无删除、无迁移）**。

### 两个实施形态

| | 3a. mpirun dsmcInitialise+ 分布式写 | 3b. in-solver 填充（推荐主形态） |
|---|---|---|
| 流程 | 24 ranks 各填自己 owned cells → 各写 processorN/lagrangian → dsmcFoam+ 各读自己的 | dsmcFoam+ 内直接填充（检测无 positions 文件时触发） |
| 文件 I/O | 写 0.7GB/rank + 读 0.7GB/rank | **零** |
| 用途 | 离线预生成（一次生成多次运行）、初始场可检查 | 常态运行（每次运行都重新初始化时） |
| 工程 | dsmcInitialise+ 并行化 + processor 写 | ctor-2 填充逻辑接入 ctor-1 路径 + OMP |
| 依赖 | 3b 的填充内核 | — |

**推荐**：3b 为主（消除文件环节），3a 作为可选副产品（复用同一填充内核
+ processor 写，服务于"生成一次、多案例复用"场景）。

## 4. 详细设计

### 4.1 触发与判定（dsmcFoam+ 构造器 1 路径）

```
if (replicatedMesh
    && startFrom == startTime
    && dsmcInitialiseDict 含 configurations        // 填充配置存在
    && 0/lagrangian/dsmc/positions 不存在           // 无初始 cloud 文件
    && 开关 dsmcInSolverFill（默认 true）)
→ 填充路径（跳过 readCloudFiltered/readFieldsFiltered）
else
→ 现有路径（读文件 / Tier-2 过滤读 / checkpoint 恢复）
```

顺序依赖：**replicatedMesh_->initialize() 必须先于填充**（cellOwner 就绪）
——Tier-2 已完成此前移 ✓。

### 4.2 填充内核（重构 ctor-2 逻辑为可复用方法）

```
dsmcCloud::initialFillFromDict(dsmcInitialiseDict, cellOwner, filterRank)
{
    // 1. 逐配置（configurations 列表顺序与串行版一致）
    for each config in dsmcAllConfigurations:
        // 2. 配置内按 owned cells 并行
        const labelList myCells = ownedCells(cellOwner, filterRank);
        #pragma omp parallel for schedule(dynamic)   // cell 填充成本不均（zone/密度）
        for cellI in myCells:
            // 3. per-cell RNG（确定性核心，见 4.3）
            dsmcFastRng rng = perCellSeed(cellI, configId);
            // 4. 采样并构造（新轻量 ctor 后由配置赋值物理量）
            parcels = config->fillCell(cellI, rng);   // 位置/速度/能级采样
            // 5. tet 定位（cell 内已知 → findTetFacePt，cell-local 便宜）
            for p in parcels: p->initCellFacePtFromCell(cellI);
            keptLocal[thread].append(parcels);
    // 6. 串行 merge（chunk 序）→ buildCellOccupancyFromScratch
}
```

要点：

- **每 rank 只遍历自己的 owned cells**（`myCells_`/owned 列表已有）——
  非 owned cell 连采样都不发生（÷24 的来源）
- 配置内的 cell 遍历 OMP `schedule(dynamic)`（zone/密度不均）
- 粒子 append：线程本地列表 + 串行 merge（与并行读同模式，避免 IDLList
  并发写）
- `addNewParcel` 的 owner 门控在填充路径**不再需要**（只遍历 owned
  cells），但保留作防线
- sigmaTcRMax 初值：**在求解器内由 dsmcInitialiseDict 密度/组分计算**
  （复用 dsmcInitialise+ 写 0/dsmcSigmaTcRMax 的同一公式——该计算代码
  在 utility 路径已存在，提取复用）；**替代 0/ 文件依赖**（3b 下该文件
  自然不存在，现有 FatalError 守卫需加"填充路径豁免"分支）

### 4.3 RNG 确定性设计（关键）

现状：填充用 `cloud_.rndGen_`（单一全局顺序流）——串行确定，但 OMP 化
后**调度顺序改变采样值**。

设计：**per-cell 种子流**

```
seed(cellI, configId) = hash(cellI, configId, globalSeed)
globalSeed = dsmcInitialiseDict.lookupOrDefault<label>("fillSeed", 42)
```

- 使用 `dsmcFastRng`（已有 splitmix64 混合，seed 并发安全——每 cell
  独立实例，无共享状态）
- **强确定性性质**：per-cell 状态只依赖 (cellI, configId, globalSeed)——
  **与 rank 数、线程数、线程调度完全无关**；分解变化时逐 cell 初始状态
  不变（只有归属 rank 变化）
- 全局 cloud.rndGen_ 不再用于填充（保留给碰撞/移动等运行期用途——运行
  期的 RNG 语义不变）
- **初始实现与串行版统计等价但非逐位一致**（不同 RNG 流）——验证标准
  相应调整为统计等价（见 §7）

### 4.4 OMP 填充的线程安全清单

| 共享资源 | 处理 |
|---|---|
| cloud IDLList append | 线程本地 DynamicList + 串行 merge（与并行读同模式） |
| `particleCount_`（origId） | 已原子化（`#pragma omp atomic capture`，本会话已改） |
| `cloud_.rndGen_` | 填充不再使用（per-cell RNG 替代） |
| mesh 几何查询（findTetFacePt/pointInCell） | const 查询 + tetBasePtIs 预热（单线程触发后并行） |
| `cellOccupancy` | 填充后由 buildCellOccupancyFromScratch 统一构建 ✓ |
| 字段（rhoNMean 等） | 填充后由 fields_.createFields 初始化 ✓ |

### 4.5 sigmaTcRMax 与 NTC 冷启动

零场初值是 NTC 的吸收态（碰撞永不发生——sigmaTcRMax 修复的教训）：

- 3b 的初值在求解器内计算：由 dsmcInitialiseDict 的 numberDensities、
  组分质量、Ttra 估计 → **sigmaTcR ≈ 期望相对速率 × σ**（与
  dsmicInitialise+ 写 0/dsmcSigmaTcRMax 的公式一致，代码已存在于
  utility 路径，提取为 dsmcCloud 静态/成员方法）
- 守卫调整：现有 FatalError 守卫增加豁免分支——`dsmcInSolverFill on`
  时不要求文件存在（初值改为内算）
- 运行期包络自增长机制不变（sigmaTcRMax 只增不减 ✓ 非零初值下无吸收
  态风险）

### 4.6 重启与 0/ 物化语义

- **首写点之后**：processorN/<time>/ 完整场+cloud ✓ checkpoint 路径
  原样工作
- **0/ 物化（可选）**：填充完成后按 processor 格式异步写一份
  （服务"中断后从 t≈0 重启"场景）——开关 `dsmcFillWriteInitialCloud`
  （默认 false——1B 级写 150GB 无必要）
- **startFrom latestTime + 无 checkpoint**：现有行为（空 cloud + 守卫
  FATAL）不变——3b 不改变 restart 语义

### 4.7 与 Tier-2 并行读取的关系

- **3b 存在时**：过滤读取路径保留（有初始 cloud 文件的旧 case 仍走
  Tier-2）；新 case 走填充路径（零读）
- **3a 产物**（processorN 分布式初始 cloud）：读取走 checkpoint 恢复
  同款路径（restoreProcessorCheckpoint 泛化），每 rank 只读自己的
  0.7GB——Tier-2 的并行读取在该格式下不必要（文件已是 owned 子集）
- 二者是**写入侧的两种形态**，读取侧统一为"per-rank owned 文件直读"

## 5. 文件改动清单

| 文件 | 改动 |
|---|---|
| `dsmcCloud.H/.C` | 新增 `initialFillFromDict(...)`（填充编排：per-config OMP + merge）；sigmaTcRMax 内算方法；守卫豁免分支 |
| `dsmcAllConfigurations/dsmcConfiguration/dsmcMeshFill/...` | 填充接口重构：接受 (cellI, rng) 的 per-cell 填充（从"遍历全部 cells"改为"填给定 cell"）；线程本地 parcel 列表 |
| `dsmcFastRng.H` | per-cell 种子构造 helper（hash(cellI, configId, globalSeed)） |
| `dsmcFoam+.C` | 填充路径触发判定 + 接线 |
| `dsmcInitialise+.C`（3a） | mpirun 化：owned-cells 填充 + processor 写 |
| sigmaTcRMax 初值计算 | 从 dsmcInitialise+ 路径提取复用 |

## 6. 验证方案

1. **守恒门**：Σ rank 填充数 == dsmcInitialiseDict 密度积分期望
   （按构造精确，容差 0）
2. **确定性门**：同 mesh + 同 dict + 同 seed，两次运行（不同线程数）
   → 逐 cell 初始状态一致（抽查 cellI 的粒子位置/速度/能级逐位一致）
3. **统计等价门**：填充 vs 串行 dsmcInitialise+ 基线（100wcell-5000w
   已有）——iteration 300 的 particles/energy/碰撞数在 RNG 噪声带内
4. **规模门**：超算 230wcell-1bparticle 重跑——初始化（提交→iteration
   10）6.5h → **< 5 min**；1B 粒子级冒烟（内存 ≤ 1TB 总量）
5. **重启门**：填充→首写→checkpoint 重启→继续运行一致
6. **旧路径回归**：开关关闭 + 有初始 cloud 文件的 case（zb、100wcell）
   行为逐位不变

## 7. 风险与缓解

| 风险 | 缓解 |
|---|---|
| per-cell RNG 改变初始实现（与历史基线非逐位一致） | 统计等价验证（§6.3）；确定性性质反而更强（跨 rank/线程不变） |
| zone/密度配置跨 rank 边界的采样偏差 | 按 cell 归属精确划分（不相交完备）+ 守恒门 |
| initCellFacePt 并发的 mesh 懒缓存 | tetBasePtIs 单线程预热（Tier-2 已验证该模式） |
| sigmaTcRMax 初值估计不足 | 运行期包络自增机制（只增不减）+ 与 dsmcInitialise+ 同公式 |
| dsmcInitialiseDict 配置类型覆盖不全（zoneFill/lineFill/laserHeating 的 per-cell 化） | 首版只支持 meshFill（覆盖全部现有 case）；其余配置遇填充路径显式 FATAL 提示走 3a/串行 |
| 填充路径的 0/ 物化缺失影响用户工作流 | 可选异步物化开关（§4.6） |

## 8. 工作量与排期

| 阶段 | 内容 | 工时 |
|---|---|---|
| P1 | 填充内核重构（per-cell 接口 + OMP）+ RNG per-cell 种子 | 2-3 天 |
| P2 | dsmcFoam+ 接线 + sigmaTcRMax 内算 + 守卫调整 | 1 天 |
| P3 | 本地验证（§6.1-6.3、6.6） | 0.5 天 |
| P4 | 超算 230wcell-1b 重跑 + 1B 冒烟（§6.4） | 0.5 天 |
| P5（可选） | 3a 分布式写 + checkpoint 泛化 | 2-3 天 |

## 9. 与优化 1/2 的关系（决策记录）

- 1（字段并行+手动解析）**独立生效**，先做可立即收益（6.5h→~30-40 min），
  且其手动解析/分块技术在 3 的读取路径（3a）中复用
- 2（Tier-2 完全体分块流式）**被 3 作废**（分布式格式下无全量文件可流式）
  ——若 3 排期确定，2 跳过
- 3 落地后：**初始化问题类别消失**（无文件往返），Tier-2 保留用于
  "有全局初始 cloud 文件"的存量 case
