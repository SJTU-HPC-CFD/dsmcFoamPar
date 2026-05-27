# DSMC 负载平衡：数据结构选型分析

日期：2026-05-28

## 1. 问题背景

hyStrath_dlb 当前的 DSMC 负载平衡方案存在两个层面的问题：

1. **宏观层面**：进程间负载平衡采用 reconstruct → re-decompose → restart 模式，开销巨大
2. **微观层面**：粒子存储采用侵入式双向链表（IDLList），缓存性能差，无法 OpenMP 并行化

本文档记录数据结构选型的分析过程和最终决策。

## 2. 当前 hyStrath 粒子存储分析

### 继承链

```
dsmcParcel → particle → IDLList<particle>::link
```

每个粒子对象内嵌 prev_/next_ 指针（16B），通过 `new` 独立分配在堆上。
Cloud<dsmcParcel> 继承自 IDLList<dsmcParcel>，粒子云本身就是一条双向链表。

### 辅助索引

```cpp
DynamicList<DynamicList<dsmcParcel*>> cellOccupancy_;
```

每个时间步需要完整重建（遍历全部粒子）。

### 链表方案的致命缺陷

- 缓存命中率极低：粒子散布在堆上，顺序遍历每次 `p = p->next_` 几乎必定 cache miss
- 无法 OpenMP 并行：链表迭代器不支持随机访问，无法被 `#pragma omp for` 分割
- 无法 SIMD 向量化：间接寻址阻止编译器自动向量化
- 每粒子 16B 额外开销（prev_/next_ 指针）
- 频繁 new/delete 导致堆碎片化

## 3. 候选方案对比

### 3.1 SPARTA 方案（AoS 连续数组）

Sandia 国家实验室的 DSMC 求解器，生产级代码，万核验证。

```cpp
struct alignas(16) DsmcParticle {
    double x[3];      // 24B  位置
    double v[3];      // 24B  速度
    double erot;      //  8B  转动能
    double evib;      //  8B  振动能
    double dtremain;  //  8B  剩余时间步
    double weight;    //  8B  权重
    int    cellId;    //  4B  所在单元
    int    typeId;    //  4B  组分
    int    flag;      //  4B  状态标志
    int    id;        //  4B  粒子ID
};  // 96B per particle
```

网格-粒子映射：cellStart[] + cellCount[] 数组（排序后）或 first[] + next[] 链表索引。

### 3.2 Cabana 方案（AoSoA）

ECP-CoPA 项目，面向百亿亿次计算的粒子方法框架。

```cpp
using DataTypes = Cabana::MemberTypes<
    double[3],  // position
    double[3],  // velocity
    double,     // erot
    double,     // evib
    double,     // dtremain
    double,     // weight
    int,        // cellId
    int,        // typeId
    int,        // flag
    int         // id
>;
using ParticleList = Cabana::AoSoA<DataTypes, MemorySpace, VecLen>;
```

数据按 VecLen（通常 8 或 16）分块，块内 SoA 布局，块间 AoS 语义。

## 4. 关键性能分析

### 4.1 碰撞阶段（占总时间 40-70%）

碰撞时每个粒子实际需要：v[3](24B) + erot(8B) + evib(8B) + typeId(4B) = 44B
碰撞时不需要：x[3](24B) + dtremain(8B) + weight(8B) + cellId(4B) + flag(4B) + id(4B) = 52B

| 指标 | SPARTA (AoS 96B) | Cabana (AoSoA slices) |
|------|-------------------|----------------------|
| 每 cell 加载量（20粒子） | 1920B (30 cache lines) | 880B (14 cache lines) |
| 碰撞带宽利用率 | 46% | 100% |
| L1 能容纳的 cell 数 | ~16 | ~36 |

但关键事实：碰撞逐 cell 处理，单个 cell 数据（无论 1920B 还是 880B）都远小于 L1 的 32KB。
排序后同 cell 粒子连续，首次加载后全部在 L1，后续随机配对都是 L1 hit。

**碰撞阶段实际性能差异：约 5-15%（Cabana 略优）。**

### 4.2 碰撞的计算密度分析

```
每次碰撞尝试：~40 FLOPs / 88B ≈ 0.45 FLOP/Byte
现代 CPU ridge point：4-8 FLOP/Byte
```

DSMC 碰撞是纯粹的访存受限（memory-bound）问题。SIMD 向量化对碰撞帮助有限。

### 4.3 碰撞的随机访问特性

碰撞配对是随机选取，不是规则的 stride-1 遍历。AoSoA 的向量化优势在随机访问模式下无法发挥。
排序后两者的 cache 行为差异很小（都是 L1 hit）。

### 4.4 粒子运动阶段（占总时间 10-20%）

| 指标 | SPARTA | Cabana |
|------|--------|--------|
| 需要数据 | x[3] + v[3] = 48B | 同 |
| 实际加载 | 96B（全结构体） | 48B（两个 slice） |
| 带宽利用率 | 50% | 100% |
| SIMD 向量化 | 否（stride-96B） | 是（stride-1） |

Cabana 在运动阶段有 2-3x 优势，但运动只占 15% 总时间，整体贡献约 2-5%。

### 4.5 OpenMP 负载平衡适配性

两者完全等价：

```cpp
#pragma omp parallel for schedule(dynamic, chunk)
for (int icell = 0; icell < nCells; icell++)
{
    int start = cellStart[icell];
    int count = cellCount[icell];
    // ... 碰撞逻辑（与数据结构无关）
}
```

数据结构不影响线程间负载分配，只影响单线程内的 cell 处理速度。

## 5. 工程复杂度对比

| 维度 | SPARTA 风格 | Cabana |
|------|------------|--------|
| 外部依赖 | 无（纯 C++ struct） | Kokkos + Cabana |
| 编译时间 | 秒级 | 分钟级（重模板） |
| 调试难度 | gdb 直接看 struct | 模板符号难读 |
| 与 OpenFOAM 集成 | 容易 | 困难（内存管理冲突） |
| 学习曲线 | 几乎为零 | 需理解 Kokkos 执行空间 |
| 改造代码量 | ~500 行 | ~2000+ 行 |
| 改造周期 | 1-2 周 | 1-2 月 |

## 6. 决策

**选择 SPARTA 风格的 AoS 连续数组方案。**

理由：

1. DSMC 碰撞是随机访问模式，AoSoA 的向量化优势无法发挥，两者性能差异 < 15%
2. 工程投入产出比：SPARTA 方案 500 行 / 1-2 周，Cabana 方案 2000+ 行 / 1-2 月
3. OpenMP 负载平衡与数据结构正交，不影响核心目标
4. SPARTA 已在 Sandia 万核级生产环境验证十年
5. 与 OpenFOAM 生态兼容性好，集成改动小
6. 如未来需要 GPU，可在 SPARTA 基础上加 Kokkos 包装，不需要从头用 Cabana

## 7. 粒子结构体设计

```cpp
struct DsmcParticle
{
    double x[3];      // 24B  位置
    double v[3];      // 24B  速度
    double erot;      //  8B  转动能
    double evib;      //  8B  振动能（标量，多模态通过辅助数组）
    double dtremain;  //  8B  剩余时间步
    double weight;    //  8B  权重（RWF）
    int    cellId;    //  4B  所在单元索引
    int    typeId;    //  4B  组分索引
    int    flag;      //  4B  状态标志（自由/吸附/待迁移）
    int    id;        //  4B  粒子ID
};  // 96B，自然 8B 对齐
```

不 padding 到 128B。理由：
- DSMC 瓶颈是粒子数量（百万到亿级），内存带宽比对齐更重要
- 排序后顺序访问占主导，96B 的带宽效率更高
- 20 粒子/cell × 96B = 1920B，仍远小于 L1 cache

### 振动能处理

与 hyStrath 的 labelList vibLevel_（每粒子变长）不同，采用 SPARTA 策略：
- 基础结构体中 evib 为标量（总振动能）
- 多模态振动量子数通过辅助数组存储：`int* vibLevels; // [nParticles × nVibModes]`
- 最大振动模态数编译期确定（MAXVIBMODE = 4，覆盖绝大多数分子）

### 网格-粒子映射

```cpp
labelList cellStart_;   // cellStart_[i] = cell i 的第一个粒子在数组中的位置
labelList cellCount_;   // cellCount_[i] = cell i 中的粒子数
```

通过按 cellId 排序粒子数组实现。排序后同 cell 粒子在内存中连续，碰撞时缓存友好。

## 8. 深度分析：SPARTA vs Cabana 各阶段性能细节

### 8.1 重新审视碰撞阶段的 cache 行为

碰撞时每个粒子实际需要的字段：

```
碰撞需要：v[3](24B) + erot(8B) + evib(8B) + typeId(4B) = 44B
碰撞不需要：x[3](24B) + dtremain(8B) + weight(8B) + cellId(4B) + flag(4B) + id(4B) = 52B
```

碰撞阶段只用到粒子数据的 46%。

假设一个 cell 有 20 个粒子，排序后连续存储：

| 指标 | SPARTA (AoS 96B) | Cabana (AoSoA slices) |
|------|-------------------|----------------------|
| 加载到 cache 的数据 | 20 × 96B = 1920B | 20 × 44B = 880B |
| 实际有用数据 | 20 × 44B = 880B | 20 × 44B = 880B |
| 带宽利用率 | 46% | 100% |
| cache line 数 | 30 | 14 |

Cabana 在碰撞阶段的 cache 占用只有 SPARTA 的一半。同样大小的 L1 cache 能容纳更多 cell 的碰撞数据。

**但这个优势在实际中有多大？** 一个 cell 的数据（无论 1920B 还是 880B）都远小于 L1 的 32KB。碰撞是逐 cell 处理的，两种方案都能把当前 cell 完整放入 L1。**带宽利用率的差异只在 L1 miss 时才体现——即首次加载 cell 数据时。**

### 8.2 碰撞的真实瓶颈：计算还是访存？

DSMC 碰撞的计算密度（arithmetic intensity）：

```
每次碰撞尝试：
  - 读取 2 个粒子的 v, erot, evib, typeId ≈ 88B
  - 计算相对速度、碰撞截面 ≈ 20-30 FLOPs
  - 接受率约 30-50%，接受后修改 v, erot, evib ≈ 50 FLOPs
  
  平均 AI ≈ 40 FLOPs / 88B ≈ 0.45 FLOP/Byte
```

现代 CPU 的 ridge point（roofline 拐点）约 4-8 FLOP/Byte。DSMC 碰撞的 AI = 0.45，**远低于 ridge point，是纯粹的访存受限（memory-bound）问题**。

这意味着：
- SIMD 向量化对碰撞帮助有限（瓶颈不在计算）
- 减少无效数据加载（Cabana 的优势）确实有价值
- 但排序后 L1 命中率已经很高，差异被掩盖

### 8.3 粒子运动阶段的对比

运动阶段：`x[i] += v[i] * dt`，只需 position 和 velocity。

| 指标 | SPARTA | Cabana |
|------|--------|--------|
| 需要的数据 | x[3] + v[3] = 48B/粒子 | 同 |
| 实际加载 | 96B/粒子（全结构体） | 48B/粒子（两个 slice） |
| 带宽利用率 | 50% | 100% |
| 可向量化 | 否（stride-96B） | 是（stride-1 within SoA block） |

粒子运动是顺序遍历所有粒子，这里 Cabana 有明确优势：
- 带宽利用率翻倍
- 内层循环可以 SIMD 向量化（AVX-512 一次处理 8 个粒子的 x += v*dt）

但粒子运动通常只占 DSMC 总时间的 10-20%（大部分时间在碰撞和边界处理）。

### 8.4 宏观量采样阶段

采样需要遍历所有粒子，累加到 cell 场：

```cpp
// 对每个粒子：读 v, erot, evib, typeId, cellId → 累加到 cell 数组
rhoN[cellId] += nParticles;
momentum[cellId] += mass * v;
linearKE[cellId] += 0.5 * mass * (v & v);
```

| 指标 | SPARTA | Cabana |
|------|--------|--------|
| 访问模式 | 顺序遍历，读全结构体 | 顺序遍历多个 slice |
| 向量化 | 困难（scatter to cell arrays） | 同样困难 |
| 排序后 | 同 cell 粒子连续，reduce 可优化 | 同 |

两者差异不大。排序后两者都能利用同 cell 粒子连续的特性做局部 reduce。

### 8.5 OpenMP 负载平衡的适配性

两者对 OpenMP 的适配完全等价：

```cpp
// 两种方案的碰撞循环结构完全相同
#pragma omp parallel for schedule(dynamic, chunk)
for (int icell = 0; icell < nCells; icell++)
{
    int start = cellStart[icell];
    int count = cellCount[icell];
    // ... 碰撞逻辑
}
```

数据结构不影响 OpenMP 负载平衡的实现方式。差异只在每个 cell 内的碰撞计算效率。

### 8.6 工程复杂度对比

| 维度 | SPARTA 风格 | Cabana |
|------|------------|--------|
| 依赖 | 无（纯 C++ struct） | Kokkos + Cabana（两个大型库） |
| 编译 | 秒级 | 分钟级（重模板实例化） |
| 调试 | gdb 直接看 struct 成员 | 模板展开后符号难读 |
| 与 OpenFOAM 集成 | 容易（替换粒子存储即可） | 困难（两套内存管理体系冲突） |
| 学习曲线 | 几乎为零 | 需要理解 Kokkos 执行空间、内存空间、policy |
| 代码行数（改造） | ~500 行 | ~2000+ 行 |
| 维护成本 | 低 | 高（Cabana/Kokkos 版本升级） |

### 8.7 GPU 前景对比

| 指标 | SPARTA | Cabana |
|------|--------|--------|
| GPU 内存合并访问 | 差（AoS 布局，warp 内线程访问不同字段） | 优（SoA 块内 warp 访问同一字段连续地址） |
| GPU 移植难度 | 需要 Kokkos 包装 + 可能重构数据布局 | 原生支持，零额外工作 |
| GPU 碰撞性能 | 需要 AoS→SoA 转换或接受性能损失 | 天然适配 |

如果未来目标包含 GPU，Cabana 的价值显著提升。但 DSMC 的碰撞随机性使得 GPU 并行化本身就很困难（warp divergence、原子操作），不是换数据结构就能解决的。

### 8.8 定量估算：碰撞阶段性能差异

假设条件：
- 100 万粒子，5000 个 cell，平均 200 粒子/cell
- 碰撞占 60% 运行时间
- L1 cache 32KB，L2 256KB
- 排序后同 cell 粒子连续

**首次加载一个 cell 的碰撞数据到 L1：**
- SPARTA：200 × 96B = 19.2KB → 1 次 L2→L1 传输，~10ns
- Cabana：200 × 44B = 8.8KB → 1 次 L2→L1 传输，~5ns

**cell 内随机碰撞配对（数据已在 L1）：**
- SPARTA：`particles[p].v[0]` → L1 hit，1-4 cycles
- Cabana：`velocities(p, 0)` → L1 hit，1-4 cycles
- 差异：几乎为零（都是 L1 命中）

**总碰撞阶段差异估算：**
- 首次加载节省：(19.2 - 8.8) / 19.2 ≈ 54% 带宽节省
- 但首次加载只占碰撞时间的一小部分（大部分时间在计算和 L1 访问）
- 实际碰撞阶段加速：约 5-15%

**粒子运动阶段差异：**
- Cabana 向量化 + 带宽翻倍 → 约 2-3x 加速
- 但运动只占 15% 总时间 → 整体贡献约 15-30% 的 15% = 2-5%

**总体性能差异估算：Cabana 比 SPARTA 快约 5-15%，但工程成本高 4-5 倍。**

### 8.9 最终结论

**对当前目标（CPU + OpenMP 负载平衡），SPARTA 方案更有价值。**

1. DSMC 碰撞是随机访问模式，不是规则的 stride-1 遍历。AoSoA 的向量化优势在碰撞阶段几乎无法发挥。排序后两者的 cache 行为差异很小。

2. 工程投入产出比：SPARTA 风格改造 hyStrath 约 500 行代码，1-2 周完成。Cabana 方案需要引入两个大型依赖库，重写粒子管理的全部代码，1-2 个月。性能差异不超过 15%。

3. OpenMP 负载平衡与数据结构正交：无论用哪种存储，`#pragma omp for schedule(dynamic)` 的效果完全相同。数据结构影响的是每个线程内的单 cell 处理速度，不影响线程间的负载分配。

4. SPARTA 已经在生产环境验证：Sandia 国家实验室用它跑过万核级 DSMC，性能和稳定性经过十年验证。

5. Cabana 的价值在于未来：如果后续需要 GPU 加速或异构计算，Cabana 的 AoSoA + Kokkos 是更好的起点。但作为第一步，先用 SPARTA 风格拿到 OpenMP 的收益，再考虑是否值得为 GPU 重构。

## 9. 下一步工作

1. 实现 DsmcParticle 结构体和粒子数组管理类
2. 实现按 cellId 的基数排序（radix sort）
3. 将碰撞循环改为基于 cellStart/cellCount 的索引访问
4. 添加 OpenMP parallel for schedule(dynamic) 到碰撞循环
5. 处理线程安全问题（RNG、反应产生新粒子）
6. 性能测试对比
