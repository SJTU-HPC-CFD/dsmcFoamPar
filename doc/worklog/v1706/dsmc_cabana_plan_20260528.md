# DSMC 粒子存储：自实现 AoSoA 方案深度分析与工作计划

日期：2026-05-28

## 1. 目标

参考 Cabana 的 AoSoA 数据结构设计和内存管理方案，**不引入 Kokkos/Cabana 外部依赖**，在 hyStrath_dlb 中自实现一个轻量级 AoSoA 粒子容器，用于 DSMC 碰撞阶段的 OpenMP 并行化。

## 2. Cabana 核心设计分析

### 2.1 AoSoA 内存布局

Cabana 的核心思想是将粒子数据组织为 **Tile（SoA 块）** 的数组：

```
AoSoA = [ Tile_0 | Tile_1 | Tile_2 | ... | Tile_{N-1} ]

每个 Tile 内部是 SoA 布局（VecLen 个粒子的同一字段连续）：
Tile_k = {
    x[3][VecLen],      // 所有粒子的位置连续
    v[3][VecLen],      // 所有粒子的速度连续
    erot[VecLen],      // 所有粒子的转动能连续
    ...
}
```

关键参数 **VecLen**（Cabana 称为 VectorLength）：
- 必须是 2 的幂（用位运算加速索引）
- CPU/OpenMP：16（匹配 AVX-512 的 16×double 或 AVX2 的 8×double）
- GPU CUDA：32（warp size）
- GPU HIP：64（wavefront size）

### 2.2 二级索引机制

源码：`Cabana/core/src/impl/Cabana_Index.hpp`

```cpp
// 粒子全局索引 i → (tile 索引 s, tile 内偏移 a)
s = i >> log2(VecLen)           // 等价于 i / VecLen
a = i & (VecLen - 1)           // 等价于 i % VecLen

// 反向：(s, a) → i
i = (s << log2(VecLen)) + a    // 等价于 s * VecLen + a
```

全部用位运算，零除法开销。

### 2.3 SoA Tile 内部结构

源码：`Cabana/core/src/Cabana_SoA.hpp`

每个字段在 Tile 内的存储规则：
- 标量字段 `T`：存储为 `T[VecLen]`
- 一维数组 `T[D0]`：存储为 `T[D0][VecLen]`（注意 VecLen 在最内层）
- 二维数组 `T[D0][D1]`：存储为 `T[D0][D1][VecLen]`

**VecLen 在最内层**确保同一字段的连续粒子在内存中相邻 → stride-1 访问 → SIMD 友好。

### 2.4 内存管理策略

源码：`Cabana/core/src/Cabana_AoSoA.hpp`

- **容量以 Tile 为单位**：capacity 总是 VecLen 的整数倍
- **精确分配**：reserve(n) 分配 ceil(n/VecLen) 个 Tile，不做指数增长
- **resize(n)**：先 reserve(n)，再更新 size 和 numSoA
- **shrinkToFit()**：释放多余 Tile，减少内存占用
- **最后一个 Tile 可能不满**：arraySize(lastTile) = size % VecLen

### 2.5 粒子删除：Swap-with-Last 压缩

源码：`Cabana/core/src/Cabana_Remove.hpp`

```
步骤 1（FindEmpty）：parallel_scan 找出需要删除的粒子位置 → indices[]
步骤 2（RemoveEmpty）：parallel_scan 将尾部的存活粒子搬到空洞位置
步骤 3：resize(new_size)，可选 shrinkToFit()
```

不使用懒删除/标记删除。删除后数组紧凑，无空洞。

### 2.6 排序/Binning 机制

源码：`Cabana/core/src/Cabana_Sort.hpp`, `Cabana_LinkedCellList.hpp`

Cabana 的空间排序流程：
1. 计算每个粒子所属 cell（CartesianGrid 映射）
2. 统计每个 cell 的粒子数（parallel count）
3. 前缀和得到 bin offsets
4. 生成 permutation vector（原始索引 → 排序后索引）
5. 两阶段重排：原数据 → scratch buffer → 按新顺序写回

**与 SPARTA 的关键区别**：
- SPARTA 用 first[]/next[] 链表索引，不移动粒子数据
- Cabana 用 permutation + 物理重排，排序后粒子在内存中按 cell 连续

### 2.7 Slice：字段级访问器

Slice 提供对单个字段的视图，跨越所有 Tile：

```cpp
auto velocities = slice<1>(aosoa);  // 获取速度字段的 slice
velocities(i, dim);                  // 访问粒子 i 的速度分量 dim
// 内部：tiles_[i/VecLen].v[dim][i%VecLen]
```

Slice 的价值：碰撞时只需要 v/erot/evib 的 slice，不加载 x/weight 等无关字段。

## 3. 自实现方案设计

### 3.1 设计原则

1. **借鉴 Cabana 的内存布局**（AoSoA + VecLen 在最内层）
2. **借鉴 Cabana 的索引机制**（位运算二级索引）
3. **借鉴 SPARTA 的 cell 索引**（first[]/next[] 而非物理重排）
4. **不用模板元编程**（Cabana 的 MemberTypes/StructMember 太重）
5. **不用 Kokkos**（直接 malloc/realloc + OpenMP）
6. **针对 DSMC 碰撞场景特化**（不做通用粒子框架）

### 3.2 Tile 结构体

```cpp
static constexpr int VECLEN = 8;  // AVX2; 可改为 16 for AVX-512
static constexpr int MAXVIBMODE = 4;

struct DsmcTile
{
    // 碰撞核心字段（碰撞时全部需要）
    double v[3][VECLEN];                    // 192B  速度
    double erot[VECLEN];                    //  64B  转动能
    double evib[VECLEN];                    //  64B  振动能

    // 碰撞辅助字段（碰撞判断/反应需要）
    label  typeId[VECLEN];                  //  64B  组分索引
    label  ELevel[VECLEN];                  //  64B  电子能级
    label  vibLevel[MAXVIBMODE][VECLEN];    // 256B  振动量子数

    // 位置/权重（反应创建新粒子、subCell 划分需要）
    double x[3][VECLEN];                    // 192B  位置
    double RWF[VECLEN];                     //  64B  径向权重因子

    // 管理字段
    label  cellI[VECLEN];                   //  64B  所在单元
    label  origIdx[VECLEN];                 //  64B  原始粒子索引（回写用）
    label  nVibModes[VECLEN];              //  64B  振动模态数
    label  flag[VECLEN];                    //  64B  状态标志
};
// 总计约 1216B per tile = 152B per particle（与 AoS 方案相同）
```

**字段排列顺序的设计意图**：碰撞最频繁访问的 v/erot/evib/typeId 放在 Tile 头部，确保碰撞时加载的第一批 cache line 就是有用数据。

### 3.3 容器类

```cpp
class DsmcParticlePool
{
    DsmcTile* tiles_;           // Tile 数组（连续分配）
    label     nParticles_;      // 实际粒子数
    label     nTiles_;          // 已用 Tile 数 = ceil(nParticles_/VECLEN)
    label     capacity_;        // 容量（Tile 数 × VECLEN）

    // SPARTA 风格 cell 索引（不移动粒子，只建链表）
    labelList cellFirst_;       // cellFirst_[cellI] = cell 中第一个粒子的全局索引
    labelList cellNext_;        // cellNext_[i] = 同 cell 下一个粒子的全局索引

    // 回写映射
    DynamicList<dsmcParcel*> parcelPtrs_;  // parcelPtrs_[i] = IDLList 中对应粒子指针

public:
    // 二级索引访问
    inline double& vx(label i)
    {
        return tiles_[i >> VECLEN_BITS].v[0][i & VECLEN_MASK];
    }
    inline double& vy(label i)
    {
        return tiles_[i >> VECLEN_BITS].v[1][i & VECLEN_MASK];
    }
    inline double& vz(label i)
    {
        return tiles_[i >> VECLEN_BITS].v[2][i & VECLEN_MASK];
    }
    inline double& erot(label i)
    {
        return tiles_[i >> VECLEN_BITS].erot[i & VECLEN_MASK];
    }
    inline label& typeId(label i)
    {
        return tiles_[i >> VECLEN_BITS].typeId[i & VECLEN_MASK];
    }
    // ... 其他字段类似

    // 内存管理
    void reserve(label n);      // 按 Tile 粒度扩容
    void resize(label n);       // 更新 nParticles_/nTiles_
    void buildFromCloud(dsmcCloud& cloud);   // 从 IDLList 构建
    void buildCellIndex();                    // 建立 first[]/next[]
    void writeBackToCloud(dsmcCloud& cloud); // 回写修改
};
```

### 3.4 与 SPARTA first[]/next[] 的结合

**为什么不用 Cabana 的物理重排（permutation）？**

Cabana 排序后粒子按 cell 物理连续，碰撞时可以用 offset+count 直接索引。但物理重排有代价：
- 需要 scratch buffer（额外 N × tupleSize 内存）
- 重排本身是 O(N) 的全量数据搬移
- 每个时间步粒子移动后需要重新排序

SPARTA 的 first[]/next[] 方案：
- 不移动粒子数据，只建立链表索引
- O(N) 建立，O(1) 额外内存（只需 cellFirst + cellNext 两个数组）
- 碰撞时通过链表遍历同 cell 粒子

**混合方案**：用 AoSoA 布局存储数据（获得 cache/SIMD 优势），用 first[]/next[] 做 cell 索引（避免物理重排开销）。

碰撞时的访问模式：
```cpp
// 遍历 cell 中的粒子，构建 plist
label np = 0;
for (label ip = cellFirst_[cellI]; ip >= 0; ip = cellNext_[ip])
{
    plist[np++] = ip;
}

// 随机选两个粒子碰撞
label candidateP = plist[rng.position(0, np-1)];
label candidateQ = plist[rng.position(0, np-1)];

// 访问速度（AoSoA 布局）
double vPx = pool.vx(candidateP);  // tiles_[candidateP/8].v[0][candidateP%8]
double vQx = pool.vx(candidateQ);
```

### 3.5 碰撞阶段的 Cache 行为分析

假设 VecLen=8，一个 cell 有 20 个粒子，排序后连续存储（占 2.5 个 Tile）：

**碰撞需要的字段**：v[3] + erot + evib + typeId = 每粒子 44B

AoSoA 布局下，碰撞字段在 Tile 头部连续：
- v[3][8] = 192B → 3 cache lines
- erot[8] = 64B → 1 cache line
- evib[8] = 64B → 1 cache line
- typeId[8] = 64B → 1 cache line
- **每 Tile 碰撞数据 = 384B = 6 cache lines**

20 个粒子 = 2.5 Tiles → 加载碰撞数据 = 2.5 × 384B = 960B = 15 cache lines

对比：
| 方案 | 20 粒子碰撞数据加载量 | Cache lines |
|------|----------------------|-------------|
| AoS (96B/粒子) | 1920B（含 52B/粒子无用数据） | 30 |
| AoSoA (Tile 头部) | 960B（全部有用） | 15 |
| 理论最优 | 880B（44B × 20） | 14 |

**AoSoA 接近理论最优**，因为碰撞字段在 Tile 头部连续，加载时几乎不带入无关数据。

### 3.6 运动阶段的 SIMD 向量化

```cpp
// AoSoA 布局下，粒子运动天然可向量化
#pragma omp simd
for (label a = 0; a < VECLEN; a++)
{
    tile.x[0][a] += tile.v[0][a] * dt;
    tile.x[1][a] += tile.v[1][a] * dt;
    tile.x[2][a] += tile.v[2][a] * dt;
}
```

x[dim][VECLEN] 和 v[dim][VECLEN] 都是 stride-1 连续内存，编译器可以直接生成 AVX 指令。

### 3.7 内存管理策略

借鉴 Cabana 但简化：

```cpp
void DsmcParticlePool::reserve(label n)
{
    label newNTiles = (n + VECLEN - 1) >> VECLEN_BITS;
    if (newNTiles <= capacity_ >> VECLEN_BITS) return;

    // 增长策略：至少翻倍（Cabana 精确分配，但我们避免频繁 realloc）
    label allocTiles = max(newNTiles, (capacity_ >> VECLEN_BITS) * 2);

    DsmcTile* newTiles = (DsmcTile*)aligned_alloc(64, allocTiles * sizeof(DsmcTile));
    if (tiles_)
    {
        memcpy(newTiles, tiles_, nTiles_ * sizeof(DsmcTile));
        free(tiles_);
    }
    tiles_ = newTiles;
    capacity_ = allocTiles << VECLEN_BITS;
}
```

与 Cabana 的区别：
- Cabana 精确分配（不过度分配）→ 我们用翻倍策略减少 realloc 次数
- Cabana 用 Kokkos View → 我们用 aligned_alloc（64B 对齐，匹配 cache line）
- Cabana 用 deep_copy → 我们用 memcpy（POD 结构体，安全）

### 3.8 粒子删除策略

碰撞中反应可能删除粒子。采用**延迟标记 + 批量压缩**：

```cpp
// 碰撞循环中：标记待删除
pool.flag(ip) = FLAG_DELETE;

// 碰撞结束后：批量压缩（Cabana swap-with-last 思路）
void DsmcParticlePool::compact()
{
    label dst = 0;
    for (label src = 0; src < nParticles_; src++)
    {
        if (flag(src) != FLAG_DELETE)
        {
            if (dst != src) copyParticle(dst, src);
            dst++;
        }
    }
    nParticles_ = dst;
    nTiles_ = (dst + VECLEN - 1) >> VECLEN_BITS;
}
```

注意：compact 后 cellFirst_/cellNext_ 失效，需要重建。但碰撞结束后不再需要 cell 索引，所以无影响。

## 4. 与原 SPARTA AoS 方案的对比

| 维度 | SPARTA AoS | 自实现 AoSoA |
|------|-----------|-------------|
| 碰撞带宽利用率 | 46%（加载 96B 用 44B） | ~92%（碰撞字段在 Tile 头部） |
| 碰撞 cache lines/20粒子 | 30 | 15 |
| 运动阶段 SIMD | 不可能（stride-96B） | 天然支持（stride-1） |
| 容器代码量 | ~50 行 | ~300 行 |
| 访问语法 | `p[i].v[0]` | `pool.vx(i)` |
| 调试 | gdb 直接看 struct | 需算 tile/slot |
| 排序 | swap 整个 struct | 不排序，用 first/next |
| 外部依赖 | 无 | 无 |
| 未来 GPU 适配 | 需要重构布局 | 天然适配（SoA 布局） |

**结论**：自实现 AoSoA 的工程成本比完整 Cabana 低一个数量级（300 行 vs 2000+ 行），但保留了 AoSoA 的核心性能优势（带宽利用率翻倍 + SIMD 向量化）。

## 5. 完整工作计划

### Phase 0：基础设施（1-2 周）

**目标**：OpenMP 编译支持 + DsmcTile/DsmcParticlePool 定义

1. 修改 `src/lagrangian/dsmc/Make/options`：添加 `-fopenmp`
2. 新建 `src/lagrangian/dsmc/particles/DsmcTile.H`：Tile 结构体定义
3. 新建 `src/lagrangian/dsmc/particles/DsmcParticlePool.H`：容器类声明
4. 新建 `src/lagrangian/dsmc/particles/DsmcParticlePool.C`：容器类实现
   - reserve() / resize() / aligned 内存分配
   - 内联访问器（vx, vy, vz, erot, evib, typeId, ...）
   - buildFromCloud()：从 IDLList 拷贝到 AoSoA
   - buildCellIndex()：建立 first[]/next[]
   - writeBackToCloud()：回写修改到 IDLList
   - compact()：删除标记粒子并压缩
5. 在 dsmcCloud 中添加 `PtrList<Random> threadRng_`

**修改文件**：
- `src/lagrangian/dsmc/Make/options`
- `src/lagrangian/dsmc/Make/files`
- 新建 `src/lagrangian/dsmc/particles/DsmcTile.H`
- 新建 `src/lagrangian/dsmc/particles/DsmcParticlePool.H`
- 新建 `src/lagrangian/dsmc/particles/DsmcParticlePool.C`
- `src/lagrangian/dsmc/clouds/dsmcCloud.H`
- `src/lagrangian/dsmc/clouds/dsmcCloud.C`

### Phase 1：集成与验证（2-3 周）

**目标**：将 DsmcParticlePool 集成到 dsmcCloud 的 evolve 流程中，验证正确性

1. 在 `dsmcCloud::evolve_moveAndCollide()` 中：
   - move() 之后调用 `particlePool_.buildFromCloud(*this)`
   - 调用 `particlePool_.buildCellIndex()`
   - 碰撞后调用 `particlePool_.writeBackToCloud(*this)`

2. 验证：对比 buildCellIndex() 结果与 cellOccupancy_ 的每 cell 粒子数

3. 单元测试：
   - 构造已知粒子分布，验证 buildFromCloud → buildCellIndex → writeBack 的往返一致性
   - 验证 compact() 后粒子数据完整性
   - 验证 reserve() 扩容后数据不丢失

**修改文件**：
- `src/lagrangian/dsmc/clouds/dsmcCloud.H`（添加 DsmcParticlePool 成员）
- `src/lagrangian/dsmc/clouds/dsmcCloud.C`（集成到 evolve 流程）

### Phase 2：OpenMP 碰撞循环（3-4 周）

**目标**：实现 `noTimeCounterOMP`，在 AoSoA 上用 OpenMP 并行执行碰撞

1. **新建 noTimeCounterOMP 类**（继承 collisionPartnerSelection）

```cpp
void noTimeCounterOMP::collide()
{
    DsmcParticlePool& pool = cloud_.particlePool();
    pool.buildFromCloud(cloud_);
    pool.buildCellIndex();

    const label nCells = mesh_.nCells();
    label totalCollisions = 0;

    #pragma omp parallel reduction(+:totalCollisions)
    {
        const int tid = omp_get_thread_num();
        Random& rng = cloud_.threadRng(tid);
        DynamicList<label> plist(256);
        List<DynamicList<label>> subCells(8);
        DynamicList<DsmcTileParticle> localNewParticles;

        #pragma omp for schedule(dynamic, 64)
        for (label cellI = 0; cellI < nCells; cellI++)
        {
            // 构建 plist
            plist.clear();
            for (label ip = pool.cellFirst(cellI); ip >= 0;
                 ip = pool.cellNext(ip))
            {
                plist.append(ip);
            }
            const label nC = plist.size();
            if (nC < 2) continue;

            // subCell 划分
            // NTC 碰撞选择
            // 碰撞执行（直接操作 pool.vx(ip) 等）
            // 反应 → localNewParticles / flag(ip) = FLAG_DELETE
        }
    }

    // 串行：处理新粒子和删除
    pool.compact();
    pool.writeBackToCloud(cloud_);
    // 将 localNewParticles 添加到 cloud
}
```

2. **碰撞模型适配**
   - 新增 `BinaryCollisionModel::collide(DsmcParticlePool&, label, label, label, Random&)` 重载
   - 主要适配 `LarsenBorgnakkeVariableHardSphere`
   - 能量分配函数接受 thread RNG

3. **反应模型适配**
   - 反应产生新粒子 → thread-local buffer
   - 反应删除粒子 → pool.flag(ip) = FLAG_DELETE
   - 碰撞结束后串行合并

4. **线程安全保证**
   - 每个 cell 只被一个线程处理 → cell 内粒子数据无竞争
   - sigmaTcRMax_[cellI]：per-cell 写入，无竞争
   - RNG：per-thread 实例
   - 新粒子：thread-local buffer
   - constProps()：只读

5. **运行时选择**
   - `dsmcProperties` 中：`collisionPartnerSelectionModel noTimeCounterOMP;`
   - 保留 noTimeCounter 作为串行回退

**修改文件**：
- 新建 `src/lagrangian/dsmc/collisionPartnerSelection/derived/noTimeCounterOMP/noTimeCounterOMP.H`
- 新建 `src/lagrangian/dsmc/collisionPartnerSelection/derived/noTimeCounterOMP/noTimeCounterOMP.C`
- `src/lagrangian/dsmc/collisions/derived/LarsenBorgnakkeVariableHardSphere/` — 添加重载
- `src/lagrangian/dsmc/reactions/derived/` — 延迟创建接口
- `src/lagrangian/dsmc/Make/files`

### Phase 3（可选）：性能优化（2-3 周）

**目标**：消除双重存储，优化热点

1. **消除 buildFromCloud/writeBack 的拷贝开销**
   - 如果碰撞是每步最大开销，拷贝开销可接受（O(N) 拷贝 vs O(N×碰撞次数) 碰撞）
   - 如果需要消除：将 DsmcParticlePool 作为主存储，Cloud 只在 move/IO 时使用
   - 这需要更大的架构改动，视 Phase 2 的性能测试结果决定

2. **碰撞阶段 SIMD 优化**
   - 对于同 Tile 内的粒子对碰撞，可以利用 SoA 布局做向量化
   - 但碰撞配对是随机的，跨 Tile 概率高 → SIMD 收益有限
   - 重点优化：碰撞后能量分配的数学计算（可向量化的部分）

3. **运动阶段向量化**
   - 粒子运动 `x += v*dt` 天然适合 AoSoA 向量化
   - 但运动阶段在 Cloud::move() 中，与 tet 追踪耦合
   - 仅对自由飞行部分（无边界交叉）可以向量化

### Phase 4（可选）：消除 cellOccupancy_ 双重索引（4-6 周）

与之前 AoS 方案的 Phase 3 相同：将所有 cellOccupancy 消费者迁移到 cellFirst/cellNext。

## 6. 验证策略

1. **Phase 0 单元测试**：
   - 构造 100 个粒子，验证 AoSoA 存储/读取一致性
   - 验证 reserve 扩容后数据完整
   - 验证 compact 后数据正确

2. **Phase 1 集成验证**：
   - buildCellIndex() 后 assert 每 cell 粒子数 == cellOccupancy_[cellI].size()
   - buildFromCloud → writeBack 往返后，所有粒子字段不变

3. **Phase 2 单线程验证**：
   - OMP_NUM_THREADS=1 运行 noTimeCounterOMP
   - 对比 noTimeCounter 的碰撞计数（相同 RNG 种子应完全一致）

4. **Phase 2 多线程统计验证**：
   - heatBath-5species 算例，2/4/8 线程
   - 稳态温度/密度与串行结果在统计误差内一致（< 1%）

5. **性能测试**：
   - hypersonicCorner 或 orion107kmNR 算例
   - 测量碰撞阶段 wall time：1/2/4/8 线程
   - 对比 AoSoA vs AoS 的单线程碰撞性能

## 7. 风险与缓解

| 风险 | 缓解措施 |
|------|----------|
| AoSoA 访问语法比 AoS 复杂 | 提供内联访问器，编译器优化后零开销 |
| 碰撞配对跨 Tile 时 cache 行为 | 排序后同 cell 粒子连续，大部分在同一/相邻 Tile |
| buildFromCloud 拷贝开销 | O(N) 线性拷贝，远小于碰撞 O(N×M) 开销 |
| 反应模型在并行区内创建粒子 | thread-local buffer + 串行合并 |
| vibLevel 超过 MAXVIBMODE=4 | 编译期 static_assert + 运行时截断警告 |
| 调试困难（tile/slot 计算） | 提供 debugPrint(label i) 辅助函数 |

## 8. 工作量估算

| 阶段 | 工作量 | 风险等级 | 核心产出 |
|------|--------|----------|----------|
| Phase 0 | 1-2 周 | 低 | DsmcTile + DsmcParticlePool 容器 |
| Phase 1 | 2-3 周 | 低-中 | 集成到 evolve 流程，验证正确性 |
| Phase 2 | 3-4 周 | 中 | noTimeCounterOMP，OpenMP 碰撞并行 |
| Phase 3 | 2-3 周 | 中 | 性能优化，消除拷贝 |
| Phase 4 | 4-6 周 | 中-高 | 消除 cellOccupancy 双重索引 |
| **总计（Phase 0-2）** | **6-9 周** | | **核心目标完成** |

## 9. 与之前 AoS 方案的关系

本方案**替代**之前的 SPARTA AoS 方案（`dsmc_sparta_plan_20260528.md`）。

核心区别：
- 数据布局从 AoS（struct 数组）改为 AoSoA（Tile 数组，Tile 内 SoA）
- 碰撞带宽利用率从 46% 提升到 ~92%
- 运动阶段获得 SIMD 向量化能力
- 未来 GPU 移植无需重构数据布局
- 工程成本增加约 250 行（容器类），但 OpenMP 碰撞循环的写法几乎不变

OpenMP 并行化的核心设计（schedule(dynamic) + per-cell 分配 + thread-local buffer）完全不变。

## 10. 参考

- Cabana AoSoA 实现：`refcode/Cabana/core/src/Cabana_AoSoA.hpp`
- Cabana 索引机制：`refcode/Cabana/core/src/impl/Cabana_Index.hpp`
- Cabana SoA Tile：`refcode/Cabana/core/src/Cabana_SoA.hpp`
- Cabana 粒子删除：`refcode/Cabana/core/src/Cabana_Remove.hpp`
- Cabana 排序/Binning：`refcode/Cabana/core/src/Cabana_Sort.hpp`
- SPARTA 粒子管理：`refcode/sparta-24Sep2025/src/particle.h`
- SPARTA 碰撞循环：`refcode/sparta-24Sep2025/src/collide.cpp`
- 数据结构选型分析：`doc/worklog/v1706/dsmc_load_balancing_data_structure_20260528.md`

