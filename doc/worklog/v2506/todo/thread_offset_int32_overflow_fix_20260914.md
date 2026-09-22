# 线程偏移 int32 溢出修复（800w-rcj 构造期崩溃）— 2026-09-14

## 症状

超算 `ocm-moss3d/800w-rcj`（7,958,528 cells，初始 75,745,943 parcels）在
dsmcCloud 构造期崩溃，新旧二进制（Jul 30 与 Sep 14 重建）相同报错：

```
FOAM FATAL ERROR: bad size -87900444
    From function Foam::List<T>::setSize(signed int)  [List.C:287]
#2 Foam::DynamicList<int>::setCapacity(int)
#3 Foam::dsmcCloud::buildCellOccupancy()
#4 Foam::dsmcCloud::buildCellOccupancyFromScratch()
#5 Foam::dsmcCloud::dsmcCloud(...)
```

作业：4616438 (m36o16)、4616456 (m12o48)、4616685 (m12o48，新二进制)。
全部 rank 同点崩溃；`openmpThreads 48` 在 controlDict 中优先于
`OMP_NUM_THREADS`，故 m36o16 实际也按 48 线程运行。

## 根因

`dsmcCloud.C` buildCellOccupancy() gathered 路径的线程偏移公式：

```cpp
generatedThreadOffsets[threadI] = threadI*gatheredParcels.size()/ompNumThreads_;
```

WM_LABEL_SIZE=32（int32 label）。构造期 replicated mesh 每 rank 持有全部
初始粒子 N=75,745,943，`threadI*N` 在 threadI≥29 时超过 2^31-1 并回绕为负。

精确验算（与报错逐位一致）：

```
off[28] = floor(28*75,745,943/48) = 44,185,133
29*75,745,943 = 2,196,632,347 ≥ 2^31 → int32 回绕 = -2,098,334,949
off[29] = -2,098,334,949/48     = -43,715,311   (C++ 向零截断)
span[28] = off[29] - off[28]    = -87,900,444   ← bad size
min(nCells, span) 取负 → localActiveCells.setCapacity(-87900444) → FATAL
```

溢出阈值：threads > 2^31/N ≈ 28.35 且构造期全量 cloud。该 case 同时满足
（48 线程、75.7M 初始粒子）；130wcell-5000wparticle（50M×48）未触发是因
29×50M < 2^31。

## 修复（8 处同款 `threadI*N/threads` 模式，统一 64 位中间量）

| 文件                 | 位置                                               | 说明                                             |
| -------------------- | -------------------------------------------------- | ------------------------------------------------ |
| dsmcCloud.C          | buildCellOccupancy appendedThreadOffsets           | append 区间                                      |
| dsmcCloud.C          | buildCellOccupancy generatedThreadOffsets          | **崩溃点**                                 |
| dsmcCloud.C          | rebuildMoveOrderedParcels (`threadI*i/nThreads`) | `threadI<=nThreads` 连末位偏移也错（静默损坏） |
| dsmcCloud.C          | 两处 moveOrderedParcels_.size() 偏移块             | appendBatch/setMoveOrdered 系列                  |
| dsmcCloudI.H         | infoMeasurements 逐粒子回退 begin/end              | 含`(threadI+1)*N` 变体                         |
| dsmcReplicatedMesh.C | 迁移 owner 分类 + pack 两处 begin/end 回退         | 含`(threadI+1)*N` 变体                         |

改法示例：

```cpp
generatedThreadOffsets[threadI] = label
(
    (static_cast<long long>(threadI)
   * static_cast<long long>(gatheredParcels.size()))
  / static_cast<long long>(ompNumThreads_)
);
```

行为等价性：非溢出输入下 64 位与 32 位整算术逐位相同，小 case 无行为变化。

## 验证

- 本地编译通过（linux64IccDPInt32Opt，dsmcFoam+ / dsmcInitialise+）
- 本地回归 `zb-cylinder-react-validate/mpi4omp2`（mpirun -np 4 × OMP 2，
  300 步）：`Total Iterations = 300`、`End main`、stuck=0、末态
  particles=2,010,878、energy=1.5698e-3，日志
  `log.threadoffset-fix-20260914`

## 待办

- 同步超算 dlb-new 重建后重跑 4616685（m12o48）；建议同时验证 m36o16
- 该修复与 8/21 chunk transfer（`large_case_chunk_transfer_20260821.md`）
  同属 Int32 溢出家族；800w 级 case 上线前建议再做一轮大 case 冒烟清单
  （构造、首步迁移、首次输出、首次 DLB）
