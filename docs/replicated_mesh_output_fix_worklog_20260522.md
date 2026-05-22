# Replicated Mesh 输出修复工作日志（2026-05-22）

## 背景

Replicated mesh 模式下所有 MPI rank 持有完整 mesh，不使用 `-parallel` 标志运行。`Pstream::parRun()` 为 false，`reduce()` 为 no-op。每个 rank 只对 `isMyCell(celli)` 的 cell 采样粒子数据。输出时需要 MPI_Allreduce 汇总所有 rank 的累积数据。

## 修复内容

### 1. 派生场 recomputation

**文件**: `src/lagrangian/dsmc/macroscopicProperties/derived/combined/dsmcVolFields/dsmcVolFields.C`

**问题**: MPI_Allreduce 汇总了所有 cumulative 数组（dsmcNCum_, nCum_, mCum_ 等），但派生场（Ttra_, rhoN_, UMean_ 等）是在 reduce 之前从局部数据计算的，导致输出只有 rank 0 拥有的 cell 有正确值。

**修复**: 在 MPI_Allreduce 之后、字段写出之前，添加 recomputation block，从 reduce 后的累积数据重新计算所有派生场：

```cpp
if (cloud_.replicatedMeshActive() && !Pstream::parRun())
{
    const scalar kBLocal = physicoChemical::k.value();
    forAll(dsmcNCum_, celli)
    {
        if (dsmcNCum_[celli] > 1e-3)
        {
            const scalar cellVolume = mesh_.cellVolumes()[celli];
            dsmcNMean_[celli] = dsmcNCum_[celli]/nAvTimeSteps;
            const scalar rhoNMean = nCum_[celli]/(nAvTimeSteps*cellVolume);
            const scalar rhoMMean = mCum_[celli]/(nAvTimeSteps*cellVolume);
            rhoN_[celli] = rhoNMean;
            rhoM_[celli] = rhoMMean;
            UMean_[celli] = momentumCum_[celli]/mCum_[celli];
            const scalar linearKEMean = 0.5*linearKECum_[celli]/(cellVolume*nAvTimeSteps);
            Ttra_[celli] = 2.0/(3.0*kBLocal*rhoNMean)
                *(linearKEMean - 0.5*rhoMMean*(UMean_[celli] & UMean_[celli]));
            p_[celli] = rhoNMean*kBLocal*Ttra_[celli];
        }
        else
        {
            dsmcNMean_[celli] = 0.001;
            rhoN_[celli] = 0.0; rhoM_[celli] = 0.0;
            UMean_[celli] = vector::zero;
            Ttra_[celli] = 0.0; p_[celli] = 0.0;
        }
    }
}
```

### 2. dsmcN_ 瞬时场 MPI_Allreduce

**问题**: `dsmcN_`（瞬时 parcel 计数，volScalarField）每步只对本 rank 的 cell 累加，输出时只有 rank 0 拥有的 ~1/4 cell 有值。

**修复**: 在 MPI_Allreduce block 中添加：
```cpp
MPI_Allreduce(MPI_IN_PLACE, dsmcN_.primitiveFieldRef().data(), nCells, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
```

### 3. Tvib 振动温度修复（isMyCell 过滤）

**问题**: DLB（Phase C auto）每 30 步改变 cellOwner。vibration 累加循环（`dsmcSpeciesEvibModCum_`）没有 `isMyCell` 过滤，而主循环（`dsmcNSpeciesCum_`）有。导致 evib/N ratio 不一致，Tvib 输出异常高（来流处 6000K，应为 ~94K）。

**诊断过程**:
1. 无 DLB 时 ratio=0.5（正确），有 DLB 时 ratio=1.8（错误）
2. PRE_REDUCE 诊断显示：某些 rank 有 evibCum 但 N=0
3. 确认 vibration 循环（line ~2825）遍历所有 cell 无过滤，主循环（line 2581）有 `isMyCell` 过滤
4. `calculateField()` 在 `evolve()` 中每步调用（line 4035），DLB 在同一步稍后执行（line 4153）

**修复**: vibration 循环中添加 `isMyCell` 过滤：
```cpp
if (cloud_.replicatedMeshActive())
{
    const auto& repMesh = cloud_.replicatedMesh();
    forAll(evibCum, celli)
    {
        if (repMesh.isMyCell(celli))
        {
            evibCum[celli] += evibCache[celli];
        }
    }
}
else
{
    forAll(evibCum, celli)
    {
        evibCum[celli] += evibCache[celli];
    }
}
```

### 4. Make/options 链接路径修正

**问题**: 编译链接时使用 `miniconda3/envs/cap/lib` 的 ParMETIS，但运行时 `run.sh` 加载 `parmetis-install/lib`。两者 ABI 不完全兼容（miniconda 版本在 MPI_Allreduce 中报 "Invalid datatype"）。

**修复**:
- `src/lagrangian/dsmc/Make/options`: `-L` 改为 `parmetis-install/lib`
- `applications/utilities/preProcessing/dsmc/dsmcInitialise+/Make/options`: 同上

### 5. 旧库清理

删除 `/home/superxcx/OpenFOAM/superxcx-v2506/platforms/linux64IcxDPInt32Opt/lib/libdsmcFoam+.so`（OpenFOAM bashrc 默认路径下的旧编译残留）。

## 验证结果

### fixoutput case（omp2×mpi4, DLB, 868步, 500步采样）
- Ma_mixture: 60000 nonzero, 0 零 ✓
- Ttra_mixture: 60000 nonzero ✓
- dsmcN_mixture: 0 零 ✓
- Tvib_mixture: 59712 nonzero, 平均 212K ✓

### fixTvib case（omp4×mpi4, DLB, 2604步, 1500步采样）
驻点线温度分布物理合理：
- 来流：Ttra ≈ Trot ≈ 88K, Tvib ≈ 94K
- 激波内：Ttra 急升至 4500K, Trot/Tvib 滞后（非平衡）
- 激波后：Tvib 从 ~100K 缓慢上升至 ~650K（振动松弛）
- 壁面附近：三温度趋向平衡
- 分布平滑，无异常跳变 ✓

### 性能
- 2604 步总运行，rank wall time max/min = 1.008（DLB 有效）
- Phase C auto DLB rebalances = 85
- 无性能回退（isMyCell 判断为内联 bool 比较）

## 关键代码位置

| 修改 | 文件 | 位置 |
|------|------|------|
| MPI_Allreduce block | dsmcVolFields.C | line ~3245 |
| Recomputation block | dsmcVolFields.C | line ~3378 |
| isMyCell vibration fix | dsmcVolFields.C | line ~2825 |
| isOutputRank guard | dsmcVolFields.C | line ~3373 |
| Make/options | src/lagrangian/dsmc/Make/options | line 33 |

## 注意事项

1. `calculateField()` 在 `evolve()` 中每步调用一次（line 4035），在 `info()` 中也调用（每 nTerminalOutputs 步）。
2. DLB 在 `evolve()` 末尾执行（line 4153），在 `calculateField()` 之后。
3. sharedSampleCache 是 static 的，6 个 dsmcVolFields 实例共享，每步只构建一次。
4. 任何新增的 per-cell 累积数组都需要：(1) 采样时加 `isMyCell` 过滤，(2) 输出时加 MPI_Allreduce。
5. Tvib 的统计平滑度取决于采样步数。1500 步给出平滑结果，200 步较粗糙（统计噪声）。

## 编译与运行

```bash
# 编译
source docs/compile.sh

# 运行（参考 case/cylinder_react/mixparallel/finalcheck/omp2_mpi4_replicatedmesh_fixTvib/run.sh）
source /home/superxcx/intel/oneapi/setvars.sh --force 2>/dev/null
source /home/superxcx/code/OpenFoam/OF-2506/OpenFOAM-v2506/etc/bashrc
export WM_PROJECT_USER_DIR=/home/superxcx/code/DSMC/dsmcFoam++/hyStrath_xcx
export FOAM_USER_APPBIN=$WM_PROJECT_USER_DIR/platforms/$WM_OPTIONS/bin
export FOAM_USER_LIBBIN=$WM_PROJECT_USER_DIR/platforms/$WM_OPTIONS/lib
export PATH=$FOAM_USER_APPBIN:$PATH
export LD_LIBRARY_PATH=$FOAM_USER_LIBBIN:/home/superxcx/code/DSMC/dsmcFoam++/parmetis-install/lib:$LD_LIBRARY_PATH
OMP_NUM_THREADS=4 OMP_PROC_BIND=close OMP_PLACES=cores mpirun -n 4 dsmcFoam+
```

---

# OpenMP 碰撞阶段反应线程安全修复（2026-05-23）

## 背景

Ma=11 高超声速反应流 case（`pal-phd-ch3.3.1-cylinder-react/omp8_mpi2`）在 OpenMP 模式下 step 43-45 崩溃（segfault），串行模式正常。Ma=7 case 无此问题。

**关键线索**：两个 case 唯一区别是 Ma=11 有大量解离反应（N2→N+N, O2→O+O），Ma=7 碰撞能量低于解离阈值，几乎不触发反应。

## 根因分析

### 调用链

```
#pragma omp parallel  (noTimeCounter.C:369/403)
  → processCell() lambda
    → cloud_.reactions().reactions()[rMId]->reaction(parcelP, parcelQ)
      → dissociationQK::dissociateParticleByPartner()
        → cloud_.addNewParcel(...)                    (dissociationQK.C:365)
          → Cloud<dsmcParcel>::addParticle(pPtr)      (dsmcCloud.C:3430)
            → this->append(pPtr)                      (Cloud.C:747) ← IDLList 链表操作，无保护
```

### 问题

`noTimeCounter::collide()` 的碰撞循环在 `#pragma omp parallel` 中执行。当解离反应发生时，`dissociationQK` 调用 `cloud_.addNewParcel()` 创建新粒子。`addNewParcel()` 内部调用 `Cloud::addParticle()` → `IDLList::append()`，这是一个链表修改操作，**完全没有线程保护**。

多线程同时 `append()` 导致链表 next/prev 指针损坏。下一步 `Cloud::move()` 遍历链表时跟随损坏指针 → segfault。

### 为什么 Ma=7 不崩溃

Ma=7 碰撞能量低于解离阈值，`addNewParcel` 几乎不被调用，race condition 极难触发。Ma=11 大量解离反应使多线程频繁同时调用 `addNewParcel`，链表损坏概率极高。

### 次要数据竞争

- `nTotDissociationReactions_[nReac]++`（计数器，不会 crash 但结果不准）
- `relax_` 成员变量（多线程共享同一 reaction 对象，影响物理正确性但不 crash）

## 修复方案

### 核心思路：Per-thread deferred buffer

碰撞阶段将新粒子存入线程局部缓冲区，碰撞结束后顺序 flush 到 Cloud。与 `dsmcFreeStreamInflowPatch.C` 中 `GeneratedParcel` 的模式相同。

### 1. 添加碰撞阶段缓冲区（dsmcCloud.H）

```cpp
bool collisionPhaseActive_;
List<DynamicList<dsmcParcel*>> collisionNewParcels_;

void beginCollisionPhase();
void endCollisionPhase();
```

### 2. 实现 begin/end（dsmcCloud.C）

```cpp
void Foam::dsmcCloud::beginCollisionPhase()
{
    if (!openmpEnabled_) return;
    collisionPhaseActive_ = true;
    const label nThreads = max(ompNumThreads_, label(1));
    collisionNewParcels_.setSize(nThreads);
    forAll(collisionNewParcels_, i)
    {
        collisionNewParcels_[i].clear();
    }
}

void Foam::dsmcCloud::endCollisionPhase()
{
    if (!collisionPhaseActive_) return;
    collisionPhaseActive_ = false;

    label totalNew = 0;
    forAll(collisionNewParcels_, threadI)
    {
        totalNew += collisionNewParcels_[threadI].size();
    }
    if (totalNew == 0) return;

    forAll(collisionNewParcels_, threadI)
    {
        forAll(collisionNewParcels_[threadI], i)
        {
            dsmcParcel* pPtr = collisionNewParcels_[threadI][i];
            Cloud<dsmcParcel>::addParticle(pPtr);
        }
        collisionNewParcels_[threadI].clear();
    }

    buildCellOccupancy();  // 重建占位数组，确保新粒子被下一步 move 正确处理
}
```

### 3. 修改 addNewParcel（dsmcCloud.C）

```cpp
void Foam::dsmcCloud::addNewParcel(...)
{
    dsmcParcel* pPtr = new dsmcParcel(...);

    if (collisionPhaseActive_)
    {
        const label threadI = currentThreadId();
        collisionNewParcels_[threadI].append(pPtr);
        return;
    }

    Cloud<dsmcParcel>::addParticle(pPtr);
    recordMoveAppendedParcel(pPtr);
}
```

两个 `addNewParcel` 重载都做了相同修改。

### 4. 在 evolve() 中包裹 collisions()

```cpp
beginCollisionPhase();
collisions();
endCollisionPhase();
```

### 5. 反应计数器加 `#pragma omp atomic`

所有反应类型的计数器增量前加 `#pragma omp atomic`：
- `dissociationQK.C`: `nTotDissociationReactions_[nReac]++`, `nDissociationReactionsPerTimeStep_[nReac]++`
- `ionisationQK.C`: `nTotIonisationReactions_[nReac]++`, `nIonisationReactionsPerTimeStep_[nReac]++`
- `exchangeQK.C`: `nTotExchangeReactions_++`, `nExchangeReactionsPerTimeStep_++`
- `chargeExchangeQK.C`: `nTotChargeExchangeReactions_++`, `nChargeExchangeReactionsPerTimeStep_++`
- `associativeIonisationQK.C`: `nTotAssociativeIonisationReactions_++`, `nAssociativeIonisationReactionsPerTimeStep_++`（两处：forward 和 reverse）

## 验证结果

### Ma=11 react case（omp8×mpi2, DLB, 2061步, 1001步采样）

- **修复前**：step 43-45 segfault（100% 复现）
- **修复后**：2061 步全部完成，正常退出 ✓

运行统计：
- 总粒子数：2,333,037（大量解离产物）
- 碰撞：129,853 collisions/step
- Phase C DLB rebalances：78
- rank wall time max/min = 1.000456（负载均衡良好）
- 总运行时间：2591s（wall clock）

### 驻点线温度分布

- 来流区（x < -1.4）：Ttra ≈ 187-191K，平滑 ✓
- 激波内（x ≈ -1.39 到 -1.30）：波动较大（统计噪声，非 bug）
- 激波后（x > -1.15）：Ttra ≈ 10000-20000K，趋势正确

激波区域波动大是因为 1001 步采样对 Ma=11 强激波不够（每 cell 仅 ~11 粒子，激波位置微振荡）。需要 5000-10000 步采样才能平滑。

### Ma=7 react case（omp8×mpi2, 对照）

正常运行无崩溃（解离反应极少，不触发 race condition）。

## 关键代码位置

| 修改 | 文件 | 说明 |
|------|------|------|
| collisionPhaseActive_ + buffer | dsmcCloud.H | line ~179 |
| beginCollisionPhase/endCollisionPhase | dsmcCloud.C | line ~857 |
| addNewParcel buffer 逻辑 | dsmcCloud.C | line ~3447, ~3469 |
| evolve() 包裹 collisions | dsmcCloud.C | line ~4060 |
| dissociationQK atomic | dissociationQK.C | line ~266 |
| ionisationQK atomic | ionisationQK.C | line ~227 |
| exchangeQK atomic | exchangeQK.C | line ~278 |
| chargeExchangeQK atomic | chargeExchangeQK.C | line ~317 |
| associativeIonisationQK atomic | associativeIonisationQK.C | line ~547, ~636 |

## 性能影响

- `collisionPhaseActive_` 检查：内联 bool 比较，零开销
- Per-thread `DynamicList::append`：无锁，零同步开销
- `endCollisionPhase()` 中 `buildCellOccupancy()`：仅当有新粒子时调用（~11ms/call），对总运行时间影响 < 1%
- `#pragma omp atomic`：单条原子指令，对碰撞热路径影响可忽略

## 已知遗留问题

1. **`relax_` data race**：多线程共享同一 reaction 对象的 `relax_` 成员。Thread A 的 `reaction()` 设 `relax_=false`（解离发生），Thread B 可能同时设 `relax_=true`（新调用开始）。影响：是否执行额外弹性碰撞的判断可能错误。不会 crash，但影响物理正确性。修复需要重构 `reaction()` 接口（返回 relax 值而非写成员变量）。

2. **`useMoveParticlePartition` 仍禁用**：该优化在 react case 中有独立的 bug（粒子 tetFace 无效），与碰撞修复无关。保持 `false` 不影响正确性，仅影响 move 阶段性能。

3. **`useOpenMPInflow` 仍禁用**：inflow 的 OpenMP 并行在 react case 中可能有独立问题，暂保持禁用。

## 后续建议

1. 增加 Ma=11 case 的采样步数（endTime 改为 0.001，有效采样 ~5900 步）
2. 调查 `useMoveParticlePartition` 在 react case 中的 tetFace 无效问题
3. 考虑重构 `relax_` 为返回值模式，彻底消除 reaction 对象的线程安全问题
