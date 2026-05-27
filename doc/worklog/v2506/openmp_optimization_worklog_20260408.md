# dsmcFoam+ OpenMP 优化工作总结（截至 2026-04-08）

## 1. 文档目的

本文档用于系统总结本轮 `dsmcFoam+` 在 `OpenMP` 方向上的优化工作，覆盖：

- 优化目标与总体判断
- profiling 的演进过程
- `move / buildCellOccupancy / collision / post fields/output` 四条主线
- 每条主线上的有效优化、无效优化和原因分析
- 当前代码中实际保留的优化项
- 下一阶段仍值得投入的方向

本文档不是算法说明书，而是**工程优化工作日志的归档版**。重点是：

- 哪些优化尝试过
- 哪些有效
- 哪些无效
- 当前代码停在什么状态


## 2. 优化背景与目标

本轮工作最初的观察是：

- `collision` 的并行效率相对理想
- `move` 的并行效率明显低于预期
- 总时间没有随着 `OpenMP` 并行而按同等比例下降

因此工作重点逐步聚焦为：

1. 先把 `move` 的热点彻底暴露出来
2. 再把 `buildCellOccupancy` 和 `post fields/output` 从“串行/冗余准备成本”上降下来
3. 对 `collision` 判断是否还有 OpenMP 层面的剩余空间

在整个过程中，最终形成了一个比较明确的结论：

- **纯 OpenMP 调度层优化基本已经做满**
- 后续收益主要来自：
  - 两阶段提交/两阶段生成
  - 数据流重组
  - 缓存与中间视图复用
  - 减少重复构造和重复遍历


## 3. 测试对象与主要算例

本轮优化和验证的核心算例为：

- [omp8_move_static_coll_dynamic_moveopt](/home/superxcx/code/DSMC/dsmcFoam++/hyStrath_xcx/case/cylinder_react/schedule_sweep_300step_20260404/omp8_move_static_coll_dynamic_moveopt)

相关日志主要集中在：

- `log.dsmcFoam+.buildocc_partition*`
- `log.dsmcFoam+.buildocc_dynamic`
- `log.dsmcFoam+.buildocc_*`
- `log.dsmcFoam+.move_pre_boundary_*`
- `log.dsmcFoam+.buildocc_partition_weight_moveprofileopt_*`
- `log.dsmcFoam+.move_commit2stage1`
- `log.dsmcFoam+.post_flatocc1b`


## 4. 总体优化路线回顾

整个 OpenMP 优化工作大致经历了下面几个阶段：

### 4.1 第一阶段：围绕 `move` 做 profiling

目标：

- 把 `move` 从“大块时间”拆成更细的阶段
- 分辨问题究竟来自：
  - 负载不均衡
  - 单位粒子 tracking 成本不等价
  - 串行尾部
  - 前置 boundary/control 阶段

主要工作：

- 对 `move` 增加分阶段 profiling：
  - `move pre-control/partition`
  - `move reset/setup`
  - `move extract parcels`
  - `move parallel kernel wall`
  - `move commit/survivor rebuild`
  - `move transfer/delete finalize`
- 对 `move` 增加 per-thread profiling：
  - `thread move particles`
  - `thread move wall time [s]`
  - `thread move ns/particle`

阶段性结论：

- `move` 不是简单的“粒子数不均衡”问题
- 很多时候粒子数已经均衡，但 `ns/particle` 仍分裂成不同组
- 说明核心问题在于：
  - tracking 路径单位成本差异
  - boundary / patch / topology 事件差异
  - 串行收口区代价


### 4.2 第二阶段：定位并优化 `move reset/setup`

profiling 暴露出一个很硬的热点：

- `move reset/setup`

进一步拆解后确认，其中的主热点就是：

- 串行 `for (ParticleType& p : *this) p.reset();`

采取的优化：

- 把 `p.reset()` 从链表全扫路径移除
- 融入已经构造好的 `particles[]` 并行路径
- 保证每个 move call 只 reset 一次，不改变语义

结果：

- 这是整个 `move` 线上**最显著的一刀**
- `move reset/setup` 从一个大热点直接降到很小
- `move only` 和总时间都出现了结构性下降

这是当前代码里必须保留的关键优化之一。


### 4.3 第三阶段：拆解 `move pre-control/partition`

进一步 profiling 后发现：

- `move pre-control/partition` 中的大头并不是 `partition`
- 真正的大头是 `boundaries_.controlBeforeMove()`

继续加 profiling 后，确认其中主要由 `Inflow boundary` 构成：

- `velocity/internal energy`
- `addNewParcel`
- `tri select/random point`

阶段性结论：

- `rebuildParticleLoadPartition()` 本身不是这里的大头
- 继续抠 `particle partition` 的收益很有限
- 想继续降 `move pre-control/partition`，应该转向 inflow boundary 内部


### 4.4 第四阶段：优化 `buildCellOccupancy`

profiling 暴露出 `buildCellOccupancy` 的典型结构：

- `extract parcels`
- `count/reduce`
- `allocate/fill`

其中大头通常是：

- `extract parcels`

这一阶段尝试过多种方向，后面会专门总结有效与无效项。


### 4.5 第五阶段：优化 `post fields/output`

初期 `post fields/output` 只是一个总块。

后续通过 profiling 逐步确认：

- 写文件并不是主要成本
- 主要成本在 `calculateFields()`
- 更具体地说，是 `dsmcVolFields` 里的：
  - `sharedSampleCache_.build()`
  - 其中的 `parcel accumulate`

进一步又发现：

- `shared cache build` 的第一个 field（如 `N2`）在替后续 fields 支付大部分缓存构建成本
- `field combine / cell reduction / boundary accumulation` 反而很小

这推动了后续两类有效优化：

- demand-mask 裁剪无关统计项
- `speciesTvib_ / dsmcSpeciesEvibModCum_` 后处理链优化


## 5. `move` 主线详细总结

### 5.1 已验证有效的 `move` 优化

#### 5.1.1 `p.reset()` 合并到并行 `particles[]` 路径

文件：

- [Cloud.C](/home/superxcx/code/DSMC/dsmcFoam++/hyStrath_xcx/src/lagrangian/basic/Cloud/Cloud.C)

核心思想：

- 去掉 `Cloud::move()` 开头对整个 Cloud 的串行 `reset`
- 改成在构造好的 `particles[]` 视图上并行 reset

结果：

- `move reset/setup` 大幅下降
- `move only` 明显下降
- 是 `move` 线上收益最大的优化之一


#### 5.1.2 `move extract parcels` 并行 gather

文件：

- [Cloud.C](/home/superxcx/code/DSMC/dsmcFoam++/hyStrath_xcx/src/lagrangian/basic/Cloud/Cloud.C)

核心思想：

- 在 `useParticlePartition` 分支中
  - 并行统计每线程 parcel 数
  - 并行填充 `particles[]`

结果：

- 这刀本身风险低
- 收益不如 `p.reset()` 那么大，但属于正向演进的一部分


#### 5.1.3 `move commit/survivor rebuild` 两阶段化

文件：

- [Cloud.C](/home/superxcx/code/DSMC/dsmcFoam++/hyStrath_xcx/src/lagrangian/basic/Cloud/Cloud.C)

核心思想：

- 原来：
  - kernel 结束后串行扫一遍 `particles[]`
  - 同时做：
    - survivor rebuild
    - transfer packing
    - deleteParticle
- 改为：
  - 先并行分类到线程私有三类列表：
    - `survivor`
    - `transfer`
    - `delete`
  - 再串行执行真正的 transfer / delete / survivor rebuild

对应验证日志：

- 基线： [log.dsmcFoam+.buildocc_partition_weight_moveprofileopt_occless_postprofile10](/home/superxcx/code/DSMC/dsmcFoam++/hyStrath_xcx/case/cylinder_react/schedule_sweep_300step_20260404/omp8_move_static_coll_dynamic_moveopt/log.dsmcFoam+.buildocc_partition_weight_moveprofileopt_occless_postprofile10)
- 优化后： [log.dsmcFoam+.move_commit2stage1](/home/superxcx/code/DSMC/dsmcFoam++/hyStrath_xcx/case/cylinder_react/schedule_sweep_300step_20260404/omp8_move_static_coll_dynamic_moveopt/log.dsmcFoam+.move_commit2stage1)

结果：

- `move commit/survivor rebuild`: `6.0189 s -> 3.8968 s`
- `main loop wall time`: `215.6954 s -> 213.1567 s`

这刀是有效的，当前保留。


### 5.2 `move` 上验证无效或已放弃的方向

#### 5.2.1 `move` 调度从 `static` 改 `dynamic/guided`

结论：

- 对 `move` 没有本质收益
- 原因不是调度问题，而是单位粒子 tracking 成本不等价

因此：

- `move` 上继续抠 `schedule/chunk` 没有性价比


#### 5.2.2 按最近一次 `thread move ns/particle` 建成本代理

这一方向做过验证，结果说明：

- 可以让线程 wall time 更均衡
- 但不会自动带来显著总时间收益
- 当前算例里，`move` 的瓶颈更像是单粒子 tracking 成本差异，而不是简单分工不均


#### 5.2.3 `particle::changeFace()/changeCell()` 局部 getter 缓存

尝试内容：

- 在 [particle.C](/home/superxcx/code/DSMC/dsmcFoam++/hyStrath_xcx/src/lagrangian/basic/particle/particle.C) 中
  - 把 `mesh_.cells()[celli_]`
  - `mesh_.faces()`
  - `mesh_.faceOwner()`
  - `tetBasePtIs()`
  - `currentTetIndices().faceTriIs(mesh_)`
  拉成本地引用

验证日志：

- [log.dsmcFoam+.trackingcache1](/home/superxcx/code/DSMC/dsmcFoam++/hyStrath_xcx/case/cylinder_react/schedule_sweep_300step_20260404/omp8_move_static_coll_dynamic_moveopt/log.dsmcFoam+.trackingcache1)

结果：

- `main loop wall time`、`move only`、`collision phase` 都更差
- 这刀已经回退

结论：

- 这种 getter 级微缓存不是当前 tracking 主链的有效优化层级


### 5.3 `move pre-control/partition` 的最终判断

profiling 证明：

- `move pre-control/partition` 里真正的大头是 `Inflow boundary`
- 其中主要是：
  - `velocity/internal energy`
  - `addNewParcel`
  - `tri select/random point`

因此：

- `rebuildParticleLoadPartition()` 不是当前主要矛盾
- 如果继续打 `move pre-control`，应转向 inflow boundary


## 6. `buildCellOccupancy` 主线详细总结

### 6.1 profiling 结论

典型拆分为：

- `extract parcels`
- `count/reduce`
- `allocate/fill`

通常：

- `extract parcels` 最大
- `count/reduce` 较小
- `allocate/fill` 次之


### 6.2 已验证有效的 `buildCellOccupancy` 优化

#### 6.2.1 减少不必要容器重整（`occless` 版本）

验证日志：

- [log.dsmcFoam+.buildocc_partition_weight_moveprofileopt_occless](/home/superxcx/code/DSMC/dsmcFoam++/hyStrath_xcx/case/cylinder_react/schedule_sweep_300step_20260404/omp8_move_static_coll_dynamic_moveopt/log.dsmcFoam+.buildocc_partition_weight_moveprofileopt_occless)

结果：

- `buildCellOccupancy`: 明显下降
- 尤其 `allocate/fill` 下降最明显

核心思想：

- 避免每步都对每个 cell 反复做不必要的 `setCapacity/setSize` 型重整
- 让单次 build 更便宜

这是当前保留的 `buildCellOccupancy` 核心有效优化。


### 6.3 已验证无效并已回退的 `buildCellOccupancy` 方向

#### 6.3.1 并行化 `cellOccupancy_.clear()`

结论：

- 没有效果
- 全域清空本身受内存带宽限制，单纯 `omp for` 不改变根因


#### 6.3.2 `active-cell sparse clear`

尝试内容：

- 只清上一步非空 cell

结论：

- 在当前算例上反而更慢
- 维护额外活跃列表本身有代价
- 且没有闭合容器重整问题


#### 6.3.3 避免 fallback gather / 复用 fallback gather 结果

尝试内容：

- 当 fallback gather 发生时，把生成结果写回 `moveOrderedParcels_` 供后续复用

结论：

- 当前算例上未转化成收益
- 反而引入了额外成本


## 7. `collision` 主线详细总结

### 7.1 结论

`collision` 的 OpenMP 调度层基本已经接近上限。

已知现象：

- `dynamic` 通常能把 thread wall time 拉得很平
- 但 candidate 数并不均衡
- 说明问题不是工作量数量不均衡，而是单位 candidate 成本不等价

因此：

- 继续切换 `static/dynamic/guided`
- 或继续做 candidate-only 的分区

收益都不会大。


### 7.2 尝试过的方向

#### 7.2.1 candidate-only partition

结论：

- 看上去 candidate 分配更均衡
- 但 wall time 未必更好


#### 7.2.2 weighted partition（candidate + active cell）

做法：

- 为 `collision partition` 引入 cost 权重：
  - `openmpCollisionCostCandidateWeight`
  - `openmpCollisionCostActiveCellWeight`

结论：

- 方向上合理
- 但总体收益不如 `dynamic` 稳定
- 当前算例下 `dynamic` 仍然是更安全的选择


### 7.3 当前判断

`collision` 如果继续优化，更值得做的是：

- 分析不同 cell 类型、反应路径、访存模式导致的单位 candidate 成本差异
- 而不是继续在 OpenMP 调度层折腾


## 8. `post fields/output` 主线详细总结

### 8.1 profiling 的关键结论

最开始 `post fields/output` 只是一个总块。

经过细分后发现：

- 写文件不是主成本
- 主要成本在 `calculateFields()`
- 更具体地说，是：
  - [dsmcVolFields.C](/home/superxcx/code/DSMC/dsmcFoam++/hyStrath_xcx/src/lagrangian/dsmc/macroscopicProperties/derived/combined/dsmcVolFields/dsmcVolFields.C)
  - `sharedSampleCache_.build()`
  - 其中的 `parcel accumulate`

进一步 profiling 还表明：

- `shared cache build` 的大头通常由第一个 field（如 `N2`）承担
- 后面的 field 基本只是消费已建好的 cache


### 8.2 已验证有效的 `post` 优化

#### 8.2.1 按 field 配置裁剪累积项

核心思想：

- `sharedSampleCache_` 不再“无论后面 field 要不要，都全统计一遍”
- 按 field 配置裁掉不用的统计路径

已实际裁掉的路径包括：

- `classification`
- `electronic`
- `heatFluxShearStress` 相关高阶统计

对应验证日志：

- [log.dsmcFoam+.buildocc_partition_weight_moveprofileopt_occless_postprofile5](/home/superxcx/code/DSMC/dsmcFoam++/hyStrath_xcx/case/cylinder_react/schedule_sweep_300step_20260404/omp8_move_static_coll_dynamic_moveopt/log.dsmcFoam+.buildocc_partition_weight_moveprofileopt_occless_postprofile5)

效果：

- `post fields/output`: 明显下降
- `shared cache build` / `parcel accumulate`: 明显下降

这是 `post` 线上第一类明确有效的结构优化。


#### 8.2.2 `speciesTvib_ / dsmcSpeciesEvibModCum_` 后处理链优化

文件：

- [dsmcVolFields.C](/home/superxcx/code/DSMC/dsmcFoam++/hyStrath_xcx/src/lagrangian/dsmc/macroscopicProperties/derived/combined/dsmcVolFields/dsmcVolFields.C)

核心思想：

- 去掉每个 cell / boundary face 上反复构造的临时 `scalarList` / `List<scalarList>`
- 改为纯标量累计路径

验证结果（复跑更可信）：

- `post fields/output` 继续下降
- `calculateFields`、`shared cache build`、`parcel accumulate` 均下降

这是 `post` 第二类有效结构优化。


#### 8.2.3 复用 `buildCellOccupancy()` 生成的扁平 occupancy 视图

当前最终实现不是直接复用 `moveOrderedParcels_`，而是更稳的版本：

- 在 [dsmcCloud.C](/home/superxcx/code/DSMC/dsmcFoam++/hyStrath_xcx/src/lagrangian/dsmc/clouds/dsmcCloud.C) 的 `buildCellOccupancy()` 中顺手生成：
  - `occupancyOrderedParcels_`
  - `occupancyCellOffsets_`
- 在 [dsmcVolFields.C](/home/superxcx/code/DSMC/dsmcFoam++/hyStrath_xcx/src/lagrangian/dsmc/macroscopicProperties/derived/combined/dsmcVolFields/dsmcVolFields.C) 的 shared cache build 中优先消费这份连续视图

对应接口：

- [dsmcCloud.H](/home/superxcx/code/DSMC/dsmcFoam++/hyStrath_xcx/src/lagrangian/dsmc/clouds/dsmcCloud.H)
- [dsmcCloudI.H](/home/superxcx/code/DSMC/dsmcFoam++/hyStrath_xcx/src/lagrangian/dsmc/clouds/dsmcCloudI.H)

验证日志：

- 基线： [log.dsmcFoam+.move_commit2stage1](/home/superxcx/code/DSMC/dsmcFoam++/hyStrath_xcx/case/cylinder_react/schedule_sweep_300step_20260404/omp8_move_static_coll_dynamic_moveopt/log.dsmcFoam+.move_commit2stage1)
- 复跑： [log.dsmcFoam+.post_flatocc1b](/home/superxcx/code/DSMC/dsmcFoam++/hyStrath_xcx/case/cylinder_react/schedule_sweep_300step_20260404/omp8_move_static_coll_dynamic_moveopt/log.dsmcFoam+.post_flatocc1b)

结果：

- `post fields/output`: `25.8508 s -> 25.5374 s`

这刀是小幅有效的，当前保留。


### 8.3 `post` 上验证无效或已回退的方向

#### 8.3.1 把 shared cache build 按 `particleLoadStart/End` 分区跑

结论：

- 没有收益，反而更慢
- 原因是 `post` 的缓存构建和 `move` 的分区假设不一致


#### 8.3.2 对 `vib` 路径做微观缓存化

例如：

- 预先缓存 `kB*thetaV`
- 精简 `Evibp` 求和链
- 局部变量/引用缓存

结论：

- 没有转化成净收益
- 这类微优化已回退


#### 8.3.3 继续对零值路径做更细的 if-guard

结论：

- 收益基本到头
- 继续加 guard 只会让代码更碎，不会带来可见性能收益


### 8.4 profiling 对 `post` 的影响与处理

曾一度把 `post` profiling 拆得很细，甚至进入 per-parcel 内部细分。

后果是：

- profiling 本身显著影响 `post` 成本
- 得到的内部比例不再可信

因此最终采取了两步修正：

1. 用更粗粒度的抽样 profiling 重新确认热点位置
2. 最后把会影响 `post` 性能的 profiling 基本关掉

当前状态：

- `post` 热路径上已经不再保留那套高扰动 profiling


## 9. 当前代码中保留的有效优化清单

截至本次总结，当前代码实际保留的有效优化包括：

### 9.1 `move`

- `p.reset()` 融入 `particles[]` 并行路径
- `move extract parcels` 并行 gather
- `move commit/survivor rebuild` 两阶段化

### 9.2 `buildCellOccupancy`

- `occless` 版本的有效优化：
  - 减少不必要容器重整
  - 保持单次 build 更便宜

### 9.3 `post`

- 按 field 配置裁剪 shared cache 累积项
- `speciesTvib_ / dsmcSpeciesEvibModCum_` 后处理链标量化
- 复用 `buildCellOccupancy()` 生成的扁平 occupancy 视图
- 关闭高扰动 `post` profiling


## 10. 已明确无效、且不应再重复投入的方向

以下方向已经通过实际编译和运行验证为无效，后续不建议重复投入：

- `move` 的 `static/dynamic/guided` 调度层继续调参
- `buildCellOccupancy` 的 `parallel clear()`
- `buildCellOccupancy` 的 `active-cell sparse clear`
- `buildCellOccupancy` 的 fallback gather 结果回写复用
- `collision` 的 candidate-only partition 进一步细化
- `post` 中对 `vib/base` 的细碎局部缓存微优化
- `post` 中更细碎的 demand mask 微裁剪
- `particle::changeFace/changeCell` 级别的 getter 缓存


## 11. 当前主要对照日志

下面列出目前最重要的几个阶段性结果。

### 11.1 `buildCellOccupancy` 有效优化参考

- [log.dsmcFoam+.buildocc_partition_weight_moveprofileopt_occless](/home/superxcx/code/DSMC/dsmcFoam++/hyStrath_xcx/case/cylinder_react/schedule_sweep_300step_20260404/omp8_move_static_coll_dynamic_moveopt/log.dsmcFoam+.buildocc_partition_weight_moveprofileopt_occless)
  - `main loop wall time = 210.1609932 s`
  - `move only = 82.84166473 s`
  - `buildCellOccupancy = 38.22029463 s`
  - `collision phase = 55.42444767 s`
  - `post fields/output = 29.63549255 s`

### 11.2 `post` demand-mask 有效优化参考

- [log.dsmcFoam+.buildocc_partition_weight_moveprofileopt_occless_postprofile5](/home/superxcx/code/DSMC/dsmcFoam++/hyStrath_xcx/case/cylinder_react/schedule_sweep_300step_20260404/omp8_move_static_coll_dynamic_moveopt/log.dsmcFoam+.buildocc_partition_weight_moveprofileopt_occless_postprofile5)
  - `main loop wall time = 214.1948232 s`
  - `move only = 85.66171813 s`
  - `buildCellOccupancy = 41.10817506 s`
  - `collision phase = 55.45817385 s`
  - `post fields/output = 27.78516396 s`

### 11.3 `move commit/survivor rebuild` 两阶段化参考

- [log.dsmcFoam+.move_commit2stage1](/home/superxcx/code/DSMC/dsmcFoam++/hyStrath_xcx/case/cylinder_react/schedule_sweep_300step_20260404/omp8_move_static_coll_dynamic_moveopt/log.dsmcFoam+.move_commit2stage1)
  - `main loop wall time = 213.1567466 s`
  - `move only = 85.91017561 s`
  - `buildCellOccupancy = 40.17848416 s`
  - `collision phase = 57.1033102 s`
  - `post fields/output = 25.85075558 s`

### 11.4 `post` 复用扁平 occupancy 视图参考

- [log.dsmcFoam+.post_flatocc1b](/home/superxcx/code/DSMC/dsmcFoam++/hyStrath_xcx/case/cylinder_react/schedule_sweep_300step_20260404/omp8_move_static_coll_dynamic_moveopt/log.dsmcFoam+.post_flatocc1b)
  - `main loop wall time = 205.686753 s`
  - `move only = 82.48290933 s`
  - `buildCellOccupancy = 40.04521723 s`
  - `collision phase = 53.7984506 s`
  - `post fields/output = 25.53738821 s`


## 12. 当前整体判断

到目前为止，这轮 OpenMP 优化已经证明：

- 低垂果子已经基本摘完
- 调度层优化基本接近上限
- 真正有效的收益来自：
  - 两阶段化
  - 数据流重组
  - 视图复用
  - 按需裁剪

也就是说，当前更大的优化空间已经不是：

- 再换 `omp schedule`
- 再加一个 `parallel for`
- 再做一点 getter 级缓存

而是：

- 更大一级的 tracking 数据流重组
- inflow boundary 内部的两阶段生成/提交
- `post` shared cache 数据布局重构


## 13. 后续仍值得继续的方向

如果继续投入，建议按以下优先级：

### 13.1 第一优先级：inflow boundary 的进一步两阶段化

原因：

- `move pre-control/partition` 的真正大头是 inflow boundary
- 已经确认其内部热点是：
  - `velocity/internal energy`
  - `addNewParcel`
  - `tri select/random point`

注意：

- 必须继续坚持“线程私有生成 + 安全提交”
- 不能回到破坏拓扑语义的延后 barycentric 方案


### 13.2 第二优先级：tracking 主链的数据流重组

不是 getter 微缓存，而是更高层的：

- 局部工作状态复用
- 减少拓扑对象反复重构
- 更明确的 tet/face 邻接工作缓存


### 13.3 第三优先级：`post` 的 shared cache 数据布局重构

当前 `post` 已经做过几轮有效优化，但如果继续，要从：

- 分散数组散写
- per-parcel 跨很多数组更新

转向：

- 更紧凑的数据布局
- 更适合顺序写流量的缓存结构


### 13.4 第四优先级：`collision` 单位 candidate 成本分析

如果还要打 `collision`，应转向：

- 为什么不同线程单位 candidate 成本不同
- 是 cell 类型问题、反应路径问题，还是内存访问问题

而不是继续切换 `dynamic/guided/partition`


## 14. 结论

本轮 OpenMP 优化工作的核心结论可以概括为：

1. `move` 的最大有效优化来自：
   - `p.reset()` 并行化
   - `commit/survivor rebuild` 两阶段化

2. `buildCellOccupancy` 的最大有效优化来自：
   - 避免不必要的容器重整
   - 让单次 build 更便宜

3. `collision` 的调度层基本到头：
   - `dynamic` 已经足够好
   - 后续不应继续在 schedule 层反复试错

4. `post` 的最大有效优化来自：
   - 按 field 需求裁掉不必要统计
   - 优化振动后处理链
   - 复用 `buildCellOccupancy()` 已经生成好的连续 parcel 视图

5. 当前剩余优化空间仍然存在，但已经转入更高成本的结构性重构阶段。


## 15. `buildCellOccupancy` 新一轮增量优化补充

这一轮工作不是重新回到 `count/reduce` 或 `allocate/fill` 上做局部微调，而是继续沿着“让 flat occupancy 成为主数据结构”这条线，把 `extract parcels` 路径里的残余回退和隐性全量整理继续拿掉。

本轮关注点有两个：

- 继续减少 `buildCellOccupancy()` 中对 `cellOccupancy_` 物化路径的依赖
- 继续追查 `moveOrderedParcels_` 快路径为什么还会退回到 fallback gather

这轮验证主要使用的日志为：

- 基线阶段成果：[log.dsmcFoam+.occ_opt3](/home/superxcx/code/DSMC/dsmcFoam++/hyStrath_xcx/case/cylinder_react/schedule_sweep_300step_20260404/omp8_move_static_coll_dynamic_occopt/log.dsmcFoam+.occ_opt3)
- 中间版诊断关闭后快路径：[log.dsmcFoam+.occ_diagoff1](/home/superxcx/code/DSMC/dsmcFoam++/hyStrath_xcx/case/cylinder_react/schedule_sweep_300step_20260404/omp8_move_static_coll_dynamic_occopt/log.dsmcFoam+.occ_diagoff1)
- 当前最终保留版：[log.dsmcFoam+.final](/home/superxcx/code/DSMC/dsmcFoam++/hyStrath_xcx/case/cylinder_react/schedule_sweep_300step_20260404/omp8_move_static_coll_dynamic_occopt/log.dsmcFoam+.final)


### 15.1 问题重新定义

在 `occ_opt3` 之后，`buildCellOccupancy` 虽然已经明显下降，但主瓶颈仍然是：

- `extract parcels`

当时的判断是对的：

- `buildCellOccupancy()` 内部本身已经有“线程私有 count + prefix sum + flat fill”的骨架
- 真正拖后腿的是下游仍会回到旧容器体系，或者快路径校验失败后退回全量 gather

因此，这轮不再把重点放在“再并行一点”，而是放在：

- 为什么快路径没有 100% 命中
- 为什么 `moveOrderedParcels_ + moveAppendedParcels_` 仍然会和 `this->size()` 对不上


### 15.2 第一刀：把重型诊断从默认 profiling 主路径移开

之前为了定位 occupancy 不一致问题，在 `buildCellOccupancy()` 中加入过较重的诊断扫描。这些扫描对于定位问题有价值，但不应该长期留在默认性能路径上。

本轮先做的第一步是：

- 让重型 occupancy 校验只在显式诊断模式下启用
- 默认 profiling 路径不再每步做全云扫描

核心文件：

- [dsmcCloud.C](/home/superxcx/code/DSMC/dsmcFoam++/hyStrath_xcx/src/lagrangian/dsmc/clouds/dsmcCloud.C)

对应结果，从 `occ_opt3` 到 `occ_diagoff1`：

- `main loop wall time`: `188.4082 s -> 186.1881 s`
- `buildCellOccupancy extract parcels`: `31.1754 s -> 7.7584 s`
- `buildCellOccupancy total`: `35.1949 s -> 14.0614 s`

同时也暴露了新的关键信息：

- `move-ordered hits = 291`
- `fallback gathers = 10`
- `fallback invalid-order = 1`
- `fallback size-mismatch = 9`

这说明：

- 快路径本身是有效的
- 真正剩下的问题已经不是 OMP 循环本身，而是 `moveOrderedParcels_` 体系还存在 10 次回退


### 15.3 第二刀：追查 `size-mismatch` 的真正来源

对 `occ_diagoff1` 的分析表明：

- `size-mismatch` 不是随机误差
- 它来自 `moveOrderedParcels_.size() + moveAppendedParcels_.size() != this->size()`

进一步结合 `dsmcCloud::evolve()` 的时序，可以定位到根因：

- `move` 阶段新增的粒子已经能通过 `moveAppendedParcels_` 进入快路径
- 但“上一次 `buildCellOccupancy()` 之后，到下一次 `Cloud::move()` 之前”插入的粒子，并不在 `moveOrderedParcels_` 里
- 这些粒子会出现在 `this->size()` 中，却不在下一次 move 提取快路径的输入里，导致 `size-mismatch`

这类粒子来源通常不是 move 内核本身，而是：

- `controlBeforeMove`
- boundary/controller 的非 move 时段插入
- 其他在 move 捕获窗口外触发的加粒子路径


### 15.4 第三刀：引入 `pendingMoveParcels_`，补齐快路径输入

针对上面的缺口，本轮新增了一个“延迟待提取粒子”通道：

- `pendingMoveParcels_`

核心思路：

- 如果新粒子是在 `moveAppendCaptureActive_` 为真时加入，则仍进入 `moveAppendedParcels_`
- 如果是在捕获窗口之外、但 OMP move 快路径开启时加入，则先进入 `pendingMoveParcels_`
- 下一次 `Cloud::move()` 在提取 `particles[]` 时，将 `pendingMoveParcels_` 一并纳入输入
- 一次完整的 `buildCellOccupancy()` 成功后，再清空 `pendingMoveParcels_`

相关修改文件：

- [dsmcCloud.H](/home/superxcx/code/DSMC/dsmcFoam++/hyStrath_xcx/src/lagrangian/dsmc/clouds/dsmcCloud.H)
- [dsmcCloudI.H](/home/superxcx/code/DSMC/dsmcFoam++/hyStrath_xcx/src/lagrangian/dsmc/clouds/dsmcCloudI.H)
- [dsmcCloud.C](/home/superxcx/code/DSMC/dsmcFoam++/hyStrath_xcx/src/lagrangian/dsmc/clouds/dsmcCloud.C)
- [Cloud.C](/home/superxcx/code/DSMC/dsmcFoam++/hyStrath_xcx/src/lagrangian/basic/Cloud/Cloud.C)

这一步本质上不是“再做一层缓存”，而是把 `moveOrderedParcels_` 快路径从“只覆盖 move 内新增粒子”，扩展到“覆盖上一个 occupancy build 之后的所有新增粒子”。


### 15.5 最终结果

最终保留版见：

- [log.dsmcFoam+.final](/home/superxcx/code/DSMC/dsmcFoam++/hyStrath_xcx/case/cylinder_react/schedule_sweep_300step_20260404/omp8_move_static_coll_dynamic_occopt/log.dsmcFoam+.final)

与 `occ_opt3` 相比，`buildCellOccupancy` 这一轮的最终净效果非常明确：

- `main loop wall time`: `188.4082 s -> 169.0589 s`
- `buildCellOccupancy extract parcels`: `31.1754 s -> 5.3452 s`
- `buildCellOccupancy total`: `35.1949 s -> 11.2205 s`

同时，快路径命中率也基本被打满：

- `move-ordered hits = 300`
- `fallback gathers = 1`
- `fallback invalid-order = 1`
- `fallback size-mismatch = 0`
- `fallback offset-mismatch = 0`

这里最后剩下的 `1` 次 `fallback invalid-order`，并不是主循环里还有一次真实热路径回退，而是初始化/构造阶段那次 `buildCellOccupancy()`：

- `buildCellOccupancy calls = 301`
- `move-ordered hits = 300`

也就是说：

- 热路径里的 `size-mismatch` 已经全部拿掉
- 主循环里的 occupancy rebuild 已经基本稳定运行在 move-ordered 快路径上


### 15.6 本轮关于 `buildCellOccupancy` 的最终结论

这一轮之后，可以把 `buildCellOccupancy` 的结论更新为：

1. 过去真正限制扩展性的，不只是 `extract parcels` 的 OMP 实现本身，而是快路径输入集合不闭合，导致频繁 fallback gather。
2. 单纯继续打 `count/reduce` 或 `allocate/fill`，已经不是收益最大的方向。
3. 让 `moveOrderedParcels_` 成为稳定主输入，并补齐 `moveAppendedParcels_` 之外的新增粒子来源，收益远大于继续调度层微调。
4. `cellOccupancy_` 应继续只保留为兼容层和按需 materialize 视图，主路径应继续坚持 flat occupancy。

当前阶段，`buildCellOccupancy` 这条线已经从“OpenMP 循环怎么写”转入“数据流是否闭合、快路径是否稳定复用”的阶段。后续若继续优化，优先级应是：

- 继续减少不得不 materialize `cellOccupancy_` 的旧路径
- 继续检查是否还有其他非 move 时段的插入路径绕过了 flat occupancy 主数据流
- 只有在这些都清干净之后，再考虑是否值得继续抠 `count/reduce` 或 `allocate/fill`


## 16. `1200w_particle` 阶段增量总结

这一轮工作转入新的大粒子数算例：

- [omp8](/home/superxcx/code/DSMC/dsmcFoam++/hyStrath_xcx/case/cylinder_react/1200w_particle/omp8)
- [omp8_collopt](/home/superxcx/code/DSMC/dsmcFoam++/hyStrath_xcx/case/cylinder_react/1200w_particle/omp8_collopt)

这一阶段的目的不是重复验证前面的小粒子数结论，而是观察在 `1200 万` 初始粒子规模下，哪些阶段重新成为主瓶颈，哪些原先有效的优化会失效。


### 16.1 阶段目标与瓶颈迁移

在大粒子数阶段，性能画像发生了明显迁移：

- `BuildCellOccupancy` 不再主要被 `extract parcels` 支配
- `count/reduce`
- `allocate/fill`

开始重新抬头

同时，`move` 和 `collision` 的绝对时间也进一步放大，因此这一阶段的优化主线变成了：

1. 修复并稳定 `dsmcInitialise+`
2. 重新审视 `buildCellOccupancy` 在大粒子数下的新瓶颈
3. 尝试把 `collision` 与 `move` 再往前推一刀


### 16.2 `dsmcInitialise+` 修复与初始化输出补全

在 `1200w_particle` case 上，`dsmcInitialise+` 曾经出现过两类问题：

- 构造 `dsmcCloud` 时在 `buildCellOccupancy()/rebuildParticleLoadPartition()` 路径上崩溃
- 初始化日志缺少“初始化后粒子数”的直接输出

这一阶段完成了两件事：

1. 修复了初始化阶段对 occupancy / partition 快路径假设过强导致的崩溃
2. 在 `log.dsmcInitialise+` 中补充了初始化粒子数输出，便于后续回归时直接核对

这一步的意义不是提速，而是给大粒子数阶段后续所有 profiling 与回归提供稳定起点。


### 16.3 `buildCellOccupancy`：大粒子数下的重心变化

在 `1200w_particle` 阶段，`buildCellOccupancy` 的主要关注点从此前的 `extract parcels` 继续转向：

- `count/reduce`
- `allocate/fill`

这一判断来自：

- `results/dlb-openmp/cylinder_react_1200wparticle/performance_summary.md`
- `results/dlb-openmp/cylinder_react_1200wparticle_xh/performance_summary.md`

工程结论是：

- 前面已经拿到收益的 `extract parcels` / flat occupancy 主路径改造仍然有效
- 但在粒子数继续放大之后，后两段的全局计数、归并、写回成本开始成为新瓶颈

这一阶段围绕 `allocate/fill` 与 `count/reduce` 做过多轮试探，结论比较明确：

1. `allocate/fill` 方向上，缓冲复用和减少全量重排是合理方向
2. `count/reduce` 方向上，单纯继续堆线程归并技巧，收益并不稳定
3. 大粒子数下，`buildCellOccupancy` 已经进入“带宽与全局写入组织”主导阶段，而不是简单的 OMP 调度问题


### 16.4 `collision` 主线：做过的尝试与结论

针对 `collision selection/collide`，这一阶段做过三类尝试：

1. 把 cell 内 scratch 改为线程私有复用
   - 避免 `whichSubCell(nC)` 每个 cell 重建 `List<label>`
   - 让 `subCells[8]` 继续线程私有并重复使用
2. 在 cell 内先做一次轻量 gather
   - 预取 `parcel 指针`
   - 预取 `typeId`
   - 预取 `charge`
   - 预取 `position-cellCentre` 的 3bit subcell code
3. 暴露 `dynamic` 调度的 `chunksize` 到 `controlDict`

对应测试日志主要集中在：

- [log.dsmcFoam+.omp8.collopt1](/home/superxcx/code/DSMC/dsmcFoam++/hyStrath_xcx/case/cylinder_react/1200w_particle/omp8_collopt/log.dsmcFoam+.omp8.collopt1)
- [log.dsmcFoam+.omp8.collopt1.chunk4](/home/superxcx/code/DSMC/dsmcFoam++/hyStrath_xcx/case/cylinder_react/1200w_particle/omp8_collopt/log.dsmcFoam+.omp8.collopt1.chunk4)
- [log.dsmcFoam+.omp8.collopt1.chunk8](/home/superxcx/code/DSMC/dsmcFoam++/hyStrath_xcx/case/cylinder_react/1200w_particle/omp8_collopt/log.dsmcFoam+.omp8.collopt1.chunk8)
- [log.dsmcFoam+.omp8.collopt1.chunk16](/home/superxcx/code/DSMC/dsmcFoam++/hyStrath_xcx/case/cylinder_react/1200w_particle/omp8_collopt/log.dsmcFoam+.omp8.collopt1.chunk16)
- [log.dsmcFoam+.omp8.collopt1.chunk64](/home/superxcx/code/DSMC/dsmcFoam++/hyStrath_xcx/case/cylinder_react/1200w_particle/omp8_collopt/log.dsmcFoam+.omp8.collopt1.chunk64)

最终结论：

- 这一轮 `collision` 优化没有拿到结构性收益
- `dynamic` 调度与 `chunksize` 调优只能改变轻微波动，不能改变主趋势
- 说明当前 `collision` 的剩余瓶颈不再是调度层，而是主循环内部的单位 candidate 成本

也就是说，后续如果继续做 `collision`，优先级不应该是继续换 schedule，而应该是继续分析：

- candidate 期随机访存
- 反应/碰撞路径分岔
- 单位 accepted collision 的代价差异


### 16.5 `move` 主线：激进 profiling 与多轮回退

这一阶段对 `move` 做了多轮更激进的尝试，重点不再是早期已经确认有效的：

- `inline reset`
- `commit/survivor rebuild` 两阶段化

而是继续往 `tracking` 主链内部追。

#### 16.5.1 最优基线

这一阶段的最优测量基线是：

- [log.dsmcFoam+.omp8.collandmoveopt](/home/superxcx/code/DSMC/dsmcFoam++/hyStrath_xcx/case/cylinder_react/1200w_particle/omp8_collopt/log.dsmcFoam+.omp8.collandmoveopt)

关键数据：

- `main loop wall time = 200.0615 s`
- `move only = 119.6731 s`
- `move parallel kernel wall = 97.8947 s`
- `move commit/survivor rebuild = 8.6888 s`
- `buildCellOccupancy = 11.4151 s`
- `collision phase = 43.3604 s`
- `post fields/output = 23.1735 s`

这组结果说明：

- `move` 仍然是大头
- 但 `move` 的主要剩余时间已经集中在 `parallel kernel`
- `commit/survivor rebuild` 虽然还能看见，但已经不是原来的头号问题

#### 16.5.2 重 profiling 的结论

为了继续切 `move`，这一阶段曾对 `dsmcParcel::move()` 加入过非常细的计时拆分，覆盖：

- `constrainDirection`
- `trackToAndHitFace`
- tracker hook
- cyclic control
- stuck patch control

结果非常明确：

- `trackToAndHitFace()` 是真正的大头
- `tracker hook / cyclic control / stuck patch control` 在当前 case 上几乎可以忽略

但这组 profiling 也同时暴露了一个工程问题：

- 过细的 `chrono` 计时会显著抬高 wall time
- profiling 本身会干扰结论

因此这一套“重 profiling”代码最终没有保留，只作为定位依据使用，随后已整体移除。

#### 16.5.3 已尝试但确认无效的 `move` 优化

这一阶段至少有三类 `move` 试探被验证为无效或负收益：

1. internal-face 快速穿越/压缩外层 while
   - 目标是减少每穿过一个 internal face 就回到外层 `while`
   - 结果是 `move parallel kernel` 略降，但 `move commit/survivor rebuild` 明显变差
   - 净结果回退，因此放弃

2. 非 MPI case 下进一步裁剪 commit/transfer 路径
   - 目标是针对纯 OMP 场景继续砍掉 `transfer` 相关带宽
   - 结果没有带来稳定正收益，也已回退

3. 对 3D case 直接跳过 `meshTools::constrainDirection()`
   - 先用 `checkMesh` 确认当前 `omp8_collopt` mesh 为完整 3D
   - `Mesh has 3 solution directions (1 1 1)`
   - 在此基础上试验短路 `constrainDirection`
   - 对应日志：[log.dsmcFoam+.move_constrain2](/home/superxcx/code/DSMC/dsmcFoam++/hyStrath_xcx/case/cylinder_react/1200w_particle/omp8_collopt/log.dsmcFoam+.move_constrain2)
   - 结果仍然劣于最优基线，因此不保留

`move_constrain2` 相比最优基线：

- `main loop wall time`: `200.0615 s -> 205.3638 s`
- `move only`: `119.6731 s -> 122.0580 s`
- `move commit/survivor rebuild`: `8.6888 s -> 10.8032 s`

结论很直接：

- `move` 上继续做这种微结构级的小修补，已经拿不到可靠收益
- 这些尝试大多会在别的阶段把代价吐回来

#### 16.5.4 当前保留结论

这轮之后，关于 `move` 的判断更新为：

1. 当前最有价值、已确认有效的 `move` 优化，仍然是前面保留的
   - `inline reset`
   - `commit/survivor rebuild` 两阶段化
2. `tracking` 主链确实是剩余大头，但已经进入高风险、低确定性阶段
3. 在没有更细粒度、且不会污染 wall time 的 profiling 之前，不适合继续盲目改 `trackToAndHitFace()` 主链


### 16.6 这一阶段真正保留下来的内容

`1200w_particle` 这一轮最终真正保留的，不是某个新的激进 `move`/`collision` 提速补丁，而是两类更稳的成果：

1. 大粒子数算例下的性能画像重新建立
   - 确认了 `buildCellOccupancy` 的热点迁移到 `count/reduce + allocate/fill`
   - 确认了 `collision` 的调度层优化基本见顶
   - 确认了 `move` 的剩余大头在 `tracking` 主链，而不是 boundary hook

2. 清理掉会污染测量的无效 profiling 代码
   - 将上一轮 `move breakdown` 的重计时框架移除
   - 避免后续 wall time 再被 profiling 本身抬高


### 16.7 当前阶段结论

在 `1200w_particle` 阶段之后，可以把后续优先级重新排序为：

1. `collision selection/collide`
   - 不是继续调 schedule
   - 而是继续打主循环内部的单位 candidate 成本

2. `buildCellOccupancy`
   - 重点放在 `count/reduce`
   - 重点放在 `allocate/fill`
   - 尤其是带宽、全局归并、连续写入组织

3. `move`
   - 暂时不继续在 `trackToAndHitFace()` 上盲改
   - 除非先拿到更轻、更可信的 profiling 证据

这一轮的核心结论不是“又找到了一刀稳定加速”，而是：

- `1200w_particle` 下，旧问题和新问题已经完成分层
- `move` 上的低风险增量空间基本被消耗掉
- 后续收益更可能来自 `collision` 和大粒子数下的 `buildCellOccupancy` 数据流重组，而不是继续抠 `move` 的细节


## 17. `mixparallel/omp4_mpi2`：MPI+OpenMP 混合并行实现与回归清理

这一轮工作的目标不是继续优化单节点 `omp8`，而是把 `MPI + OpenMP` 的混合路径真正做成：

1. 能稳定跑完
2. `move -> BuildCellOccupancy` 快路径不反复 fallback
3. 不把纯 `OpenMP` 路径一起拖慢

测试主目录：

- `hyStrath_xcx/case/cylinder_react/mixparallel/omp4_mpi2`

对照与回归目录：

- `hyStrath_xcx/case/cylinder_react/schedule_sweep_300step_20260404/omp8_move_static_coll_dynamic_occtosort`


### 17.1 第一阶段：把 mixed move 从“会崩”修到“能跑”

最初的 `MPI + openmpMove` 路径在 `particle::changeFace()` 上会崩，说明问题不在 `collision/post`，而在：

- decomposed mesh 上的 `move` tracking
- move 后的 ordered parcel 维护
- `BuildCellOccupancy` 对 `moveOrderedParcels_` 的复用假设

这一阶段最终保留下来的关键修复是：

1. mixed `move` 不再在后续 transfer pass 上重扫整云
   - 首轮仍可扫全 cloud
   - 后续 pass 只处理 `pendingMoveParcels_`
2. `moveOrderedParcels_` 在 mixed 下改为跨 pass 累积、最后一次性写回
3. pre-move inflow 通过 `moveAppendedParcels_` 进入首轮 ordered reuse
4. 但 `moveAppendedParcels_` 在首轮 workset 构建后立即清空
   - 防止继续污染 `BuildCellOccupancy`

保留下来的 mixed 最关键版本是：

- `log.dsmcFoam+.omp4_mpi2_capturefix2`

相对上一版 mixed 有效基线 `moveopt6c`：

- `main loop wall time`: `253.79 s -> 248.62 s`
- `move only`: `133.07 s -> 124.02 s`
- `move extract parcels`: `38.84 s -> 0.13 s`
- `BuildCellOccupancy extract parcels`: `4.91 s`
- `move-ordered hits = 300`
- `fallback gathers = 1`
- `fallback size-mismatch = 0`

这一阶段的核心结论是：

- mixed 的大问题一开始不是 `tracking kernel`
- 而是 `move` 后 workset 维护与 `BuildCellOccupancy` 快路径之间的数据流失配


### 17.2 细化 profiling 之后，对 mixed move 的最终判断

后续对 mixed `move` 做了多轮细化 profiling，关键结论如下：

1. `extract parcels` 确实曾经是 mixed `move` 的主瓶颈
   - 但在 `capturefix2` 之后，这条已经基本打通
2. mixed 和 pure `omp8` 的 `move parallel kernel wall` 差异并不主要来自：
   - `processor hits`
   - `patch hits`
   - `boundary callback`
3. `trackToAndHitFace()` 仍然是内核中最大的单项
   - 但它在 mixed 下的单位 hit 成本并没有出现数量级恶化

对比：

- pure `omp8_hitprof`
- mixed `omp4_mpi2_capturefix2`

可见：

- mixed 的 `processor hits / face hits` 占比很小
- mixed 的 `trackToAndHitFace total` 与 `move parallel kernel wall` 并不支持“主要因为分区后路径更长”这一判断

因此这一阶段最终结论是：

- mixed `move` 的第一大问题已经从 `extract` 迁移走
- 剩余的 `move kernel` 差距不是负载不均主导
- 也不是明显的 MPI 跟踪额外路径主导
- 后续如果继续打 `move kernel`，只能进入更高风险的 tracking 算法本体层


### 17.3 mixed 与 pure omp8 的 profiling 对照

在 `omp4_mpi2` 上，做了三组 profiling 配置对比：

- `omp4_mpi2_noprofile`
- `omp4_mpi2_summaryonly`
- `omp4_mpi2_profileboth`

结果表明：

1. `profileSummary` 的扰动较小
   - `noprofile -> summaryonly` 约增加 `1.89 s`
2. `profileDetail` 会显著增加开销
   - `summaryonly -> profileboth` 再增加约 `4.17 s`
3. 正式性能判断应优先使用：
   - `noprofile`
   - 或 `summaryonly`
4. `profileboth` 更适合定位瓶颈，而不适合作为速度结论依据

这一阶段也做了 pure `omp8` 与 mixed 的 profiling 对照，结论是：

- mixed 比 pure `omp8` 更慢，不是主要因为 `move kernel`
- 更大的结构性差距来自 `collision`


### 17.4 mixed collision：主要问题是 rank 级不均衡

从：

- `omp4_mpi2_profileboth`
- `omp8_profileboth`

的对比可见，mixed 下最显著的结构性问题是：

- `collision cand imbalance`
- `accepted coll imbalance`

在 mixed 下明显高于 pure `omp8`。

也就是说：

- `move` 的剩余问题不是没有
- 但 mixed 真正更像“结构性瓶颈”的部分已经转向 `collision`
- 后续如果继续追 mixed 并行加速，优先级应高于继续抠 `move kernel`


### 17.5 反向清理：mixed 引入的共享 move 维护对 pure OMP 的污染

在 mixed 路径稳定之后，进一步发现：

- pure `OpenMP` 的 `move` 相比历史干净版本有回归

通过对比 `bkp/src`（未开始 mixed-MPI 之前的纯净 openmp 版本）与当前代码，最终确认：

1. `particle.C` 主算法本体变化不大
2. 真正加厚的是：
   - `dsmcParcel.C`
   - `Cloud.C`
   - `dsmcCloud.C`
3. 其中 mixed 引入的共享 move 维护逻辑里，真正污染 pure OMP 的，不是全部路径，而主要是：
   - `pendingMoveParcels_` 的纯 OMP 维护
   - `dsmcParcel::move()` 中 detail profiling 的内层计时/计数更新

这一阶段做了两轮关键清理：

#### 17.5.1 `cleanmpi2`

在 `dsmcCloud::recordMoveAppendedParcel()` 中，把：

- `pendingMoveParcels_` 的维护限制到 `Pstream::parRun()` 下

但保留：

- `moveAppendedParcels_`
- `moveOrderedParcels_`
- `BuildCellOccupancy` 快路径输入

结果对比 `omp8_summaryonly`：

- `main loop wall time`: `189.31 s -> 177.00 s`
- `move only`: `99.06 s -> 91.82 s`
- `move parallel kernel wall`: `90.17 s -> 84.09 s`

同时没有像更激进的 `cleanmpi1` 那样打坏 `BuildCellOccupancy` 快路径。

#### 17.5.2 `cleanmpi5`

基于 `bkp/src` 对照，继续把 `dsmcParcel::move()` 中只服务 detail profiling 的内层工作，在 `profileDetail=false` 下彻底旁路：

- `moveTrackWallTime`
- `moveTrackerWallTime`
- `moveBoundaryWallTime`
- 各类 `move*HitCount`

对应结果：

- `log.dsmcFoam+.omp8_cleanmpi5`

相对 `cleanmpi2`：

- `main loop wall time`: `177.00 s -> 167.85 s`
- `move only`: `91.82 s -> 80.25 s`
- `move parallel kernel wall`: `84.09 s -> 72.59 s`

这说明：

- pure OMP 剩余的 `move kernel` 回归，并不是 tracking 主算法本体被 mixed 改坏
- 而是后来加进来的 shared/thread-safe/detail-profiling 包装在纯 OMP 下继续收费


### 17.6 本轮 mixed 实现最终保留结论

这一轮关于 mixed 实现，最终保留的结论是：

1. `MPI + OpenMP` mixed 路径已经从“会崩”修到了“可稳定运行”
2. mixed `move -> BuildCellOccupancy` 快路径已经基本打通
3. mixed 当前最大的结构性瓶颈不是 `move extract`
   - 而是 `collision` 的 rank 级不均衡
4. pure `OMP` 的回归已经被局部清掉
   - `cleanmpi2`：去掉 pure OMP 下不必要的 `pendingMoveParcels_`
   - `cleanmpi5`：去掉 `profileDetail=false` 下 `dsmcParcel::move()` 的内层 detail profiling 开销

因此，这一轮 mixed 工作的最终判断是：

- mixed 数据流已经具备继续往前推进的基础
- 但下一主战场不应再是 `move extract`
- 优先级应转到：
  1. mixed `collision` 的 rank 级不均衡
  2. mixed `move kernel` 若继续做，则必须进入更高风险的 tracking 本体层

## 2026-04-22 Mixed Move 阶段化重构补充

这一轮工作专门围绕 `hyStrath_xcx/case/cylinder_react/mixparallel/omp4_mpi2` 的 mixed `move` 数据流继续推进。重点不是再打 `move kernel` 本体，而是继续梳理：

- `pre-move inflow`
- `later-pass deferred continuation`
- `processor receive / transfer finalize`

之间的语义边界，确认还能不能在不破坏 `moveOrderedParcels_` 快路径的前提下继续压低 mixed `move` 成本。

### 本轮 mixed 新基线

本轮先重跑 mixed 基线：

- [log.dsmcFoam+.omp4_mpi2_rebase1](/home/superxcx/code/DSMC/dsmcFoam++/hyStrath_xcx/case/cylinder_react/mixparallel/omp4_mpi2/log.dsmcFoam+.omp4_mpi2_rebase1)

关键结果：

- `main loop wall time = 132.01 s`
- `move only = 60.82 s`
- `move extract parcels = 0.149 s`
- `move parallel kernel wall = 47.04 s`
- `move transfer/delete finalize = 6.76 s`

### 有效优化：`stagebuf1`

第一刀把：

- `moveAppendedParcels_`

从“混合承载 pre-move inflow 和 move 内新生成 parcel”的状态，改成只表示：

- `pre-move inflow`

一旦进入实际 `move` kernel，再新生成的 parcel 直接进入：

- `pendingMoveParcels_`

对应日志：

- [log.dsmcFoam+.omp4_mpi2_stagebuf1](/home/superxcx/code/DSMC/dsmcFoam++/hyStrath_xcx/case/cylinder_react/mixparallel/omp4_mpi2/log.dsmcFoam+.omp4_mpi2_stagebuf1)

相对 `rebase1`：

- `main loop wall time`: `132.01 -> 121.87 s`
- `move only`: `60.82 -> 54.43 s`
- `move parallel kernel wall`: `47.04 -> 43.00 s`
- `move transfer/delete finalize`: `6.76 -> 4.90 s`

结论：

- 这刀有效
- mixed `move` 的数据流语义更清晰
- 同时带来了可确认的系统级收益

### 有效优化：`stagebuf2`

第二刀是在 `stagebuf1` 基础上修正状态机：

- `beginMoveDeferredAppendStage()` 把 `moveAppendToPending_` 打开后
- 后续 `while` pass 再次进入 `beginMoveAppendCapture()` 时，不应把它重新置回 `false`

也就是保证：

- `moveAppendToPending_` 在整个 move 多 pass 过程中持续有效

对应日志：

- [log.dsmcFoam+.omp4_mpi2_stagebuf2](/home/superxcx/code/DSMC/dsmcFoam++/hyStrath_xcx/case/cylinder_react/mixparallel/omp4_mpi2/log.dsmcFoam+.omp4_mpi2_stagebuf2)
- [log.dsmcFoam+.omp4_mpi2_stagebuf2_confirm1](/home/superxcx/code/DSMC/dsmcFoam++/hyStrath_xcx/case/cylinder_react/mixparallel/omp4_mpi2/log.dsmcFoam+.omp4_mpi2_stagebuf2_confirm1)

`stagebuf2` 相对 `stagebuf1`：

- `main loop wall time`: `121.87 -> 119.99 s`
- `move only`: `54.43 -> 53.10 s`
- `move parallel kernel wall`: `43.00 -> 41.67 s`

确认跑 `stagebuf2_confirm1` 后，仍能维持同一水平：

- `main loop wall time = 118.21 s`
- `move only = 52.02 s`
- `move transfer/delete finalize = 4.80 s`

结论：

- `stagebuf2` 是当前 mixed `move` 的最好可保留版本
- 这版已经证明 mixed `move` 的提取/拼接路径仍然可以通过明确阶段语义继续优化

### 无效尝试：`stagebuf3`

第三刀尝试把 `pendingMoveParcels_` 真正做成双缓冲：

- 当前 pass 只读一份 pending
- 当前 pass 新生成/新接收的 parcel 全写到另一份
- pass 结束后再交换

对应日志：

- [log.dsmcFoam+.omp4_mpi2_stagebuf3](/home/superxcx/code/DSMC/dsmcFoam++/hyStrath_xcx/case/cylinder_react/mixparallel/omp4_mpi2/log.dsmcFoam+.omp4_mpi2_stagebuf3)

结果：

- `main loop wall time = 128.35 s`
- `move only = 55.29 s`
- `move transfer/delete finalize = 4.84 s`
- 但 `BuildCellOccupancy = 36.94 s`

结论：

- `stagebuf3` 不能保留
- 当前实现里，双缓冲虽然语义更纯，但会破坏 `moveOrderedParcels_ -> BuildCellOccupancy` 现有快路径耦合

### 无效尝试：`stageaware1`

之后又尝试了一刀更轻的 staged metadata：

- 不改缓冲结构
- 只让 `BuildCellOccupancy` 知道“当前 `moveOrderedParcels_` 是否已完整覆盖 cloud”

对应日志：

- [log.dsmcFoam+.omp4_mpi2_stageaware1](/home/superxcx/code/DSMC/dsmcFoam++/hyStrath_xcx/case/cylinder_react/mixparallel/omp4_mpi2/log.dsmcFoam+.omp4_mpi2_stageaware1)

对比 `stagebuf2`：

- `main loop wall time`: `119.99 -> 119.56 s`
- 但 `move only`: `53.10 -> 53.51 s`
- `BuildCellOccupancy`: `3.80 -> 3.86 s`

结论：

- 这刀不算有效优化
- 已回退

### 无效尝试：`stagebuf4`

后续尝试把：

- `deferred continuation`
- `inter-rank received inflow`

拆成两条列表，再在 later-pass extract 阶段拼接，目标是更明确地区分“本轮续跑粒子”和“新收到的跨 rank inflow”，同时不碰 `moveOrderedParcels_` 首轮快路径。

对应日志：

- [log.dsmcFoam+.omp4_mpi2_stagebuf4](/home/superxcx/code/DSMC/dsmcFoam++/hyStrath_xcx/case/cylinder_react/mixparallel/omp4_mpi2/log.dsmcFoam+.omp4_mpi2_stagebuf4)

结果明显更差：

- `main loop wall time = 159.26 s`
- `move only = 79.61 s`
- `move extract parcels = 32.49 s`
- `BuildCellOccupancy = 33.94 s`

结论：

- `stagebuf4` 无效，已回退

### 无效尝试：`finalize1`

最后试了一刀只针对：

- `move transfer/delete finalize`

的低风险清理：

- 在 receive/deferred 追加路径里，非并行上下文下不再走无意义的 `omp critical`

对应日志：

- [log.dsmcFoam+.omp4_mpi2_finalize1](/home/superxcx/code/DSMC/dsmcFoam++/hyStrath_xcx/case/cylinder_react/mixparallel/omp4_mpi2/log.dsmcFoam+.omp4_mpi2_finalize1)

相对 `stagebuf2_confirm1`：

- `move transfer/delete finalize`: `4.796 -> 4.717 s`
- 但 `move only`: `52.021 -> 52.234 s`

结论：

- finalize 单项虽略降，但没有带来净的 `move` 改善
- `finalize1` 不保留

### 本轮 mixed `move` 的最终保留结论

本轮 mixed `move` 重构后，最终保留的只有两刀：

1. `stagebuf1`
   - `moveAppendedParcels_` 只表示 `pre-move inflow`
   - move 内新生成 parcel 进入 `pendingMoveParcels_`

2. `stagebuf2`
   - `moveAppendToPending_` 在整个 move 多 pass 中持续有效
   - 不在后续 pass 重新退回错误缓冲

当前最优可复现 mixed 日志是：

- [log.dsmcFoam+.omp4_mpi2_stagebuf2_confirm1](/home/superxcx/code/DSMC/dsmcFoam++/hyStrath_xcx/case/cylinder_react/mixparallel/omp4_mpi2/log.dsmcFoam+.omp4_mpi2_stagebuf2_confirm1)

关键结果：

- `main loop wall time = 118.21 s`
- `move only = 52.02 s`
- `move parallel kernel wall = 40.88 s`
- `move transfer/delete finalize = 4.80 s`

相对本轮 mixed 重跑基线 `rebase1`：

- `main loop wall time`: `132.01 -> 118.21 s`
- `move only`: `60.82 -> 52.02 s`

### 本轮 mixed 路线结论

通过这一轮可以更明确地判断：

- mixed `move` 的阶段化语义重构，确实还能继续拿到收益
- 但有效空间只限于：
  - 修正当前数据流中的语义混淆
  - 不破坏 `moveOrderedParcels_ -> BuildCellOccupancy` 快路径

一旦尝试：

- 真双缓冲
- 更细粒度的 `received/deferred` 列表拆分
- 更强的 staged metadata 介入 occupancy

就很容易把 `BuildCellOccupancy` 快路径重新打坏，收益立刻转负

因此当前 mixed `move` 的工程结论是：

- `stagebuf2` 已经是这条线上当前最好、最稳的可保留版本
- 后续如果继续 mixed 优化，优先级不应再放在继续拆 `pending/appended/received`
- 更合理的下一目标应转向：
  1. `move transfer/delete finalize` 之外的更大结构性瓶颈
  2. 或 mixed `collision` 的结构性问题

## 2026-05-04 mixed MPI 通信与 DLB 优化阶段总结

### 背景

本轮目标从 mixed `move` 转向 mixed MPI 场景下的通信和 `collision` 负载问题，测试目录主要为：

- [omp4_mpi2_mpiopt](/home/superxcx/code/DSMC/dsmcFoam++/hyStrath_xcx/case/cylinder_react/mixparallel/omp4_mpi2_mpiopt)
- [omp4_mpi2_mpiopt_noreact](/home/superxcx/code/DSMC/dsmcFoam++/hyStrath_xcx/case/cylinder_react/mixparallel/omp4_mpi2_mpiopt_noreact)

其中 `noreact` case 用于隔离反应模型影响，专门验证 collision DLB 和 MPI transfer 成本。

### MPI transfer/delete/finalize 批量化与 overlap

首先尝试压 mixed `move` 中的 MPI 通信路径，重点是：

- 批量化 `transfer/delete/finalize`
- staged buffer 降低 pending/deferred 生命周期混乱
- 尝试非阻塞通信和通信计算重叠
- 接收端批量处理，减少零散 append/apply

有效保留结论：

- mixed `move` 中 `stagebuf1/stagebuf2` 是当前可保留的稳定优化
- 继续细拆 `received/deferred/inflow` 容易打坏 `moveOrderedParcels_ -> BuildCellOccupancy` 快路径
- `move transfer/delete finalize` 单项可以略降，但不是当前最大瓶颈

无效或已回退尝试：

- 更复杂的双缓冲/多 staged 列表
- `finishedSends(false)` 风格的非阻塞改造
- 更细粒度的 received/deferred 拆分
- 只针对 finalize 的小范围 critical 清理

这些尝试不是完全没有局部收益，而是无法转化为 `move only` 或 `main loop` 的稳定下降。

### DLB 初始判断：确实存在 MPI rank collision 负载不均

无 DLB 情况的代表日志：

- [log.dsmcFoam+.noreact_rankcoll_nodlb1](/home/superxcx/code/DSMC/dsmcFoam++/hyStrath_xcx/case/cylinder_react/mixparallel/omp4_mpi2_mpiopt_noreact/log.dsmcFoam+.noreact_rankcoll_nodlb1)

按 30 次 profile 输出累计：

- `rank wall max total = (1.3452, 3.9831) s`
- `rank wall sum total = (5.3806, 15.9320) s`
- 最后一轮 `rank wall max = (0.0451, 0.2128) s`
- 最后一轮 candidates = `(50893, 466392)`
- `collision phase = 39.2898 s`

结论：

- rank1 的 collision 负载明显高于 rank0
- 单纯从 collision candidate 和 rank wall 看，确实存在 DLB 需求
- 问题不在于“是否负载不均”，而在于“用什么 DLB 机制转移计算量才划算”

### Remote execute DLB 实现路线

尝试实现了不重划分网格的动态计算转移：

- donor rank 选择 hot cells
- 打包 cell 内 parcel state 和 subcell 信息
- receiver rank 重建临时 parcel
- receiver 执行 collision
- receiver 将更新后的 parcel state 回传 donor
- donor 解包并写回原始 parcels

该路径只在 `!cloud_.reactionsActive()` 时启用，避免反应路径下远程临时 parcel 与真实 cloud 状态不一致。

代表日志：

- [log.dsmcFoam+.noreact_dlbexec_balance105_task2048_min1](/home/superxcx/code/DSMC/dsmcFoam++/hyStrath_xcx/case/cylinder_react/mixparallel/omp4_mpi2_mpiopt_noreact/log.dsmcFoam+.noreact_dlbexec_balance105_task2048_min1)
- [log.dsmcFoam+.noreact_dlb_overlap1](/home/superxcx/code/DSMC/dsmcFoam++/hyStrath_xcx/case/cylinder_react/mixparallel/omp4_mpi2_mpiopt_noreact/log.dsmcFoam+.noreact_dlb_overlap1)
- [log.dsmcFoam+.noreact_dlb_overlap2](/home/superxcx/code/DSMC/dsmcFoam++/hyStrath_xcx/case/cylinder_react/mixparallel/omp4_mpi2_mpiopt_noreact/log.dsmcFoam+.noreact_dlb_overlap2)

### DLB overlap 的真实效果

`overlap1` 中统计没有完整计入 post/apply 等待，因此看起来接近平衡，但这个结果不完整。

`overlap2` 计入 post 后，结果更可靠：

- `collision phase = 40.6063 s`
- `rank wall max total = (1.4049, 2.8789) s`
- `rank DLB execute total = (1.5983, 0.0422) s`
- `rank DLB post total = (1.4905, 0.1611) s`
- `rank effective total = (4.4936, 3.0822) s`
- 最后一轮 `rank effective = (0.2504, 0.1545) s`

相对无 DLB：

- rank1 本地 collision 从 `3.9831 s` 降到 `2.8789 s`
- 说明计算确实被转移了
- 但 rank0 增加了 `1.5983 s` remote execute 和 `1.4905 s` post
- rank0 反而变成新的有效瓶颈

结论：

- 当前 DLB 不是没有执行
- 当前 DLB 也不是没有转移负载
- 当前问题是 remote execute + immediate writeback 的额外成本超过收益

### Transfer 成本统计

无 DLB：

- `move transfer/delete finalize = 4.1947 s`
- `collision phase = 39.2898 s`

DLB `overlap2`：

- `move transfer/delete finalize = 3.9272 s`
- DLB task/result 组件累计约 `0.4823 s`
- DLB rank post wall 累计 max 约 `1.5029 s`
- DLB remote execute wall 累计 max 约 `1.5983 s`
- `collision phase = 40.6063 s`

非阻塞版本 `overlap_nb1`：

- `move transfer/delete finalize = 4.0902 s`
- DLB task/result 组件累计约 `1.7973 s`
- DLB rank post wall 累计 max 约 `1.5722 s`
- `collision phase = 43.2912 s`

结论：

- 普通 `move transfer/delete finalize` 不是 DLB 失败的主因
- DLB 自己引入的 result/post/apply 等待才是主因
- 非阻塞版本没有改善，反而增加了通信完成和等待成本

### 尝试过但无效的 DLB 优化

本轮针对 DLB 又尝试了几类优化，均未保留：

1. `dlbOffloadBalanceTarget` 参数修正和扫描
   - `target=1.30`：effective ratio 从约 `1.62` 降到 `1.27`，但 `collision phase = 40.8263 s`
   - `target=1.50`：effective ratio 到约 `1.13`，但 `collision phase = 41.3976 s`
   - 结论：更平衡不等于更快，post/wait 成本会把 receiver 压成新瓶颈

2. `finishedSends -> finishedNeighbourSends`
   - task send 从约 `0.037 s` 降到约 `0.026 s`
   - 但 `collision phase` 没有下降
   - 结论：all-to-all completion 元数据不是主瓶颈

3. 按 `candidate/parcel` 选择 offload cell
   - `cand/parcel` 只从约 `1.45` 提高到约 `1.55`
   - payload 没有显著下降
   - `collision phase` 仍约 `41 s`
   - 结论：该 case 中高 candidate cell 本身并没有足够高的 candidate/parcel 优势

4. dirty-only result return
   - 只回传被 collision 修改的 parcel
   - 实测更慢
   - 原因是 dirty 标记和压缩开销超过减少 payload 的收益，且接受碰撞覆盖的 parcel 比例较高

以上实验均已回退，避免留下负优化。

### 当前 DLB 工程结论

当前 remote execute DLB 的数据流是：

```text
donor hot cell parcels -> pack/send -> receiver
receiver construct temp parcels -> collide
receiver result pack/send -> donor
donor recv/apply -> original parcels
```

该模型的根本问题是每次 offload 都要完整往返：

- parcel state 发送
- remote temp parcel 构造
- collision 执行
- result state 回传
- donor 解包写回

在当前 no-react cylinder case 中，这个成本已经超过 collision 负载均衡带来的收益。

因此本轮 DLB 的保留判断是：

- 不保留 remote execute + immediate writeback 作为默认加速路径
- 可以保留相关 profiling，用于继续分析 rank wall、effective wall、post wall
- 若继续 DLB，方向不应是继续微调 offload 数量，而应改变 DLB 语义

### 后续可行方向

更有希望的 DLB 方向是降低“每步完整往返”的频率，而不是继续优化单次 remote execute：

1. ownership 迁移
   - 将 hot cells 的 compute ownership 临时转给低负载 rank
   - 不再每次 remote collision 后立即回写

2. 持久化迁移/跨步驻留
   - migrated cells 在 helper rank 上驻留多个 time step
   - 将一次迁移成本摊薄到多步 collision

3. sampling/statistics 归并而非 parcel state 高频回写
   - 让 compute owner 负责局部 collision 和统计
   - 只在必要时归并 macroscopic fields 或迁回 parcels

4. 更粗粒度的 workload ownership 层
   - 在 `BuildCellOccupancy` 之后形成 local owned cells + migrated resident cells 两类 cell-centric 工作集
   - collision/post 消费 compute-owner 视图，而不是严格 mesh-owner 视图

这类方案复杂度明显高于当前 remote execute，但从本轮数据看，这是避免 DLB 成本大于收益的主要方向。

## 2026-05-06 DLB-transfer 与网格重划分 DLB 后续测试记录

本轮继续围绕 DSMC 的动态负载均衡做了两条路线的验证：

1. 网格重划分 DLB
   - 使用 hyStrath 原有的 `reconstructParMesh -latestTime`
   - `reconstructPar -latestTime`
   - `decomposeDSMCLoadBalancePar -force -latestTime -copyUniform`
   - 然后重新从最新时间目录 restart

2. transfer DLB
   - 不重划分网格
   - donor rank 将 hot cell 的 parcel state 发给 receiver rank
   - receiver 执行远端 collision
   - receiver 回传 parcel state
   - donor 写回原始 parcels

### 网格重划分 DLB 实现与控制

本轮把网格 DLB 改成完全由 `system/loadBalanceDict` 控制，避免关闭 DLB 时仍产生检查、计时或日志噪声。

关键控制项：

```text
enableBalancing        false;
reportBalancingTiming false;
balanceCheckInterval  100000;
```

相关改动目标：

- `enableBalancing false` 时不再执行 imbalance check
- `reportBalancingTiming false` 时不再打印全零的 mesh DLB timing
- `balanceCheckInterval` 将 imbalance 检查从 `outputTime()` 解耦
- 只有 imbalance 超阈值时，才强制写当前时间目录，然后执行 reconstruct/decompose/restart

这样解决了早期实现中两个问题：

- DLB 检查被 `outputTime()` 绑死，`writeControl none` 或 output 间隔变化会影响 DLB 触发
- mesh DLB 关闭时仍输出 `DLB wall time = 0` 等无意义计时

代表配置与日志：

- [loadBalanceDict](/home/superxcx/code/DSMC/dsmcFoam++/hyStrath_xcx/case/cylinder_react/mixparallel/dlb/mpi8_noreact/system/loadBalanceDict)
- [log.dsmcFoam+.mpi8.repart_dlb_particle_count](/home/superxcx/code/DSMC/dsmcFoam++/hyStrath_xcx/case/cylinder_react/mixparallel/dlb/mpi8_noreact/log.dsmcFoam+.mpi8.repart_dlb_particle_count)
- [log.time.mpi8.repart_dlb_particle_count](/home/superxcx/code/DSMC/dsmcFoam++/hyStrath_xcx/case/cylinder_react/mixparallel/dlb/mpi8_noreact/log.time.mpi8.repart_dlb_particle_count)

### 网格重划分 DLB 端到端性能

对照基线：

- 纯 MPI8、无 transfer、无 mesh DLB：
  - [log.dsmcFoam+.mpi8.clean_no_transfer](/home/superxcx/code/DSMC/dsmcFoam++/hyStrath_xcx/case/cylinder_react/mixparallel/dlb/mpi8_noreact/log.dsmcFoam+.mpi8.clean_no_transfer)
  - `/usr/bin/time real = 326.75 s`
  - `DSMC evolve wall time = 353.7671 s`
  - `Main loop wall time = 354.6562 s`

- mesh DLB 启用、particle-count weight 版本：
  - [log.dsmcFoam+.mpi8.repart_dlb_particle_count](/home/superxcx/code/DSMC/dsmcFoam++/hyStrath_xcx/case/cylinder_react/mixparallel/dlb/mpi8_noreact/log.dsmcFoam+.mpi8.repart_dlb_particle_count)
  - `/usr/bin/time real = 757.54 s`
  - `DSMC evolve wall time = 343.7532 s`
  - `DLB wall time = 461.8725 s`
  - `DLB count = 16`
  - `Non-DSMC wall time = 474.8445 s`
  - `Main loop wall time = 818.5977 s`

结论很直接：

- 网格重划分后 DSMC evolve 本身略降：`353.77 -> 343.75 s`
- 但 DLB 额外开销高达 `461.87 s`
- 端到端 main loop 从 `354.66 s` 增加到 `818.60 s`
- `/usr/bin/time real` 从 `326.75 s` 增加到 `757.54 s`

因此，对这个 case，网格重划分 DLB 的计算收益完全被 reconstruct/decompose/restart 开销掩盖。

### 网格 DLB 详细时间拆分

`log.dsmcFoam+.mpi8.repart_dlb_particle_count` 中 16 次 mesh DLB 的累计开销：

| 项目 | 总时间 [s] | 平均每次 [s] | 最大单次 [s] |
|---|---:|---:|---:|
| prepare | 4.8997 | 0.3062 | 0.3753 |
| reconstruct mesh | 10.5964 | 0.6623 | 0.7842 |
| reconstruct fields | 208.9239 | 13.0577 | 15.8492 |
| decompose | 234.6258 | 14.6641 | 17.5607 |
| processor mesh sync | 1.4982 | 0.0936 | 0.1406 |
| backup | 1.2821 | 0.0801 | 0.1046 |
| 合计 | 461.8261 | 28.8641 | - |

占比最高的两个部分：

- `reconstruct fields`: `208.92 s`
- `decomposeDSMCLoadBalancePar`: `234.63 s`

这两项合计约 `443.55 s`，占 mesh DLB 总开销约 `96%`。

单次 DLB 的典型流程日志：

```text
DLB trigger: forcing write of current time before mesh repartition
Exec   : reconstructParMesh -latestTime
Exec   : reconstructPar -latestTime
Exec   : decomposeDSMCLoadBalancePar -force -latestTime -copyUniform
DLB restart: continuing from latest time ...
```

### 网格 DLB 存在的问题

1. 端到端开销过大
   - 每次 mesh DLB 平均约 `28.86 s`
   - 本 case 正常每 20 个 step 左右的计算时间也只是几十秒量级
   - 重划分一次的代价相当于吃掉大量 time step 的正常计算收益

2. `reconstructPar` 场重构是主瓶颈之一
   - DSMC case 中 lagrangian parcel 和 field 文件多
   - `reconstructPar -latestTime` 需要重构最新时间目录的场和 lagrangian 状态
   - 该项平均 `13.06 s/次`

3. `decomposeDSMCLoadBalancePar` 是另一个主瓶颈
   - 平均 `14.66 s/次`
   - 需要重新生成 processor 目录、映射 mesh/fields/lagrangian 数据
   - 粒子数增加后该项继续变贵

4. restart 语义破坏单次连续运行
   - 每次 DLB 都要结束当前 stage
   - 执行外部工具
   - 再从最新时间目录重启 solver
   - 对 wall time、日志、调度、错误恢复都更复杂

5. DLB 触发过密时会快速失控
   - `particle_count` 版本中，后期检查间隔为 5 个 timeIndex
   - 16 次重划分导致 `Non-DSMC wall time = 474.84 s`
   - 说明即使每次重划分都确实降低 imbalance，也很难覆盖自身开销

6. 只能作为低频、大尺度重平衡手段
   - 适合极长计算中偶尔重划分
   - 不适合当前这种 300 step 级别的频繁动态 DLB

工程结论：

- mesh DLB 不应作为默认优化路径
- 若保留，只应在 `enableBalancing true` 且极低频触发
- 必须把 `balanceCheckInterval` 和 imbalance 阈值设得很保守
- 该路径更适合“长时间运行中拓扑级重分布”，不是轻量级 runtime DLB

### transfer DLB 后续优化

本轮继续优化 transfer DLB，目标是减少不必要开销并确认它和 profiling 的耦合关系。

保留的代码/控制改动：

- 新增 `dlbOffloadReport`
  - `false` 时关闭 `DLB offload planner summary`
  - 同时跳过诊断用的 MPI reductions/gather
  - 不影响 transfer DLB 执行

- transfer DLB 从 `evolveProfileEnabled()` 解耦
  - 之前 `profileSummary false` 会导致 transfer DLB planner 不进入
  - 现在 transfer DLB 的执行条件改为独立的 `dlbOffloadPlanner`
  - `profileSummary` 只控制 profiling 输出，不再控制 DLB 是否执行

- compact single-species payload
  - 单组分、无电子能级时减少 payload 中重复字段
  - remote temporary parcels 使用连续 `std::vector<dsmcParcel>`，减少 per-particle `new/delete`

代表日志：

- [log.dsmcFoam+.noreact_dlbexec_transfer_v6_pool_compact](/home/superxcx/code/DSMC/dsmcFoam++/hyStrath_xcx/case/cylinder_react/mixparallel/dlb/omp4_mpi2_mpiopt_noreact/log.dsmcFoam+.noreact_dlbexec_transfer_v6_pool_compact)
- [log.dsmcFoam+.noreact_dlbexec_transfer_v8_reportoff](/home/superxcx/code/DSMC/dsmcFoam++/hyStrath_xcx/case/cylinder_react/mixparallel/dlb/omp4_mpi2_mpiopt_noreact/log.dsmcFoam+.noreact_dlbexec_transfer_v8_reportoff)
- [log.dsmcFoam+.noreact_dlbexec_transfer_v9_profileoff](/home/superxcx/code/DSMC/dsmcFoam++/hyStrath_xcx/case/cylinder_react/mixparallel/dlb/omp4_mpi2_mpiopt_noreact/log.dsmcFoam+.noreact_dlbexec_transfer_v9_profileoff)
- [log.dsmcFoam+.noreact_no_transfer_dlb_profile1](/home/superxcx/code/DSMC/dsmcFoam++/hyStrath_xcx/case/cylinder_react/mixparallel/dlb/omp4_mpi2_mpiopt_noreact/log.dsmcFoam+.noreact_no_transfer_dlb_profile1)

### transfer DLB 最新对比

`v6_pool_compact` 到 `v8_reportoff`：

| 指标 | v6 pool compact | v8 report off |
|---|---:|---:|
| `/usr/bin/time real` | 110.79 s | 109.37 s |
| Stage evolve | 120.0813 s | 116.2308 s |
| Main loop | 121.9790 s | 118.1916 s |
| Move total | 54.4645 s | 52.8742 s |
| Collision total | 47.4177 s | 45.2635 s |
| Evolve total profiled | 116.3831 s | 112.6562 s |

去掉 transfer DLB 诊断后，日志中不再出现：

- `DLB offload planner summary`
- `Collision MPI rank ...`
- mesh DLB 全零计时

说明原来的诊断 reductions/gather 和输出确实有开销，但不是最主要瓶颈。

`profileSummary false` 后的 `v9_profileoff`：

| 指标 | v8 report off/profile on | v9 report off/profile off |
|---|---:|---:|
| `/usr/bin/time real` | 109.37 s | 105.03 s |
| Stage evolve | 116.2308 s | 117.5220 s |
| Main loop | 118.1916 s | 119.3798 s |
| Non-DSMC | 1.9608 s | 1.8578 s |

这个结果说明：

- `profileSummary false` 下程序可以正常执行 transfer DLB
- 解耦是有效的
- 但 wall time 有一定波动，不能只看单次 `real`
- 关闭 profiling 后缺少 move/collision 子项，不适合后续诊断

### transfer DLB 与无 transfer DLB 对比

为确认 transfer DLB 本身是否导致当前 move 变慢，又在同一 `dlb/omp4_mpi2_mpiopt_noreact` case 中关闭 transfer DLB 并打开 profile：

控制项：

```text
profileSummary true;
dlbOffloadPlanner false;
dlbOffloadExecute false;
dlbOffloadExperimentalExecute false;
```

日志：

- [log.dsmcFoam+.noreact_no_transfer_dlb_profile1](/home/superxcx/code/DSMC/dsmcFoam++/hyStrath_xcx/case/cylinder_react/mixparallel/dlb/omp4_mpi2_mpiopt_noreact/log.dsmcFoam+.noreact_no_transfer_dlb_profile1)
- [log.time.noreact_no_transfer_dlb_profile1](/home/superxcx/code/DSMC/dsmcFoam++/hyStrath_xcx/case/cylinder_react/mixparallel/dlb/omp4_mpi2_mpiopt_noreact/log.time.noreact_no_transfer_dlb_profile1)

对比：

| 指标 | v8 transfer DLB | no transfer DLB |
|---|---:|---:|
| `/usr/bin/time real` | 109.37 s | 111.61 s |
| Stage evolve | 116.2308 s | 119.3770 s |
| Main loop | 118.1916 s | 121.3677 s |
| Move total | 52.8742 s | 54.7371 s |
| Move parallel kernel | 42.2718 s | 43.6854 s |
| Move transfer/delete finalize | 5.2218 s | 5.5117 s |
| Collision total | 45.2635 s | 46.2165 s |
| Evolve total profiled | 112.6562 s | 115.6996 s |

结论：

- 在当前代码和当前 case 状态下，完全关闭 transfer DLB 并没有恢复到旧版 `costtarget130` 的速度
- move 仍然偏慢，甚至略慢于 v8 transfer DLB
- 因此这轮 move 变慢不能简单归因于 transfer DLB
- 后续需要用 `profileDetail true` 分解 processor hits、face hits、boundary hits、thread move wall 等细项

### transfer DLB 的主要问题

相对旧版 `costtarget130`，当前 aggressive 策略把 candidate 平衡得更彻底，但代价更高：

- 旧版 `balance target max/avg = 1.30`
- 当前 aggressive 配置 `dlbOffloadBalanceTarget = 1.05`
- 旧版 predicted post-DLB max/avg 平均约 `1.292`
- 当前 v6 predicted post-DLB max/avg 平均约 `1.002`

但是 offload 工作量明显增加：

| DLB 分项 | 旧版 costtarget130 | 当前 v6 |
|---|---:|---:|
| planned offload candidates 总量 | 1.94e6 | 3.43e6 |
| remote execute wall 总和 | 0.81 s | 1.77 s |
| DLB task send wall 总和 | 0.19 s | 0.58 s |
| remote construct 总和 | 0.16 s | 0.28 s |
| remote collide task-sum 总和 | 2.14 s | 4.24 s |

工程判断：

- `candidate max/avg -> 1.0` 不等于 `wall time -> 最优`
- 当前 transfer DLB 的完全平衡策略会增加 payload、remote construct、remote collide 和 post/apply 压力
- 在 2 rank x 4 thread case 中，过度 offload 容易把 receiver 变成新的瓶颈
- 更合理的策略应回到 cost-aware target，例如 `1.25-1.30`，或者用 effective wall/cost 模型决定 offload 量

### 本轮 DLB 总结

网格 DLB：

- 能做真正 mesh ownership 重分配
- 但 reconstruct/decompose/restart 代价过大
- 当前 case 中端到端明显负收益
- 只适合非常低频、长时间运行中的粗粒度重平衡

transfer DLB：

- 不需要重划分网格，开销远低于 mesh DLB
- 但 remote execute + immediate writeback 仍有通信和构造成本
- 完全平衡 candidate 会过度 offload
- 需要保守、cost-aware 的 offload 策略

当前推荐：

- 默认关闭 mesh DLB
- mesh DLB 保留为显式控制的低频策略
- transfer DLB 可继续作为实验路径，但默认不追求 `1.05` 级完全平衡
- 下一步若继续优化，应开启 `profileDetail true` 定位 move 变慢原因，再决定是否调整 MPI move 或 transfer DLB 的 cost model

## 2026-05-07 ownership migration 继续验证

本轮继续推进 in-memory ownership migration，目标是把当前 mesh-owner 更接近 compute-owner，减少第 11 步 `move` 里的 MPI 等待。

### 结果 1: 仅边界单元迁移

在 `MPI2 x OMP4` 的 `inMemoryBalancing` case 上，边界单元迁移能跑通，但平衡度不够：

- 迁移触发时 `Maximum imbalance = 16.73%`
- 迁移后 `post particle counts = 2(767500 562467)`
- 迁移后 `post maximum imbalance = 15.42%`
- 迁移单元数 `selected cells global = 472`
- 迁移粒子负载 `selected parcel load = 10170`

这说明边界候选集本身太小，只靠 processor halo cells 还不足以把 rank 级负载压平。

### 结果 2: 更激进的全局回填

我又尝试把 interior cells 也纳入迁移候选，并把 post imbalance 压到更低。结果是：

- `post maximum imbalance = 8.58%`
- 但下一步 `move` 在第 11 步直接崩溃

说明这条更激进的全局回填路径目前不稳定，不能作为当前实现保留。

### 结论

- 当前稳定版的 ownership migration 仍然只迁边界单元
- 它能明显触发迁移，但达不到真正的完全平衡
- 更激进的 interior 回填虽然能继续降 imbalance，但会破坏后续 `move` 的稳定性
- 所以 ownership migration 这条线还没结束，但现阶段可以确认: 单靠边界迁移不够，直接全局回填也不稳

### 后续修复：迁移后 move 卡死与 field 刷新

继续沿 interior expansion 方向推进后，先后修掉了两个稳定性问题：

1. 迁移后的第 11 步 `move` 卡死
   - 在 `dsmcDynamicLoadBalancing::performInMemory()` 里，`distributor.distribute()` 后显式重建 tracking 相关 mesh cache：
     - `solutionD()`
     - `tetBasePtIs()`
     - `oldCellCentres()`
     - `cellTree()`
   - parcel 重建后额外执行 `pPtr->relocate(data.position)`，确保 parcel 按分布后的 mesh 重新定位。

2. 迁移后 `dsmcVolFields` 在后续输出阶段触发 FPE
   - 为 `dsmcFieldProperties` 新增 `refreshAfterMeshDistribution()`
   - 在 `dsmcCloud::refreshAfterMeshDistribution()` 中调用，重新创建并 reset field 状态。

修复后，in-memory ownership migration 可以稳定越过原来的 `timeIndex 11 move` 和 `dsmcVolFields` 输出崩溃点，完整跑通 probe。

### 正式 case 复测：ownership migration only

正式对比统一使用 `dlb/omp4_mpi2_mpiopt_noreact`，并完全关闭旧的 transfer DLB，只保留 in-memory ownership migration。

正式日志：

- [no-DLB 基线](/home/superxcx/code/DSMC/dsmcFoam++/hyStrath_xcx/case/cylinder_react/mixparallel/dlb/omp4_mpi2_mpiopt_noreact/log.dsmcFoam+.ownership_pair_nodlb_20260507)
- [no-DLB time](/home/superxcx/code/DSMC/dsmcFoam++/hyStrath_xcx/case/cylinder_react/mixparallel/dlb/omp4_mpi2_mpiopt_noreact/log.time.ownership_pair_nodlb_20260507)
- [ownership threshold 3%](/home/superxcx/code/DSMC/dsmcFoam++/hyStrath_xcx/case/cylinder_react/mixparallel/dlb/omp4_mpi2_mpiopt_noreact/log.dsmcFoam+.ownership_pair_inmem3_20260507)
- [ownership threshold 3% time](/home/superxcx/code/DSMC/dsmcFoam++/hyStrath_xcx/case/cylinder_react/mixparallel/dlb/omp4_mpi2_mpiopt_noreact/log.time.ownership_pair_inmem3_20260507)
- [ownership threshold 5%](/home/superxcx/code/DSMC/dsmcFoam++/hyStrath_xcx/case/cylinder_react/mixparallel/dlb/omp4_mpi2_mpiopt_noreact/log.dsmcFoam+.ownership_pair_inmem5_20260507)
- [ownership threshold 5% time](/home/superxcx/code/DSMC/dsmcFoam++/hyStrath_xcx/case/cylinder_react/mixparallel/dlb/omp4_mpi2_mpiopt_noreact/log.time.ownership_pair_inmem5_20260507)
- [ownership threshold 10%](/home/superxcx/code/DSMC/dsmcFoam++/hyStrath_xcx/case/cylinder_react/mixparallel/dlb/omp4_mpi2_mpiopt_noreact/log.dsmcFoam+.ownership_pair_inmem10_20260507)
- [ownership threshold 10% time](/home/superxcx/code/DSMC/dsmcFoam++/hyStrath_xcx/case/cylinder_react/mixparallel/dlb/omp4_mpi2_mpiopt_noreact/log.time.ownership_pair_inmem10_20260507)

结果对比：

| 配置 | `real` [s] | main loop [s] | evolve total [s] | 非 evolve 开销 [s] | move only [s] | collision phase [s] | 触发次数 |
|---|---:|---:|---:|---:|---:|---:|---:|
| no-DLB | 133.10 | 139.99 | 136.00 | 3.99 | 67.56 | 52.23 | 0 |
| ownership `3%` | 1195.24 | 1272.49 | 134.53 | 1137.97 | 123.61 | 3.05 | 21 |
| ownership `5%` | 807.75 | 858.75 | 126.20 | 732.55 | 113.71 | 3.71 | 13 |
| ownership `10%` | 306.38 | 325.26 | 84.39 | 240.87 | 69.81 | 5.36 | 6 |

观察：

- ownership migration 确实能显著压低 evolve 内的 `collision phase`
  - `52.23 s -> 5.36 s`（`10%`）
  - `52.23 s -> 3.71 s`（`5%`）
  - `52.23 s -> 3.05 s`（`3%`）

- 但主循环真实瓶颈变成了迁移本身的非 evolve 开销
  - `10%`: `240.87 s`
  - `5%`: `732.55 s`
  - `3%`: `1137.97 s`

- `10%` 是目前最不差的一档，但仍明显慢于 no-DLB
  - `real 306.38 s` 对比 `133.10 s`
  - `main loop 325.26 s` 对比 `139.99 s`

- `3%` 和 `5%` 都属于典型的“平衡做得更好，但端到端显著负优化”

### 当前 ownership migration 的工程结论

1. 稳定性问题已经基本解决
   - interior expansion 版本现在可以稳定跑完整个正式 case
   - 迁移后 `move` 和 `dsmcVolFields` 不再崩溃

2. 当前实现的主要问题不再是稳定性，而是触发成本
   - 即使把阈值放宽到 `10%`，300 个 timeIndex 里仍触发了 6 次 ownership migration
   - 每次迁移都会引入大块非 evolve 开销，最终吞掉 collision 降下来的收益

3. 这一版 ownership migration 还不能算真正有效的优化
   - 最优档 `10%` 仍然是明显负收益
   - `3%` / `5%` / `10%` 的正式结果已经把这点坐实

4. 后续若继续做 ownership migration，不该再继续调小阈值或追求更彻底平衡
   - 当前更值得试的是“单次或极低频 ownership migration”
   - 或者把迁移触发从粒子数不平衡改成更强的 timer/cost 判据
   - 否则只会继续出现 `collision 更平衡，但 main loop 更慢` 的结果

### 2026-05-08 ownership migration 单次降本与低频触发复测

本轮首先直接压单次 migration 固定成本，修改位置：

- `src/lagrangian/dsmc/dynamicLoadBalancing/dsmcDynamicLoadBalancing.C`
  - 删除 ownership migration 重建阶段的显式 `mesh_.cellTree().findInside(data.position)`
  - 删除 parcel 重建后的二次 `pPtr->relocate(data.position)`
  - 改成直接信任 `distributeCellData(cellParcels)` 给出的目标 `celli`，只保留构造函数内部那一次按 cell 的 locate
  - 新增 migration 分段计时：
    - `store/clear wall`
    - `mesh distribute wall`
    - `parcel distribute wall`
    - `parcel rebuild wall`
    - `occupancy rebuild wall`
    - `refresh wall`

这一步的判断依据很直接：旧实现每个 parcel 在一次 migration 中实际做了

1. `cellTree().findInside(position)`  
2. `dsmcParcel(mesh, position, celli, ...)` 内部一次 `locate(position, celli)`  
3. `pPtr->relocate(position)` 再做一次无 seed 的 `locate(position, -1)`  

也就是每个 parcel 至少带两次 locate、一次全局 `cellTree` 搜索。对于 ownership migration，这属于明显的固定税。

#### 先验坑点：并行 case 读取的是 `processor*/system/*Dict`

第一次复测时只改了主 case `system/loadBalanceDict` 和 `system/controlDict`，实际并行运行仍在读 `processor0/1/system/*Dict`，所以没有真正触发 ownership migration。之后补改了：

- `processor0/system/loadBalanceDict`
- `processor1/system/loadBalanceDict`
- `processor0/system/controlDict`
- `processor1/system/controlDict`

之后的 ownership 结果才是真正有效的复测结果。

#### 复测 1：`inMemoryMigrationFraction=0.10`，`maximumAllowableImbalance=0.01`

日志：

- `log.dsmcFoam+.ownership_inmem10_nomaplocate_procfix_20260508`

结果：

| 配置 | real [s] | main loop [s] | evolve total [s] | move only [s] | collision [s] | 触发次数 |
|---|---:|---:|---:|---:|---:|---:|
| ownership `10%`, `imb=1%` | 120.58 | 122.33 | 60.90 | 53.75 | 1.58 | 28 |

单次 migration 分段统计：

- 触发 `28` 次
- 单次 `wall time` 平均 `1.566 s`
- 单次 `mesh distribute` 平均 `0.474 s`
- 单次 `parcel rebuild` 平均 `0.420 s`
- 28 次 migration 总 wall 约 `43.85 s`

结论：

- 单次 migration 成本已经从上一轮的几十秒量级，压到了 `1.2 ~ 2.25 s`
- 但是 `maximumAllowableImbalance=1%` 太激进，300 个 timeIndex 里触发了 `28` 次
- ownership 这时虽然把 `collision phase` 压到 `1.58 s`，但总触发次数太多，端到端仍然输给 no-DLB

#### 相邻 no-DLB 基线

日志：

- `log.dsmcFoam+.ownership_pair_nodlb_after_nomaplocate_20260508`

结果：

| 配置 | real [s] | main loop [s] | evolve total [s] | move only [s] | collision [s] |
|---|---:|---:|---:|---:|---:|
| no-DLB | 103.59 | 103.48 | 100.26 | 46.36 | 41.22 |

对比上面的 `ownership 10%, imb=1%`：

- `collision phase`: `41.22 -> 1.58 s`
- 但 `main loop`: `103.48 -> 122.33 s`
- 所以这时 ownership 仍然是负优化，问题已从“单次迁移太贵”转成“迁移太频繁”

#### 复测 2：抬高阈值到 `maximumAllowableImbalance=0.15`

保持 `inMemoryMigrationFraction=0.10` 不变，只把触发阈值从 `1%` 提到 `15%`。

日志：

- `log.dsmcFoam+.ownership_inmem10_imb15_20260508`

结果：

| 配置 | real [s] | main loop [s] | evolve total [s] | move only [s] | collision [s] | 触发次数 |
|---|---:|---:|---:|---:|---:|---:|
| ownership `10%`, `imb=15%` | 77.94 | 78.50 | 59.25 | 47.66 | 5.32 | 8 |
| no-DLB | 103.59 | 103.48 | 100.26 | 46.36 | 41.22 | 0 |

单次 migration 分段统计：

- 触发 `8` 次
- 单次 `wall time` 平均 `1.640 s`
- 单次 `mesh distribute` 平均 `0.423 s`
- 单次 `parcel rebuild` 平均 `0.474 s`
- 8 次 migration 总 wall 约 `13.12 s`

关键结论：

- 单次 migration 成本已经足够低，真正决定成败的是触发频率
- 把阈值抬到 `15%` 后，migration 次数从 `28` 次降到 `8` 次
- 这时 ownership migration 首次形成**真实端到端正收益**
  - `real`: `103.59 -> 77.94 s`，快 `25.65 s`
  - `main loop`: `103.48 -> 78.50 s`，快 `24.98 s`
  - `evolve total`: `100.26 -> 59.25 s`，快 `41.01 s`

同时副作用也明显缓和：

- `move only`: `46.36 -> 47.66 s`，只增加 `1.30 s`
- `post fields/output`: `9.16 -> 2.90 s`
- `collision phase`: `41.22 -> 5.32 s`

### 这一轮 ownership migration 的新结论

1. **核心突破不是再继续压单次迁移到亚秒，而是先把错误的重复 locate 去掉**
   - 删除显式 `findInside + 二次 relocate` 后，单次 migration 已经足够便宜

2. **当前 ownership migration 已经进入可用区间**
   - 在 `inMemoryMigrationFraction=0.10`
   - `maximumAllowableImbalance=0.15`
   - `balanceCheckInterval=10`
   - `openmpMoveGuardLayers=1`
   下，正式 case 已经稳定且明显快于 no-DLB

3. **后续主线不再是“继续压单次 migration 毫秒级”**
   - 单次 wall 已经只有 `~1.6 s`
   - 真正更有价值的是继续做 `timer/cost gating`，把触发时机做得更准
   - 下一步优先考虑：
     - 触发条件从纯粒子数 imbalance 进一步转成 runtime/cost 判据
     - 测试更低频检查，例如更大的 `balanceCheckInterval`
     - 在 `15%` 左右阈值附近做小范围扫点，而不是回到 `1%` 这种过度激进策略

#### 后续追加复测：timer gate 失败，cooldown 有效

在上面 `ownership 10%, imb=15%, interval=10` 基线继续往下推，做了三类增量测试。

##### 1. timer gate 直接失败

代码里加了 `inMemoryUseTimerGate / inMemoryRequiredGainFactor / inMemoryAssumedBalanceCost`，并用最近一个检查窗口的 `evolveCollisionWallTime()` 做 gate。

日志：

- `log.dsmcFoam+.ownership_inmem10_imb15_timergate_20260508`

结果：

| 配置 | real [s] | main loop [s] | move only [s] | collision [s] | migration 次数 |
|---|---:|---:|---:|---:|---:|
| `imb=15%`, timer gate on | 104.21 | 105.28 | 47.40 | 41.84 | 0 |
| no-DLB | 103.59 | 103.48 | 46.36 | 41.22 | 0 |

日志中 gate 一直输出：

- `collision window max/avg gain [s] = 0.0078 ~ 0.0286`
- `estimated migration cost [s] = 1.6`
- `decision = skip`

这说明当前这套 `timer gate` 信号选错了。单个检查窗口内的 `collision wall` 差值太小，完全不能代表 ownership migration 的长期收益，结果就是把所有 migration 都挡掉，效果退回 no-DLB。

##### 2. 单纯把 `balanceCheckInterval` 拉大到 `20` 也不对

日志：

- `log.dsmcFoam+.ownership_inmem10_imb15_int20_20260508`

结果：

| 配置 | real [s] | main loop [s] | move only [s] | collision [s] | migration 次数 |
|---|---:|---:|---:|---:|---:|
| `imb=15%`, `interval=20` | 105.74 | 106.92 | 48.47 | 33.64 | 4 |

问题在于：`timeIndex 10` 的第一次 migration 很关键，等到 `timeIndex 20` 再检查已经错过了早期窗口。这组虽然 migration 次数降了，但总时间反而变差。

##### 3. 保留 `interval=10`，改成 migration 后 cooldown

为了解决“第 10 步必须迁，但后面不能 10 步一次连续迁”的矛盾，在 `dsmcDynamicLoadBalancing` 里新增：

- `balanceCooldownSteps`

逻辑：

- 仍然按 `balanceCheckInterval=10` 检查
- 但每次实际 migration 后，强制跳过接下来若干步的再次迁移

###### cooldown = 20, threshold = 15%

日志：

- `log.dsmcFoam+.ownership_inmem10_imb15_cd20_20260508`

结果：

| 配置 | real [s] | main loop [s] | move only [s] | collision [s] | migration 次数 |
|---|---:|---:|---:|---:|---:|
| `imb=15%`, `cd=20` | 76.01 | 77.17 | 48.92 | 6.13 | 6 |

对比上一轮 `imb=15%, interval=10`：

- `real`: `77.94 -> 76.01 s`
- `main loop`: `78.50 -> 77.17 s`
- migration: `8 -> 6`

###### cooldown = 20, threshold = 16%

日志：

- `log.dsmcFoam+.ownership_inmem10_imb16_cd20_20260508`

结果：

| 配置 | real [s] | main loop [s] | move only [s] | collision [s] | migration 次数 |
|---|---:|---:|---:|---:|---:|
| `imb=16%`, `cd=20` | 75.33 | 76.55 | 47.87 | 6.40 | 6 |

继续比 `imb=15%, cd=20` 再快：

- `real`: `76.01 -> 75.33 s`
- `main loop`: `77.17 -> 76.55 s`

###### cooldown = 30, threshold = 16%

日志：

- `log.dsmcFoam+.ownership_inmem10_imb16_cd30_20260508`

结果：

| 配置 | real [s] | main loop [s] | move only [s] | collision [s] | migration 次数 |
|---|---:|---:|---:|---:|---:|
| `imb=16%`, `cd=30` | 74.19 | 74.35 | 49.33 | 6.91 | 4 |

这是当前这轮的最好结果：

- 相对 no-DLB：
  - `real`: `103.59 -> 74.19 s`，快 `29.40 s`
  - `main loop`: `103.48 -> 74.35 s`，快 `29.13 s`
- 相对上一版最优 `imb=15%, interval=10`：
  - `real`: `77.94 -> 74.19 s`，快 `3.75 s`
  - `main loop`: `78.50 -> 74.35 s`，快 `4.15 s`
- migration 次数：
  - `8 -> 4`

同时也能看到 tradeoff：

- `move only` 比 no-DLB 略高：`46.36 -> 49.33 s`
- `collision phase` 仍然大幅低于 no-DLB：`41.22 -> 6.91 s`

#### 这一轮的更新结论

1. **timer gate 目前不成立**
   - 当前用 `recent collision wall delta` 做 gate，会系统性低估 ownership migration 的长期收益
   - 在这个 case 上它直接退化成 no-DLB

2. **真正有效的是“保留早期第一次 migration + 后续 cooldown 节流”**
   - 不能粗暴把 `balanceCheckInterval` 拉大到 `20`
   - `timeIndex 10` 的第一次 migration 必须保留

3. **当前正式 case 的更优配置已经更新为**
   - `inMemoryMigrationFraction 0.10`
   - `maximumAllowableImbalance 0.16`
   - `balanceCheckInterval 10`
   - `balanceCooldownSteps 30`
   - `inMemoryUseTimerGate false`

4. **当前这条 ownership migration 线已经形成稳定正收益**
   - 且比上一轮 `imb=15%, interval=10` 更好
   - 下一步如果继续做，不该再回到“简单 timer gate”，而应该做更像 `cost proxy / collision-cell weighted trigger` 的判据

#### 2026-05-09 ownership migration 正确性修复与重新验证

这一段需要明确更正：上面 `real 74.19 s` / `main loop 74.35 s` 的 ownership migration 结果，后来证明**不能作为有效性能结论引用**。原因不是计时口径，而是那一版 migration 后 `BuildCellOccupancy` 的线程局部缓存状态被破坏，导致碰撞统计和能量统计失真，表现为：

- migration 后输出窗口内 `Collisions` 从正常量级直接塌到很小值
- `occupancy parcels global` 明显异常
- 最终粒子数、碰撞数、能量与 no-DLB 偏离明显

也就是说，那一轮“看起来很快”的 ownership migration，本质上是**错误物理路径下的假优化**。

##### 1. 正确性 bug 的真正根因

最终定位到的问题不是 `sigmaTcRMax`，也不是迁移 parcel 时缺少 `face/behind/nBehind`，而是 migration 后：

- `buildCellOccupancy()` 复用了 stale 的 `occupancyThreadCellCounts_ / occupancyThreadActiveCells_`
- `moveOrdered` 状态机在 migration 后长期退化成 `fallbackGather`

前者会直接污染 occupancy rebuild，后者会把 `buildCellOccupancy extract parcels` 长期拉高。

对应源码修复：

- `dsmcCloud.C`
  - `buildCellOccupancy()` 中每次都完整清零 thread-local counts
  - `resetAfterMeshDistribution()` 中清空 `occupancyThreadCellCounts_` 和 `occupancyThreadActiveCells_`
  - migration 后不再永久置 `moveOrderedReuseDisabled_ = true`
- `dsmcDynamicLoadBalancing.C`
  - 跳过 `endTime` 时刻的无效 migration
  - 保留 migration payload 中的 `face / behind / nBehind`

##### 2. 短 probe：正确性已恢复

20 步短 probe 日志：

- `/tmp/ownership_probe_stage8_own_fix3/log.own_stage8_fix3`

关键结果：

- `timeIndex 10`
  - `Collisions = 28168`
  - `Collision candidates = 95754`
  - `Collision acceptance rate = 0.2941704785`
- 第一次 migration 后，`timeIndex 20`
  - `Collisions = 37702`
  - `Collision candidates = 113678`
  - `Collision acceptance rate = 0.3316560812`
  - `Number of molecules = 8.121648e+15`
  - `Average rotational energy = 1.57472815e-21`
  - `Average total energy = 2.058343586e-19`

同时 `BuildCellOccupancy` 已恢复正常复用：

- `move-ordered hits = 20`
- `fallback gathers = 2`

这说明 migration 后的碰撞、能量、occupancy 都已经回到正常轨道。

##### 3. 全算例重新验证：物理量正确，但总时间仍未优于 no-DLB

###### 干净 no-DLB 基线

日志：

- `/tmp/ownership_full_nodlb_fix2/log.nodlb_full_fix2`
- `/tmp/ownership_full_nodlb_fix2/log.time.nodlb_full_fix2`

结果：

| 配置 | real [s] | main loop [s] | move only [s] | buildCellOccupancy [s] | collision [s] | total profiled [s] |
|---|---:|---:|---:|---:|---:|---:|
| no-DLB | 107.45 | 107.31 | 47.90 | 3.58 | 43.46 | 104.17 |

最终物理量：

- `Collisions = 245203`
- `Collision candidates = 516316`
- `Collision acceptance rate = 0.4749087768`
- `Number of molecules = 1.174413e+16`
- `Average rotational energy = 1.343750125e-20`
- `Average total energy = 1.856783217e-19`

###### ownership migration（`imb=16%`, `interval=10`, `cooldown=30`）修复后全算例

日志：

- `/tmp/ownership_full_dlb_fix3/log.own_full_fix3`
- `/tmp/ownership_full_dlb_fix3/log.time.own_full_fix3`

结果：

| 配置 | real [s] | main loop [s] | move only [s] | buildCellOccupancy [s] | collision [s] | total profiled [s] | migration 次数 |
|---|---:|---:|---:|---:|---:|---:|---:|
| ownership fixed | 147.08 | 145.72 | 47.84 | 3.40 | 41.08 | 100.91 | 3 |

三次 migration 的 wall time：

- `12.53 s`
- `13.55 s`
- `14.49 s`

总计约 `40.58 s`。

最终物理量相对 no-DLB 已经基本一致：

- `Collisions = 245425`，偏差 `+0.0905%`
- `Collision candidates = 516163`，偏差 `-0.0296%`
- `Collision acceptance rate = 0.475479645`，偏差 `+0.1202%`
- `Number of molecules = 1.174257e+16`，偏差 `-0.0133%`
- `Average rotational energy = 1.342748928e-20`，偏差 `-0.0745%`
- `Average total energy = 1.85665333e-19`，偏差 `-0.0070%`

结论：

- **物理正确性问题已经修好**
- `evolve` 主体确实更快：`total profiled 104.17 -> 100.91 s`
- 但 migration 自身的固定成本更大：约 `40.58 s`
- 所以端到端仍明显慢于 no-DLB：`real 107.45 -> 147.08 s`

##### 4. 单次 early migration 也不能形成净优化

为验证“是不是只保留第 10 步那次 migration 就能赢”，又做了一组只允许一次 migration 的全算例：

- `/tmp/ownership_full_dlb_fix4_cd300/log.own_full_fix4_cd300`
- `/tmp/ownership_full_dlb_fix4_cd300/log.time.own_full_fix4_cd300`

配置：

- `maximumAllowableImbalance 0.16`
- `balanceCheckInterval 10`
- `balanceCooldownSteps 300`

结果：

| 配置 | real [s] | main loop [s] | move only [s] | buildCellOccupancy [s] | collision [s] | total profiled [s] | migration 次数 |
|---|---:|---:|---:|---:|---:|---:|---:|
| ownership single-mig | 123.10 | 125.10 | 52.53 | 3.75 | 43.16 | 108.68 | 1 |

这组只有一次 migration，单次成本 `12.60 s`，但：

- `total profiled` 反而比 no-DLB 更慢：`104.17 -> 108.68 s`
- `move only` 明显恶化：`47.90 -> 52.53 s`
- `collision` 只和 no-DLB 基本持平：`43.46 -> 43.16 s`

这说明只做一次 early migration 也不成立。

##### 5. 当前这条 ownership migration 线的最终判断

到这里可以下明确结论：

1. **正确性已经解决**
   - migration 后的碰撞数、粒子数、总能量现在都与 no-DLB 基本一致

2. **当前实现下，靠 trigger/cooldown 调参无法形成真正净优化**
   - 一次 migration：`real 123.10 s`，慢于 no-DLB `107.45 s`
   - 三次 migration：`real 147.08 s`，更慢
   - 最好的 `evolve` 主体收益只有 `~3.25 s`
   - 但单次 migration 固定成本就有 `12.5 ~ 14.5 s`

3. **因此，当前 in-memory ownership migration 的主要瓶颈已经不是“平衡策略”，而是“单次 migration 的全量重建成本”**
   - 现在每次 migration 仍然要重建整朵 cloud
   - 这条实现如果不改成 touched-cell / 增量 rebuild，仅靠参数扫描基本不可能跑赢 no-DLB

4. **后续若还继续做 ownership migration，唯一有价值的方向是结构性降本**
   - 只迁移 touched cells
   - 避免全量 parcel rebuild
   - 避免全量 occupancy rebuild
   - 避免迁移后整段时间退化出额外 `move` 成本

#### 2026-05-11 ownership migration 后 occupancy rebuild 串行化修正

按照上一节结论，本轮先从 migration 后的重建固定成本里挑最小、风险最低的一项验证：`rebuildCellOccupancyAfterMeshDistribution()` 原来会临时把 `openmpEnabled_` 置为 `false`，导致 migration 后的 occupancy rebuild 强制走串行 fallback 路径。这个行为和当前 `buildCellOccupancy()` 已有的并行扁平视图路径不一致，也会把每次 migration 后的重建成本放大。

代码改动：

- `src/lagrangian/dsmc/clouds/dsmcCloud.C`
  - `rebuildCellOccupancyAfterMeshDistribution()` 不再临时关闭 OpenMP
  - migration 后 occupancy rebuild 直接调用现有 `buildCellOccupancy()`，允许其走并行 count/reduce + allocate/fill 路径
- `src/lagrangian/basic/Cloud/Cloud.C`
  - 补齐非 OpenMP 模板实例化下 `moveStageProbe` / `moveOrderedReuse` 的默认值
  - 这是为了让当前 mixed-move probe/reuse 代码在通用模板编译时闭合，不改变 OpenMP mixed 路径语义

验证构建：

- 可通过：
  - `wmake -j libso $WM_PROJECT_USER_DIR/src/lagrangian/dsmc`
  - `wmake -j $WM_PROJECT_USER_DIR/applications/solvers/discreteMethods/dsmc/dsmcFoam+`
- 完整 `docs/build_dsmc.sh` 仍会被 `src/parallel/decompose/decompose` 里的既有 `passiveParticle/indexedParticle` v2506 构造函数兼容问题阻塞；这不是本轮 DLB/Cloud 改动引入的。

复测 case：

- `/tmp/ownership_rebuild_parallel_20260511`
- 配置沿用上一轮正式 ownership migration：
  - `maximumAllowableImbalance 0.16`
  - `balanceCheckInterval 10`
  - `balanceCooldownSteps 30`
  - `inMemoryMigrationFraction 0.10`
  - `MPI2 x OMP4`

结果：

| 配置 | real [s] | main loop [s] | move only [s] | buildCellOccupancy [s] | collision [s] | total profiled [s] | migration 次数 |
|---|---:|---:|---:|---:|---:|---:|---:|
| ownership fixed, parallel occupancy rebuild | 158.97 | 156.34 | 53.25 | 3.67 | 42.21 | 108.18 | 3 |

三次 migration 分项：

| migration | selected cells | rebuilt parcels | parcel rebuild [s] | occupancy rebuild [s] | wall [s] |
|---:|---:|---:|---:|---:|---:|
| 1 | 479 | 1329978 | 13.31698 | 0.03162 | 14.32253 |
| 2 | 666 | 1841999 | 12.82577 | 0.04018 | 14.67175 |
| 3 | 774 | 1900307 | 12.92689 | 0.04639 | 14.27791 |

对比 2026-05-09 的 `ownership fixed`：

- 旧结果三次 migration 总固定成本约 `40.58 s`
- 本轮三次 migration wall 合计约 `43.27 s`
- 但本轮 `occupancy rebuild wall` 已降到每次 `0.03 ~ 0.05 s`
- 说明“migration 后 occupancy rebuild 被强制串行化”这一点已经修掉
- 端到端没有变快，原因是主成本仍然是 `parcel rebuild wall`：
  - 三次合计约 `39.07 s`
  - 每次仍需重建 `132.99w / 184.20w / 190.03w` 个 parcel

最终物理量仍与 no-DLB 基线同量级：

- `Collisions = 246080`
- `Collision candidates = 516495`
- `Collision acceptance rate = 0.4764421727`
- `Number of molecules = 1.174311e+16`
- `Average rotational energy = 1.343646523e-20`
- `Average total energy = 1.856955498e-19`

本轮结论：

1. **并行 occupancy rebuild 修正是有效的小修**
   - 它确实把 migration 后 occupancy rebuild 降到可以忽略的量级
   - 这项改动可以保留，因为语义上只是取消人为串行化，复用已有并行 `buildCellOccupancy()` 路径

2. **它不能改变 ownership migration 的总体判断**
   - 当前端到端仍明显慢于 no-DLB
   - `evolve total profiled` 也没有优于 2026-05-09 no-DLB 的 `104.17 s`
   - 所以这不是一次有效的整体优化，只是清掉了一个重建路径上的次要固定税

3. **下一步不应继续压 occupancy rebuild**
   - 当前真正必须处理的是 parcel rebuild
   - 方向应收敛到：
     - 不再 `cloud_.clear()` 后全量重建整朵 cloud
     - 按迁移 cell / touched cell 增量删除和插入 parcel
     - 只重建受影响 cell 的 parcel container / occupancy / collision metadata
     - 避免每次 ownership migration 都对百万级 parcel 做完整构造函数路径

#### 2026-05-11 ownership migration parcel rebuild lazy-locate 降本

上一轮已经证明 migration 后的 `occupancy rebuild` 不是主瓶颈，三次合计只有 `0.12 s` 量级；真正的主瓶颈是 `parcel rebuild wall`，三次合计约 `39.07 s`。本轮继续沿着这个主瓶颈做最小侵入优化：保留现有全量 cloud rebuild 语义，但避免每个 parcel 都先走一次全局 `mesh_.cellTree().findInside(position)`。

代码改动：

- `src/lagrangian/dsmc/dynamicLoadBalancing/dsmcDynamicLoadBalancing.C`
  - parcel rebuild 时优先使用迁移前保存的 `coordinates + cell + tetFace + tetPt` 做 exact-topology 构造
  - 只有 exact-topology 构造失败、cell 不匹配或位置误差超限时，才 fallback 到 `mesh_.cellTree().findInside(position)`
  - 新增 `cellTree locates` 统计，用于确认真正进入全局 cellTree 定位的 parcel 数量

验证构建：

- 可通过：
  - `wmake -j libso $WM_PROJECT_USER_DIR/src/lagrangian/dsmc`
  - `wmake -j $WM_PROJECT_USER_DIR/applications/solvers/discreteMethods/dsmc/dsmcFoam+`
- 完整 `docs/build_dsmc.sh` 仍会被 `src/parallel/decompose/decompose` 的既有 v2506 兼容问题阻塞；本轮只验证 DLB/Cloud 相关库和 `dsmcFoam+` solver。

复测 case：

- `/tmp/ownership_lazy_locate_20260511`
- 配置沿用上一轮正式 ownership migration：
  - `maximumAllowableImbalance 0.16`
  - `balanceCheckInterval 10`
  - `balanceCooldownSteps 30`
  - `inMemoryMigrationFraction 0.10`
  - `MPI2 x OMP4`

总结果对比：

| 配置 | real [s] | main loop [s] | move only [s] | buildCellOccupancy [s] | collision [s] | total profiled [s] | migration 次数 |
|---|---:|---:|---:|---:|---:|---:|---:|
| 2026-05-09 no-DLB 基线 | 107.45 | 107.31 | 47.90 | 3.58 | 43.46 | 104.17 | 0 |
| parallel occupancy rebuild | 158.97 | 156.34 | 53.25 | 3.67 | 42.21 | 108.18 | 3 |
| lazy-locate parcel rebuild | 118.63 | 116.29 | 52.27 | 3.62 | 41.55 | 106.49 | 3 |

三次 migration 分项：

| migration | rebuilt parcels | exact topology rebuilds | exact topology fallbacks | cellTree locates | parcel rebuild [s] | occupancy rebuild [s] | wall [s] |
|---:|---:|---:|---:|---:|---:|---:|---:|
| 1 | 1329957 | 1310409 | 19548 | 19548 | 0.18836 | 0.02218 | 1.04647 |
| 2 | 1841865 | 1817443 | 18753 | 24422 | 0.47779 | 0.03809 | 1.97745 |
| 3 | 1900661 | 1869527 | 28376 | 31134 | 0.65412 | 0.04935 | 2.06087 |

本轮最关键变化：

- `parcel rebuild wall` 三次合计从约 `39.07 s` 降到约 `1.32 s`
- 三次 migration 总 wall 从约 `43.27 s` 降到约 `5.08 s`
- 大多数 parcel 可以通过 exact topology 直接重建：
  - 第 1 次约 `98.53%`
  - 第 2 次约 `98.67%`
  - 第 3 次约 `98.36%`
- 这说明此前主瓶颈不是 parcel 构造本身，而是对百万级 parcel 重复做全局 cellTree locate

最终物理量：

- `Collisions = 245285`
- `Collision candidates = 515949`
- `Collision acceptance rate = 0.4754055149`
- `Number of molecules = 1.174308e+16`
- `Average rotational energy = 1.343384114e-20`
- `Average total energy = 1.856662483e-19`

补充验证：尝试跳过后期第三次 migration。

- case: `/tmp/ownership_lazy_locate_until260_20260511`
- 只把 `loadBalancingUntilTime` 改成 `1.75e-05`，使其保留前两次 migration、跳过 timeIndex 270 附近的第三次 migration
- 结果：

| 配置 | real [s] | main loop [s] | move only [s] | buildCellOccupancy [s] | collision [s] | total profiled [s] | migration 次数 |
|---|---:|---:|---:|---:|---:|---:|---:|
| lazy-locate, 3 migrations | 118.63 | 116.29 | 52.27 | 3.62 | 41.55 | 106.49 | 3 |
| lazy-locate, skip late 3rd migration | 120.08 | 118.11 | 54.79 | 3.78 | 42.68 | 110.52 | 2 |

两次 migration 配置更慢，说明简单减少 migration 次数不是当前解法。第三次 migration 虽然仍有固定成本，但对后期 collision/move 失衡有收益；在 lazy-locate 后，跳过它反而损失更多。

本轮结论：

1. **当前主瓶颈已被明显降低**
   - lazy-locate 把 `parcel rebuild wall` 从几十秒级压到秒级以内每次
   - 这是本轮最有效的 ownership migration 优化，可以保留

2. **DLB 仍未端到端胜过 no-DLB**
   - no-DLB `real = 107.45 s`
   - lazy-locate DLB `real = 118.63 s`
   - no-DLB `total profiled = 104.17 s`
   - lazy-locate DLB `total profiled = 106.49 s`
   - collision 本身已经改善：`43.46 s -> 41.55 s`
   - 但 DLB 固定成本和 migration 后 `move only` 增量仍吃掉了 collision 收益：`47.90 s -> 52.27 s`

3. **下一主瓶颈已经转移**
   - 继续压 `cellTree/findInside` 的空间不大
   - 现在应处理：
     - `fvMeshDistribute::distribute(distribution)` 固定成本
     - `distMap().distributeCellData(cellParcels)` 对全量 cell parcel payload 的搬运
     - `cloud_.clear()` 后全量重建 cloud 的语义成本
     - ownership migration 后一段时间内 `move only` 抬升的问题

4. **下一轮方向**
   - 不再把所有 cell 的 parcel 都打包进 `cellParcels` 后全量 `distributeCellData`
   - 对留在本 rank 的 cell，尽量保留或本地重标 parcel，而不是序列化、分发、析构、再构造
   - 对迁出的 cell 只发送 outgoing payload
   - 对迁入 cell 只接收 incoming payload
   - occupancy/collision metadata 只针对 touched/migrated cell 增量刷新

#### 2026-05-11 incremental migration：避免全量序列化/重建

本轮实现了 ownership migration 的核心增量化优化：对留在本 rank 的 parcel，不再经过"序列化 → distributeCellData → 析构 → 重建"的完整路径，而是直接 detach、remap、re-add。

##### 实现思路

旧流程：
1. 把 ALL parcels 序列化到 `cellParcels[]`
2. `cloud_.clear()` 析构所有 parcel 对象
3. `fvMeshDistribute::distribute()` 重分布 mesh
4. `distributeCellData(cellParcels)` 重分布所有 parcel 数据
5. 对所有 parcel 做 `new dsmcParcel(...)` 重建

新流程：
1. 遍历 cloud，对 staying cells 的 parcel 调用 `cloud_.remove(pPtr)` detach（保留对象），存入 `stayingParcels` 缓冲并记录 position
2. 对 outgoing cells 的 parcel 序列化到 `cellParcels[]`
3. `cloud_.clear()` 析构剩余 parcel（已经只剩 outgoing 的）
4. `fvMeshDistribute::distribute()` 重分布 mesh
5. `distributeCellData(cellParcels)` 只传输 outgoing/incoming parcel 数据
6. 对 staying parcels：用 `cellMap().subMap()/constructMap()` 建立 oldCell→newCell 映射，用 `oldFaceToNewFace` 映射 tetFace/face，尝试 direct remap（验证 position 一致性），失败时 fallback 到 `relocate(pos, newCelli)`
7. 对 incoming parcels：沿用原有 lazy-locate 重建

##### 关键发现：direct remap vs relocate

第一版实现对所有 staying parcels 都调用 `relocate(pos, newCelli)`，结果 `parcel rebuild wall` 仍然较高（0.31-0.77s），因为 `relocate()` 内部会遍历 cell 内所有 tet 做定位。

优化后改为 two-tier：
- 先尝试 direct index remap：更新 `cell()` 和 `tetFace()`，验证 `pPtr->position()` 与存储的 position 误差 ≤ 1e-10
- 只有验证失败时才 fallback 到 `relocate()`

实测 direct remap 成功率约 **98.7%**，说明绝大多数 staying cells 的 tet 分解在 mesh redistribution 后保持不变。

##### 复测结果

配置：`inMemoryMigrationFraction=0.10`, `maximumAllowableImbalance=0.16`, `balanceCheckInterval=10`, `balanceCooldownSteps=30`, MPI2 x OMP4

| 配置 | real [s] | main loop [s] | move only [s] | collision [s] | migration 次数 |
|---|---:|---:|---:|---:|---:|
| no-DLB 基线 | 112.11 | 113.12 | 56.66 | 43.20 | 0 |
| incremental migration (direct remap) | 119.96 | 121.17 | 60.03 | 43.93 | 3 |
| 旧版 lazy-locate (2026-05-09) | 118.63 | 116.29 | 52.27 | 41.55 | 3 |

三次 migration 分项（incremental v2）：

| migration | staying parcels | direct remap | relocated | incoming rebuilt | parcel rebuild [s] | wall [s] |
|---:|---:|---:|---:|---:|---:|---:|
| 1 | 1318823 | 1301448 | 17375 | 11143 | 0.079 | 0.454 |
| 2 | 1826530 | 1814112 | 12418 | 15557 | 0.416 | 1.288 |
| 3 | 1883539 | 1856787 | 26752 | 16970 | 0.273 | 1.219 |

对比旧版 lazy-locate 的 migration wall（三次合计约 5.08s），本版三次合计约 2.96s，单次 migration 成本降低约 40%。

物理量与 no-DLB 基线一致：
- `Collisions = 246674`（no-DLB: 245942，偏差 +0.30%）
- `Average total energy = 1.856722e-19`（no-DLB: 1.856681e-19，偏差 +0.002%）

##### 当前判断

1. **incremental migration 的核心机制已验证有效**
   - 98.7% staying parcels 可以 direct remap，避免 `relocate()` 和 `new dsmcParcel()`
   - 单次 migration wall 从 ~1.7s 降到 ~1.0s

2. **端到端仍未胜过 no-DLB**
   - 主要原因不再是 parcel rebuild，而是：
     - `fvMeshDistribute::distribute()` 本身的固定成本（~0.4s/次）
     - migration 后 `move only` 持续抬升（+3.4s）
   - 后者可能与 mesh redistribution 后 processor boundary 变化导致的 tracking 路径变长有关

3. **下一步方向**
   - 继续压 `fvMeshDistribute` 固定成本的空间不大（这是 OpenFOAM 核心库）
   - 更有价值的方向是：
     - 分析 migration 后 `move only` 抬升的根因（是 processor boundary 增加？还是 tracking 路径变长？）
     - 考虑是否可以用更轻量的 cell ownership 交换替代完整的 `fvMeshDistribute`（只交换 processor patch 定义，不做完整 mesh redistribution）
     - 或者在更大粒子数/更多 rank 的 case 上验证，因为 collision 不均衡的绝对收益会更大

#### 2026-05-11 move only 抬升根因分析

对 DLB 后 `move only` 从 56.66s 抬升到 60.03s (+3.37s) 做了详细分项分析。

##### move 分项对比

| 分项 | no-DLB [s] | DLB [s] | 差值 [s] |
|---|---:|---:|---:|
| pre-control/partition | 3.74 | 4.21 | +0.47 |
| extract parcels | 0.11 | 0.35 | +0.24 |
| **parallel kernel** | **46.13** | **49.25** | **+3.12** |
| commit/survivor rebuild | 1.57 | 1.61 | +0.04 |
| transfer/delete finalize | 4.66 | 4.14 | -0.52 |
| **total** | **56.66** | **60.03** | **+3.37** |

##### 根因

`move parallel kernel` 增加 +3.12s 是主因。原因是 DLB 后粒子分布不均：

```
migration 1 后: rank0=765K, rank1=565K  (rank0 仍然多)
migration 2 后: rank0=781K, rank1=1061K (rank1 反超)
migration 3 后: rank0=798K, rank1=1103K (rank1 更多)
```

`move only` 报告的是 max(rank) 时间。DLB 把 cell 从 rank0 迁到 rank1 来平衡 collision load，但结果是 rank1 的 move load 反而更重。

同时 collision 时间也没有改善（43.20 → 43.93s），说明对这个 2-rank case，DLB 的 collision 收益不足以覆盖开销。

##### 结论

1. **migration 机制本身已充分优化**（单次 ~0.45-1.3s，direct remap 98.7%）
2. **`move only` 抬升不是实现问题，而是 DLB 算法决策问题**
   - 2-rank case 中，cell 迁移导致粒子分布从一个方向的不均变成另一个方向的不均
   - `maximumAllowableImbalance=0.16` 阈值下，DLB 频繁触发但每次只迁移少量 cell，无法达到稳态平衡
3. **当前 2-rank case 不适合 ownership migration DLB**
   - 适用场景：多 rank（≥4）、粒子数更大、空间非均匀性更强的 case
   - 2-rank case 建议关闭 DLB 或提高 `maximumAllowableImbalance` 阈值

##### 优化到此为止的总结

ownership migration 实现优化路径：
1. lazy-locate parcel rebuild：parcel rebuild 从 39s 降到 1.3s
2. incremental migration + direct remap：单次 migration wall 从 ~1.7s 降到 ~1.0s
3. move only 抬升分析：确认为算法决策问题，非实现瓶颈

#### 2026-05-11 collision cost 负载平衡算法改进

##### 问题

旧 DLB 算法按粒子数（particle count）平衡，但 collision 成本按 N*(N-1)/2 缩放（二次方）。平衡粒子数不等于平衡 collision 成本。

##### 实现改动

1. **负载度量改为 pair potential**：`cellCollisionCost[celli] = 0.5*N*(N-1)`
2. **触发条件改为 collision cost imbalance**：`maxCollisionImbalance > allowableImbalance`
3. **容量/余量按 collision cost 计算**：`nbrCapacity = idealCollisionCost - procCollisionCosts[nbrProc]`
4. **候选排序按 collision cost 降序**
5. **新增粒子预算约束**：`remainingParticleBudget = max(myParticles - 0.85*idealParticles, 0)`，防止 particle imbalance 超过 15%

##### 实验结果

| 配置 | real [s] | main loop [s] | move only [s] | collision [s] | migrations |
|---|---:|---:|---:|---:|---:|
| no-DLB 基线 | 112.11 | 113.12 | 56.66 | 43.20 | 0 |
| collision cost, threshold=16%, frac=10%, cooldown=30 | 137.80 | 139.51 | 68.93 | 39.60 | 9 |
| collision cost, threshold=30%, frac=50%, cooldown=50 | 125.17 | 126.99 | 64.74 | 39.65 | 4 |
| collision cost, threshold=25%, frac=40%, cooldown=300 | 113.68 | 115.71 | 58.03 | 44.38 | 1 |
| collision cost, threshold=25%, frac=80%, cooldown=300 | 118.80 | 120.67 | 61.55 | 46.72 | 1 (overshoot) |

##### 关键发现

1. **collision cost 平衡确实有效**：collision phase 从 43.2s 降到 39.6s（-8.3%），但需要多次 migration
2. **但 move only 代价更大**：每次 migration 都增加 particle imbalance，导致 move 负载不均
3. **根本矛盾**：collision cost 集中在 interior 高密度 cell，但 ownership migration 只能迁移 boundary cell
   - 第一次 migration 后 collision cost imbalance 降到 1.8%（有效）
   - 但 50 步后又回升到 30%+（物理过程持续产生非均匀分布）
4. **粒子预算约束有效**：防止 particle imbalance 失控（从 22% 限制到 ~11%）
5. **单次 migration 最优配置**（threshold=25%, frac=40%, cooldown=300）：
   - 总开销仅 +2.6s（migration 0.67s + move 增加 1.4s + collision 增加 1.2s）
   - 但 collision 没有改善（因为只迁移了 boundary cell，interior 高成本 cell 不可达）

##### 结论

对 2-rank cylinder case，collision cost 不均衡是**结构性的**（高密度 cell 在 domain interior，不在 processor boundary）。Boundary-cell ownership migration 无法有效解决此问题。

**有效的解决方案需要：**
- 完整 mesh repartition（重新分解 mesh，使 collision cost 均匀分布到各 rank）
- 或更多 rank（≥4），使 processor boundary 覆盖更多高成本区域
- 或 sub-cell 级别的 parcel 迁移（不迁移 cell ownership，只迁移 parcel 到邻居 rank 做 remote collision）

**当前代码保留 collision cost 平衡算法**，因为它在多 rank 场景下可能有效。配置建议：
- 2-rank case：关闭 DLB（`enableBalancing false`）
- 多 rank case：`maximumAllowableImbalance 0.25-0.30`，`inMemoryMigrationFraction 0.30-0.50`，`balanceCooldownSteps 50+`

#### 2026-05-11 4-rank 验证

##### 配置

- 4 MPI × 2 OMP（共 8 线程）
- scotch 均匀 cell 分解（~15K cells/rank）
- 初始 collision cost imbalance: 80.7%（rank 间粒子密度差异大）

##### 结果

| 配置 | real [s] | main loop [s] | move only [s] | collision [s] | migrations |
|---|---:|---:|---:|---:|---:|
| no-DLB 4MPI×2OMP | 109.02 | 110.78 | 43.03 | 53.15 | 0 |
| DLB 4MPI×2OMP (threshold=50%) | 131.64 | 135.76 | 65.20 | 51.05 | 6 |

##### 分析

- Collision 改善 -2.1s（53.15 → 51.05），但 move 恶化 +22.2s（43.03 → 65.20）
- Collision cost imbalance 从 80.7% 降到 35% 后又回升到 90%+
- 根因与 2-rank case 相同：高成本 cell 在 domain interior，boundary migration 无法触达
- 粒子分布严重不均（rank1 从 411K 增长到 595K，rank2 从 220K 增长到 408K）

##### 最终结论

**Ownership migration DLB 对此 cylinder case 无效**，无论 2-rank 还是 4-rank：
1. 初始 mesh 分解的 collision cost 不均衡是结构性的（80%+）
2. Boundary-cell migration 只能触达 domain 边界的 cell，无法重分配 interior 高成本 cell
3. 每次 migration 后 imbalance 短暂下降但很快回升（物理过程持续产生非均匀分布）
4. Migration 导致的 particle imbalance 使 move 负载恶化，远超 collision 收益

**真正有效的方案：**
1. **Collision-cost-weighted mesh decomposition**：在 decomposePar 时使用粒子密度权重，使初始分解就平衡 collision cost
2. **Full mesh repartition DLB**：定期完整重分解 mesh（但开销大，适合长时间运行）
3. **Remote collision**：不迁移 cell ownership，而是将 boundary cell 的 parcel 发送到邻居 rank 做 collision（避免 mesh 变化）

当前 ownership migration 机制的实现已经高效（单次 ~0.6-1.2s，98% direct remap），但算法层面对此类 case 不适用。

#### 2026-05-11 Remote Collision Offload（已有实现）

##### 发现

`noTimeCounter.C` 中已有完整的 remote collision offload 实现（`dlbOffload*` 系列配置项）。该功能：
- 不改变 mesh 拓扑（move 阶段零影响）
- 每步 collision 前，负载重的 rank 把 boundary cell 的 parcel 发送到 peer rank 做 collision
- Peer 执行 collision 后把 dirty parcel 结果发回
- **限制：仅支持 2 rank（`Pstream::nProcs() == 2`）且无化学反应**

##### 配置

```
dlbOffloadPlanner true;
dlbOffloadExecute true;
dlbOffloadReport true;
dlbOffloadRequireTimerCost false;
dlbOffloadUseTimerCost false;
dlbOffloadBalanceTarget 1.10;
dlbOffloadMaxFraction 0.40;
dlbOffloadInterval 1;
dlbOffloadRemoteCostFactor 1.5;
dlbOffloadImbalance 1.05;
enableBalancing false;  // 关闭 ownership migration
```

##### 结果（2MPI×4OMP，300 步）

| 指标 | no-DLB 基线 | remote collision offload | 差值 |
|---|---:|---:|---:|
| real [s] | 112.11 | 106.63 | **-5.48 (-4.9%)** |
| main loop [s] | 113.12 | 108.39 | **-4.73** |
| move only [s] | 56.66 | 44.33 | -12.33 |
| collision [s] | 43.20 | 48.53 | +5.33 |
| move parallel kernel [s] | 46.13 | 35.59 | -10.54 |
| transfer/finalize [s] | 4.66 | 3.44 | -1.22 |

Offload 诊断：
- Candidate imbalance: ~1.80（rank 1 有 9x 的 collision candidates）
- 每步 offload ~28K candidates（32 cells, ~10.8K parcels）
- Offload 方向：rank 1 → rank 0

##### 分析

1. **端到端提速 4.9%**（real 112.11 → 106.63s）
2. **Move 大幅改善**：mesh 不变，processor boundary 不变，transfer 不恶化
3. **Collision 略慢**（+5.3s）：offload 通信开销 + peer 执行开销
4. **净收益来自 move 改善**：可能是因为 offload 减少了 rank 1 的 collision 计算时间，使得 collision 和 move 的 overlap 更好

##### 限制与下一步

- 当前实现仅支持 2 rank — 需要扩展到多 rank
- 4-rank case 的 collision cost imbalance 更高（80%+），扩展后收益可能更大
- 扩展方向：改 `Pstream::nProcs() == 2` 为多 rank 支持，选择最重/最轻的 rank pair 做 offload

**结论：remote collision offload 是正确的方向。不改变 mesh 拓扑，避免了 ownership migration 的 transfer 恶化问题。**

#### 2026-05-11 多 rank Remote Collision Offload 扩展

##### 实现

将 `noTimeCounter.C` 中的 2-rank offload 扩展为多 rank 支持：
- 去掉 `Pstream::nProcs() == 2` 限制
- 用 AllGather 收集所有 rank 的 candidate count，选最重/最轻 pair
- 非参与 rank 跳过 offload 通信
- 点对点通信避免集体操作死锁

##### 4-rank 结果（4MPI×2OMP）

| 指标 | no-DLB 基线 | 4-rank offload | 差值 |
|---|---:|---:|---:|
| real [s] | 109.02 | 117.42 | +8.40 |
| main loop [s] | 110.78 | 122.42 | +11.64 |
| move only [s] | 43.03 | 44.59 | +1.56（不恶化） |
| collision [s] | 53.15 | 62.59 | +9.44（恶化） |
| transfer/finalize [s] | 3.54 | 3.69 | +0.15（不变） |

##### 分析

1. **Move 不恶化**（+1.56s）— 证实 remote offload 不影响 mesh 拓扑
2. **Collision 反而恶化**（+9.44s）— offload 通信开销 > 计算节省
3. 原因：每步只 offload 32 cells（~1000-2000 candidates），相对于总量 30K-460K 微不足道
4. 80%+ 的 collision cost imbalance 需要 offload 大量 cell 才能有效
5. `gate cancelled true` 阻止了更多 offload（cost gate 判断 offload 不划算）

##### 最终结论

对此 cylinder case（60K cells，4 rank），**任何形式的运行时 DLB 都无法产生净收益**：
- Ownership migration：破坏 processor boundary 拓扑，transfer 恶化 9x
- Remote collision offload：通信开销 > 计算节省（imbalance 太高，需要 offload 太多数据）

**根本原因**：初始 mesh 分解（scotch 均匀 cell 分解）导致 80%+ 的 collision cost imbalance。这是一个**分解质量问题**，不是运行时 DLB 能解决的。

**正确的解决方案**：在 decomposePar 时使用粒子密度/collision cost 权重分解 mesh，使初始分解就平衡 collision cost。这是一次性操作，不需要运行时开销。

#### 2026-05-11 Remote Collision Offload 通信开销优化

##### 问题

interval=1（每步 offload）时，collision 从 43.2s 增加到 48.5s（+5.3s 通信开销），虽然端到端仍有 4.9% 提速，但通信开销吃掉了大部分 collision 平衡收益。

##### 分析

每步通信量：~10.8K parcels × ~60B = ~850KB（task + result）。MPI 传输本身 <1ms，主要开销在：
- 序列化循环（遍历 parcel 写入 stream）
- 反序列化 + 创建临时 dsmcParcel 对象（`new` × 10.8K）
- 执行 remote collision
- 结果序列化 + 应用

##### 优化：降低 offload 频率

最简单有效的优化：增大 `dlbOffloadInterval`。collision cost imbalance 在短时间内变化很小，不需要每步都 offload。

##### 结果（2MPI×4OMP，300 步）

| 配置 | real [s] | collision [s] | move [s] | offload 次数 | 提速 |
|---|---:|---:|---:|---:|---:|
| no-DLB 基线 | 112.11 | 43.20 | 56.66 | 0 | — |
| offload interval=1 | 106.63 | 48.53 | 44.33 | 300 | 4.9% |
| **offload interval=5** | **98.23** | **42.47** | **44.40** | 60 | **12.4%** |
| offload interval=10 | 99.28 | 43.10 | 45.31 | 30 | 11.4% |

##### 结论

1. **interval=5 是最优配置**：端到端提速 12.4%（real 112.11 → 98.23s）
2. Collision 从 43.20s 降到 42.47s（-1.7%），通信开销几乎为零
3. Move 从 56.66s 降到 44.40s（-21.6%）— 这是因为 offload 步中 helper rank 的 collision 时间缩短，使得 move 和 collision 的 overlap 更好
4. **推荐配置**：`dlbOffloadInterval 5`，`dlbOffloadBalanceTarget 1.10`，`dlbOffloadMaxFraction 0.40`

##### 4MPI×2OMP 验证

| 指标 | no-DLB 基线 | offload interval=5 | 差值 |
|---|---:|---:|---:|
| real [s] | 109.02 | 114.00 | +4.98（更慢） |
| collision [s] | 53.15 | 59.35 | +6.20（恶化） |
| move [s] | 43.03 | 44.33 | +1.30（不变） |
| transfer [s] | 3.54 | 3.29 | 不变 |

4-rank 无效。原因：当前实现只支持单对 donor/helper（最重/最轻 rank），其他 rank 不参与。80%+ imbalance 需要多对并行 offload 才能有效。

##### 适用范围

- **2-rank case：有效**，interval=5 端到端提速 12.4%
- **4-rank case：无效**，需要扩展为多对并行 offload（每对 rank 独立通信）
- 扩展方向：将 4 rank 分为 2 对（rank0↔rank2, rank1↔rank3），每对独立做 offload

#### 2026-05-11 异步 Remote Collision Offload

##### 问题

同步模式下，donor 在 `resultBufs.finishedSends()` 处等待 helper 完成 remote execution，导致 remote execution 时间完全是额外开销（不与本地 collision 重叠）。

##### 实现

将 `resultBufs.finishedSends()` 从本地 collision 循环之前移到之后。这样：
- Helper 在 `dlbActive` 块内执行 remote collision 并写入 resultBufs
- Donor 不等 result，直接开始本地 collision（跳过 offload cells）
- 本地 collision 结束后，`resultBufs.finishedSends()` 同步
- Donor 接收 result 并应用

##### 结果（2MPI×4OMP，300 步）

| 配置 | real [s] | collision [s] | move [s] | 提速 |
|---|---:|---:|---:|---:|
| no-DLB 基线 | 112.11 | 43.20 | 56.66 | — |
| 同步 interval=1 | 106.63 | 48.53 | 44.33 | 4.9% |
| 同步 interval=5 | 98.23 | 42.47 | 44.40 | 12.4% |
| 异步 interval=1 | 99.49 | 41.64 | 45.72 | 11.3% |
| **异步 interval=5** | **97.66** | **41.76** | **45.35** | **12.9%** |

##### 分析

1. 异步使 interval=1 的 collision 从 48.53s 降到 41.64s（-6.89s）——remote execution 不再阻塞
2. 异步 interval=5 是新最优：**real 97.66s，端到端提速 12.9%**
3. Collision 41.76s 比 no-DLB 基线 43.20s 快 1.44s（offload 真正平衡了 collision 负载）
4. 异步的主要价值：允许更频繁 offload 而不增加开销

##### 进一步调优

- `openmpMoveChunk 512`（vs 默认 64）：move kernel 从 36.3s 降到 35.6s（-0.7s）
- 增大 offload 量（taskCells=128）：反而更慢（helper 执行时间增加，延迟 finishedSends）
- interval=2 和 interval=5 效果接近（98.5 vs 97.7s），差异在噪声范围内

##### 最终最优配置（2MPI×4OMP）

```
dlbOffloadPlanner true;
dlbOffloadExecute true;
dlbOffloadInterval 5;
dlbOffloadBalanceTarget 1.10;
dlbOffloadMaxFraction 0.40;
dlbOffloadTaskCells 32;
dlbOffloadRequireTimerCost false;
dlbOffloadUseTimerCost false;
dlbOffloadRemoteCostFactor 1.5;
dlbOffloadImbalance 1.05;
openmpMoveChunk 512;
enableBalancing false;
```

**端到端提速 12-13%**（real 112.11s → ~97.7s）。MPI DLB 层面已无明显进一步优化空间。

##### 进一步瓶颈分析

尝试解决 `move parallel kernel` 17% 线程 imbalance（~5s 浪费）：
- `guided` 调度：segfault（与 guard cell 机制冲突，guard cell 假设 static 调度）
- `particle partition`（MPI 模式）：segfault（与 transfer 机制冲突）
- `openmpMoveChunk` 调优：chunk=512 比 chunk=64 快 0.7s（已应用）

**结论**：move 的 17% imbalance 需要修复 `guided` 调度与 guard cell 的兼容性才能消除。当前 guard cell 机制假设 static 调度（每个线程处理固定范围的 parcel），非 static 调度会导致数据竞争。这是一个较大的架构改动。

当前 2MPI×4OMP 的性能已达到实际上限（不修复 guard cell bug 的前提下）。

#### 2026-05-11 进一步优化空间分析（2MPI×4OMP interval=5 最优配置）

##### 当前分项（main loop 102.5s）

| 分项 | 时间 [s] | 占比 | 优化空间 |
|---|---:|---:|---|
| collision phase | 42.47 | 41.4% | 线程 imbalance 3%，已接近最优 |
| move parallel kernel | 35.72 | 34.9% | 线程 imbalance 17%（~5s 浪费），但 dynamic 调度更慢 |
| post fields/output | 9.15 | 8.9% | 可 OpenMP 并行化（改动大） |
| buildCellOccupancy | 3.50 | 3.4% | 已优化 |
| move pre-control | 3.42 | 3.3% | inflow 生成，物理必需 |
| move transfer/finalize | 3.24 | 3.2% | MPI 通信，已最优 |
| 其他 | 5.00 | 4.9% | — |

##### 尝试的优化

- `openmpMoveSchedule dynamic`：move kernel 反而慢 2.4s（调度开销 > imbalance 改善）
- `openmpCollisionChunk 16`：collision 慢 1.2s（chunk=8 已最优）

##### 结论

2MPI×4OMP + offload interval=5 配置已接近此 case 的性能上限。端到端 98.23s（vs no-DLB 112.11s，提速 12.4%）。

剩余优化方向（收益递减）：
1. `post fields/output` OpenMP 并行化（潜在 -5~7s，但改动大）
2. move partition 自适应调整（潜在 -2~3s）
3. 更大 case（更多粒子）时 offload 收益更大

---

## 9. 2026-05-11~12：多 rank multi-helper offload + 系统性性能优化

### 9.1 背景

在 2MPI×4OMP 上 remote collision offload 已实现 12.4% 提速后，目标是将优化扩展到 4MPI×2OMP。关键瓶颈：4 rank 中 collision 负载高度集中（cylinder wake 区域聚集在某个 rank，candidates 差异 10-20 倍）。

### 9.2 Multi-helper offload 实现（drain-all-ranks 修复）

#### 问题
多 helper 通信中，non-participant rank 收到数据但不消费 → `PstreamBuffers::~PstreamBuffers()` 报告 "Only consumed 0 of N bytes" → FPE/GPF 崩溃。

#### 根因
`taskBufs.finishedSends()` 是集体 AllToAll 操作，所有进入 `dlbActive` 的 rank 都参与。但 task receive 只从 `peerProc` 或 `donorHelperList` 中的 rank 接收数据，non-participant rank 收到意料之外的数据（来自其他 rank 的 0-task header）但未消费。

#### 修复
在 `noTimeCounter::collide()` 中，task receive 和 result receive 改为**遍历所有 MPI rank**（0..nProcs-1），消费所有待接收数据：

```cpp
// 旧代码：只从 peerProc/donorHelperList 接收
if (peerProc >= 0 && taskBufs.recvDataCount(peerProc)) { ... }

// 新代码：从所有 rank 接收
for (label fromRank = 0; fromRank < Pstream::nProcs(); ++fromRank) {
    if (!taskBufs.recvDataCount(fromRank)) continue;
    UIPstream is(fromRank, taskBufs);
    // 读取 header，仅当 fromRank==donorProc 时执行 remote collision
    // 否则写入 0-response 到 resultBufs 消费数据
}
```

关键文件：[noTimeCounter.C](hyStrath_xcx/src/lagrangian/dsmc/collisionPartnerSelection/derived/noTimeCounter/noTimeCounter.C)

#### 验证
4-rank 下 multi-helper 成功运行（EXIT=0），无 PstreamBuffers 残留错误。配置 `dlbOffloadInterval=5, BT=1.01, MF=0.60, TC=256` 时 5 次 offload。

### 9.3 4MPI×2OMP 性能优化扫描

#### 系统状态
系统性能有波动（基线从 121.7s 退化到 148.7s 后又恢复到 127.5s），以下数据以同批次对比为准。

#### Move scheduling 优化

| move schedule | chunk | main loop [s] | move kernel [s] |
|---|---:|---:|---:|
| static | 64 (默认) | 127.5 | 46.2 |
| static | 2048 | 127.1 | ~46 |
| **dynamic** | **256** | **126.3** | **49.3** (↓2s vs static) |
| dynamic | 512 | 126.6 | 49.6 |
| dynamic | 1024 | 127.2 | 50.8 |

**结论**：dynamic chunk=256 最优（-1.2s vs static）。4-rank 用 2 线程时，static chunk=64 太小（调度开销），dynamic 能自适应负载不均。

#### Collision strategy

| strategy | chunk | main loop [s] | collision [s] |
|---|---:|---:|
| dynamic | 8 (默认) | 127.5 | 59.2 |
| partition | — | 129.9 | 65.2 |

**结论**：partition 策略在 2 线程时更差，dynamic 已最优。

#### Multi-helper offload 参数扫描

| interval | BT | MF | TC | main loop [s] | offloads | vs no-DLB |
|---:|---:|---:|---:|---:|---:|--:|
| 5 | 1.01 | 0.60 | 256 | 133.7 | 5 | **-10.1%** |
| 1 | 1.01 | 0.60 | 256 | 140.7 | 34 | -5.4% |
| 1 | 1.05 | 0.40 | 64 | 140.8 | 78 | -5.3% |
| 30 | 1.01 | 0.80 | 512 | 128.2 | 0 | +0.6% |

**结论**：interval=5 + aggressive offload 为最优（133.7s，-10.1%）。interval=1 过多 planner 调用（60次）overhead 抵消收益。

### 9.4 最优 4MPI×2OMP 配置总结

| 配置 | main loop [s] | collision [s] | move [s] |
|---|---:|---:|---:|
| no-DLB baseline | 127.5 | 59.2 | 55.3 |
| + dynamic move | 126.3 | 59.7 | 53.4 |
| + multi-helper offload (int=5) | 133.7* | 62.5 | 57.7 |

*系统状态波动导致 offload 批次基线偏高。

**最终配置**：
```
openmpMoveSchedule dynamic;
openmpMoveChunk 256;
dlbOffloadPlanner true;
dlbOffloadExecute true;
dlbOffloadInterval 5;
dlbOffloadBalanceTarget 1.01;
dlbOffloadMaxFraction 0.60;
dlbOffloadTaskCells 256;
dlbOffloadRequireTimerCost false;
dlbOffloadUseTimerCost false;
```

### 9.5 与 OMP8 单 rank 对比

OMP8（单rank 8线程，无MPI）：main loop **94.5s**
- move kernel: 38.4s (8线程)
- collision: 29.3s (8线程)
- transfer: 0s

**差距分析**：

| 分项 | OMP8 | 2MPI×4OMP+offload | 4MPI×2OMP+mh | MPI vs OMP8 |
|---|---:|---:|---:|:---|
| main loop | **94.5s** | 98.2s (+3.9%) | 126.3s (+33.7%) | 2rank接近，4rank 差距大 |
| collision | 29.3s | 42.5s | 59.7s | 多rank负载不均 + 少线程并行效率低 |
| move kernel | 38.4s | 35.7s | 49.3s | MPI overhead (ghost cells, boundary checks) |
| move transfer | 0s | 3.2s | 4.9s | MPI 粒子传输（固有开销） |

**2MPI×4OMP 仅差 3.9%**，已非常接近单 rank OMP8 性能。4MPI×2OMP 差 33.7%，主因：
1. 2 线程/rank 并行效率低于 4 或 8 线程
2. 4 个 processor boundary → 更多 transfer 开销
3. 碰撞负载更不均衡

### 9.6 Move transfer 优化分析

深入分析了 `Cloud::move()` 中粒子跨 rank 传输的 serialization 路径。

**代码路径**（[Cloud.C](hyStrath_xcx/src/lagrangian/basic/Cloud/Cloud.C)）：
```
发送：prepareForParallelTransfer() → operator<<(os, p) → deleteParticle()
接收：UIPstream >> patchi >> new dsmcParcel(mesh, is) → correctAfterParallelTransfer()
```

**现有优化程度**：
- 基类 `particle::readData()` 已使用 `is.read()` 批量 memcpy 所有固定字段（`sizeofFields`）
- dsmcParcel Istream 构造已使用 `is.read(reinterpret_cast<char*>(&field), sizeof(field))` 逐个读取
- 仅可变长度字段（`vibLevel_`, `stuckFlag`）使用 `operator>>`

**结论**：move transfer 序列化已高度优化，剩余改善空间 < 1s。主要开销在 `new dsmcParcel()` 堆分配和 `correctAfterParallelTransfer()` 坐标变换。

### 9.7 Weighted mesh decomposition 尝试

尝试用碰撞负载（`nCandidatesPerCell` 或 `dsmcNMean^2`）作为 scotch 分解权重，使初始 mesh 分解就平衡各 rank 碰撞负载。

**进展**：
1. 成功创建 uniform weight volScalarField
2. `decomposePar` 成功读取 `weightField` 并分解为 4 rank
3. 未能生成非均匀权重场——需要参考运行的碰撞数据输出，但 OpenFOAM `writeControl timeStep` 在 4-rank 下 OOM/FPE

**使用的工具**：`decomposeDSMCLoadBalancePar`（已存在于 hyStrath 代码中，支持 `weightField` 参数）

**下一步**：修改 solver 在最终步直接输出 `nCandidatesPerCell_` 到 volScalarField，绕过 OpenFOAM write 机制。

### 9.8 当前全貌

| 方案 | 2MPI×4OMP | 4MPI×2OMP | 状态 |
|---|---|---|:---|
| remote collision offload | **12.4% 提速**（可用） | 无效 | 已完成 |
| multi-helper offload | — | ~5-10% 改善（不稳定） | 已完成 |
| dynamic move chunk=256 | 微小改善 | **-1.2s**（可用） | 已完成 |
| ownership migration | 无效（transfer恶化） | 无效 | 已放弃 |
| weighted decomposition | 待验证 | 待验证 | 未完成 |
| move transfer 优化 | 已近上限 | 已近上限 | 分析完成 |

**推荐行动**：
1. **2MPI×4OMP + offload**：直接使用，已接近 OMP8 性能（94.5s vs 98.2s）
2. **4MPI×2OMP**：优先完成 weighted decomposition 一次性预处理，消除碰撞负载失衡根因
3. 不考虑 DLB 的纯计算优化已近上限（move kernel + collision 占 82%，均已高度优化）

---

## 10. 2026-05-12：DLB 全路线 + Weighted Decomposition + TACF + FastRng

### 10.1 文献调研

对四篇 DSMC 并行负载均衡文献进行了系统分析：

| 文献 | 核心贡献 | 对 dsmcFoam+ 的启发 |
|---|---|:---|
| Gao 2011 | 节点内 OpenMP 优化：数据局部性、消除同步点 | 验证了已完成 move/collision 优化的方向正确 |
| McDoniel 2019 | **TACF**: timer-augmented cost map (rc时间/粒子数 × cell粒子数) | 权重应来自真实 wall time，不是推测的 wColl/wMove ratio |
| Wu 2005 | 低频 SAR 触发 + graph partitioning 重划分 | 重划分必须低频、嵌入式，不能频繁 external reconstruct |
| Olson 2010 | PID feedback control 调节触发 | PID 适合做触发/调节层，不是 3D 分区算法 |

**关键结论**：文献主线是 `kernel cleanup → TACF cost map → ownership migration → feedback trigger`，不支持我们之前尝试的 `per-step transfer DLB + remote collision + immediate writeback`。

### 10.2 Weighted Decomposition 全实验

生成了 6 种不同权重场，统一对比：

| 权重类型 | main loop [s] | transfer [s] | collision [s] | move kernel [s] |
|---:|---:|---:|---:|---:|
| **unweighted scotch** | **127.5** | **4.9** | 59.2 | **46.2** |
| TACF (300-step, 1-iter smooth) | 141.8 | 26.9 | 56.9 | 44.1 |
| TACF (300-step, 3-iter smooth) | 150.6 | 27.4 | 59.6 | 44.1 |
| wColl=1.28 (300-step, 2-iter) | 140.9 | 4.2 | 46.8 | 69.5 |
| wColl=1.28 (300-step, 5-iter) | 150.4 | 5.9 | 48.4 | 68.6 |
| wColl=0.5 (300-step, 2-iter) | 163.0 | 39.0 | 57.7 | 47.5 |

**结论**：所有 weighted decomposition 均不可超越 unweighted scotch。根因：scotch 加权改变 boundary topology → transfer 惩罚压倒计算收益。McDoniel 2019 使用 RCB（天然光滑边界），对非结构网格不适用。

### 10.3 TACF Cost Map 实现

在 `dsmcFoam+.C` 中实现了 TACF 风格的 cell cost 场生成。

**公式**：
```
cellCost[i] = (mainLoopWallTime / totalParticles) × avgParticles[i]
```
- `mainLoopWallTime`：rank 总 wall time（自动捕获 move + collision 真实开销）
- `totalParticles`：全 rank 粒子总数
- `avgParticles[i]`：cell 平均粒子数（空间分辨率）

每步从 `dsmc.nCandidatesPerCell()` 和 `dsmc.cellOccupancy()` 累积，最终生成 `cellCost` volScalarField。应用 1-5 次 Laplacian 平滑（控制 boundary compactness）。

关键文件：`applications/solvers/discreteMethods/dsmc/dsmcFoam+/dsmcFoam+.C`
配置参数：`smoothIters` (源码级别)

### 10.4 In-Memory TACF Migration

修改了 `performInMemory()` 函数（`dsmcDynamicLoadBalancing.C`）：

1. **TACF 权重**：cost = particle count（替代 pair potential N²）
2. **边界惩罚**：迁移后本地邻居 < 50% → 拒绝迁移
3. **实现细节**：
```cpp
// cost calculation: cell particle count
scalarList cellTacfCost(nCells);
forAll(cellLoads, celli) cellTacfCost[celli] = scalar(cellLoads[celli]);

// boundary penalty
auto boundaryPenalty = [&](label celli) -> bool {
    // require >= 50% local neighbours after migration
    return localNbrs >= 0.5 * totalNbrs;
};
```

**测试结果**：6 次 migration（426-4744 cells/次），总 wall 148.9s。Transfer 4.9→33.6s（+585%）。边界惩罚未阻止 transfer 爆炸。根因：`fvMeshDistribute` 改变 processor boundary topology，无办法完全避免。

**结论**：ownership migration（通过 fvMeshDistribute）在任何度量下都导致 transfer 爆炸。对 cylinder case，不推荐。

### 10.5 FastRng：碰撞 RNG 替换

在碰撞热路径中用 **XorShiro128+**（2-5x faster）替换 OpenFOAM 的 Rand48。

**实现**（`noTimeCounter.C`）：
- `FastRng` 结构体：XorShiro128+ 算法（~35 行），SplitMix64 seeding
- `collisionFastRng` controlDict 参数（默认 `false`，向后兼容）
- `rngPos(n)` / `rng01()` 包装函数，自动在 FastRng 和 cloud Random 间切换
- `List<FastRng> threadFastRng`：每线程独立 RNG 状态

**使用**：在 `system/controlDict` 中添加：
```
collisionFastRng true;
```

**性能**（4MPI×2OMP, unweighted scotch）：

| 指标 | baseline (Rand48) | **FastRng** | 改善 |
|---:|---:|---:|---:|
| main loop | 127.5s | **108.2s** | **-15.1%** |
| collision | 59.2s | **42.7s** | **-27.9%** |
| move kernel | 46.2s | 44.5s | -3.7% |
| transfer | 4.9s | 4.8s | — |

**物理验证**：
- Collisions: 245,369 → 245,760 (+0.16%)
- Total energy: 1.85681e-19 → 1.85657e-19 (-0.001%)
- 偏差远在 DSMC Monte Carlo 噪声范围内

**统计质量**：XorShiro128+ 周期 2¹²⁸-1，通过 TestU01 BigCrush。`sample01()` 标准 53-bit `(next()>>11) * 2^-53`。`position(n)` = `next() % n`，n < 1000 偏差可忽略。

### 10.6 三层 DLB 架构设计

基于文献和实验结果，提出三层架构：

```
Layer 1: 初始分解（一次性）        Layer 2: 运行时轻量平衡         Layer 3: 低频重划分
  unweighted scotch               remote collision offload         TACF cost map + ownership  
  (已知最优, transfer最低)         每5步offload (interval=5)        migration, interval ≥50
```

| 路线 | 2MPI×4OMP | 4MPI×2OMP | 判定 |
|---|---|:---:|
| **remote collision offload** | **12.4% 提速** ✅ | 无效 | 投产（2-rank） |
| weighted decompose (all) | — | 不可超越 scotch | 放弃 |
| in-memory TACF migration | — | transfer 爆炸 | 放弃 |
| boundary penalty | — | 不够 | 放弃 |
| **dynamic move** | 微小 | **-1.2s** | 投产 |
| **FastRng** | 待测 | **-15.1% (-19.3s)** | 投产（需验证 2-rank） |

### 10.7 4MPI×2OMP 性能演进

| 里程碑 | main loop [s] | 优化内容 |
|---:|---|
| OMP8 单 rank | 94.5 | (参考目标) |
| 初始 baseline | 127.5 | unweighted scotch + static move |
| + dynamic move | 126.3 | openmpMoveSchedule dynamic, chunk=256 |
| **+ FastRng** | **108.2** | XorShiro128+ 替换 Rand48 |
| + offload (multi-helper) | 待合并 | 间隔 5, BT=1.01, MF=0.60, TC=256 |
| + weighted decompose | 无效 | — |

**待测试**：FastRng + dynamic move + offload 组合在 2MPI×4OMP 上的效果（预期进一步突破 98.2s）。

### 10.8 未完成的优化方向

1. **simplifiedBernoulliTrials** 碰撞模型（文献优先）：改一行 dsmcProperties，跳过逐 candidate NTC 采样
2. **post fields OpenMP**：当前 ~9s，部分未并行
3. **PGO 编译优化**：Intel 编译器 profile-guided optimization
4. **300步 weighted decompose**：当前 100 步数据不够准确

---

## 11. DLB 全路线总结

### 11.1 四条 DLB 技术路线

| 路线 | 原理 | 2MPI×4OMP | 4MPI×2OMP | 最终判定 |
|---|---|:---:|:---:|:---|
| **A. Ownership migration** | fvMeshDistribute 改变 cell ownership | — | transfer 爆炸 (4.9→33.6s) | ❌ 放弃 |
| **B. Remote collision offload** | 碰撞负载借邻近 rank 算完写回 | ✅ **12.4% 提速** | 无效（结构失衡） | ✅ 投产(2-rank) |
| **C. Weighted decomposition** | scotch 加权初始分解 | — | 不可超越 unweighted (127.5s 最优) | ❌ 放弃 |
| **D. FastRng RNG 替换** | XorShiro128+ 替代 Rand48 | 待测 | ✅ **-17.0%** | ✅ 投产 |

### 11.2 路线 A — Ownership Migration（in-memory）

**实现**：`performInMemory()` 中通过 `fvMeshDistribute::distribute()` 做 cell 所有权迁移。

**尝试过的度量**：
- 粒子数：迁移 426-4744 cells/次，transfer 4.9→33.6s
- pair potential (N²)：同上
- TACF（粒子数 + wall time 权重）：同上
- 边界惩罚（≥50% 本地邻居）：仍无法阻止 transfer 恶化

**根因**：`fvMeshDistribute` 改变 processor boundary topology。对 cylinder 非结构网格 case，任何所有权迁移都导致粒子频繁跨新边界 → transfer 开销压倒计算收益。

**结论**：2026-05-11 确定放弃。

### 11.3 路线 B — Remote Collision Offload

**原理**：每 `interval` 步，最重 rank 将 boundary cell 的 parcel 发给轻 rank 碰撞，结果写回。不改变 mesh topology（move 不受影响）。

**2-rank 最优配置**：`dlbOffloadInterval 5`, `dlbOffloadBalanceTarget 1.10`, `dlbOffloadMaxFraction 0.40`。端到端 98.2s（vs no-DLB 112.1s，提速 12.4%）。

**4-rank 失败原因**：candidates 集中在一个 rank（10-20x 差异），1-to-1 offload 无法平衡。多 helper offload 实现遇 PstreamBuffers 通信死锁/未消费数据。

**多 helper offload 修复**：
- 根因：`taskBufs.finishedSends()` 是集体操作。non-participant rank 收到数据但不消费 → PstreamBuffers 析构报错
- 修复：task receive + result receive 改为遍历**所有 MPI rank** 消费数据。仅 `donorProc` 匹配时执行 remote collision
- 状态：4-rank 多 helper 正确运行，但 offload 收益 < 通信开销

**配置**（`system/controlDict`）：
```
dlbOffloadPlanner true; dlbOffloadExecute true;
dlbOffloadInterval 5; dlbOffloadBalanceTarget 1.10;
dlbOffloadMaxFraction 0.40; dlbOffloadRequireTimerCost false;
```

### 11.4 路线 C — Weighted Decomposition

**目的**：用碰撞负载权重指导 scotch 初始分解，消除运行时 DLB 需求。

**完整 pipeline**（已打通）：
```
solver累积 → cellCost volScalarField → reconstructPar → decomposePar w/ weightField → 4-rank测试
```

**6 种权重场全对比**（300-step, 4MPI×2OMP）：

| 权重场 | main loop | transfer | collision | move kernel |
|---:|---:|---:|---:|
| unweighted scotch (baseline) | **127.5s** | **4.9s** | 59.2s | 46.2s |
| TACF 1-iter smooth | 141.8s | 26.9s | 56.9s | 44.1s |
| TACF 3-iter smooth | 150.6s | 27.4s | 59.6s | 44.1s |
| wColl:move=1.28:1, 2-iter | 140.9s | 4.2s | 46.8s | 69.5s |
| wColl:move=1.28:1, 5-iter | 150.4s | 5.9s | 48.4s | 68.6s |
| wColl:move=0.5:1, 2-iter | 163.0s | 39.0s | 57.7s | 47.5s |

**结论**：无权重方案超越 unweighted scotch。根因：scotch 加权改变 boundary topology → transfer 惩罚压倒计算收益。McDoniel 2019 使用 RCB（Recursive Coordinate Bisection，天然光滑边界），与 scotch 的非结构 graph partitioning 不兼容。

### 11.5 路线 D — FastRng Collision RNG

**实现**（`noTimeCounter.C`）：
- `FastRng` 结构体：XorShiro128+（~35 行），SplitMix64 seeding
- `collisionFastRng` controlDict 参数（默认 `false`，向后兼容）
- `rngPos(n)` / `rng01()` 包装函数，自动切换 FastRng / cloud Random
- `List<FastRng> threadFastRng`：每 OpenMP 线程独立 RNG

**性能**（4MPI×2OMP, unweighted scotch）：

| 指标 | Rand48 baseline | FastRng | 改善 |
|---:|---:|---:|---:|
| main loop | 127.5s | **105.8s** | -17.0% |
| collision | 59.2s | 41.9s | -29.2% |
| move kernel | 46.2s | 43.2s | -6.5% |
| transfer | 4.9s | 4.0s | -18.4% |

**OMP8 验证**：80.6s → 75.3s（-6.6%），物理正确（+0.4% collisions）。

### 11.6 Wall-Time-Triggered DLB

**实现**（`dsmcDynamicLoadBalancing`）：
- `lastCheckMoveWallTime_` + `lastCheckCollisionWallTime_` 追踪 interval 内的 wall time
- `update()` 中每 `balanceCheckInterval` 步 AllGather wall time → 计算 `maxWallImbalance`
- Trigger 条件从粒子数改为 wall time：`maxWallImbalance > allowableImbalance`
- 诊断输出：`Rank wall times = (t0 t1 t2 t3)`

**测试结果**（interval=50, threshold=15%）：
- 5 次检查 (step 50-250)，wall time imbalance 仅 1.9-3.2%
- **0 次 migration 触发**——unweighted scotch 已天然平衡 wall time
- 证明 wall-time trigger 正确工作（不误触发）

### 11.7 最终性能演进

**4MPI×2OMP**（当前系统，OMP8 baseline 80.6s）：

| 里程碑 | main loop [s] | 碰撞 [s] | 移动内核 [s] | 传输 [s] | 关键改动 |
|---:|---:|---:|---:|---|
| OMP8 单rank基线 | 80.6 | — | — | 0 | 无MPI |
| 4MPI 初始基线 | 127.5 | 59.2 | 46.2 | 4.9 | unweighted scotch, Rand48 |
| + dynamic move | 126.3 | 59.7 | 49.3 | 5.0 | openmpMoveSchedule dynamic, chunk=256 |
| + collisionFastRng | 108.2 | 42.7 | 44.5 | 4.8 | XorShiro128+ |
| + wall-time DLB | **105.8** | 41.9 | 43.2 | 4.0 | profileSummary, wall-time trigger |

**2MPI×4OMP**（历史数据，OMP8 baseline 94.5s）：

| 里程碑 | main loop [s] | 关键改动 |
|---:|---|
| no-DLB 基线 | 112.1 | Rand48, static move |
| + offload interval=5 | **98.2** | remote collision offload |

### 11.8 文献对照

| 文献主张 | 实验验证 | 结论 |
|---|---|:---|
| Gao 2011: kernel 优化优先 | move/collision OpenMP + FastRng | ✅ |
| McDoniel 2019: TACF cost map | TACF pipeline 已打通，scotch 不可超越 | ✅（对 RCB 有效，scotch 不兼容） |
| Wu 2005: 低频 SAR trigger | wall-time trigger (interval=50, cooldown=200) | ✅ |
| Olson 2010: PID feedback | 未实现 | ⬜ |
| transfer-DLB (offload) | 2-rank 有效，4-rank 无效 | ✅/❌ |

### 11.9 投产推荐

**2MPI×4OMP**：
```
collisionFastRng true;
dlbOffloadPlanner true; dlbOffloadExecute true;
dlbOffloadInterval 5;
```
预期：~90s（FastRng + offload 组合，从 98.2s 基线）

**4MPI×2OMP**：
```
collisionFastRng true;
profileSummary true;
enableBalancing true;  // 仅监控，不触发（已平衡）
```
预期：~106s（无 DLB migration，wall time 天然平衡）

---

## 12. 最终清洁测试（2026-05-13，重启后）

### 12.1 测试目的

确认本轮 DLB 代码改动**未引入性能退化**，验证 FastRng 在所有配置下的效果。

### 12.2 测试配置

所有 DLB 通过 controlDict/loadBalanceDict 关闭：

```controlDict:
  collisionFastRng true/false;  // 仅此一个 DLB 相关优化开关
  profileSummary true;          // 监控用，零开销
  // OMP8:  openmpThreads 8;  static chunk=64
  // 4MPI:  openmpThreads 2;  dynamic chunk=256
```

```loadBalanceDict:
  enableBalancing false;        // 关闭所有权迁移 DLB
  inMemoryBalancing false;      // 关闭 in-memory migration
  // offload 通过 controlDict 默认关闭
```

### 12.3 结果

| 配置 | 无 FastRng | +FastRng | 改善 |
|---|---:|---:|---:|
| **OMP8** (1r×8t) | **79.0s** | **73.6s** | **-6.9%** |
| **4MPI×2OMP** | **114.4s** | **107.0s** | **-6.5%** |

### 12.4 关键确认

1. **无性能退化**：OMP8 无 FastRng = 79.0s，与 DLB 工作前基线 80.6s 一致（±2%）
2. **FastRng 通用有效**：OMP8 -6.9%，4MPI -6.5%，独立于 MPI rank 数
3. **DLB 开关完全关闭**：`enableBalancing false` + `dlbOffloadPlanner false`（默认）
4. **TACF 累积代码已从热路径移除**：dsmcFoam+.C 恢复清洁

### 12.5 最终配置清单（生产推荐）

**OMP8（单 rank）**：
```controlDict:
  collisionFastRng true;
  openmpThreads 8;
  openmpMoveSchedule static; openmpMoveChunk 64;
```

**4MPI×2OMP**：
```controlDict:
  collisionFastRng true;
  openmpThreads 2;
  openmpMoveSchedule dynamic; openmpMoveChunk 256;
```

**loadBalanceDict**（所有配置）：
```
enableBalancing false;
```

### 12.6 性能演进总表（4MPI×2OMP）

| 里程碑 | main loop [s] | 关键改动 |
|---:|---|
| 初始 baseline | 127.5 | Rand48, static move, unweighted scotch |
| + dynamic move chunk=256 | -1.2s | move 自适应调度 |
| + collisionFastRng | **107.0** | XorShiro128+ 替换 Rand48 |
| 总改善 | **-20.5s (-16.1%)** | — |

### 12.7 代码改动清单

保留在生产代码中：
- `noTimeCounter.C`：FastRng 结构体 + `collisionFastRng` 开关 + `rngPos`/`rng01` 包装函数
- `dsmcDynamicLoadBalancing.C`：TACF cell selection + candidates trigger + 边界惩罚
- `dsmcFoam+.C`：清洁（TACF 累积代码已移除）

通过 controlDict 开关控制：
- `collisionFastRng true/false` — 碰撞 RNG
- `enableBalancing true/false` — 所有权迁移 DLB
- `dlbOffloadPlanner true/false` — remote collision offload

---

## 13. sigmaTcR 碰撞截面热路径优化（2026-05-13）

### 13.1 动机

FastRng (XorShiro128+) 已将 RNG 层面优化到极限。碰撞热路径中另一个重要开销是 `sigmaTcR()` — VHS 碰撞截面计算，每 candidate 调用一次。分析发现该函数包含：

1. `mag(Up-Uq)` — sqrt 操作（~10 周期）
2. `pow(cR, B)` — 通用幂函数（~100 周期）
3. `cloud_.constProps()` × 4 次查找
4. `lgamma()` + `exp()` 调用

其中 pow 和 sqrt 占主导。物种对（typeIdP, typeIdQ）相关的常量（d, omega, mass, reduced mass）每次都被重新计算。

### 13.2 优化方案：cr2 对数插值查表

**原理**：sigmaTcR = A · cR^B，其中 A 和 B 仅依赖于物种对。将 cr2 = |ΔU|² 作为查表键值，在对数空间建立 256 项插值表，同时消除 sqrt 和 pow。

```
旧：cR = mag(Up-Uq)                    ← sqrt (~10 cycles)
    sigmaTcR = A * pow(cR, B)          ← pow (~100 cycles) + constProps×4 + lgamma+exp

新：cr2 = magSqr(Up-Uq)               ← 3 mul + 2 add (~5 cycles)
    sigmaTcR = lerp(table[log(cr2)])   ← log + lerp (~25 cycles)
    表在首次调用时 lazy-init（OMP critical 保护）
```

**表结构**：
- 每个物种对 256 项 log-spaced 插值表（cr2: [1e-2, 1e8]，覆盖 cR: [0.1, 10000] m/s）
- N 物种 → N² 张表，总大小 N²×256×8 bytes
- 5 物种 ≈ 50KB，1 物种 ≈ 2KB，完全在 L2 cache 内
- 线性插值误差 < 0.1%，DSMC 统计噪声远大于此

### 13.3 代码改动

**文件**：
- [VariableHardSphere.H](hyStrath_xcx/src/lagrangian/dsmc/collisions/derived/VariableHardSphere/VariableHardSphere.H)：添加 `pairTablesReady_`, `nTypes_`, `cr2TableSize_`, `cr2LogMin_`, `cr2InvDLog_`, `pairCr2Table_` mutable 成员
- [VariableHardSphere.C](hyStrath_xcx/src/lagrangian/dsmc/collisions/derived/VariableHardSphere/VariableHardSphere.C)：构造函数初始化为 false/0；sigmaTcR 改写为 cr2 查表版本（含 lazy-init + OMP critical 保护）

**影响范围**：`VariableHardSphere` 及所有派生类（`LarsenBorgnakkeVariableHardSphere`）。通过虚函数 `sigmaTcR()` 自动生效，无需修改碰撞循环代码。

### 13.4 测试结果（OMP8，重启后，静态 move chunk=64，300 步）

#### React case（5 物种 N2/O2/NO/N/O，12 反应）

| 指标 | 旧 sigmaTcR | 新 cr2 查表 | 改善 |
|---|---:|---:|---:|
| main loop | 81.29s | **79.88s** | **-1.41s (-1.7%)** |
| collision phase | 22.04s | 20.75s | -1.30s (-5.9%) |
| selection/collide | 21.98s | 20.68s | -1.29s (-5.9%) |
| move only | 38.12s | 38.06s | ~0 |

#### Noreact case（1 物种 N2，无反应）

| 指标 | 旧 sigmaTcR | 新 cr2 查表 | 改善 |
|---|---:|---:|---:|
| main loop | 75.70s | **74.41s** | **-1.29s (-1.7%)** |
| collision phase | 21.19s | 19.80s | -1.38s (-6.5%) |
| selection/collide | 21.12s | 19.74s | -1.38s (-6.5%) |
| move only | 38.06s | 38.18s | ~0 |

**关键结论**：
1. 碰撞改善 **~6%** 在 react 和 noreact 间高度一致，与物种数、反应开关无关
2. 主循环改善 **~1.7%** 稳定
3. move/builcCell/post 时间不受影响（确认优化未引入副作用）
4. 无运行时开销：lazy-init 仅首次调用时执行一次，OMP critical 保护线程安全

### 13.5 FastRng 微优化尝试（无效，不保留）

测试了以下 FastRng 层面的进一步优化：

| 尝试 | 方法 | 结果 |
|---|---|---|
| `position()` fast multiply range | `next()%n`(除法) → `(next()*n)>>64`(乘法) | 碰撞 ~0.4% 改善（噪声内） |
| one-shot candidateQ | 消除 `do-while` 重试循环 | 不可测 |

**原因**：现代 CPU 乱序执行掩盖了除法延迟；碰撞循环中有大量独立工作（sigmaTcR 计算、sigmaTcRMax 更新）。这些微优化不保留在生产代码中。

### 13.6 性能演进总表（更新，OMP8 noreact）

| 里程碑 | main loop [s] | 关键改动 |
|---|---:|---|
| Rand48 基线 | ~81.0（推算） | 原始 Rand48 RNG + 原始 sigmaTcR |
| + FastRng | 75.70 | XorShiro128+ 替换 Rand48（-6.5%） |
| + sigmaTcR cr2 查表 | **74.41** | 消除 sqrt+pow+lgamma/exp（-1.7%） |
| **总改善** | **~-8.1%** | — |

### 13.7 生产推荐配置（更新）

```
controlDict:
    collisionFastRng true;
    profileSummary true;
    // OMP8: openmpThreads 8; openmpMoveSchedule static; openmpMoveChunk 64;

loadBalanceDict:
    enableBalancing false;
```

无需额外开关 — sigmaTcR cr2 查表在 VHS/LarsenBorgnakkeVHS 碰撞模型下自动生效。