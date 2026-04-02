
# dsmcFoam+ 从 OpenFOAM-v1706 到 OpenFOAM-v2506 移植说明文档（草案）

## 1. 文档目的

本文档用于说明 `hyStrath` 中 `dsmcFoam+` 求解器及其相关库从 `OpenFOAM-v1706` 迁移到 `OpenFOAM-v2506` 的背景、实施范围、关键兼容性问题、当前验证结果以及已知遗留问题。  
本文档依据你提供的移植过程记录整理而成，重点反映**实际做过且在记录中能够确认**的工作，不对未在记录中出现的细节做臆测。

---

## 2. 原始版本与目标版本

### 2.1 原始版本
- 宿主平台：`OpenFOAM-v1706`
- 原始求解器来源：`hyStrath` 中的 `dsmcFoam+`
- 原始实现特征：
  - `dsmcFoam+` 不是简单基于原版 `dsmcFoam` 的少量补丁，而是在 `OpenFOAM-v1706` 原生 DSMC 框架之上进行了较大规模扩展。
  - `hyStrath` 自带定制的：
    - `lagrangian/basic`
    - `lagrangian/dsmc`
    - `dsmcCloud`
    - `dsmcParcel`
    - solver 入口 `dsmcFoam+`
  - 原始版本中还集成了动态负载均衡、AMR、扩展边界、反应模型、变量时间步、额外坐标系等功能。

### 2.2 目标版本
- 宿主平台：`OpenFOAM-v2506`
- 当前移植目标：
  - 先在新目录 `hyStrath_xcx` 中构建后续移植基线；
  - 第一阶段仅解决：
    - `liblagrangian+`
    - `libdsmcFoam+`
    - `dsmcFoam+`
    的编译与基础运行；
  - `AMR` 与 `dynamicLoadBalancing` 暂不纳入主迁移链路，留待后续单独处理。

---

## 3. 移植范围与阶段划分

### 3.1 总体原则
本次迁移并不是简单把 `dsmcFoam+.C` 从旧版本拷贝到新版本，而是将 `hyStrath dsmcFoam+` 所依赖的整套定制 `lagrangian/DSMC` 栈重新挂接到 `OpenFOAM-v2506` 的宿主 API 上。

### 3.2 第一阶段实际边界
第一阶段以“**最小可编译、最小可运行**”为目标，优先保留核心链路，暂时剥离高级扩展特性。主要边界如下：

#### 保留
- `liblagrangian+`
- `libdsmcFoam+`
- `dsmcFoam+`
- 最小 `dsmcCloud`
- 最小 `dsmcParcel`
- 基础碰撞与统计必需状态
- 基础初始化与基础场输出链路

#### 暂缓或剥离
- `dynamicLoadBalancing`
- `AMR`
- `axisymmetric` / `spherical` 等扩展坐标系
- `variableTimeStepModel`
- 大量反应模块
- 高级碰撞伙伴选择模块
- 大量扩展边界模块
- 宏观性质统计扩展
- 电子/QK 等进一步物理模型

---

## 4. 编译环境

根据记录，迁移与验证至少涉及以下环境信息：

### 4.1 OpenFOAM 环境
- 目标环境使用 `OpenFOAM-v2506`
- 实际编译时并未直接依赖源码树内未完成构建的 `wmake` 工具链，而是切换到已经可用的外部环境：
  - 环境命令：`of2506env`
  - 已存在可执行参考：`/home/superxcx/code/OpenFoam/OF-2506/OpenFOAM-v2506/platforms/linux64IcxDPInt32Opt/bin/dsmcFoam`

### 4.2 编译器/平台信息
从已知路径可以确认目标平台至少包含如下构建目标信息：
- 平台目录：`platforms/linux64IcxDPInt32Opt`
- 说明当前可用目标环境为：
  - Linux 64 位
  - Intel oneAPI/`icx` 系列编译环境
  - Double Precision
  - 32-bit label/int
  - Opt 优化构建

### 4.3 环境问题记录
迁移早期曾遇到宿主环境不可直接编译的问题：
- 当前 `OpenFOAM-v2506` 源码树缺少 `wmake` 生成的 `platforms/tools/.../wmkdepend`
- 因此后续改为使用已经预编译好的 `of2506env` 环境继续移植

> 注：具体 `MPI` 版本、`wmake` 变量、系统发行版、完整编译器版本号在当前记录中未完整给出，正式版文档可在后续补充。

---

## 5. 主要修改文件

根据移植记录，当前阶段可确认的核心修改文件如下。

### 5.1 目录基线
- 新工作目录：`hyStrath_xcx`
- 该目录作为 `dsmcFoam+` 向 `OpenFOAM-v2506` 迁移的后续开发目录

### 5.2 核心云与求解器骨架
- `hyStrath_xcx/src/lagrangian/dsmc/clouds/dsmcCloud.H`
- `hyStrath_xcx/src/lagrangian/dsmc/clouds/dsmcCloud.C`
- `hyStrath_xcx/applications/solvers/discreteMethods/dsmc/dsmcFoam+/dsmcFoam+.C`

### 5.3 最小 parcel 适配
- `hyStrath_xcx/src/lagrangian/dsmc/parcels/dsmcParcel.H`
- `hyStrath_xcx/src/lagrangian/dsmc/parcels/dsmcParcel.C`
- `hyStrath_xcx/src/lagrangian/dsmc/parcels/dsmcParcelIO.C`

### 5.4 类型注册与链接修复
- `hyStrath_xcx/src/lagrangian/dsmc/dsmcTypes.C`

### 5.5 构建边界与库组织
- `hyStrath_xcx/src/lagrangian/basic/Make/files`
- `hyStrath_xcx/src/lagrangian/dsmc/Make/files`
- 以及对应的 `Make/options`

### 5.6 第二阶段最小可运行链路中出现的补充文件
在进一步打通最小运行算例时，还修改/补充了以下文件：
- `hyStrath_xcx/src/lagrangian/dsmc/collisions/derived/LarsenBorgnakkeVariableHardSphere/LarsenBorgnakkeVariableHardSphere.H`
- `hyStrath_xcx/src/lagrangian/dsmc/collisions/derived/LarsenBorgnakkeVariableHardSphere/LarsenBorgnakkeVariableHardSphere.C`
- `hyStrath_xcx/applications/utilities/preProcessing/dsmc/dsmcInitialise+/dsmcInitialise+.C`

> 注：如果后续要形成正式的“修改文件清单”，建议再补一份 `git diff --name-only` 导出表，作为附录。

---

## 6. 关键兼容性问题

本次迁移的核心难点不是求解器主程序，而是 `OpenFOAM-v1706` 到 `OpenFOAM-v2506` 宿主 `lagrangian` API 的跨版本变化。已确认的关键兼容性问题如下。

### 6.1 原 `hyStrath` 的 `lagrangian/basic` 旧接口无法直接复用
早期编译已经表明，阻塞点首先出现在 `hyStrath` 自带的旧版 `lagrangian/basic`，尤其是：
- `Cloud.C`
- `particleI.H`
- `particleTemplates.C`

这些代码仍按 `v1706` 风格使用：
- `globalMeshData`
- `tetIndices`
- 旧的粒子/网格映射接口

而 `OpenFOAM-v2506` 的宿主接口已显著变化，导致旧底座直接移植成本过高。

### 6.2 `particle/DSMCParcel` 底层追踪机制变化
`OpenFOAM-v2506` 的粒子追踪已经转向 **barycentric** 机制。  
这意味着：
- 粒子内部状态表示与 `v1706` 不同；
- `particle` / `Cloud` / `autoMap()` / 追踪相关调用方式发生变化；
- `hyStrath` 旧版 `particle`/`Cloud` 底座不适合原样覆盖 `v2506`。

这也是迁移策略转向“采用 `v2506` 原生 `lagrangian/basic` 底座，只保留上层 `dsmcFoam+` 特有逻辑”的关键原因。

### 6.3 `Cloud::autoMap()` / 动网格映射机制变化
记录中已提前识别出高风险点：
- `v2506` 的 `Cloud::autoMap()` 依赖预先调用 `storeGlobalPositions()`
- `hyStrath` 原有 `dsmcCloud::autoMap(const mapPolyMesh&)` 基于旧流程实现

因此：
- 主迁移阶段不应直接绑定 `AMR`
- `AMR` 与动态重映射必须在后续独立阶段重新适配

### 6.4 `dsmcCloud` 对高级模块存在硬依赖
第一阶段开始时已确认：
- 现有 `dsmcCloud.H` 对第二阶段模块存在硬依赖
- 因此无法通过简单复制进入 `v2506`
- 必须先做“瘦身版 `dsmcCloud`”或提供 `stub/no-op` 级替代实现

### 6.5 类型注册、模板实例化与链接顺序问题
在最小版构建过程中，还出现了典型的 OpenFOAM 风格链接问题：
- 对象文件中缺失 `typeName/debug` 定义
- 需要采用更贴近 `v2506` 的注册方式或显式模板实例化
- 求解器链接顺序需要调整为先看到 `libdsmcFoam+`，再链接 `liblagrangian+`，否则 `Cloud<dsmcParcel>` 等符号无法正确收束

### 6.6 并行归约接口兼容问题
在进一步打通最小运行链路时，曾遇到：
- `v2506` 对 `Field` 的 `reduce(sumOp<Field>)` 不接受
- 处理方式改为逐单元归约，以保证并行兼容性

### 6.7 初始粒子坐标写盘格式问题
算例运行时还暴露出一个与 `v2506` barycentric 追踪直接相关的问题：
- `0/lagrangian/dsmc/coordinates` 中写出的值为 `(-1e+300 ...)`
- 这会导致 `v2506` 的 barycentric 追踪在第一步直接崩溃

因此在第二阶段最小可运行链路中，初始写盘格式被纳入必修复问题。

---

## 7. 移植策略与实现思路

### 7.1 总体策略
本次迁移采用了“**先保底座，后恢复功能**”的策略，而不是试图一次性把全部 `hyStrath` 特性平移到 `v2506`。具体思路为：

1. 先在 `hyStrath_xcx` 中建立新工作副本；
2. 明确第一阶段只追求最小编译和最小运行；
3. 不继续硬补 `v1706` 风格的 `lagrangian/basic`；
4. 改为使用 `v2506` 原生 `particle/Cloud` 底座；
5. 将 `dsmcParcel`、`dsmcCloud` 缩减为最小适配层；
6. 后续再逐步恢复碰撞、边界、初始化与统计链路；
7. `AMR` 与动态负载均衡单独另开阶段处理。

### 7.2 第一阶段实施结果
第一阶段得到的不是完整功能版 `dsmcFoam+`，而是一个面向 `OpenFOAM-v2506` 的**最小核心版**：
- `liblagrangian+` 可编译；
- `libdsmcFoam+` 可编译；
- `dsmcFoam+` 求解器可编译；
- 求解器可正常执行 `-help`。

### 7.3 第二阶段最小运行链路
在后续进一步恢复基础功能时，打通了：
- `dsmcInitialise+` 初始化
- `dsmcFoam+` 真实时间推进
- `noTimeCounter`
- `LarsenBorgnakkeVariableHardSphere`
- `dsmcDiffuseWallPatch`
- `dsmcFreeStreamInflowPatch`
- `dsmcMeshFill`
- `dsmcVolFields` 所需基础输出链路

---

## 8. 编译方法

以下编译方法依据记录整理，适合写入当前草案。正式版可在你本地再核对后替换为最终命令。

### 8.1 环境准备
先进入已经可用的 `OpenFOAM-v2506` 环境：
```bash
of2506env
```

### 8.2 进入迁移目录
```bash
cd /home/superxcx/code/DSMC/dsmcFoam++
```

### 8.3 以 `hyStrath_xcx` 作为迁移代码目录
后续所有编译均以 `hyStrath_xcx` 为主目录，而不是继续直接修改原始 `hyStrath`。

### 8.4 编译顺序建议
根据记录，推荐优先按以下顺序编译：
1. `liblagrangian+`
2. `libdsmcFoam+`
3. `dsmcFoam+`
4. 需要时再编译 `dsmcInitialise+`

### 8.5 编译注意事项
- 不建议直接依赖未构建完成的 `OpenFOAM-v2506` 源码树工具链；
- 应使用已可工作的 `of2506env`；
- 若出现链接异常，应优先检查：
  - `FOAM_USER_LIBBIN`
  - `Make/options`
  - 库链接顺序
  - 类型注册与模板实例化

> 注：记录中未保留每一步完整 `wmake libso`/`wmake` 命令串，正式版建议你再补一节“精确编译命令”。

---

## 9. 运行方法

### 9.1 第一阶段最小验证
第一阶段完成后，至少可执行：
```bash
dsmcFoam+ -help
```
用于验证：
- 可执行文件已正确生成；
- 运行时链接无缺库；
- 基本求解器入口工作正常。

### 9.2 第二阶段最小算例运行
后续记录表明，已经在 `orion107kmNR` 算例上打通最小运行链路，流程至少包括：

1. 使用 `dsmcInitialise+` 进行初始化；
2. 使用 `dsmcFoam+` 进入真实时间推进；
3. 输出基础 `dsmcVolFields` 相关场。

### 9.3 运行注意事项
- 当前阶段建议优先使用已生成好的网格；
- 不建议把 `snappyHexMesh` 问题与求解器迁移问题混在一起处理；
- 若算例在第一步即崩溃，应优先检查：
  - 初始粒子坐标写盘格式
  - barycentric 追踪兼容性
  - parcel 初始状态是否完整

---

## 10. 验证结果

### 10.1 第一阶段验证结果
可确认结果如下：
- `hyStrath_xcx` 已建立为迁移目录；
- `liblagrangian+` 已在 `v2506` 环境下成功出库；
- `libdsmcFoam+` 已成功编译；
- `dsmcFoam+` 已成功编译；
- `dsmcFoam+ -help` 能正常启动。

这说明：
- `OpenFOAM-v2506` 宿主 API 已被基本接通；
- 最小版 `particle` / `Cloud` / `parcel` / `cloud` / solver 链路成立；
- 第一阶段“最小可编译”目标完成。

### 10.2 第二阶段验证结果
可确认结果如下：
- `orion107kmNR` 算例可使用 `dsmcInitialise+` 完成初始化；
- `dsmcFoam+` 可真实进入时间推进；
- 程序不再在第一步因坐标/追踪问题立即崩溃。

这说明：
- 基础物理链路已不止停留在编译通过；
- 最小可运行路径已初步建立。

### 10.3 尚未完成的验证
当前记录中**尚不能确认**以下项目已经系统完成：
- 与 `OpenFOAM-v1706` 原版/旧版 `dsmcFoam+` 的逐项数值对比
- 长时间积分稳定性验证
- 并行一致性验证
- 多算例回归验证
- 性能对比验证
- `AMR` / `dynamicLoadBalancing` 回归验证

因此，当前“验证结果”更准确的表述应为：
> 已完成“可编译 + 最小可运行”验证，尚未完成“完整数值一致性与功能完备性”验证。

---

## 11. 已知遗留问题

### 11.1 AMR 尚未纳入主迁移
由于 `autoMap()`、全局位置存储和动网格映射逻辑变化明显，`AMR` 需要后续单独适配，暂不能认为已经完成。

### 11.2 dynamicLoadBalancing 尚未恢复
动态负载均衡已被明确从主迁移链路中拆出，后续应单独处理数据迁移、重分布和粒子重挂接逻辑。

### 11.3 高级物理模型尚未完整恢复
当前记录只确认恢复了最小碰撞、初始化、入口/壁面和基础场输出链路；大量高级模型仍未回迁。

### 11.4 数值一致性尚未形成正式结论
虽然程序已能跑通最小链路，但尚不能据此直接宣称与 `v1706` 版本数值完全一致。

### 11.5 构建命令与环境信息仍需补齐
正式技术文档建议补充：
- 完整 `wmake` 命令
- 环境变量
- 编译器版本号
- MPI 版本号
- 依赖库版本
- 运行节点/系统环境

### 11.6 网格生成工具链问题与求解器问题需分离
记录中提到 `snappyHexMesh` 存在 `ptscotch` 动态库符号冲突。  
该问题属于外部工具链/环境问题，不应与 `dsmcFoam+` 主迁移的成败混为一谈，但仍建议在“已知问题”中保留说明。

---

## 12. 当前结论

基于现有记录，可以给出如下阶段性结论：

1. `dsmcFoam+` 从 `OpenFOAM-v1706` 向 `OpenFOAM-v2506` 的迁移已经完成了最关键的一步：  
   **在 `hyStrath_xcx` 中建立了面向 `v2506` 的最小核心版实现。**

2. 当前已完成的层级是：  
   **编译打通 + 最小运行打通。**

3. 当前尚未完成的层级是：  
   **完整功能恢复、完整物理模型恢复、数值一致性验证、AMR/动态负载均衡恢复。**

因此，现阶段最合适的技术表述是：

> 已完成 `dsmcFoam+` 从 `OpenFOAM-v1706` 到 `OpenFOAM-v2506` 的第一阶段核心迁移，在 `hyStrath_xcx` 中构建了基于 `v2506` 原生 `lagrangian` 底座的最小版 `liblagrangian+`、`libdsmcFoam+` 与 `dsmcFoam+`，实现了编译通过和最小算例运行，为后续恢复高级物理模型、AMR 与动态负载均衡奠定了代码基础。

---

## 13. 建议补充的附录

正式版文档建议再增加以下附录：

### 附录 A：完整修改文件清单
建议由 `git diff --name-only` 生成。

### 附录 B：完整编译命令
逐条列出 `wmake libso`、`wmake` 及环境变量。

### 附录 C：算例验证记录
建议包含：
- 算例名
- 初始化是否通过
- 是否进入时间推进
- 是否输出场
- 是否并行运行
- 备注

### 附录 D：后续迁移路线
建议分为：
1. 恢复基础物理功能
2. 恢复高级边界/反应模型
3. 恢复 AMR
4. 恢复 dynamicLoadBalancing
5. 做数值一致性验证
6. 做并行与性能验证


################################################

十二、实施顺序

基础设施
dsmcCloud 增加 OpenMP 开关、线程数、线程私有 RNG
collisionSelection 去掉缓存 rndGen_
第一版碰撞并行
只改 noTimeCounter::collide()
先用线程私有 RNG
先用局部计数 + 末端 reduction
第一版验证
串行结果统计对照
8 线程/16 线程短跑
检查守恒、碰撞数、反应数、壁面通量是否量级一致
第二版真正 dual-decomposition
实现 collision weights -> multiSegmentPartition
用连续 cell 段替代简单 dynamic 调度
第三版再看 Move/Index
这时才处理 buildCellOccupancy 和 relocation