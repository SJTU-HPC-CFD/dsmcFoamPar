# dsmcFoam+ OFv1706 hyStrath_dlb 工作入口清单

更新时间：2026-07-21

本文是本仓库后续工作的第一入口。开始任何代码修改、算例测试或性能分析前，先按
这里确认环境、case、运行方式和记录规则。

## 1. 工作目标

本轮工作是在 OFv1706 版本上复现并继续优化此前 OF-v2506 版本中的
`dsmcFoam+` 性能优化路线，重点包括：

1. OpenMP 路径和 OMP 负载均衡优化；
2. MPI replicated mesh 路径和 DLB 优化；
3. MPI+OpenMP mixed 路径优化；
4. 完整性能剖析、正确性验证和可复现实验记录。

最终目标不是单次跑得更快，而是在正确性可控、日志完整、配置可复现的前提下，
量化 move、collision、build occupancy、migration、DLB、I/O 等阶段的优化收益。

## 2. 固定路径

工作目录：

```bash
cd /home/superxcx/code/OpenFoam/OF-1706/hyStrath_dlb
```

参考代码：

```text
/home/superxcx/code/DSMC/dsmcFoam++/hyStrath_xcx
```

参考代码只读，不要修改。历史工作日志可以用于理解思路，但后续判断以当前仓库
源码、参考代码实际实现、当前 `controlDict` 和新测试日志为准。

主要工作日志目录：

```text
doc/worklog/v2506
doc/worklog/v2506/detail_omp
doc/worklog/v2506/detail_mpi
doc/worklog/v2506/detail_mix
```

## 3. 环境和编译

激活环境：

```bash
source /home/superxcx/code/OpenFoam/OF-1706/hyStrath_dlb/doc/scripts/env.sh
```

编译：

```bash
source /home/superxcx/code/OpenFoam/OF-1706/hyStrath_dlb/doc/scripts/build-dsmcFoam.sh
```

编译脚本内部使用 `wmake -j`，输出应落在本仓库本地路径：

```text
platforms/linux64IccDPInt32Opt/bin/dsmcFoam+
platforms/linux64IccDPInt32Opt/bin/dsmcInitialise+
```

运行或测试前检查二进制来源：

```bash
which dsmcFoam+
which dsmcInitialise+
which mpirun
```

`dsmcFoam+` 必须解析到本仓库 `platforms/.../bin`，`mpirun` 应来自当前 oneAPI/OpenFOAM
环境，不要混用系统 MPI。

## 4. 当前控制项约定

当前 FastRng 只保留统一开关：

```text
fastRng true;
```

不要再使用旧控制项：

```text
collisionFastRng
moveFastRng
```

大规模或性能测试默认使用低开销 profile：

```text
profileSummary true;
profileDetail false;
```

`profileDetail true` 只用于迁移、通信或 DLB 的逐步调试，不能作为大规模默认设置。

当前 production correctness 边界优先使用：

```text
writeControl runTime;
writeInterval 1.e-3;
replicatedMeshMigrateInterval 1;
replicatedMeshWriteMode processor;
```

历史 `replicatedMeshMigrateInterval 10` 只作为诊断和旧结果对照，不作为当前最新
默认边界。

## 5. 运行模式

### 5.1 OpenMP 8 核

```bash
source /home/superxcx/code/OpenFoam/OF-1706/hyStrath_dlb/doc/scripts/env.sh
export OMP_NUM_THREADS=8
export OMP_DYNAMIC=false
dsmcFoam+ > log.omp8 2>&1
```

### 5.2 MPI replicated mesh 8 核

replicated mesh 是 raw-MPI replicated-mesh 路径：

- 不执行 `decomposePar`；
- 不加 `-parallel`；
- 粒子迁移由 replicated mesh owner/migration 路径处理。

```bash
source /home/superxcx/code/OpenFoam/OF-1706/hyStrath_dlb/doc/scripts/env.sh
export OMP_NUM_THREADS=1
export OMP_DYNAMIC=false
mpirun -np 8 dsmcFoam+ > log.mpi8 2>&1
```

### 5.3 MPI+OpenMP mixed

`MPI4xOMP2` 是当前常用 mixed 验证模式：

```bash
source /home/superxcx/code/OpenFoam/OF-1706/hyStrath_dlb/doc/scripts/env.sh
export OMP_NUM_THREADS=2
export OMP_DYNAMIC=false
export OMP_PROC_BIND=close
export OMP_PLACES=cores
mpirun -np 4 dsmcFoam+ > log.mpi4omp2 2>&1
```

`MPI2xOMP4`：

```bash
source /home/superxcx/code/OpenFoam/OF-1706/hyStrath_dlb/doc/scripts/env.sh
export OMP_NUM_THREADS=4
export OMP_DYNAMIC=false
export OMP_PROC_BIND=close
export OMP_PLACES=cores
mpirun -np 2 dsmcFoam+ > log.mpi2omp4 2>&1
```

mixed replicated mesh 同样不要加 `-parallel`。

### 5.4 标准 OpenFOAM decomposed MPI 对照

只有 `mpi8origin` 这种标准 decomposed 对照 case 使用 `decomposePar` 和 `-parallel`：

```bash
source /home/superxcx/code/OpenFoam/OF-1706/hyStrath_dlb/doc/scripts/env.sh
export OMP_NUM_THREADS=1
export OMP_DYNAMIC=false
decomposePar -force > log.decomposePar 2>&1
mpirun -np 8 dsmcFoam+ -parallel > log.mpi8origin 2>&1
```

## 6. 主要算例

历史主基准：

```text
run/hyStrath/dsmcFoam+/xcx_test/pal-phd-3.3.1-cylinder-react/timecompare/ourmesh
```

历史备份和原始日志：

```text
run/hyStrath/dsmcFoam+/xcx_test/pal-phd-3.3.1-cylinder-react/timecompare/ourmeshbkp
```

当前 `zb-cylinder-react` 验证集：

```text
run/hyStrath/dsmcFoam+/xcx_test/zb-cylinder-react-validate
```

该验证集已整理为 5 个干净子 case：

| mode           | 含义                           | 运行方式                                                    |
| -------------- | ------------------------------ | ----------------------------------------------------------- |
| `omp8`       | 当前 8-core 默认最快基线       | `OMP_NUM_THREADS=8 dsmcFoam+`                             |
| `mpi2omp4`   | replicated mesh mixed          | `OMP_NUM_THREADS=4 mpirun -np 2 dsmcFoam+`                |
| `mpi4omp2`   | replicated mesh mixed 常用验证 | `OMP_NUM_THREADS=2 mpirun -np 4 dsmcFoam+`                |
| `mpi8`       | replicated mesh pure MPI       | `OMP_NUM_THREADS=1 mpirun -np 8 dsmcFoam+`                |
| `mpi8origin` | 标准 decomposed MPI 对照       | `decomposePar -force && mpirun -np 8 dsmcFoam+ -parallel` |

`zb-cylinder-react-validate` 中已去掉旧 `*.slurm`、`initial.sh`、旧 log、旧备份、
schedule sweep、processor 目录和非零时间目录，只保留 case 本体。

## 7. 正确性检查

任何性能结论前必须检查：

1. 求解器正常结束：`Total Iterations = ...` 和 `End main`；
2. 无 fatal、segfault、MPI abort、NaN、floating point exception；
3. `stuck particles = 0`；
4. 粒子数非零，并与对照模式在可接受范围内；
5. 碰撞数、总能量、温度、密度等统计量与对照模式一致或差异可解释；
6. mixed/replicatedMesh 日志必须确认实际 rank/thread 数。

`MPI4xOMP2` 日志至少应看到：

```text
Replicated mesh: initialized with 4 MPI ranks
OpenMP max threads = 2
Total Iterations = 300
End main
```

不要把后台运行、错误 MPI 栈、错误二进制、错误 `OMP_NUM_THREADS` 或错误 `-parallel`
得到的时间当作有效结果。

## 8. 性能记录规则

保留完整 `dsmcFoam+` stdout log。批量测试必须建立独立结果目录，并保存：

```text
log
controlDict
decomposeParDict
运行命令
环境说明
结果 CSV/Markdown 汇总
```

记录重点包括：

- 外部 `real/user/sys`；
- `full evolve`；
- `move`；
- `build_occupancy`；
- `coll`；
- `post`；
- `comm total`；
- `migration wall`；
- DLB checks/rebalances；
- particles per rank；
- rank wall max/min；
- final particles、stuck particles、collisions、total energy。

性能对比至少区分：

1. 源码状态；
2. `controlDict` 状态；
3. 运行命令和 MPI/OMP 实际启动状态；
4. 单次结果和重复统计；
5. 功能验证和性能归因。

单次慢/快不能直接作为优化结论。需要量化 `fastRng`、DLB、migration interval、
pinning 等效果时，使用固定配置的多重复 A/B 测试。

## 9. 工作顺序

默认推进顺序：

1. OMP 路径和 OMP DLB；
2. MPI replicated mesh DLB；
3. MPI+OpenMP mixed；
4. 完整性能剖析、文档整理和论文/报告对比。

但实际修改顺序以当前 bug、正确性风险和性能剖析热点为准。修改共享路径前先确认
是否会影响 OMP、replicatedMesh、mixed 和 `mpi8origin` 四类运行方式。

## 10. 禁止事项

- 不要修改参考代码 `/home/superxcx/code/DSMC/dsmcFoam++/hyStrath_xcx`；
- 不要在 replicated mesh 运行中使用 `-parallel`；
- 不要把 `replicatedMeshMigrateInterval` 当成 DLB 检查周期；
- 不要只看 `real` 时间而忽略正确性指标；
- 不要用旧 `collisionFastRng` / `moveFastRng` 写新配置；
- 不要用 `profileDetail true` 做大规模默认性能测试；
- 不要删除或覆盖用户已有的有效工作日志和手动下载/整理文件。
