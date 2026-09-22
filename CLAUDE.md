# hyStrath_dlb 工作规范（CLAUDE）

更新时间：2026-09-08

本文件用于统一本仓库后续工作约束，优先级高于临时记忆/旧 worklog，任何新增工作都以此为准。

## 0. 范围与关键硬约束

1. 本项目工作目录固定为：
   `/home/superxcx/code/OpenFoam/OF-1706/hyStrath_dlb`
2. 参考代码仓库只读，不得修改：
   `/home/superxcx/code/DSMC/dsmcFoam++/hyStrath_xcx`
3. 优先以当前仓库源码、当前 `controlDict` 和当前日志为判定依据；历史日志仅作经验参考。
4. 各次性能/正确性实验必须使用独立结果目录，不得覆盖已有有效日志、结果与备份。
5. 任何新结论必须在源码、输入、网格、参数、命令、运行环境、进程/线程设置、日志齐全时才可提交。

## 1. 固定路径与工作流入口

```bash
cd /home/superxcx/code/OpenFoam/OF-1706/hyStrath_dlb
```

主要工作日志目录：

- `doc/worklog/v2506/`
- `doc/worklog/v2506/detail_omp/`
- `doc/worklog/v2506/detail_mpi/`
- `doc/worklog/v2506/detail_mix/`

核心目标：优先按以下顺序推进（如有正确性或环境问题可临时打断）：

1. OMP 路径与 OMP 负载均衡优化
2. MPI replicated mesh 与 DLB
3. MPI+OpenMP mixed
4. 完整性能剖析 + 正确性验证 + 文档归档

## 2. 环境与编译

进入环境：

```bash
source /home/superxcx/code/OpenFoam/OF-1706/hyStrath_dlb/doc/scripts/env.sh
```

编译：

```bash
bash doc/scripts/build-dsmcFoam.sh
```

构建产物应落在：

- `platforms/linux64IccDPInt32Opt/bin/dsmcFoam+`
- `platforms/linux64IccDPInt32Opt/bin/dsmcInitialise+`

运行或测试前必须确认二进制来源：

```bash
which dsmcFoam+
which dsmcInitialise+
which mpirun
```

`dsmcFoam+` 必须解析到本仓库 `platforms/.../bin`；`mpirun` 必须来自当前 oneAPI/OpenFOAM 环境，禁止混用系统 MPI。

## 3. 配置与默认运行边界

- FastRng 统一使用：
  - `fastRng true;`
  - 不再使用：
    - `collisionFastRng`
    - `moveFastRng`
- 默认性能测试 profile：
  - `profileSummary true;`
  - `profileDetail false;`
- `profileDetail true` 仅用于迁移/通信/DLB 的小范围逐步调试，不作为大规模默认设置。
- 当前推荐生产正确性边界：
  - `writeControl runTime;`
  - `writeInterval 1.e-3;`
  - `replicatedMeshMigrateInterval 1;`
  - `replicatedMeshWriteMode processor;`
- `replicatedMeshMigrateInterval 10` 仅作旧对照或诊断配置，不作为当前默认。

## 4. 运行方式（必须严格）

### 4.1 OpenMP 8 核（默认基线）

```bash
source /home/superxcx/code/OpenFoam/OF-1706/hyStrath_dlb/doc/scripts/env.sh
export OMP_NUM_THREADS=8
export OMP_DYNAMIC=false
dsmcFoam+ > log.omp8 2>&1
```

### 4.2 MPI replicated mesh 8 核

- 不执行 `decomposePar`
- 不加 `-parallel`
- 粒子迁移使用 replicated mesh owner/migration 路径

```bash
source /home/superxcx/code/OpenFoam/OF-1706/hyStrath_dlb/doc/scripts/env.sh
export OMP_NUM_THREADS=1
export OMP_DYNAMIC=false
mpirun -np 8 dsmcFoam+ > log.mpi8 2>&1
```

### 4.3 MPI+OpenMP mixed（常用）

```bash
source /home/superxcx/code/OpenFoam/OF-1706/hyStrath_dlb/doc/scripts/env.sh
export OMP_NUM_THREADS=2
export OMP_DYNAMIC=false
export OMP_PROC_BIND=close
export OMP_PLACES=cores
mpirun -np 4 dsmcFoam+ > log.mpi4omp2 2>&1
```

```bash
source /home/superxcx/code/OpenFoam/OF-1706/hyStrath_dlb/doc/scripts/env.sh
export OMP_NUM_THREADS=4
export OMP_DYNAMIC=false
export OMP_PROC_BIND=close
export OMP_PLACES=cores
mpirun -np 2 dsmcFoam+ > log.mpi2omp4 2>&1
```

### 4.4 标准 OpenFOAM decomposed MPI 对照（仅 `mpi8origin`）

```bash
source /home/superxcx/code/OpenFoam/OF-1706/hyStrath_dlb/doc/scripts/env.sh
export OMP_NUM_THREADS=1
export OMP_DYNAMIC=false
decomposePar -force > log.decomposePar 2>&1
mpirun -np 8 dsmcFoam+ -parallel > log.mpi8origin 2>&1
```

### 4.5 mixed / replicated mesh 日志必检

`MPI4xOMP2`/`MPI2xOMP4` 的结果必须出现：

```
Replicated mesh: initialized with ... MPI ranks
OpenMP max threads = ...
Total Iterations = 300
End main
```

## 5. 核心算例与目录

- 历史主基准：
  - `run/hyStrath/dsmcFoam+/xcx_test/pal-phd-3.3.1-cylinder-react/timecompare/ourmesh`
- 历史备份：
  - `run/hyStrath/dsmcFoam+/xcx_test/pal-phd-3.3.1-cylinder-react/timecompare/ourmeshbkp`
- 当前验证集：
  - `run/hyStrath/dsmcFoam+/xcx_test/zb-cylinder-react-validate`

`zb-cylinder-react-validate` 保持 5 个干净子 case：

- `omp8`：8 核 OMP
- `mpi2omp4`：replicated mesh mixed
- `mpi4omp2`：replicated mesh mixed（常用）
- `mpi8`：replicated mesh pure MPI
- `mpi8origin`：标准 decomposed MPI 对照

验证集已去除旧 slurm、初始脚本、历史日志、备份、schedule sweep、旧 processor 目录和异常时间目录，仅保留 case 本体。

## 6. 正确性检查（性能结论前必做）

1. 求解器正常结束：出现 `Total Iterations = ...` 与 `End main`
2. 无 `fatal`、`segfault`、`MPI abort`、NaN、`floating point exception`
3. `stuck particles = 0`
4. 粒子数非零且与对照模式在可接受范围
5. 碰撞数、总能量、温度、密度需与对照可比；有偏差必须可解释
6. mixed/replicatedMesh 日志必须确认实际 `rank` 与 `thread` 数

不得将后台运行成功、错误二进制、错误 `OMP_NUM_THREADS` 或错误 `-parallel` 下的计时作为有效结果。

## 7. 性能记录与归因

每次性能实验保留完整 `dsmcFoam+` 标准输出，并在独立结果目录保存：

- `log`
- `controlDict`
- `decomposeParDict`
- 运行命令
- 环境说明
- 结果 `CSV/Markdown`

重点关注并完整记录：

- 外部 `real/user/sys`
- `full evolve`
- `move`
- `build_occupancy`
- `coll`
- `post`
- `comm total`
- `migration wall`
- DLB checks/rebalances
- `particles per rank`
- `rank wall max/min`
- `final particles / stuck particles / collisions / total energy`

性能比较必须区分：

1. 源码版本
2. `controlDict` 配置
3. 运行命令与 MPI/OMP 启动状态
4. 复现实验（至少三次）统计
5. 正确性和性能归因链路

单次快/慢不构成优化结论；必须使用固定配置的 A/B 重复统计来验证 fastRng、DLB、migration interval、pinning 等改动效果。

## 8. 实验归档规范（建议目录）

结果目录建议包含：

- `run.log`（启动与命令）
- `environment.txt`（模块、`env`、`mpirun` 来源）
- `controlDict`、`decomposeParDict`（如有）
- `solver.stdout`
- `timing.csv`（或等价）
- `correctness.md`

发现问题、失败、超时、无收益时保留原始日志，不删除、不覆盖、不得“清理即证明修复”。

## 9. 禁止事项

- 不要修改 `/home/superxcx/code/DSMC/dsmcFoam++/hyStrath_xcx`
- 不要在 replicated mesh 运行中加 `-parallel`
- 不要把 `replicatedMeshMigrateInterval` 当作 DLB 检查周期解释
- 不要只看 `real` 时间做结论，必须通过正确性检查
- 不要在新配置中使用 `collisionFastRng`、`moveFastRng`
- 不要用 `profileDetail true` 做大规模默认性能测试
- 不要删除或覆盖用户已有有效日志与手工整理文件

---

以上规则按“先正确、后计时、再归档”执行；默认用 2026-07-21 的 `dsmcFoam++ OFv1706 hyStrath_dlb` 工作入口清单更新到本项目 AGENTS。
