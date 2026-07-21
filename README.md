# hyStrath_dlb 仓库阅读指南

## 项目简介

本项目基于 OpenFOAM v1706 和 hyStrath 的 `dsmcFoam+`，面向非结构网格 DSMC
求解器开展并行化与性能优化工作。主要内容包括 OpenMP 路径优化、MPI replicated
mesh 与动态负载均衡（DLB）、MPI+OpenMP 混合并行、随机数和粒子处理热点优化，
以及大规模算例的性能剖析与扩展性测试。

项目同时重视数值正确性和实验可复现性。所有优化均需要结合标准验证算例，对
粒子数、宏观场、碰撞统计和输出结果进行检查，并分别量化 move、collision、
occupancy、migration、DLB 和 I/O 等阶段的性能变化。

## 注意

任何代码修改，都要记录在 `doc/worklog/v2506` 下的日志中，并在日志中说明修改目的、修改内容、验证结果和性能变化。

## 阅读入口

首次阅读本代码仓库时，请先从项目入口开始：

```text
doc/worklog/v2506/dsmcFoam++tips.md
```

该文件记录当前工作环境、编译命令、运行方式、验证规则，以及 replicated mesh
等关键运行约束。

## 主要阅读顺序

1. 项目入口和当前工作流程：

   ```text
   doc/worklog/v2506/dsmcFoam++tips.md
   ```
2. 论文内容和原始材料：

   ```text
   paper/raw
   paper/raw/完整技术报告.md
   paper/raw/论文内容安排.docx
   paper/raw/引言逻辑详细说明.md
   ```
3. 工作日志和技术报告：

   ```text
   doc/worklog/v2506
   doc/worklog/v2506/dsmcFoam_plus_omp_technical_report_20260606.md
   doc/worklog/v2506/dsmcFoam_plus_mpi_replicated_mesh_dlb_technical_report_20260608.md
   doc/worklog/v2506/dsmcFoam_plus_mpi_replicated_mesh_further_optimization_technical_report_20260617.md
   doc/worklog/v2506/dsmcFoam_plus_mpi_omp_mixed_technical_report_20260610.md
   doc/worklog/v2506/dsmcFoam_plus_chaosuan_correctness_performance_technical_report_20260624.md
   ```
4. 当前验证算例：

   ```text
   run/hyStrath/dsmcFoam+/xcx_test/zb-cylinder-react-validate
   ```

## 使用说明

- 将 `doc/worklog/v2506/dsmcFoam++tips.md` 作为本分支编译、运行、性能剖析和
  正确性验证的操作依据。
- 通过 `doc/worklog/v2506` 了解完整的优化流程、测试过程和历史决策。
- 通过 `paper/raw` 查阅论文内容及论文撰写所需的原始材料。
- 使用 `zb-cylinder-react-validate` 进行当前版本的标准验证测试。

## 首次使用

以下流程适用于当前仓库的 Linux、Intel oneAPI、OpenFOAM v1706 环境。脚本默认
使用以下基础环境路径：

```text
~/intel/oneapi
~/code/OpenFoam/OF-1706/OpenFOAM-v1706
```

如果本机路径不同，需要先修改 `doc/scripts/env.sh` 和
`doc/scripts/build-dsmcFoam.sh` 中对应的路径。

### 1. 获取并安装 ParMETIS

`dsmcFoam+` 的 MPI 分解和相关负载均衡功能依赖 ParMETIS。仓库已经在
`package/parmetis` 中提供 GKlib、METIS 和 ParMETIS 源码压缩包，三个文件都必须
存在：

```text
package/parmetis/gklib-master.tar.gz
package/parmetis/metis-master.tar.gz
package/parmetis/parmetis-main.tar.gz
```

先加载 OpenFOAM 和编译器环境，然后执行仓库提供的安装脚本：

```bash
cd /home/superxcx/code/OpenFoam/OF-1706/hyStrath_dlb
source doc/scripts/env.sh

cd package/parmetis
bash makeParMETIS.sh
```

默认安装位置为：

```text
/home/superxcx/code/OpenFoam/OF-1706/hyStrath_dlb/package/parmetis/parmetis-install
```

安装完成后，将当前 shell 的 ParMETIS 路径指向该目录。由于现有环境脚本包含
默认路径，必须在 `source env.sh` 之后补充下面的覆盖命令：

```bash
export PARMETIS_DIR=/home/superxcx/code/OpenFoam/OF-1706/hyStrath_dlb/package/parmetis/parmetis-install
export LD_LIBRARY_PATH=$FOAM_USER_LIBBIN:$PARMETIS_DIR/lib:$FOAM_PROJECT_LIBBIN:$FOAM_PROJECT_MPI_LIBBIN:$LD_LIBRARY_PATH
```

检查安装结果：

```bash
test -f "$PARMETIS_DIR/include/parmetis.h"
ls "$PARMETIS_DIR/lib"/*metis* "$PARMETIS_DIR/lib"/*GKlib*
```

### 2. 编译 dsmcFoam+

在同一个已经配置好 ParMETIS 的 shell 中执行：

```bash
cd /home/superxcx/code/OpenFoam/OF-1706/hyStrath_dlb
source doc/scripts/build-dsmcFoam.sh
```

编译脚本会依次构建 `liblagrangian+`、`libdecompose`、`libgeneralMolecule`、
`libdsmcFoam+`、`dsmcFoam+` 和 `dsmcInitialise+`。输出文件应位于当前代码位置的：

```text
platforms/linux64IccDPInt32Opt/bin/dsmcFoam+
platforms/linux64IccDPInt32Opt/bin/dsmcInitialise+
```

### 3. 准备验证算例

以当前标准验证算例的 `mpi4omp2` 模式为例：

```bash
cd run/hyStrath/dsmcFoam+/xcx_test/zb-cylinder-react-validate/mpi4omp2
source doc/scripts/env.sh
```

### 4. 运行验证算例

`mpi4omp2` 是 replicated mesh 的 MPI+OpenMP 混合模式。该模式不执行
`decomposePar`，也不添加 `-parallel`：

```bash
export OMP_NUM_THREADS=2
export OMP_DYNAMIC=false
export OMP_PROC_BIND=close
export OMP_PLACES=cores
mpirun -np 4 dsmcFoam+ > log.dsmcFoam.mpi4omp2 2>&1
```

其他标准模式如下：

```bash
# OpenMP 8 核
cd ../omp8
export OMP_NUM_THREADS=8
export OMP_DYNAMIC=false
export OMP_PROC_BIND=close
export OMP_PLACES=cores
dsmcFoam+ > log.dsmcFoam.omp8 2>&1

# replicated mesh 8 个 MPI rank，不使用 -parallel
cd ../mpi8
export OMP_NUM_THREADS=1
export OMP_DYNAMIC=false
mpirun -np 8 dsmcFoam+ > log.dsmcFoam.mpi8 2>&1

# 标准 decomposed MPI 对照模式
cd ../mpi8origin
bash clean.sh
decomposePar > log.decomposePar 2>&1
mpirun -np 8 dsmcFoam+ -parallel > log.dsmcFoam.mpi8origin 2>&1
```

性能测试时保留完整日志，并同时记录源码版本、控制文件、MPI rank 数、
OpenMP 线程数、写出模式和实际二进制路径。

### 5. 首次运行检查

运行完成后，至少检查以下内容：

```bash
tail -n 80 log.dsmcFoam.mpi4omp2
rg -n "ExecutionTime|move phase|collision phase|DLB|Fatal|error" log.dsmcFoam.mpi4omp2
```

重点确认没有 `Fatal error`、MPI 异常或粒子数异常，并核对 `move`、`collision`、
`migration`、`DLB` 和总运行时间等性能指标。正式性能对比前，应先完成短步数的
正确性验证，再使用统一配置进行多次重复测试。
