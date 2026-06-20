# METIS MPI Minimal Validation

目的：把 replicated-mesh 初始化里的
`METIS_PartGraphKway()` 单独拉出来验证，不依赖 OpenFOAM 求解器。

## 文件

- `metis_mpi_minimal.cpp`
  - MPI 程序
  - 仅 rank 0 读取 `owner/neighbour`
  - 构造与 `dsmcReplicatedMesh::computeCellOwnerScotch()` 一致的 CSR
  - 默认自动推导 `nCells = max(owner, neighbour) + 1`
  - `--cells` 仅作为可选覆盖项保留
  - 调用 `METIS_PartGraphKway()`
  - 支持 `--fpe` 打开 `FE_DIVBYZERO | FE_INVALID | FE_OVERFLOW`

- `build_metis_mpi_minimal.sh`
  - 编译器和 MPI 路径都从当前已激活环境获取
  - 只需要显式设置 `PARMETIS_DIR`
  - 优先使用 `mpiicpc`，没有则退回 `mpiicpx`

- `run_metis_mpi_minimal_zb64.sh`
  - 对 `zb-cylinder-react/mpi8` 网格做 `np=64` 验证

## 使用

```bash
source /home/superxcx/code/OpenFoam/OF-1706/hyStrath_dlb/doc/scripts/env.sh
cd /home/superxcx/code/OpenFoam/OF-1706/hyStrath_dlb/doc/worklog/v2506/detail_mpi/metis_minimal_validation_20260618
bash build_metis_mpi_minimal.sh
bash run_metis_mpi_minimal_zb64.sh
```

如果你要切换到别的 PARMETIS 安装，只改这一项：

```bash
export PARMETIS_DIR=/path/to/parmetis-install
bash build_metis_mpi_minimal.sh
NP=64 bash run_metis_mpi_minimal_zb64.sh
```

也可以手动指定：

```bash
mpirun -np 64 ./metis_mpi_minimal \
  --owner /path/to/owner \
  --neighbour /path/to/neighbour \
  --parts 64 \
  --fpe \
  --verbose
```

如需强制覆盖自动推导的 cell 数，也可以手动加：

```bash
mpirun -np 64 ./metis_mpi_minimal \
  --owner /path/to/owner \
  --neighbour /path/to/neighbour \
  --cells 60000 \
  --parts 64
```

## 结果解释

- 如果这个最小程序也在超算上 `SIGFPE`：
  - 问题已经与 OpenFOAM 主程序无关
  - 优先检查超算上的 `libmetis/libGKlib`、MPI wrapper、运行时库路径

- 如果这个最小程序正常，而 `dsmcFoam+` 初始化仍然在
  `computeCellOwnerScotch()` 崩：
  - 问题更可能在 OpenFOAM/replicatedMesh 初始化上下文
  - 下一步应把求解器里的 `metis` 初始化路径改成 fallback / wrapped path
