# zb-cylinder-react 300-step 性能重测（2026-06-22）

## 口径

- 单轮 5 组：`OMP8`、`MPI2xOMP4`、`MPI4xOMP2`、`MPI8`、`MPI8origin`。
- 300 steps：
  - `deltaT 6.640048894e-08`
  - `endTime 1.9920146682e-05`
- no-write：
  - `writeInterval 1.e-3`
  - `writeInterval > endTime`，运行中不触发 output。
- replicatedMesh 组使用当前正确性口径：
  - `replicatedMeshMigrateInterval 1`
  - `replicatedMeshWriteMode processor`
  - `replicatedMeshProcessorWriteTimeMesh false`
  - `replicatedMeshAutoDLB true`
  - `replicatedMeshDLBAlpha 1`
  - `replicatedMeshDLBMinRemainingSteps 50`
- 每组复制到本目录 `cases/` 下运行，不修改 `run/hyStrath/.../xcx_test` 原始 case。

## 结果

| mode | np | OMP/rank | real [s] | full evolve [s] | move+collide [s] | move [s] | buildOcc [s] | collision [s] | migration [s] | DLB rebalances | ok |
|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---|
| OMP8 | 1 | 8 | 55.63 | 53.37 | 53.19 | 38.97 | 3.71 | 9.41 | - | - | yes |
| MPI2xOMP4 | 2 | 4 | 63.37 | 61.66 | 61.47 | 39.58 | 3.40 | 8.28 | 11.07 | 2 | yes |
| MPI4xOMP2 | 4 | 2 | 67.16 | 65.58 | 65.37 | 40.78 | 3.33 | 13.19 | 8.16 | 2 | yes |
| MPI8 | 8 | 1 | 101.03 | 100.56 | 100.30 | 64.51 | 8.42 | 14.37 | 17.56 | 5 | yes |
| MPI8origin | 8 | 1 | 114.59 | 122.65 | 122.53 | 91.30 | 10.73 | 26.52 | - | - | yes |

## 结论

- 5 组均正常结束，`Total Iterations = 300`，`End main`，无 fatal/segfault，`stuck=0`。
- 当前正确性口径 `replicatedMeshMigrateInterval=1` 下，MPI8 replicatedMesh 仍快于
  原始 decompose 版 MPI8origin：
  - real: `101.03 s` vs `114.59 s`
  - full evolve: `100.56 s` vs `122.65 s`
- 与 OMP8 相比，replicatedMesh MPI8 明显变慢，主要来自每步 migration 和 move/build
  成本增加：
  - MPI8 migration wall `17.56 s`
  - MPI8 move `64.51 s`，OMP8 move `38.97 s`
- 混合并行下本轮最佳是 `MPI2xOMP4`：
  - real `63.37 s`
  - full evolve `61.66 s`
  - 但仍慢于 OMP8 的 `55.63 s`
- 本轮不能直接和 2026-06-17 的 `migrateInterval=10` 历史最佳作同口径性能比较；
  那一轮口径后来已被判定不满足当前 raw-MPI replicatedMesh 的流场正确性约束。

## 文件

```text
run_zb_full_series_300step_retest_20260622.sh
parse_zb_full_series_300step_retest_20260622.py
run_manifest.tsv
results.csv
aggregate.csv
logs/
controlDicts/
cases/
```
