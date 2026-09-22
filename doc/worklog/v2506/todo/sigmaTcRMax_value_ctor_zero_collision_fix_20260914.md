# dsmcSigmaTcRMax 值构造器读失效（全程零碰撞）— 2026-09-14

## 症状

超算 `pal-phd-3.3.1-react/template/mpi16_omp4_check`（job 4617247，Sep 14 16:06，
int32 修复后的新二进制）：求解器正常运行到 0.0005，但每个 info 块打印
**"No collisions"**——全窗口全局零碰撞候选。参考正确结果为同目录
`mpi16_omp4`（job 4018927，bkp 代码）。

本地 zb-cylinder-react-validate/mpi4omp2 回归同样中招（30×"No collisions"），
即 9/14 的 int32 修复验证是假阳性（门检查无碰撞项，已补）。

## 根因链（v1706 源码级）

1. processor checkpoint restart 功能包将 `sigmaTcRMax_` 从**读构造器**改为
   **值构造器**（`dimensionedScalar("zero", ..., 0.0)`），readOption 由
   helper `dsmcSigmaReadOption()` 提供，fresh start 分支返回 `MUST_READ`
2. v1706 值构造器末尾调 `readIfPresent()`
   （GeometricField.C:213 等 6 处构造器同款）
3. `readIfPresent()` 对 `MUST_READ` **只打警告不读文件**
   （GeometricField.C:96-104：MUST_READ 分支无 read 调用，直接
   `return false`）——16 条 "suggests that a read constructor ... would be
   more appropriate" 警告即此
4. 字段全 0：`0/dsmcSigmaTcRMax`（1.82e-16）被无视
5. NTC `selectedPairs ∝ sigmaTcRMax`（noTimeCounter.C:1006-1010）→
   nCandidates = 0 → 全程零碰撞；sigmaTcRMax 的运行期更新只在候选对循环
   内（noTimeCounter.C:1119）→ **0 是吸收态，永不自愈**
6. 铁证：check 写出 `processor0/0.0005/dsmcSigmaTcRMax` = `uniform 0`
   （参考版为 nonuniform 正常场）；move/inflow/能量不受影响 → 无 FATAL
   静默错误

另注：试图用 `sigmaTcRMax_.read()` 显式补救失败——`regIOobject::readData`
默认实现 `return false` 且 GeometricField 未重载；`readFields()` 为 private
不可外部调用。

## 修复

`dsmcCloud.C`：

1. `dsmcSigmaReadOption()` fresh start 分支 `MUST_READ` →
   **`READ_IF_PRESENT`**（v1706 "可选读"标准惯用法：值构造器 +
   READ_IF_PRESENT → readIfPresent 真正执行 readFields 读盘）；
   restart 分支保持 READ_IF_PRESENT（无文件→零初值，随后
   `restoreProcessorCheckpoint()` 显式恢复，该路径本就正确）
2. 构造器体加响亮守卫：fresh start 若 `0/dsmcSigmaTcRMax` 与
   processor checkpoint 均不存在 → FatalError（保持历史 MUST_READ 的
   强约束；静默零场 = 零碰撞陷阱）

## 验证（zb-cylinder-react-validate/mpi4omp2，300 步）

| 指标 | 零碰撞 bug 版 | 修复版 | 6/25 历史基线 |
|---|---:|---:|---:|
| "No collisions" 计数 | 30 | **0** | — |
| Collisions（前三窗口） | 无 | 7943/16453/25611 | 正常 |
| final particles | 2,010,878 | **1,957,830** | 1,957,715 |
| Total energy | 1.5697e-3 | **1.9622e-3** | 1.9619e-3 |
| stuck | 0 | 0 | 0 |

末态物理量与 6/25 基线逐位吻合（差异为 RNG 轨迹噪声）。

## 流程改进

- `extract_metrics.py` 新增 `collisions` 门（rank0 碰撞计数 > 0），
  默认开启（`--no-expect-collisions` 关闭）——根除"零碰撞假阳性"
- 教训：**功能包跨版本同步时，GeometricField 构造器形态变更必须连同
  readOption 语义一起审**（读构造器↔MUST_READ / 值构造器↔READ_IF_PRESENT
  是配对的）；本次 int32 修复回归的 PASS 结论因门检查不含碰撞项而失效

## 超算同步清单（重建后重跑 mpi16_omp4_check）

```
src/lagrangian/dsmc/clouds/dsmcCloud.C      # int32 + sigmaTcRMax 两处修复
src/lagrangian/dsmc/clouds/dsmcCloudI.H     # int32 偏移
src/lagrangian/dsmc/replicatedMesh/dsmcReplicatedMesh.C  # int32 偏移
```
