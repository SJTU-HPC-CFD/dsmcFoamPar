# Tvib 最小粒子数阈值修复工作日志（2026-05-28）

## 背景

在多组分 DSMC 模拟中，部分微量组分（如 NO、N、O）在某些 cell 中的累积 parcel 数极少（1-2个），导致该组分的振动温度 Tvib 出现非物理的异常大值。

### 问题根源

Tvib 的计算公式为：
```
iMean = ΣEvib / (kB * θv * N_parcels)
Tvib  = θv / ln(1 + 1/iMean)
```

当 N_parcels 很小时：
1. iMean 的统计波动极大（离散量子态导致跳变）
2. 单个 parcel 处于非零能级即可使 Tvib 飙升至数千甚至数万 K
3. 原有判据 `dsmcNSpeciesCum_[i][celli] > SMALL`（SMALL ≈ 1e-15）本质上只防除零，无统计过滤作用

### 影响范围

- Tvib_ 仅用于后处理输出（dsmcVolFields 写场文件）
- **不影响 QK 化学反应**：QK 反应直接使用粒子自身的振动能级 `p.vibLevel()`，不依赖宏观 Tvib 场
- 问题本质是输出/可视化问题，不污染物理模拟

## 修复内容

### 新增可配置参数 `nMinParcelsTvib` 和 `iMeanMinTvib`

**文件**:
- `src/lagrangian/dsmc/macroscopicProperties/derived/combined/dsmcVolFields/dsmcVolFields.H`
- `src/lagrangian/dsmc/macroscopicProperties/derived/combined/dsmcVolFields/dsmcVolFields.C`

**改动**:

1. **dsmcVolFields.H** — 新增成员变量：
```cpp
label nMinParcelsTvib_;
scalar iMeanMinTvib_;
```

2. **dsmcVolFields.C** — 构造函数初始化列表：
```cpp
nMinParcelsTvib_(1),
iMeanMinTvib_(0.01),
```

3. **dsmcVolFields.C** — `createField()` 中从 `controlDict` 读取参数：
```cpp
nMinParcelsTvib_ =
    mesh_.time().controlDict().lookupOrDefault<label>("nMinParcelsTvib", 1);

iMeanMinTvib_ =
    mesh_.time().controlDict().lookupOrDefault<scalar>("iMeanMinTvib", 0.01);
```

4. **dsmcVolFields.C** — 内部 cell 的 Tvib 计算判据修改：
```cpp
// 原始代码
dsmcNSpeciesCum_[i][celli] > SMALL
// 修改后
dsmcNSpeciesCum_[i][celli] >= nMinParcelsTvib_

// 原始代码
if (iMean > SMALL)
// 修改后
if (iMean > iMeanMinTvib_)
```

5. **dsmcVolFields.C** — 边界面的 Tvib 计算判据修改：
```cpp
// 原始代码
if (iMean > SMALL)
// 修改后
if (iMean > iMeanMinTvib_)
```

## 使用方法

在 `system/controlDict` 中添加（可选）：

```
nMinParcelsTvib 100;       // 最小累积 parcel 数阈值，默认为 1
iMeanMinTvib    0.01;      // 最小平均振动量子数阈值，默认为 0.01
```

### nMinParcelsTvib

- 默认值 1：与原始行为基本一致（原为 > 1e-15，现为 >= 1）
- 建议值 100：过滤微量组分在低密度区域的统计噪声
- 当某组分在某 cell 的累积 parcel 数低于阈值时，该组分不参与 Tvib 加权平均

### iMeanMinTvib

- 默认值 0.01：过滤振动模态未有效激发时的虚假 Tvib
- 物理含义：至少 1% 的分子被振动激发，才认为 Tvib 有意义
- 解决的问题：主要组分（如 N2）parcel 数充足但几乎全在基态时，少量激发态粒子（如激波回散）导致的虚假高 Tvib
- 对应的 Tvib 截断：N2(θv=3371K) < 730K，O2(θv=2256K) < 489K 时被过滤
- 同时作用于内部 cell 和边界面

两个参数互补：
- `nMinParcelsTvib`：解决"样本太少"（微量组分）
- `iMeanMinTvib`：解决"激发太弱"（主要组分在低温区）

## 备注

- 边界面（boundary faces）的 Tvib 计算：`iMeanMinTvib` 已应用，`nMinParcelsTvib` 未应用（边界使用通量密度而非 parcel 计数）
- 阈值作用于累积量 `dsmcNSpeciesCum_`（跨所有 nAvTimeSteps 的总 parcel 数），而非单步 parcel 数
- SPARTA 代码中无类似阈值处理，仅判断 count==0 和 ibar==0
