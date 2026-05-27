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

### 新增可配置参数 `nMinParcelsTvib`

**文件**:
- `src/lagrangian/dsmc/macroscopicProperties/derived/combined/dsmcVolFields/dsmcVolFields.H`
- `src/lagrangian/dsmc/macroscopicProperties/derived/combined/dsmcVolFields/dsmcVolFields.C`

**改动**:

1. **dsmcVolFields.H** — 新增成员变量：
```cpp
label nMinParcelsTvib_;
```

2. **dsmcVolFields.C** — 构造函数初始化列表：
```cpp
nMinParcelsTvib_(1),
```

3. **dsmcVolFields.C** — `createField()` 中从 `controlDict` 读取参数：
```cpp
nMinParcelsTvib_ =
    mesh_.time().controlDict().lookupOrDefault<label>("nMinParcelsTvib", 1);
```

4. **dsmcVolFields.C** — 内部 cell 的 Tvib 计算判据修改：
```cpp
// 原始代码
dsmcNSpeciesCum_[i][celli] > SMALL

// 修改后
dsmcNSpeciesCum_[i][celli] >= nMinParcelsTvib_
```

## 使用方法

在 `system/controlDict` 中添加（可选）：

```
nMinParcelsTvib 10;    // 可选，默认为 1
```

- 默认值 1：与原始行为基本一致（原为 > 1e-15，现为 >= 1）
- 建议值 10：消除大部分统计噪声，适用于多组分高超声速模拟
- 当某组分在某 cell 的累积 parcel 数低于阈值时，该组分不参与 Tvib 加权平均

## 备注

- 边界面（boundary faces）的 Tvib 计算未修改，因为边界使用通量密度 `speciesRhoNBF_` 而非 parcel 计数，且边界统计问题较内部 cell 轻微
- 阈值作用于累积量 `dsmcNSpeciesCum_`（跨所有 nAvTimeSteps 的总 parcel 数），而非单步 parcel 数
