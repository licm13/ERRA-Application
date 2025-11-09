# Ensemble Rainfall-Runoff Analysis (ERRA) Python Package

# ERRA（集成降雨-径流分析）Python 包

> **A comprehensive Python implementation of the ERRA framework for rainfall-runoff analysis, including advanced features for nonlinear and non-stationary hydrologic systems.**
>
> **用于降雨-径流分析的 ERRA 框架的综合 Python 实现，包括非线性和非平稳水文系统的高级功能。**

[![License: GPL v3](https://img.shields.io/badge/License-GPLv3-blue.svg)](https://www.gnu.org/licenses/gpl-3.0)
[![Python 3.9+](https://img.shields.io/badge/python-3.9+-blue.svg)](https://www.python.org/downloads/)

---

## Table of Contents / 目录

- [Overview](#overview--概述)
- [Key Features](#key-features--主要特性)
- [Installation](#installation--安装)
- [Quick Start](#quick-start--快速开始)
- [Advanced Features](#advanced-features--高级功能)
- [Core API](#core-api--核心-api)
- [Repository Structure](#repository-structure--仓库结构)
- [Examples](#examples--示例)
- [References](#references--参考文献)
- [License](#license--许可证)

---

## Overview / 概述

**ERRA (Ensemble Rainfall-Runoff Analysis)** is a powerful framework for analyzing the relationship between precipitation and streamflow using ensemble-based deconvolution methods. This Python package provides a complete implementation of ERRA, including advanced features originally developed in R by James Kirchner.

**ERRA（集成降雨-径流分析）**是一个强大的框架，使用基于集成的反卷积方法分析降水与河川流量之间的关系。这个 Python 包提供了 ERRA 的完整实现，包括 James Kirchner 最初在 R 中开发的高级功能。

### What ERRA Does / ERRA 的功能

ERRA estimates **Runoff Response Distributions (RRDs)** - the impulse response functions that describe how a catchment transforms precipitation into streamflow over time. Key capabilities include:

ERRA 估计**径流响应分布（RRDs）**——描述流域如何随时间将降水转化为河川流量的脉冲响应函数。主要功能包括：

- **Linear RRD estimation** via ridge regression with Tikhonov regularization
- **Nonlinear response analysis** to capture intensity-dependent behavior
- **Non-stationary analysis** to detect changes in catchment response over time or under different conditions
- **Robust estimation** to handle outliers and measurement errors
- **Multiple precipitation drivers** (rain vs. snow, convective vs. stratiform, etc.)

- **线性 RRD 估计**：通过带 Tikhonov 正则化的岭回归
- **非线性响应分析**：捕捉依赖于强度的行为
- **非平稳分析**：检测流域响应随时间或不同条件的变化
- **鲁棒估计**：处理异常值和测量误差
- **多个降水驱动变量**（雨 vs. 雪、对流 vs. 层状等）

---

## Key Features / 主要特性

### 🚀 Complete ERRA Implementation

- ✅ **Linear RRD estimation** with ridge regression and Tikhonov regularization
- ✅ **Nonlinear Response Functions (NRF)** via intensity-based splitting (xknots)
- ✅ **Non-stationary analysis** via covariate-based splitting (split_params)
- ✅ **Robust estimation** using Iteratively Reweighted Least Squares (IRLS)
- ✅ **Broken-stick lag representation** for computational efficiency
- ✅ **Multiple precipitation drivers** with automatic label management
- ✅ **Quantile filtering** to remove trends and seasonal patterns
- ✅ **Time aggregation** for different temporal resolutions

### 📦 Modern Python Package

- Standard `pip install` workflow
- Clean API with comprehensive docstrings (bilingual: English & Chinese)
- Type hints for better IDE support
- Minimal dependencies (numpy, pandas, scipy, matplotlib)

### 📊 Rich Visualization

- Automatic plotting of RRDs with error bars
- Fitted vs. observed discharge comparisons
- Comprehensive residual diagnostics
- Broken-stick representations
- Support for Chinese fonts in plots

### 📚 Extensive Documentation & Examples

- Detailed docstrings with parameter explanations and recommended values
- Five complete example scripts demonstrating different use cases
- Master demonstration showcasing all advanced features
- Bilingual documentation (English & Chinese)

---

## Installation / 安装

### Prerequisites / 前提条件

- Python ≥ 3.9
- pip (Python package installer)

### Install from source / 从源码安装

```powershell
# Clone the repository / 克隆仓库
git clone https://github.com/licm13/ERRA.git
cd ERRA

# Install in editable mode / 以可编辑模式安装
python -m pip install -U pip; pip install -e .

# Or install with development dependencies / 或安装开发依赖
python -m pip install -U pip; pip install -e .[dev]
```

### Dependencies / 依赖

The package automatically installs:
包会自动安装以下依赖：

- `numpy>=1.20.0` - Numerical computations / 数值计算
- `pandas>=1.3.0` - Data structures / 数据结构
- `scipy>=1.7.0` - Scientific computing / 科学计算
- `matplotlib>=3.4.0` - Plotting / 绘图

---

## Quick Start / 快速开始

### Basic Linear Analysis / 基本线性分析

```python
import numpy as np
from erra import erra, plot_erra_results

# Load your data (or use synthetic data for testing)
# 加载您的数据（或使用合成数据进行测试）
precipitation = np.random.exponential(5, 1000)  # mm/day
discharge = np.random.gamma(2, 3, 1000)  # m³/s

# Run ERRA analysis / 运行 ERRA 分析
result = erra(
    p=precipitation,
    q=discharge,
    m=60,  # Maximum lag: 60 days
    nu=0.1,  # Regularization strength
    dt=1.0,  # Time step: 1 day
)

# Access results / 访问结果
print(result.rrd)  # Runoff Response Distribution
print(result.stderr)  # Standard errors
print(f"R² = {np.corrcoef(discharge[60:], result.fitted)[0,1]**2:.3f}")

# Plot results / 绘制结果
plot_erra_results(
    result,
    observed_q=discharge,
    output_dir="./figures",
    save_plots=True,
)
```

### Running Examples / 运行示例

```powershell
# Run the master demonstration (showcases all features)
# 运行大师级演示（展示所有功能）
python examples\master_demonstration.py

# Run specific examples
# 运行特定示例
python examples\gao2025_dynamic_linkages.py
python examples\sharif_ameli2025_functional_simplicity.py
```

---

## Advanced Features / 高级功能

### 1. Nonlinear Response Analysis with `xknots`

Analyze how runoff response varies with precipitation intensity:
分析径流响应如何随降水强度变化：

```python
result = erra(
    p=precipitation,
    q=discharge,
    m=60,
    xknots=[50, 80, 95],  # Split at 50th, 80th, 95th percentiles
    xknot_type='percentiles',  # Interpret xknots as percentiles
    show_top_xknot=False,  # Don't show the unreliable top knot
)

# Access Nonlinear Response Functions / 访问非线性响应函数
print(result.nrf)  # NRF at each intensity segment
print(result.xknot_values)  # Actual intensity thresholds used
```

**Recommended xknot values / 推荐的 xknot 值:**
- `[50, 80, 95]` - Standard split (median, high, very high)
- `[33, 67]` - Tertile split (low, medium, high)
- `[25, 50, 75]` - Quartile split

### 2. Non-stationary Analysis with `split_params`

Analyze how response changes under different antecedent conditions:
分析响应如何在不同前期条件下变化：

```python
# Example: Split by antecedent wetness
# 示例：按前期湿度分割
antecedent_q = discharge.copy()  # Use lagged discharge as wetness proxy

split_params = {
    'crit': [antecedent_q],  # Criterion variable(s)
    'crit_label': ['Wetness'],  # Label(s)
    'crit_lag': [1],  # Lag by 1 timestep
    'pct_breakpts': [True],  # Use percentiles
    'breakpts': [[50, 90]],  # Split at 50th and 90th percentiles
    'thresh': [0],  # Ignore zeros when computing percentiles
    'by_bin': [True],  # Compute breakpoints within each bin
}

result = erra(
    p=precipitation,
    q=discharge,
    m=60,
    split_params=split_params,
)

# Access split results / 访问分割结果
print(result.split_labels)  # Labels for each subset
print(result.split_criteria)  # Criteria bounds for each subset
```

### 3. Robust Estimation

Handle outliers and measurement errors:
处理异常值和测量误差：

```python
result = erra(
    p=precipitation,
    q=discharge,
    m=60,
    robust=True,  # Enable IRLS
    robust_maxiter=10,  # Maximum iterations
    robust_tolerance=1e-4,  # Convergence tolerance
)
```

### 4. Broken-stick Lag Representation

Efficient representation for long lag times:
长时滞的高效表示：

```python
result = erra(
    p=precipitation,
    q=discharge,
    m=120,  # Long maximum lag
    nk=15,  # Use 15 knots for broken-stick representation
    nu=0.1,
)

# Access broken-stick representation / 访问断棍表示
print(result.lag_knots)  # RRD at knot lags
```

---

## Core API / 核心 API

### Main Function: `erra()`

```python
erra(
    p,  # Precipitation (vector or matrix)
    q,  # Discharge (vector)
    wt=None,  # Optional weights
    m=60,  # Maximum lag
    nk=0,  # Number of broken-stick knots (0 = none)
    nu=0.0,  # Tikhonov regularization (0-1)
    fq=0.0,  # Quantile filter (0-1, 0 = none)
    dt=1.0,  # Time step
    agg=1,  # Aggregation factor
    labels=None,  # Precipitation labels
    xknots=None,  # Nonlinear intensity knots
    xknot_type='percentiles',  # How to interpret xknots
    show_top_xknot=False,  # Show top knot in output
    split_params=None,  # Non-stationary splitting parameters
    robust=False,  # Use robust estimation (IRLS)
    robust_maxiter=10,  # IRLS max iterations
    robust_tolerance=1e-4,  # IRLS convergence tolerance
) -> ERRAResult
```

### Key Parameters & Recommended Values

#### `m` (Maximum Lag)
- **Hourly data**: 60-120 (2.5-5 days)
- **Daily data**: 30-60 (1-2 months)
- **每小时数据**: 60-120（2.5-5 天）
- **每日数据**: 30-60（1-2 个月）

#### `nu` (Regularization Strength)
- **Clean data**: 0 (no regularization)
- **Noisy data**: 0.01-0.1 (light smoothing)
- **Very noisy data**: 0.1-0.5 (moderate smoothing)
- **清洁数据**: 0（无正则化）
- **噪声数据**: 0.01-0.1（轻度平滑）
- **非常噪声的数据**: 0.1-0.5（中度平滑）

#### `fq` (Quantile Filter)
- **No filtering**: 0 (default)
- **Remove baseflow**: 0.1-0.3
- **Strong detrending**: 0.5 (median)
- **无滤波**: 0（默认）
- **移除基流**: 0.1-0.3
- **强去趋势**: 0.5（中位数）

#### `nk` (Broken-stick Knots)
- **Full resolution**: 0 (default)
- **Efficient representation**: 10-20
- **Rapid screening**: 5-10
- **完全分辨率**: 0（默认）
- **高效表示**: 10-20
- **快速筛选**: 5-10

---

## Repository Structure / 仓库结构

```
ERRA/
├── src/erra/                      # Main Python package / 主 Python 包
│   ├── __init__.py                # Package initialization / 包初始化
│   ├── erra_core.py               # Core ERRA implementation / 核心 ERRA 实现
│   ├── nonlin.py                  # Nonlinear response functions / 非线性响应函数
│   ├── splitting.py               # Data splitting for non-stationarity / 非平稳数据分割
│   └── utils.py                   # Plotting and utilities / 绘图和工具
│
├── examples/                      # Example scripts / 示例脚本
│   ├── master_demonstration.py    # Comprehensive demo of all features / 所有功能的综合演示
│   ├── gao2025_dynamic_linkages.py
│   ├── sharif_ameli2025_functional_simplicity.py
│   ├── tu2025_permafrost_transition.py
│   └── complex_sensitivity_study.py
│
├── reference_materials/           # Original R code and papers / 原始 R 代码和论文
│   ├── R_implementation/          # Original R scripts / 原始 R 脚本
│   ├── papers/                    # Application papers / 应用论文
│   └── theory_pdfs/               # Theory documentation / 理论文档
│
├── figures/                       # Generated figures / 生成的图表
├── pyproject.toml                 # Package configuration / 包配置
└── README.md                      # This file / 本文件
```

---

## Examples / 示例

The `examples/` directory contains five comprehensive demonstrations:
`examples/` 目录包含五个综合演示：

### 1. `master_demonstration.py` ⭐

**A comprehensive showcase of ALL advanced features:**
**所有高级功能的综合展示：**

- Multiple precipitation drivers (convective vs. stratiform)
- Nonlinear response analysis (xknots)
- Non-stationary analysis (wetness-based splitting)
- Broken-stick lag representation
- Robust estimation with outliers

**Run it:** `python examples/master_demonstration.py`

### 2. `gao2025_dynamic_linkages.py`

Demonstrates the dynamic linkages between different precipitation types (convective bursts, stratiform rain, groundwater recharge) and streamflow.

演示不同降水类型（对流暴雨、层状降雨、地下水补给）与河川流量之间的动态联系。

### 3. `sharif_ameli2025_functional_simplicity.py`

Illustrates the concept of "functional simplicity" - how complex spatial heterogeneity can manifest as simple dominant mechanisms at the catchment scale.

说明"功能简洁性"的概念——复杂的空间异质性如何在流域尺度上表现为简单的主导机制。

### 4. `tu2025_permafrost_transition.py`

Simulates permafrost degradation impacts on runoff generation, showing declining sensitivity of discharge to precipitation over time.

模拟多年冻土退化对径流生成的影响，显示流量对降水的敏感度随时间下降。

### Observed data preparation for SCI reproductions / SCI 案例观测数据准备

Each SCI-inspired example can operate on real observations. Run the corresponding downloader once before executing the example to fetch NOAA GHCN precipitation and USGS NWIS discharge records (outputs saved under `code/examples/data/processed/`).

每个 SCI 案例均支持真实观测数据。首次运行示例前，请执行以下脚本获取 NOAA GHCN 降水与 USGS NWIS 流量数据（结果存放在 `code/examples/data/processed/` 目录）。

```bash
python code/examples/data_prep/gao2025_fetch_data.py
python code/examples/data_prep/sharif_ameli2025_fetch_data.py
python code/examples/data_prep/tu2025_fetch_data.py
```

> 数据来源 / Data sources: NOAA National Centers for Environmental Information, Global Historical Climatology Network (Daily Summaries); USGS National Water Information System (Daily Discharge). See the [Open Data Sources](#open-data-sources--开放数据来源) section for full citations.

### 5. `complex_sensitivity_study.py`

A comprehensive stress test with multiple drivers, variable weights, and non-stationary noise to probe ERRA's robustness.

具有多个驱动变量、可变权重和非平稳噪声的综合压力测试，以探测 ERRA 的鲁棒性。

---

## References / 参考文献

### Primary References / 主要参考文献

**Kirchner, J.W.** (2024). Characterizing nonlinear, nonstationary, and heterogeneous hydrologic behavior using Ensemble Rainfall-Runoff Analysis (ERRA): proof of concept. *Hydrology and Earth System Sciences*, 28, 4427-4454.
📄 https://doi.org/10.5194/hess-28-4427-2024

**Kirchner, J.W.** (2022). Impulse response functions for heterogeneous, nonstationary, and nonlinear systems, estimated by deconvolution and demixing of noisy time series. *Sensors*, 22(9), 3291.
📄 https://doi.org/10.3390/s22093291

### Application Studies / 应用研究

See `reference_materials/papers/` for three recent studies applying ERRA to real-world catchments.

查看 `reference_materials/papers/` 获取将 ERRA 应用于实际流域的三个最新研究。

### Open Data Sources / 开放数据来源

- **USGS National Water Information System (NWIS)** — Daily discharge data for USGS 11477000, 09506000, and 15515500 (accessed via `https://waterservices.usgs.gov/nwis/dv/`). 美国地质调查局 NWIS 日尺度流量数据。
- **NOAA National Centers for Environmental Information (NCEI)** — Global Historical Climatology Network Daily Summaries for stations USC00043110, USW00023160, and USW00026411 (accessed via `https://www.ncei.noaa.gov/access/services/data/v1`). 美国国家海洋与大气管理局（NOAA）NCEI GHCN 日尺度降水、温度与积雪数据。

---

## License / 许可证

This project is licensed under the GNU General Public License v3.0 (GPL-3.0).

本项目采用 GNU 通用公共许可证 v3.0 (GPL-3.0) 授权。

```
Copyright (C) 2025 ETH Zurich and James Kirchner
Copyright (C) 2025 Python Implementation Contributors
```

See the [LICENSE](LICENSE) file or https://www.gnu.org/licenses/gpl-3.0.en.html for details.

---

## Citation / 引用

If you use this package in your research, please cite:

如果您在研究中使用此包，请引用：

```bibtex
@article{kirchner2024erra,
  title={Characterizing nonlinear, nonstationary, and heterogeneous hydrologic behavior using Ensemble Rainfall-Runoff Analysis (ERRA): proof of concept},
  author={Kirchner, James W.},
  journal={Hydrology and Earth System Sciences},
  volume={28},
  pages={4427--4454},
  year={2024},
  doi={10.5194/hess-28-4427-2024}
}

@article{kirchner2022impulse,
  title={Impulse response functions for heterogeneous, nonstationary, and nonlinear systems, estimated by deconvolution and demixing of noisy time series},
  author={Kirchner, James W.},
  journal={Sensors},
  volume={22},
  number={9},
  pages={3291},
  year={2022},
  doi={10.3390/s22093291}
}
```

---

## Feedback / 反馈

Questions or suggestions can be submitted through the [GitHub issue tracker](https://github.com/licm13/ERRA/issues).

问题或建议可以通过 [GitHub issue 跟踪器](https://github.com/licm13/ERRA/issues)提交。

---

## Acknowledgments / 致谢

This Python implementation builds upon the original R code by **James Kirchner** (ETH Zurich). We thank him for making his theoretical framework and code publicly available.

此 Python 实现基于 **James Kirchner**（苏黎世联邦理工学院）的原始 R 代码。我们感谢他公开其理论框架和代码。

---

**Happy analyzing! / 祝分析愉快！** 🌊📊

---

## Examples and Demos / 示例与演示

- Advanced examples live in `examples/` and import the installed package `erra`.
    进阶示例位于 `examples/`，通过已安装的 `erra` 包运行。

- Quick demos live in `code/python-version-example/` and ship with a local `erra.py` for
    self-contained runs. Prefer the package under `src/erra` for development work.
    快速演示位于 `code/python-version-example/`，包含本地 `erra.py` 便于即开即用；开发时更推荐使用 `src/erra` 包。

Data location / 数据路径：

- Default demo data are under `reference_materials/R_implementation/demonstration-scripts/Source data/`.
    默认演示数据位于 `reference_materials/R_implementation/demonstration-scripts/Source data/`。
- You can override via environment variable `ERRA_DATA_DIR`.
    可通过环境变量 `ERRA_DATA_DIR` 自定义数据目录。

PowerShell example / PowerShell 示例：

```powershell
$env:ERRA_DATA_DIR = "${PWD}\reference_materials\R_implementation\demonstration-scripts\Source data"
python code\python-version-example\example.py
```

---

## Plotting Guide / 绘图指南

This package provides comprehensive plots via `plot_erra_results()`.
本包通过 `plot_erra_results()` 提供完整绘图功能。

Generated plots / 生成的图片：

1. RRD with error bars (`*_rrd_with_errors.png`)  带误差棒的 RRD 图
2. Fitted vs observed (`*_fitted_vs_observed.png`)  拟合对比图（含 R² 与 1:1 线）
3. Residuals analysis (`*_residuals_analysis.png`)  残差分析（时间序列、直方图、Q-Q、残差vs拟合）
4. Broken-stick (`*_broken_stick.png`)  折线（断棍）表示（当 `nk>0`）

Basic usage / 基本用法：

```python
from erra import erra, plot_erra_results
import pandas as pd
from pathlib import Path

result = erra(p=df[["p"]], q=df["q"], m=48)
script_name = Path(__file__).stem
figures_dir = Path(__file__).resolve().parent / "figures"

plot_erra_results(
        result=result,
        observed_q=df["q"],
        output_dir=figures_dir,
        filename_prefix=script_name,
        save_plots=True,
        show_plots=False,
        use_chinese=True,   # True: bilingual / False: English only
)
```

Notes / 注意：

- Figures are saved at 300 DPI; the `figures` folder is created automatically.
    图片默认 300 DPI；`figures` 目录自动创建。
- If Chinese fonts are missing, set `use_chinese=False` to avoid glyph warnings.
    系统缺少中文字体时，建议设置 `use_chinese=False` 以避免字体告警。

---

## Theory Assets / 理论资源

Key theory PDFs and the original R implementation are included under `reference_materials/`.
理论文档与原始 R 实现收录在 `reference_materials/` 目录。

```text
reference_materials/
├── theory_pdfs/                      # Intro and core theory PDFs / 理论与方法综述
├── papers/                           # Application papers / 应用论文
└── R_implementation/
        ├── erra_scripts_v1.06/           # Original R scripts / 原始 R 脚本
        └── demonstration-scripts/
                └── Source data/              # MOPEX 等演示数据
```

For a quick Python demo, see `code/python-version-example/`.
快速 Python 演示请参见 `code/python-version-example/`。

---

## Scripts overview / 脚本总览

Where to find and what each script does — with concise bilingual notes.
脚本位置与作用概览——附简要中英文说明。

- code/examples/master_demonstration.py
    - EN: Master demo covering all advanced features: multiple drivers, nonlinear (xknots), non-stationary splitting, broken-stick (nk), robust IRLS; synthetic data; saves figures under the same folder.
    - 中文：大师级综合演示，涵盖多驱动、非线性（xknots）、非平稳分割、断棍（nk）、鲁棒IRLS；使用合成数据；图片保存在同级目录。

- code/examples/gao2025_dynamic_linkages.py
    - EN: Reproduces Gao (2025)-style dynamic linkages among convective, stratiform, recharge proxy; uses weights; typical m≈45; outputs `gao2025_dynamic_*` figures.
    - 中文：复刻 Gao（2025）风格的对流/层状/补给三类驱动及其动态联系；包含观测权重；m≈45；输出 `gao2025_dynamic_*` 图件。

- code/examples/sharif_ameli2025_functional_simplicity.py
    - EN: Demonstrates “functional simplicity” with a wet/dry Markov chain driving forcings; dt=0.5 (12-hour steps); shows contrasting fast/slow responses.
    - 中文：通过干/湿两态马尔可夫链驱动示例，展示“功能简洁性”；dt=0.5（12小时步长）；体现快/慢响应差异。

- code/examples/tu2025_permafrost_transition.py
    - EN: Permafrost transition demo with mid-series kernel change (degrading sensitivity); illustrates non-stationarity in time.
    - 中文：多年冻土过渡示例，中途切换冲激响应核（敏感度下降）；体现时间上的非平稳性。

- code/examples/complex_sensitivity_study.py
    - EN: Stress-test comparing weak vs strong regularization (nu) and different fq; returns result variants and saves comparative plots.
    - 中文：综合压力测试，对比弱/强正则与不同 fq；返回多组结果并保存对比图。

- code/examples/example.py
    - EN: Minimal linear analysis using MOPEX SacoR dataset; bilingual plots; auto-resolves data folder or honors `ERRA_DATA_DIR`.
    - 中文：最小线性示例（MOPEX SacoR 数据）；中英文图；自动定位数据目录或使用 `ERRA_DATA_DIR`。

- code/examples/example_en.py
    - EN: Same as above with English-only plots; useful when Chinese fonts are unavailable.
    - 中文：纯英文版，系统无中文字体时更稳妥。

- code/python-version-example/erra.py
    - EN: Self-contained demo implementation for quick runs in the `code/python-version-example/` folder. Its `ERRAResult` differs from `src/erra/erra_core.ERRAResult`, so linters may warn; prefer `src/erra` for development.
    - 中文：`code/python-version-example/` 目录的自包含演示实现；其 `ERRAResult` 与包内类型不同，静态检查可能提示；开发时优先使用 `src/erra`。

- Package modules / 包内模块
    - `src/erra/erra_core.py`: EN: Algorithmic core with linear/nonlinear/splitting/robust/nk; returns `ERRAResult`. 中文：算法核心，包含线性/非线性/分割/鲁棒/nk；返回 `ERRAResult`。
    - `src/erra/utils.py`: EN: Plotting utilities (`plot_erra_results`): RRD+errors, fitted vs observed, residuals, broken-stick. 中文：绘图工具。
    - `src/erra/nonlin.py`: EN: Nonlinear helpers (x′ transform, NRF construction and labels). 中文：非线性辅助（x′ 变换、NRF 组装与标签）。
    - `src/erra/splitting.py`: EN: Covariate-based splitting utilities and validators. 中文：按协变量分割工具与校验。

Run tips / 运行提示：

- Use installed package for scripts in `code/examples/`:
    使用已安装包运行 `code/examples/`：

    ```powershell
    python -m pip install -U pip; pip install -e .
    python code\examples\master_demonstration.py
    ```

- Ensure data path is available (for MOPEX demos):
    确保演示数据路径可用：

    ```powershell
    $env:ERRA_DATA_DIR = "${PWD}\reference_materials\R_implementation\demonstration-scripts\Source data"
    ```
