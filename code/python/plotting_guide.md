# ERRA Python 绘图功能使用说明 / ERRA Python Plotting Guide

## 概述 / Overview

此Python模块为ERRA（Ensemble Rainfall-Runoff Analysis）分析提供了完整的绘图功能。图片会自动保存到当前脚本路径下的`figures`文件夹中，文件名以脚本名作为前缀。

This Python module provides comprehensive plotting functionality for ERRA (Ensemble Rainfall-Runoff Analysis). Figures are automatically saved to a `figures` folder in the current script's directory, with filenames prefixed by the script name.

## 功能特点 / Features

- 自动创建`figures`文件夹 / Automatically creates `figures` folder
- 基于脚本名的智能文件命名 / Intelligent filename based on script name  
- 支持中英文双语或纯英文标签 / Supports bilingual or English-only labels
- 高质量图片输出（300 DPI） / High-quality image output (300 DPI)
- 完整的残差分析 / Comprehensive residual analysis

## 生成的图片类型 / Generated Plot Types

1. **RRD系数图 / RRD Coefficients Plot** (`*_rrd_with_errors.png`)
   - 带误差棒的径流响应分布 / Runoff Response Distribution with error bars
   - 显示系数随时滞的变化 / Shows coefficient variation with lag

2. **拟合对比图 / Fitted vs Observed Plot** (`*_fitted_vs_observed.png`)
   - 时间序列对比 / Time series comparison
   - 散点图with 1:1线和R²值 / Scatter plot with 1:1 line and R² value

3. **残差分析图 / Residuals Analysis** (`*_residuals_analysis.png`)
   - 残差时间序列 / Residual time series
   - 残差分布直方图 / Residual distribution histogram
   - 正态Q-Q图 / Normal Q-Q plot
   - 残差vs拟合值 / Residuals vs fitted values

4. **折线表示图 / Broken-stick Plot** (`*_broken_stick.png`)
   - 仅当使用`nk > 0`参数时生成 / Only generated when `nk > 0` parameter is used

## 使用方法 / Usage

### 基本用法 / Basic Usage

```python
from erra import erra, plot_erra_results
import pandas as pd
from pathlib import Path

# 运行ERRA分析 / Run ERRA analysis
result = erra(
    p=precipitation_data,
    q=discharge_data,
    m=48,
    # ... 其他参数 / other parameters
)

# 生成图片 / Generate plots
script_name = Path(__file__).stem  # 获取脚本名 / Get script name
figures_dir = Path(__file__).resolve().parent / "figures"

plot_erra_results(
    result=result,
    observed_q=discharge_data,
    output_dir=figures_dir,
    filename_prefix=script_name,
    save_plots=True,
    show_plots=False,  # 设为False避免阻塞 / Set False to avoid blocking
    use_chinese=True   # True for bilingual, False for English only
)
```

### 参数说明 / Parameter Description

- `result`: ERRA分析结果 / ERRA analysis results
- `observed_q`: 原始观测流量数据（可选） / Original observed discharge (optional)
- `output_dir`: 图片保存目录（默认为./figures） / Output directory (default: ./figures)
- `filename_prefix`: 文件名前缀（建议使用脚本名） / Filename prefix (recommend script name)
- `save_plots`: 是否保存图片 / Whether to save plots
- `show_plots`: 是否显示图片 / Whether to show plots
- `use_chinese`: 是否使用中文标签 / Whether to use Chinese labels
- `figsize`: 图片尺寸 / Figure size
- `dpi`: 图片分辨率 / Image resolution

### 文件命名规则 / File Naming Convention

图片文件名格式：`{filename_prefix}_{plot_type}.png`

示例 / Examples:
- `example_rrd_with_errors.png`
- `analysis_fitted_vs_observed.png`
- `my_script_residuals_analysis.png`

## 示例脚本 / Example Scripts

1. **`example.py`** - 带中英文双语标签 / With bilingual labels
2. **`example_en.py`** - 纯英文标签 / English-only labels

## 中文字体支持 / Chinese Font Support

模块会自动尝试检测和配置中文字体：
- Windows: SimHei, SimSun, Microsoft YaHei
- Linux: WenQuanYi Zen Hei, Noto Sans CJK SC
- 如果没有找到中文字体，建议使用`use_chinese=False`

The module automatically tries to detect and configure Chinese fonts:
- Windows: SimHei, SimSun, Microsoft YaHei  
- Linux: WenQuanYi Zen Hei, Noto Sans CJK SC
- If no Chinese fonts are found, recommend using `use_chinese=False`

## 依赖包 / Dependencies

```bash
pip install matplotlib scipy pandas numpy
```

## 注意事项 / Notes

1. 图片会自动保存为PNG格式，300 DPI高分辨率 / Figures are automatically saved as PNG, 300 DPI high resolution
2. `figures`文件夹会自动创建 / `figures` folder is automatically created
3. 使用脚本名作为前缀可以避免不同脚本间的文件冲突 / Using script name as prefix avoids file conflicts between different scripts
4. 建议在批处理时设置`show_plots=False`以避免GUI阻塞 / Recommend setting `show_plots=False` for batch processing to avoid GUI blocking