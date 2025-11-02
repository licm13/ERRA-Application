# ERRA Theory Assets / ERRA 理论资源

本目录包含三份原始 PDF 文档（方法简介与两篇 SCI 论文）、原作者提供的
R 语言脚本以及本文档中新增的 Python 版本示例。通过这些资料，可以完整
复现 Ensemble Rainfall-Runoff Analysis (ERRA) 的理论背景与实践流程。

## Directory layout / 目录结构

```
ERRA-theory/
├── an-introduction-to-erra-v1.06.pdf
├── Kirchner_2022_Impulse Response Functions for Nonlinear, Nonstationary, and Heterogeneous.pdf
├── Kirchner_2024_Characterizing nonlinear, nonstationary, and heterogeneous hydrologic behavior.pdf
├── erra_scripts_v1.06/              # 原始 R 语言实现
└── python/                          # 新增 Python 版本
    ├── __init__.py
    ├── erra.py                      # 核心算法（含中英文注释）
    └── example.py                   # 使用示例（含中英文注释）
```

## Python port overview / Python 版本说明

- `erra.py` 采用 NumPy 与 Pandas 实现 ERRA 核心的加权岭回归求解流程，提
  供中英文注释和文档字符串，便于国内外开发者阅读。
- `example.py` 读取仓库随附的 MOPEX Saco River 小时数据，演示如何使用
  Python 版本计算 Runoff Response Distribution (RRD)。
- 通过 `__init__.py` 暴露 `erra` 函数与 `ERRAResult` 数据类，可在其他
  Python 项目中直接导入使用。

## Dependencies / 依赖

- Python >= 3.9
- NumPy
- Pandas

可以使用以下命令安装依赖：

```bash
pip install numpy pandas
```

## Usage / 使用方法

1. 确认当前工作目录为仓库根目录。
2. 使用模块运行示例脚本：

```bash
python ERRA-theory/python/example.py
```

运行后终端会输出设计矩阵规模、前几个滞后的 RRD 系数及其标准误差。

### Integrating into your workflow / 集成到自定义流程

```python
from pathlib import Path
import sys

import pandas as pd

sys.path.append(str(Path("ERRA-theory/python").resolve()))
from erra import erra

p = pd.DataFrame({"rain": rain_series})
q = flow_series

result = erra(p=p, q=q, m=72, agg=4, nu=0.05, fq=0.5)
print(result.rrd.head())
```

如需分段线性（broken-stick）表示，可传入 `nk` 参数：

```python
result = erra(p=p, q=q, m=72, nk=6, nu=0.1)
print(result.lag_knots)
```

## Notes / 注意事项

- Python 版本聚焦于核心回归求解流程，暂未完整覆盖 R 代码中的所有高级
  功能（例如多重过滤与鲁棒迭代求解）。若需使用这些功能，可参考原始
  R 脚本或在此基础上扩展。
- 输入数据需按固定时间步长排列，缺测值会自动剔除对应的回归行。

欢迎在 GitHub issue 中反馈问题或提交改进建议。
