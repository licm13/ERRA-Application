"""Ensemble Rainfall-Runoff Analysis (ERRA) Package.

ERRA 集成降雨-径流分析包

This package provides a comprehensive Python implementation of the Ensemble
Rainfall-Runoff Analysis (ERRA) framework, including advanced features for
nonlinear and non-stationary hydrologic analysis.

本包提供了 ERRA 框架的完整 Python 实现，包括非线性和非平稳水文分析的高级功能。

Main Functions / 主要函数:
--------------------------
- erra(): Core ERRA analysis function / 核心 ERRA 分析函数
- plot_erra_results(): Comprehensive plotting utilities / 综合绘图工具

Key Features / 主要特性:
-----------------------
- Linear rainfall-runoff analysis / 线性降雨-径流分析
- Nonlinear response functions (NRF) / 非线性响应函数
- Non-stationary analysis with data splitting / 带数据分割的非平稳分析
- Robust estimation (IRLS) / 鲁棒估计
- Broken-stick lag representation / 断棍时滞表示
- Tikhonov regularization / Tikhonov 正则化

References / 参考文献:
---------------------
Kirchner, J.W. (2024). Characterizing nonlinear, nonstationary, and
heterogeneous hydrologic behavior using Ensemble Rainfall-Runoff Analysis
(ERRA): proof of concept. Hydrology and Earth System Sciences, 28, 4427-4454.
https://doi.org/10.5194/hess-28-4427-2024

Kirchner, J.W. (2022). Impulse response functions for heterogeneous,
nonstationary, and nonlinear systems, estimated by deconvolution and
demixing of noisy time series. Sensors, 22(9), 3291.
https://doi.org/10.3390/s22093291
"""

from .erra_core import erra, ERRAResult
from .utils import plot_erra_results

__version__ = "1.1.0"
__author__ = "James Kirchner, Python Implementation Contributors"
__license__ = "GPL-3.0"

__all__ = [
    "erra",
    "ERRAResult",
    "plot_erra_results",
    "__version__",
]
