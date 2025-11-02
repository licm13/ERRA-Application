"""Synthetic reproduction of Gao et al. (2025) ERRA application.

合成复刻 Gao 等（2025）关于 ERRA 的应用研究，生成降水—地下水补给—径流的动态联
系示例。脚本通过随机种子确保可重复性，并生成中英文注释的图件。
"""
from __future__ import annotations

from pathlib import Path
from typing import Tuple

import numpy as np
import pandas as pd

# ---------------------------------------------------------------------------
# Import ERRA utilities / 导入 ERRA 工具
# ---------------------------------------------------------------------------
from erra import erra, plot_erra_results


def create_precipitation_series(n: int, seed: int = 42) -> Tuple[pd.DataFrame, pd.Series]:
    """Generate synthetic precipitation drivers.

    构建三个代表性降水驱动：
    - convective_bursts / 对流暴雨：短时高强度
    - stratiform_rain / 层状降雨：持续低强度
    - recharge_proxy / 补给指标：代表入渗—地下水补给

    Parameters
    ----------
    n : int
        Number of time steps / 时间步数量
    seed : int, optional
        Random seed for reproducibility / 随机种子
    """

    rng = np.random.default_rng(seed)
    time_index = pd.date_range("2000-01-01", periods=n, freq="D")

    convective = rng.gamma(shape=2.0, scale=5.0, size=n)
    stratiform = rng.gamma(shape=1.2, scale=3.0, size=n)

    recharge_base = 0.4 * stratiform + rng.normal(0.0, 0.3, size=n)
    recharge = np.clip(recharge_base, a_min=0.0, a_max=None)

    precip = pd.DataFrame(
        {
            "convective_bursts": convective,
            "stratiform_rain": stratiform,
            "recharge_proxy": recharge,
        },
        index=time_index,
    )

    weight = pd.Series(
        1.0 + 0.2 * np.sin(np.linspace(0, 8 * np.pi, n)),
        index=time_index,
        name="weights",
    )
    return precip, weight


def simulate_streamflow(precip: pd.DataFrame, weight: pd.Series, m: int) -> pd.Series:
    """Create a synthetic discharge record with groundwater influence.

    利用预设的冲激响应核将降水驱动卷积成流量序列，并加入权重变化代表观测
    质量差异。
    """

    rng = np.random.default_rng(123)
    n = len(precip)
    kernels = {
        "convective_bursts": 0.8 * np.exp(-np.arange(m) / 2),
        "stratiform_rain": 0.5 * np.exp(-np.arange(m) / 6),
        "recharge_proxy": 0.3 * np.exp(-np.arange(m) / 15),
    }

    discharge = np.zeros(n)
    for name, kernel in kernels.items():
        series = precip[name].to_numpy()
        convolved = np.convolve(series, kernel, mode="full")
        # Take only the first n elements to match the discharge array size
        discharge += convolved[:n]

    discharge += 0.1 * rng.standard_normal(n)
    discharge *= weight.to_numpy() / weight.mean()

    return pd.Series(discharge, index=precip.index, name="discharge")


def main() -> None:
    """Run the Gao et al. style ERRA demo / 运行 Gao 等风格示例。"""

    m = 45
    precip, weight = create_precipitation_series(n=540)
    discharge = simulate_streamflow(precip, weight, m=m)

    result = erra(
        p=precip,
        q=discharge,
        wt=weight,
        m=m,
        nu=5e-3,
        fq=0.2,
        dt=1.0,
        labels=["Convective", "Stratiform", "Recharge"],
    )

    print("设计矩阵规模 / design matrix shape:", result.design_shape)
    print("岭回归正则化系数 / ridge strength:", result.reg_strength)
    print("前五个滞后的 RRD 系数 / first five lagged RRD values:")
    print(result.rrd.head())

    figures_dir = Path(__file__).resolve().parent / "figures"
    plot_erra_results(
        result=result,
        observed_q=discharge,
        output_dir=figures_dir,
        filename_prefix="gao2025_dynamic",
        show_plots=False,
        save_plots=True,
        figsize=(11, 7),
        dpi=300,
        use_chinese=True,
    )


if __name__ == "__main__":
    main()
