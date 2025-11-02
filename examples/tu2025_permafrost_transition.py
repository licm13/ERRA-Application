"""Synthetic ERRA reproduction for Tu et al. (2025).

脚本模拟 Tu 等（2025）研究中的多年冻土退化情景，通过在时间序列中途切换冲激响应
核展示径流对降水敏感度的下降。生成的结果包含中英文注释图件。
"""
from __future__ import annotations

from pathlib import Path
from typing import Tuple

import numpy as np
import pandas as pd

from erra import erra, plot_erra_results


def create_climate_series(n: int = 600, seed: int = 314) -> Tuple[pd.DataFrame, pd.Series]:
    """Build synthetic precipitation and thaw-depth drivers.

    生成三类驱动：降雨、融雪与主动层厚度（Active Layer Thickness, ALT）。
    """

    rng = np.random.default_rng(seed)
    t = np.arange(n)
    seasonality = 0.6 + 0.4 * np.sin(2 * np.pi * t / 180)

    rainfall = rng.gamma(shape=2.5, scale=2.0, size=n) * seasonality
    snowmelt = np.clip(rng.normal(1.2 * seasonality, 0.5, size=n), a_min=0.0, a_max=None)
    alt = 0.5 + 0.4 * np.sin(2 * np.pi * (t - 45) / 365) + 0.05 * rng.standard_normal(n)
    alt = np.clip(alt, 0.1, None)

    forcings = pd.DataFrame(
        {
            "rainfall": rainfall,
            "snowmelt": snowmelt,
            "active_layer": alt,
        },
        index=pd.RangeIndex(n),
    )

    weights = pd.Series(1.0 + 0.2 * np.cos(2 * np.pi * t / 90), index=forcings.index, name="weights")
    return forcings, weights


def simulate_permafrost_response(forcings: pd.DataFrame, m: int) -> pd.Series:
    """Convert climate drivers to discharge with a mid-series transition.

    前半段（多年冻土完好）冲激响应核短而强；后半段（退化）核更加缓慢且幅度减小。
    """

    n = len(forcings)
    split = n // 2
    kernels_stage1 = {
        "rainfall": 0.9 * np.exp(-np.arange(m) / 2.5),
        "snowmelt": 0.7 * np.exp(-np.arange(m) / 4.0),
        "active_layer": 0.4 * np.exp(-np.arange(m) / 10.0),
    }
    kernels_stage2 = {
        "rainfall": 0.5 * np.exp(-np.arange(m) / 4.0),
        "snowmelt": 0.3 * np.exp(-np.arange(m) / 6.0),
        "active_layer": 0.2 * np.exp(-np.arange(m) / 14.0),
    }

    discharge = np.zeros(n)
    data = forcings.to_numpy()
    for i in range(n):
        stage_kernels = kernels_stage1 if i < split else kernels_stage2
        for driver_idx, column in enumerate(forcings.columns):
            kernel = stage_kernels[column]
            for lag in range(min(m, i + 1)):
                discharge[i] += kernel[lag] * data[i - lag, driver_idx]

    discharge += 0.12 * np.random.default_rng(777).standard_normal(n)
    return pd.Series(discharge, index=forcings.index, name="discharge")


def main() -> None:
    """Run the Tu et al. style ERRA demo / 运行 Tu 等风格示例。"""

    m = 50
    forcings, weights = create_climate_series()
    discharge = simulate_permafrost_response(forcings, m=m)

    result = erra(
        p=forcings,
        q=discharge,
        wt=weights,
        m=m,
        nu=1e-2,
        fq=0.25,
        dt=1.0,
        labels=["Rain", "Snowmelt", "ALT"],
    )

    print("设计矩阵规模 / design matrix shape:", result.design_shape)
    print("RRD 最前几行 / initial RRD rows:")
    print(result.rrd.head())

    figures_dir = Path(__file__).resolve().parent / "figures"
    plot_erra_results(
        result=result,
        observed_q=discharge,
        output_dir=figures_dir,
        filename_prefix="tu2025_permafrost",
        show_plots=False,
        save_plots=True,
        figsize=(11, 7),
        dpi=300,
        use_chinese=True,
    )


if __name__ == "__main__":
    main()
