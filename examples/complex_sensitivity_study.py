"""Advanced ERRA stress test with multilingual annotations.

综合性压力测试脚本：创建具有多尺度非平稳性的合成流域，检验 ERRA 在不同正则化强
度下的鲁棒性，并输出中英文注释图件。
"""
from __future__ import annotations

from pathlib import Path
from typing import Dict, Tuple

import numpy as np
import pandas as pd

from erra import erra, plot_erra_results


def build_multiscale_forcings(n: int = 720, seed: int = 2025) -> Tuple[pd.DataFrame, pd.Series]:
    """Create heterogeneous drivers with time-varying weights.

    构建 4 个驱动变量：对流降雨、层状降雨、雪融和蒸散抑制因子，并加入逐月变化
    的观测权重以模拟数据质量变化。
    """

    rng = np.random.default_rng(seed)
    t = np.arange(n)

    convective = rng.gamma(shape=2.2, scale=2.5, size=n)
    stratiform = rng.gamma(shape=1.5, scale=1.5, size=n) * (0.8 + 0.4 * np.sin(2 * np.pi * t / 90))
    snow = np.maximum(0.0, rng.normal(0.6 + 0.4 * np.cos(2 * np.pi * t / 180), 0.3, size=n))
    evap_suppression = 0.5 + 0.4 * np.sin(2 * np.pi * t / 365) + 0.1 * rng.standard_normal(n)

    forcings = pd.DataFrame(
        {
            "convective": convective,
            "stratiform": stratiform,
            "snowmelt": snow,
            "evap_suppression": evap_suppression,
        },
        index=pd.RangeIndex(n),
    )

    weights = pd.Series(1.0 + 0.3 * np.sin(2 * np.pi * t / 120), index=forcings.index, name="weights")
    return forcings, weights


def simulate_discharge(forcings: pd.DataFrame, m: int) -> pd.Series:
    """Generate discharge with non-stationary kernels and noise bursts.

    对每个驱动设置两段式冲激响应核，并在特定时间加入噪声峰值，以模拟极端事件。
    """

    n = len(forcings)
    split = int(n * 0.65)
    base_kernels: Dict[str, Tuple[float, float]] = {
        "convective": (1.1, 2.0),
        "stratiform": (0.6, 8.0),
        "snowmelt": (0.5, 15.0),
        "evap_suppression": (-0.4, 12.0),
    }

    discharge = np.zeros(n)
    data = forcings.to_numpy()
    cols = list(forcings.columns)
    rng = np.random.default_rng(4096)

    for i in range(n):
        stage_factor = 1.0 if i < split else 0.7
        for j, col in enumerate(cols):
            amplitude, scale = base_kernels[col]
            kernel = (amplitude * stage_factor) * np.exp(-np.arange(m) / (scale * (0.8 + 0.4 * (i / n))))
            for lag in range(min(m, i + 1)):
                discharge[i] += kernel[lag] * data[i - lag, j]

    # Add intermittent noise bursts / 加入间歇性噪声峰值
    shock_indices = rng.choice(np.arange(50, n - 50), size=5, replace=False)
    for idx in shock_indices:
        discharge[idx : idx + 10] += rng.normal(0, 0.5, size=min(10, n - idx))

    discharge += 0.2 * rng.standard_normal(n)
    return pd.Series(discharge, index=forcings.index, name="discharge")


def run_erra_variants(forcings: pd.DataFrame, discharge: pd.Series, weights: pd.Series, m: int) -> Dict[str, object]:
    """Evaluate ERRA under weak and strong regularisation.

    比较弱正则（近似最小二乘）与较强正则两种设定的差异，并返回结果字典。
    """

    configs = {
        "weak_regularisation": {"nu": 1e-3, "fq": 0.1},
        "strong_regularisation": {"nu": 4e-2, "fq": 0.2},
    }

    results = {}
    for name, params in configs.items():
        results[name] = erra(
            p=forcings,
            q=discharge,
            wt=weights,
            m=m,
            nu=params["nu"],
            fq=params["fq"],
            agg=2,
            dt=0.5,
            labels=["Conv", "Strat", "Snow", "Evap"],
        )
    return results


def main() -> None:
    """Run the complex stress-test / 运行综合压力测试。"""

    m = 40
    forcings, weights = build_multiscale_forcings()
    discharge = simulate_discharge(forcings, m=m)
    results = run_erra_variants(forcings, discharge, weights, m=m)

    weak = results["weak_regularisation"]
    strong = results["strong_regularisation"]

    print("弱正则 RRD 前三行 / Weak regularisation head:")
    print(weak.rrd.head(3))
    print("强正则 RRD 前三行 / Strong regularisation head:")
    print(strong.rrd.head(3))

    figures_dir = Path(__file__).resolve().parent / "figures"
    plot_erra_results(
        result=weak,
        observed_q=discharge,
        output_dir=figures_dir,
        filename_prefix="complex_sensitivity_weak",
        show_plots=False,
        save_plots=True,
        figsize=(12, 8),
        dpi=300,
        use_chinese=True,
    )

    plot_erra_results(
        result=strong,
        observed_q=discharge,
        output_dir=figures_dir,
        filename_prefix="complex_sensitivity_strong",
        show_plots=False,
        save_plots=True,
        figsize=(12, 8),
        dpi=300,
        use_chinese=True,
    )


if __name__ == "__main__":
    main()
