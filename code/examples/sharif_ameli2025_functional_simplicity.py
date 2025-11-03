"""Synthetic ERRA study mirroring Sharif & Ameli (2025).

通过合成的湿润/干燥状态转移模型复刻 Sharif 和 Ameli（2025）的暴雨径流分析，突显
"功能简洁性"概念。脚本包含详细中英文注释并输出带中文字体的图件。
"""
from __future__ import annotations

from pathlib import Path
from typing import Tuple

import numpy as np
import pandas as pd

from erra import erra, plot_erra_results


def generate_state_series(n: int, seed: int = 7) -> np.ndarray:
    """Create a two-state (dry/wet) Markov chain.

    构建一个简单的马尔可夫链：0 表示干燥，1 表示湿润。湿润状态下暴雨更频繁，
    以模拟研究中的“功能简洁性”——虽然机理复杂，但主要响应路径清晰。
    """

    rng = np.random.default_rng(seed)
    states = np.zeros(n, dtype=int)
    state = 0
    for i in range(n):
        if state == 0:
            state = rng.choice([0, 1], p=[0.82, 0.18])
        else:
            state = rng.choice([1, 0], p=[0.70, 0.30])
        states[i] = state
    return states


def create_forcing_series(n: int, seed: int = 99) -> Tuple[pd.DataFrame, pd.Series]:
    """Generate precipitation proxies influenced by the wet/dry states.

    - ``burst_rain`` / 暴雨：湿润状态时概率增大
    - ``antecedent_moisture`` / 前期含水量：利用指数平滑积累
    - ``fractured_bedrock`` / 裂隙产流：湿润时增强、干燥时弱化
    """

    rng = np.random.default_rng(seed)
    states = generate_state_series(n, seed=seed + 10)

    burst = rng.gamma(shape=3.0, scale=4.0, size=n) * (0.6 + 0.8 * states)
    frontal = rng.gamma(shape=1.5, scale=2.0, size=n) * (0.9 + 0.3 * states)

    antecedent = np.zeros(n)
    for i in range(1, n):
        antecedent[i] = 0.85 * antecedent[i - 1] + 0.4 * frontal[i]

    fractured = 0.2 * frontal + 0.8 * burst * (0.4 + 0.6 * states)

    precip = pd.DataFrame(
        {
            "burst_rain": burst,
            "antecedent_moisture": antecedent,
            "fractured_bedrock": fractured,
        },
        index=pd.RangeIndex(n),
    )

    weights = pd.Series(1.0 + 0.1 * states, index=precip.index, name="weights")
    return precip, weights


def simulate_discharge(precip: pd.DataFrame, m: int) -> pd.Series:
    """Convert synthetic forcings to streamflow.

    设置快速（暴雨）与慢速（裂隙-基流）两类冲激响应核，实现“功能简洁性”。
    """

    n = len(precip)
    kernels = {
        "burst_rain": 1.0 * np.exp(-np.arange(m) / 1.8),
        "antecedent_moisture": 0.3 * np.exp(-np.arange(m) / 7.0),
        "fractured_bedrock": 0.45 * np.exp(-np.arange(m) / 12.0),
    }

    discharge = np.zeros(n)
    for col, kernel in kernels.items():
        convolved = np.convolve(precip[col].to_numpy(), kernel, mode="full")
        # Take only the first n elements to match the discharge array size
        discharge += convolved[:n]

    noise = 0.15 * np.random.default_rng(2024).standard_normal(n)
    discharge += noise
    return pd.Series(discharge, index=precip.index, name="discharge")


def main() -> None:
    """Run the Sharif & Ameli style ERRA demo / 运行 Sharif & Ameli 风格示例。"""

    m = 30
    precip, weights = create_forcing_series(n=400)
    discharge = simulate_discharge(precip, m=m)

    result = erra(
        p=precip,
        q=discharge,
        wt=weights,
        m=m,
        nu=8e-3,
        fq=0.3,
        dt=0.5,  # Represent 12-hour steps / 表示 12 小时时间步
        labels=["Burst", "Antecedent", "Fractured"],
    )

    print("设计矩阵规模 / design matrix shape:", result.design_shape)
    print("RRD 前几行 / leading RRD rows:")
    print(result.rrd.head())

    figures_dir = Path(__file__).resolve().parent / "figures"
    plot_erra_results(
        result=result,
        observed_q=discharge,
        output_dir=figures_dir,
        filename_prefix="sharif_ameli2025",
        show_plots=False,
        save_plots=True,
        figsize=(11, 7),
        dpi=300,
        use_chinese=True,
    )


if __name__ == "__main__":
    main()
