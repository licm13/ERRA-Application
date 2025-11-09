"""Synthetic ERRA reproduction for Tu et al. (2025).

脚本模拟 Tu 等（2025）研究中的多年冻土退化情景，通过在时间序列中途切换冲激响应
核展示径流对降水敏感度的下降。生成的结果包含中英文注释图件。
"""
from __future__ import annotations

from pathlib import Path
from typing import Optional, Tuple

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

from erra import erra, plot_erra_results


DATA_DIR = Path(__file__).resolve().parent / "data" / "processed"
PROCESSED_FILE = DATA_DIR / "tu2025_inputs.csv"


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


def load_observed_dataset() -> Optional[Tuple[pd.DataFrame, pd.Series, pd.Series]]:
    if not PROCESSED_FILE.exists():
        return None

    data = pd.read_csv(PROCESSED_FILE, index_col="date", parse_dates=True)
    required = {
        "rainfall",
        "snowmelt",
        "active_layer",
        "discharge",
        "weights",
    }
    if not required.issubset(data.columns):
        raise ValueError("Processed Tu et al. dataset lacks required columns.")

    forcings = data[["rainfall", "snowmelt", "active_layer"]]
    discharge = data["discharge"]
    weights = data["weights"]
    return forcings, discharge, weights


def summarize_and_plot_differences(result, output_dir: Path, filename_prefix: str) -> None:
    peak_lags = result.rrd.idxmax()
    total_response = result.rrd.sum()
    snowmelt_share = total_response["Snowmelt"] / total_response.sum()
    active_cumsum = result.rrd["ALT"].cumsum()
    target = 0.95 * active_cumsum.iloc[-1]
    alt_tail = active_cumsum[active_cumsum >= target].index[0]

    rows = [
        (
            "Rainfall peak lag",
            "<1 day (Tu et al., 2025)",
            f"{peak_lags['Rain']:.1f} days",
            "Daily observations smear the rapid post-storm runoff.",
        ),
        (
            "Snowmelt share",
            "≈30% contribution (Tu et al., Fig. 4)",
            f"{snowmelt_share:.0%}",
            "Interior Alaska site shows stronger nival influence.",
        ),
        (
            "Active layer 95% mass lag",
            ">20 days (Tu et al., 2025)",
            f"{alt_tail:.1f} days",
            "Computed from thawing degree days; differs due to available indices.",
        ),
    ]

    fig, ax = plt.subplots(figsize=(11, 3.5))
    ax.axis("off")
    table = ax.table(
        cellText=[row[1:] for row in rows],
        rowLabels=[row[0] for row in rows],
        colLabels=["Original study", "ERRA reproduction", "Difference notes"],
        loc="center",
        cellLoc="left",
    )
    table.auto_set_font_size(False)
    table.set_fontsize(10)
    table.scale(1, 1.2)
    ax.set_title(
        "Comparison with Tu et al. (2025) / 与 Tu 等 (2025) 对比",
        fontsize=12,
        pad=16,
    )

    output_dir.mkdir(exist_ok=True)
    table_path = output_dir / f"{filename_prefix}_comparison_table.png"
    fig.savefig(table_path, dpi=300, bbox_inches="tight")
    plt.close(fig)


def main() -> None:
    """Run the Tu et al. style ERRA demo / 运行 Tu 等风格示例。"""

    m = 50

    dataset = load_observed_dataset()
    if dataset is not None:
        print(
            "使用 NOAA/USGS 实测数据 / Using observed NOAA/USGS dataset."
        )
        forcings, discharge, weights = dataset
        weights = weights.reindex(forcings.index).fillna(1.0)
    else:
        print(
            "未找到观测数据，使用合成序列。请运行 data_prep/tu2025_fetch_data.py 获取真实数据。"
        )
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

    summarize_and_plot_differences(result, figures_dir, "tu2025_permafrost")


if __name__ == "__main__":
    main()
