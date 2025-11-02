"""ERRA analysis example with English-only plots.

ERRA Python示例脚本：使用纯英文图表的版本，适用于没有中文字体的系统。
"""
from __future__ import annotations

from pathlib import Path
import sys

import pandas as pd

if __package__ is None or __package__ == "":
    # Allow running the script via ``python example_en.py``
    sys.path.append(str(Path(__file__).resolve().parent))
    from erra import erra, plot_erra_results  # type: ignore
else:
    from .erra import erra, plot_erra_results


DATA_DIR = Path(__file__).resolve().parents[1] / "demonstration-scripts" / "Source data"


def main() -> None:
    """Run ERRA analysis with English-only plots."""

    # Get script name for figure naming
    script_name = Path(__file__).stem  # Gets 'example_en' from 'example_en.py'
    figures_dir = Path(__file__).resolve().parent / "figures"

    dataset = DATA_DIR / "MOPEX SacoR hourly.txt"
    df = pd.read_csv(dataset, sep=r'\s+')  # Using newer pandas syntax

    result = erra(
        p=df[["p"]],
        q=df["q"],
        m=48,
        agg=6,
        fq=0.5,
        nu=1e-2,
        dt=6.0 / 24.0,  # 6 hours converted to days
        labels=["Saco_precip"],
    )

    print("Design matrix shape:", result.design_shape)
    print("First five lag coefficients:")
    print(result.rrd.head())
    print("Standard errors:")
    print(result.stderr.head())
    
    # Generate plots with English text only
    plot_erra_results(
        result=result,
        observed_q=df["q"],
        output_dir=figures_dir,
        filename_prefix=script_name,
        show_plots=True,
        save_plots=True,
        figsize=(10, 6),
        dpi=300,
        use_chinese=False  # Use English-only labels
    )


if __name__ == "__main__":
    main()