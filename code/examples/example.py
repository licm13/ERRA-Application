"""Quick demonstration of the Python ERRA port.

ERRA Python 示例脚本：读取 MOPEX Saco River 小时数据，聚合到 6 小时
尺度后估算径流响应分布，并打印前几个时滞的结果。

Data location:
- 默认从仓库根目录下的
    ``reference_materials/R_implementation/demonstration-scripts/Source data``
    查找演示数据。
- 也可以通过环境变量 ``ERRA_DATA_DIR`` 指定包含
    ``MOPEX SacoR hourly.txt`` 等文件的目录。
"""
from __future__ import annotations

from pathlib import Path
import sys
import os

import pandas as pd

if __package__ is None or __package__ == "":
    # Allow running the script via ``python example.py``
    sys.path.append(str(Path(__file__).resolve().parent))
    from erra import erra, plot_erra_results  # type: ignore
else:
    from .erra import erra, plot_erra_results


def _get_repo_root(start: Path) -> Path:
    """Best-effort detection of the repository root (via pyproject.toml or .git)."""
    cur = start
    for p in [cur, *cur.parents]:
        if (p / "pyproject.toml").exists() or (p / ".git").exists():
            return p
    return start.parents[len(start.parents) - 1]


def _resolve_data_dir() -> Path:
    """Resolve the demo data directory with env-var and fallbacks.

    Priority:
    1) ERRA_DATA_DIR env var
    2) reference_materials/R_implementation/demonstration-scripts/Source data
    3) code/demonstration-scripts/Source data (legacy)
    4) Search for any "**/Source data" containing known files
    """
    env_path = os.getenv("ERRA_DATA_DIR")
    if env_path:
        p = Path(env_path).expanduser().resolve()
        if p.exists():
            return p

    repo_root = _get_repo_root(Path(__file__).resolve())
    candidates = [
        repo_root
        / "reference_materials"
        / "R_implementation"
        / "demonstration-scripts"
        / "Source data",
        repo_root / "code" / "demonstration-scripts" / "Source data",
    ]
    for c in candidates:
        if c.exists():
            return c

    for d in repo_root.glob("**/Source data"):
        if d.is_dir() and any((d / fname).exists() for fname in [
            "MOPEX SacoR hourly.txt",
            "MOPEX NantahalaR hourly.txt",
            "Plynlimon hourly data for ERRA.txt",
        ]):
            return d

    raise FileNotFoundError(
        "Could not locate demo 'Source data' folder. Set ERRA_DATA_DIR or place files under reference_materials/R_implementation/demonstration-scripts/Source data."
    )


DATA_DIR = _resolve_data_dir()


def main() -> None:
    """Run a minimal ERRA analysis using bundled demonstration data.

    使用仓库自带的示例数据运行一次简易 ERRA 分析。"""

    # Get script name for figure naming
    script_name = Path(__file__).stem  # Gets 'example' from 'example.py'
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
        dt=6.0 / 24.0,  # 6 小时 -> 日比例 / 6 hours converted to days
        labels=["Saco_precip"],
    )

    print("设计矩阵规模 / design matrix shape:", result.design_shape)
    print("前五个滞后的 RRD 系数 / first five lag coefficients:")
    print(result.rrd.head())
    print("标准误差 / standard errors:")
    print(result.stderr.head())
    
    # Generate comprehensive plots using the new plotting function
    plot_erra_results(
        result=result,
        observed_q=df["q"],
        output_dir=figures_dir,
        filename_prefix=script_name,
        show_plots=True,
        save_plots=True,
        figsize=(10, 6),
        dpi=300
    )


if __name__ == "__main__":
    main()
