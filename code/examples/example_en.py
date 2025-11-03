"""ERRA analysis example with English-only plots.

ERRA Python示例脚本：使用纯英文图表的版本，适用于没有中文字体的系统。

Data location:
- By default this script looks for demo data under
    ``reference_materials/R_implementation/demonstration-scripts/Source data``
    at the repository root.
- You can override the location by setting the environment variable
    ``ERRA_DATA_DIR`` to the folder that contains files such as
    ``MOPEX SacoR hourly.txt``.
"""
from __future__ import annotations

from pathlib import Path
import os
import sys

import pandas as pd

if __package__ is None or __package__ == "":
    # Allow running the script via ``python example_en.py``
    sys.path.append(str(Path(__file__).resolve().parent))
    from erra import erra, plot_erra_results  # type: ignore
else:
    from .erra import erra, plot_erra_results


def _get_repo_root(start: Path) -> Path:
    """Return the repository root by looking for a marker file/dir.

    Falls back to the top-most parent if no marker is found.
    """
    cur = start
    for p in [cur, *cur.parents]:
        if (p / "pyproject.toml").exists() or (p / ".git").exists():
            return p
    return start.parents[len(start.parents) - 1]


def _resolve_data_dir() -> Path:
    """Resolve the directory containing demo source data.

    Priority:
    1) ERRA_DATA_DIR env var (if set and exists)
    2) reference_materials/R_implementation/demonstration-scripts/Source data
    3) code/demonstration-scripts/Source data (legacy/backup)
    4) Fallback: search for any '**/Source data' that contains the MOPEX files
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

    # Last resort: search
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