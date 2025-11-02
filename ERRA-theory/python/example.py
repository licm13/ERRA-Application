"""Quick demonstration of the Python ERRA port.

ERRA Python 示例脚本：读取 MOPEX Saco River 小时数据，聚合到 6 小时
尺度后估算径流响应分布，并打印前几个时滞的结果。
"""
from __future__ import annotations

from pathlib import Path
import sys

import pandas as pd

if __package__ is None or __package__ == "":
    # Allow running the script via ``python example.py``
    sys.path.append(str(Path(__file__).resolve().parent))
    from erra import erra  # type: ignore
else:
    from .erra import erra


DATA_DIR = Path(__file__).resolve().parents[1] / "demonstration-scripts" / "Source data"


def main() -> None:
    """Run a minimal ERRA analysis using bundled demonstration data.

    使用仓库自带的示例数据运行一次简易 ERRA 分析。"""

    dataset = DATA_DIR / "MOPEX SacoR hourly.txt"
    df = pd.read_csv(dataset, delim_whitespace=True)

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


if __name__ == "__main__":
    main()
