"""Utility functions for ERRA package.

工具函数模块

This module provides plotting utilities and helper functions for ERRA analysis.
本模块提供 ERRA 分析的绘图工具和辅助函数。
"""

from __future__ import annotations

from pathlib import Path
from typing import TYPE_CHECKING, Optional, Tuple

import numpy as np
import pandas as pd

if TYPE_CHECKING:
    from .erra_core import ERRAResult


def plot_erra_results(
    result: ERRAResult,
    observed_q: Optional[np.ndarray | pd.Series] = None,
    output_dir: Optional[str | Path] = None,
    filename_prefix: str = "erra",
    show_plots: bool = True,
    save_plots: bool = True,
    figsize: Tuple[int, int] = (10, 6),
    dpi: int = 300,
    use_chinese: bool = True,
) -> None:
    """Create comprehensive plots for ERRA results.

    为ERRA结果创建综合图表。

    Parameters / 参数
    ----------
    result : ERRAResult
        ERRA analysis results / ERRA分析结果
    observed_q : array-like, optional
        Original observed discharge for comparison / 原始观测流量用于对比
    output_dir : str or Path, optional
        Directory to save plots. If None, saves to './figures' / 保存图片的目录
    filename_prefix : str, optional
        Prefix for saved figure filenames / 保存图片的文件名前缀
    show_plots : bool, optional
        Whether to display plots / 是否显示图片
    save_plots : bool, optional
        Whether to save plots to files / 是否保存图片到文件
    figsize : tuple, optional
        Figure size for plots / 图片尺寸
    dpi : int, optional
        DPI for saved figures / 保存图片的DPI
    use_chinese : bool, optional
        Whether to include Chinese text in plots / 是否在图中包含中文文本
    """
    try:
        import matplotlib.pyplot as plt
        from scipy import stats

        # Configure matplotlib for Chinese fonts if needed
        if use_chinese:
            _configure_chinese_fonts()

    except ImportError as e:
        raise ImportError(
            "matplotlib and scipy are required for plotting. "
            "Install with: pip install matplotlib scipy"
        ) from e

    if save_plots:
        if output_dir is None:
            output_dir = Path.cwd() / "figures"
        else:
            output_dir = Path(output_dir)
        output_dir.mkdir(exist_ok=True)

    # 1. RRD with error bars / 带误差棒的RRD
    _plot_rrd_with_error_bars(
        result,
        output_dir,
        filename_prefix,
        figsize,
        dpi,
        show_plots,
        save_plots,
        use_chinese,
    )

    # 2. Fitted vs observed (if observed data provided)
    if observed_q is not None:
        _plot_fitted_vs_observed(
            result,
            observed_q,
            output_dir,
            filename_prefix,
            figsize,
            dpi,
            show_plots,
            save_plots,
            use_chinese,
        )

    # 3. Residuals analysis / 残差分析
    _plot_residuals_analysis(
        result, output_dir, filename_prefix, figsize, dpi, show_plots, save_plots, use_chinese
    )

    # 4. Broken-stick representation (if available)
    if result.lag_knots is not None:
        _plot_broken_stick(
            result,
            output_dir,
            filename_prefix,
            figsize,
            dpi,
            show_plots,
            save_plots,
            use_chinese,
        )

    if save_plots:
        print(f"\n图片已保存到 / Figures saved to: {output_dir}")


def _configure_chinese_fonts():
    """Configure matplotlib to support Chinese fonts.

    配置 matplotlib 以支持中文字体。
    """
    import matplotlib.pyplot as plt
    import matplotlib.font_manager as fm

    # Try to find available Chinese fonts
    chinese_fonts = [
        "SimHei",
        "SimSun",
        "Microsoft YaHei",
        "WenQuanYi Zen Hei",
        "Noto Sans CJK SC",
        "Source Han Sans SC",
    ]

    available_fonts = [f.name for f in fm.fontManager.ttflist]
    chinese_font = None

    for font in chinese_fonts:
        if font in available_fonts:
            chinese_font = font
            break

    if chinese_font:
        plt.rcParams["font.sans-serif"] = [chinese_font] + plt.rcParams[
            "font.sans-serif"
        ]
        plt.rcParams["axes.unicode_minus"] = False  # Fix negative sign display


def _plot_rrd_with_error_bars(
    result, output_dir, filename_prefix, figsize, dpi, show_plots, save_plots, use_chinese
):
    """Plot RRD with error bars.

    绘制带误差棒的RRD。
    """
    import matplotlib.pyplot as plt

    plt.figure(figsize=figsize)

    for col in result.rrd.columns:
        rrd_values = result.rrd[col].values
        stderr_values = result.stderr[col].values

        plt.errorbar(
            result.lags,
            rrd_values,
            yerr=stderr_values,
            label=f"RRD {col}",
            capsize=3,
            alpha=0.8,
        )

    if use_chinese:
        plt.xlabel("Lag (time units) / 时滞 (时间单位)")
        plt.ylabel("RRD Coefficient / RRD系数")
        plt.title("Runoff Response Distribution with Error Bars\n径流响应分布 (带误差棒)")
    else:
        plt.xlabel("Lag (time units)")
        plt.ylabel("RRD Coefficient")
        plt.title("Runoff Response Distribution with Error Bars")

    plt.legend()
    plt.grid(True, alpha=0.3)
    plt.tight_layout()

    if save_plots:
        fig_path = output_dir / f"{filename_prefix}_rrd_with_errors.png"
        plt.savefig(fig_path, dpi=dpi, bbox_inches="tight")
        if use_chinese:
            print(f"保存RRD图 / Saved RRD plot: {fig_path}")
        else:
            print(f"Saved RRD plot: {fig_path}")

    if show_plots:
        plt.show()
    else:
        plt.close()


def _plot_fitted_vs_observed(
    result, observed_q, output_dir, filename_prefix, figsize, dpi, show_plots, save_plots, use_chinese
):
    """Plot fitted vs observed discharge.

    绘制拟合值与观测值对比图。
    """
    import matplotlib.pyplot as plt

    if isinstance(observed_q, pd.Series):
        observed_q = observed_q.values

    plt.figure(figsize=(figsize[0], figsize[1] * 1.5))

    # Time series comparison (subset for visibility)
    plot_length = min(1000, len(result.fitted), len(observed_q))

    plt.subplot(2, 1, 1)
    time_idx = np.arange(plot_length)

    if use_chinese:
        fitted_label = "Fitted / 拟合值"
        observed_label = "Observed / 观测值"
        ylabel = "Discharge / 流量"
        title1 = f"Fitted vs Observed Discharge (First {plot_length} points)\n拟合值与观测值对比 (前{plot_length}个点)"
    else:
        fitted_label = "Fitted"
        observed_label = "Observed"
        ylabel = "Discharge"
        title1 = f"Fitted vs Observed Discharge (First {plot_length} points)"

    plt.plot(
        time_idx,
        result.fitted[:plot_length],
        "b-",
        label=fitted_label,
        alpha=0.8,
        linewidth=1,
    )
    plt.plot(
        time_idx,
        observed_q[:plot_length],
        "r-",
        label=observed_label,
        alpha=0.6,
        linewidth=1,
    )
    plt.ylabel(ylabel)
    plt.title(title1)
    plt.legend()
    plt.grid(True, alpha=0.3)

    # Scatter plot
    plt.subplot(2, 1, 2)
    min_len = min(len(result.fitted), len(observed_q))
    plt.scatter(observed_q[:min_len], result.fitted[:min_len], alpha=0.5, s=1)

    # Add 1:1 line
    min_val = min(np.min(observed_q[:min_len]), np.min(result.fitted[:min_len]))
    max_val = max(np.max(observed_q[:min_len]), np.max(result.fitted[:min_len]))
    plt.plot([min_val, max_val], [min_val, max_val], "k--", label="1:1 Line")

    # Calculate R²
    r2 = np.corrcoef(observed_q[:min_len], result.fitted[:min_len])[0, 1] ** 2
    plt.text(
        0.05,
        0.95,
        f"R² = {r2:.3f}",
        transform=plt.gca().transAxes,
        bbox=dict(boxstyle="round", facecolor="white", alpha=0.8),
    )

    if use_chinese:
        xlabel = "Observed Discharge / 观测流量"
        ylabel = "Fitted Discharge / 拟合流量"
        title2 = "Scatter Plot: Fitted vs Observed\n散点图：拟合值与观测值"
    else:
        xlabel = "Observed Discharge"
        ylabel = "Fitted Discharge"
        title2 = "Scatter Plot: Fitted vs Observed"

    plt.xlabel(xlabel)
    plt.ylabel(ylabel)
    plt.title(title2)
    plt.legend()
    plt.grid(True, alpha=0.3)

    plt.tight_layout()

    if save_plots:
        fig_path = output_dir / f"{filename_prefix}_fitted_vs_observed.png"
        plt.savefig(fig_path, dpi=dpi, bbox_inches="tight")
        if use_chinese:
            print(f"保存拟合对比图 / Saved fitted vs observed plot: {fig_path}")
        else:
            print(f"Saved fitted vs observed plot: {fig_path}")

    if show_plots:
        plt.show()
    else:
        plt.close()


def _plot_residuals_analysis(
    result, output_dir, filename_prefix, figsize, dpi, show_plots, save_plots, use_chinese
):
    """Plot residuals analysis.

    绘制残差分析图。
    """
    import matplotlib.pyplot as plt
    from scipy import stats

    plt.figure(figsize=(figsize[0] * 1.2, figsize[1] * 1.3))

    # Time series of residuals
    plt.subplot(2, 2, 1)
    plt.plot(result.residuals, "b-", alpha=0.7, linewidth=0.8)
    plt.axhline(y=0, color="r", linestyle="--", alpha=0.7)

    if use_chinese:
        plt.xlabel("Time Index / 时间索引")
        plt.ylabel("Residuals / 残差")
        plt.title("Residual Time Series\n残差时间序列")
    else:
        plt.xlabel("Time Index")
        plt.ylabel("Residuals")
        plt.title("Residual Time Series")
    plt.grid(True, alpha=0.3)

    # Histogram of residuals
    plt.subplot(2, 2, 2)
    plt.hist(result.residuals, bins=50, alpha=0.7, edgecolor="black", linewidth=0.5)

    if use_chinese:
        plt.xlabel("Residuals / 残差")
        plt.ylabel("Frequency / 频次")
        plt.title("Residual Distribution\n残差分布")
    else:
        plt.xlabel("Residuals")
        plt.ylabel("Frequency")
        plt.title("Residual Distribution")
    plt.grid(True, alpha=0.3)

    # Q-Q plot
    plt.subplot(2, 2, 3)
    stats.probplot(result.residuals, dist="norm", plot=plt)
    if use_chinese:
        plt.title("Q-Q Plot (Normal)\n正态Q-Q图")
    else:
        plt.title("Q-Q Plot (Normal)")
    plt.grid(True, alpha=0.3)

    # Residuals vs fitted
    plt.subplot(2, 2, 4)
    plt.scatter(result.fitted, result.residuals, alpha=0.5, s=1)
    plt.axhline(y=0, color="r", linestyle="--", alpha=0.7)

    if use_chinese:
        plt.xlabel("Fitted Values / 拟合值")
        plt.ylabel("Residuals / 残差")
        plt.title("Residuals vs Fitted\n残差与拟合值")
    else:
        plt.xlabel("Fitted Values")
        plt.ylabel("Residuals")
        plt.title("Residuals vs Fitted")
    plt.grid(True, alpha=0.3)

    plt.tight_layout()

    if save_plots:
        fig_path = output_dir / f"{filename_prefix}_residuals_analysis.png"
        plt.savefig(fig_path, dpi=dpi, bbox_inches="tight")
        if use_chinese:
            print(f"保存残差分析图 / Saved residuals analysis: {fig_path}")
        else:
            print(f"Saved residuals analysis: {fig_path}")

    if show_plots:
        plt.show()
    else:
        plt.close()


def _plot_broken_stick(
    result, output_dir, filename_prefix, figsize, dpi, show_plots, save_plots, use_chinese
):
    """Plot broken-stick representation.

    绘制断棍表示图。
    """
    import matplotlib.pyplot as plt

    plt.figure(figsize=figsize)

    # Plot original RRD
    for col in result.rrd.columns:
        plt.plot(
            result.lags, result.rrd[col].values, "-", alpha=0.6, label=f"Original RRD {col}"
        )

    # Plot broken-stick representation
    for col in result.lag_knots.columns:
        plt.plot(
            result.lag_knots.index,
            result.lag_knots[col].values,
            "o-",
            linewidth=2,
            markersize=4,
            label=f"Broken-stick {col}",
        )

    if use_chinese:
        plt.xlabel("Lag (time units) / 时滞 (时间单位)")
        plt.ylabel("RRD Coefficient / RRD系数")
        plt.title("Broken-stick Representation of RRD\nRRD的折线表示")
    else:
        plt.xlabel("Lag (time units)")
        plt.ylabel("RRD Coefficient")
        plt.title("Broken-stick Representation of RRD")

    plt.legend()
    plt.grid(True, alpha=0.3)
    plt.tight_layout()

    if save_plots:
        fig_path = output_dir / f"{filename_prefix}_broken_stick.png"
        plt.savefig(fig_path, dpi=dpi, bbox_inches="tight")
        if use_chinese:
            print(f"保存折线图 / Saved broken-stick plot: {fig_path}")
        else:
            print(f"Saved broken-stick plot: {fig_path}")

    if show_plots:
        plt.show()
    else:
        plt.close()
