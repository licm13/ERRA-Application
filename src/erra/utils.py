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

try:
    import matplotlib.pyplot as plt
    from scipy import stats

    PLOTTING_AVAILABLE = True
except ModuleNotFoundError:
    # More specific than ImportError - only catches missing packages
    PLOTTING_AVAILABLE = False

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
    if not PLOTTING_AVAILABLE:
        raise ImportError(
            "matplotlib and scipy are required for plotting. "
            "Install with: pip install matplotlib scipy"
        )

    # Configure matplotlib for Chinese fonts if needed
    if use_chinese:
        _configure_chinese_fonts()

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
        result,
        output_dir,
        filename_prefix,
        figsize,
        dpi,
        show_plots,
        save_plots,
        use_chinese,
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


def _bilingual_text(english: str, chinese: str = "", use_chinese: bool = True) -> str:
    """Create bilingual text for plots.

    创建双语图表文本。

    Parameters / 参数
    ----------
    english : str
        English text / 英文文本
    chinese : str, optional
        Chinese text. If empty, returns English only / 中文文本
    use_chinese : bool
        Whether to include Chinese text / 是否包含中文

    Returns / 返回
    -------
    str
        Formatted text / 格式化文本
    """
    if use_chinese and chinese:
        return f"{english} / {chinese}"
    return english


def _save_and_show_plot(
    fig_path: Optional[Path],
    dpi: int,
    show_plot: bool,
    save_plot: bool,
    plot_name: str,
    use_chinese: bool,
) -> None:
    """Save and/or show a plot with consistent handling.

    保存和/或显示图表的统一处理。

    Parameters / 参数
    ----------
    fig_path : Path or None
        Path to save figure / 保存图片的路径
    dpi : int
        DPI for saved figure / 保存图片的DPI
    show_plot : bool
        Whether to show the plot / 是否显示图片
    save_plot : bool
        Whether to save the plot / 是否保存图片
    plot_name : str
        Name of the plot for messages / 图片名称用于消息
    use_chinese : bool
        Whether to use bilingual messages / 是否使用双语消息
    """
    if save_plot and fig_path is not None:
        plt.savefig(fig_path, dpi=dpi, bbox_inches="tight")
        if use_chinese:
            print(f"保存{plot_name} / Saved {plot_name}: {fig_path}")
        else:
            print(f"Saved {plot_name}: {fig_path}")

    if show_plot:
        plt.show()
    else:
        plt.close()


def _plot_rrd_with_error_bars(
    result,
    output_dir,
    filename_prefix,
    figsize,
    dpi,
    show_plots,
    save_plots,
    use_chinese,
):
    """Plot RRD with error bars.

    绘制带误差棒的RRD。
    """
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

    plt.xlabel(_bilingual_text("Lag (time units)", "时滞 (时间单位)", use_chinese))
    plt.ylabel(_bilingual_text("RRD Coefficient", "RRD系数", use_chinese))
    title = "Runoff Response Distribution with Error Bars"
    if use_chinese:
        title += "\n径流响应分布 (带误差棒)"
    plt.title(title)

    plt.legend()
    plt.grid(True, alpha=0.3)
    plt.tight_layout()

    fig_path = (
        output_dir / f"{filename_prefix}_rrd_with_errors.png" if save_plots else None
    )
    _save_and_show_plot(fig_path, dpi, show_plots, save_plots, "RRD plot", use_chinese)


def _plot_fitted_vs_observed(
    result,
    observed_q,
    output_dir,
    filename_prefix,
    figsize,
    dpi,
    show_plots,
    save_plots,
    use_chinese,
):
    """Plot fitted vs observed discharge.

    绘制拟合值与观测值对比图。
    """
    if isinstance(observed_q, pd.Series):
        observed_q = observed_q.values

    plt.figure(figsize=(figsize[0], figsize[1] * 1.5))

    # Time series comparison (subset for visibility)
    plot_length = min(1000, len(result.fitted), len(observed_q))

    plt.subplot(2, 1, 1)
    time_idx = np.arange(plot_length)

    fitted_label = _bilingual_text("Fitted", "拟合值", use_chinese)
    observed_label = _bilingual_text("Observed", "观测值", use_chinese)
    ylabel1 = _bilingual_text("Discharge", "流量", use_chinese)
    title1 = f"Fitted vs Observed Discharge (First {plot_length} points)"
    if use_chinese:
        title1 += f"\n拟合值与观测值对比 (前{plot_length}个点)"

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
    plt.ylabel(ylabel1)
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

    xlabel = _bilingual_text("Observed Discharge", "观测流量", use_chinese)
    ylabel2 = _bilingual_text("Fitted Discharge", "拟合流量", use_chinese)
    title2 = "Scatter Plot: Fitted vs Observed"
    if use_chinese:
        title2 += "\n散点图：拟合值与观测值"

    plt.xlabel(xlabel)
    plt.ylabel(ylabel2)
    plt.title(title2)
    plt.legend()
    plt.grid(True, alpha=0.3)

    plt.tight_layout()

    fig_path = (
        output_dir / f"{filename_prefix}_fitted_vs_observed.png" if save_plots else None
    )
    _save_and_show_plot(
        fig_path, dpi, show_plots, save_plots, "fitted vs observed plot", use_chinese
    )


def _plot_residuals_analysis(
    result,
    output_dir,
    filename_prefix,
    figsize,
    dpi,
    show_plots,
    save_plots,
    use_chinese,
):
    """Plot residuals analysis.

    绘制残差分析图。
    """
    plt.figure(figsize=(figsize[0] * 1.2, figsize[1] * 1.3))

    # Time series of residuals
    plt.subplot(2, 2, 1)
    plt.plot(result.residuals, "b-", alpha=0.7, linewidth=0.8)
    plt.axhline(y=0, color="r", linestyle="--", alpha=0.7)
    plt.xlabel(_bilingual_text("Time Index", "时间索引", use_chinese))
    plt.ylabel(_bilingual_text("Residuals", "残差", use_chinese))
    title1 = "Residual Time Series"
    if use_chinese:
        title1 += "\n残差时间序列"
    plt.title(title1)
    plt.grid(True, alpha=0.3)

    # Histogram of residuals
    plt.subplot(2, 2, 2)
    plt.hist(result.residuals, bins=50, alpha=0.7, edgecolor="black", linewidth=0.5)
    plt.xlabel(_bilingual_text("Residuals", "残差", use_chinese))
    plt.ylabel(_bilingual_text("Frequency", "频次", use_chinese))
    title2 = "Residual Distribution"
    if use_chinese:
        title2 += "\n残差分布"
    plt.title(title2)
    plt.grid(True, alpha=0.3)

    # Q-Q plot
    plt.subplot(2, 2, 3)
    stats.probplot(result.residuals, dist="norm", plot=plt)
    title3 = "Q-Q Plot (Normal)"
    if use_chinese:
        title3 += "\n正态Q-Q图"
    plt.title(title3)
    plt.grid(True, alpha=0.3)

    # Residuals vs fitted
    plt.subplot(2, 2, 4)
    plt.scatter(result.fitted, result.residuals, alpha=0.5, s=1)
    plt.axhline(y=0, color="r", linestyle="--", alpha=0.7)
    plt.xlabel(_bilingual_text("Fitted Values", "拟合值", use_chinese))
    plt.ylabel(_bilingual_text("Residuals", "残差", use_chinese))
    title4 = "Residuals vs Fitted"
    if use_chinese:
        title4 += "\n残差与拟合值"
    plt.title(title4)
    plt.grid(True, alpha=0.3)

    plt.tight_layout()

    fig_path = (
        output_dir / f"{filename_prefix}_residuals_analysis.png" if save_plots else None
    )
    _save_and_show_plot(
        fig_path, dpi, show_plots, save_plots, "residuals analysis", use_chinese
    )


def _plot_broken_stick(
    result,
    output_dir,
    filename_prefix,
    figsize,
    dpi,
    show_plots,
    save_plots,
    use_chinese,
):
    """Plot broken-stick representation.

    绘制断棍表示图。
    """
    plt.figure(figsize=figsize)

    # Plot original RRD
    for col in result.rrd.columns:
        plt.plot(
            result.lags,
            result.rrd[col].values,
            "-",
            alpha=0.6,
            label=f"Original RRD {col}",
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

    plt.xlabel(_bilingual_text("Lag (time units)", "时滞 (时间单位)", use_chinese))
    plt.ylabel(_bilingual_text("RRD Coefficient", "RRD系数", use_chinese))
    title = "Broken-stick Representation of RRD"
    if use_chinese:
        title += "\nRRD的折线表示"
    plt.title(title)

    plt.legend()
    plt.grid(True, alpha=0.3)
    plt.tight_layout()

    fig_path = (
        output_dir / f"{filename_prefix}_broken_stick.png" if save_plots else None
    )
    _save_and_show_plot(
        fig_path, dpi, show_plots, save_plots, "broken-stick plot", use_chinese
    )
