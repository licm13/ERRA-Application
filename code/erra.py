"""Python implementation of the Ensemble Rainfall-Runoff Analysis (ERRA).

Python 版本的 ERRA 实现，提供与原始 R 代码一致的核心思路：
通过加权岭回归求解降雨-径流系统的脉冲响应（Runoff Response
Distribution, RRD）。模块提供简洁的 API，方便与 NumPy/Pandas
生态结合。
"""
from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path
from typing import Iterable, Optional, Sequence, Tuple

import numpy as np
import pandas as pd


@dataclass
class ERRAResult:
    """Container for ERRA outputs / ERRA 结果容器。"""

    lags: np.ndarray
    """Lag times corresponding to the RRD coefficients (time units follow the
    input data spacing). / RRD 系数对应的时滞（单位与输入时间步一致）。"""

    rrd: pd.DataFrame
    """Estimated runoff response distributions for each precipitation driver.
    每个降水驱动变量的径流响应分布估计值。列名称与输入降水变量相同。"""

    stderr: pd.DataFrame
    """Standard errors of the RRD coefficients (same layout as ``rrd``).
    RRD 系数的标准误差，与 ``rrd`` 具有相同结构。"""

    fitted: np.ndarray
    """Fitted discharge series. / 拟合的流量序列。"""

    residuals: np.ndarray
    """Residuals of the regression. / 回归残差。"""

    weights: np.ndarray
    """Observation weights used in the regression. / 回归中使用的观测权重。"""

    design_shape: Tuple[int, int]
    """Shape of the design matrix. / 设计矩阵的形状。"""

    reg_strength: float
    """Tikhonov regularization strength. / Tikhonov 正则化强度。"""

    lag_knots: Optional[pd.DataFrame]
    """Optional broken-stick representation when ``nk`` > 0. / 当 ``nk``>0 时的
    分段线性（折线）表示。"""

    def to_dataframe(self) -> pd.DataFrame:
        """Return a tidy DataFrame with lag-indexed RRD values.

        返回一个以时滞为索引的整理数据表。"""

        data = {
            "lag": self.lags,
        }
        for col in self.rrd.columns:
            data[f"rrd_{col}"] = self.rrd[col].to_numpy()
            data[f"stderr_{col}"] = self.stderr[col].to_numpy()
        return pd.DataFrame(data)


def erra(
    p: Iterable[Sequence[float]] | pd.DataFrame | pd.Series | np.ndarray,
    q: Sequence[float] | pd.Series | np.ndarray,
    wt: Optional[Sequence[float]] = None,
    m: int = 60,
    nk: int = 0,
    nu: float = 0.0,
    fq: float = 0.0,
    dt: float = 1.0,
    agg: int = 1,
    labels: Optional[Sequence[str]] = None,
) -> ERRAResult:
    """Estimate runoff response distributions using ERRA methodology.

    利用 ERRA 方法估算径流响应分布。该函数实现了原始 R 脚本中的核心
    思路：构建设计矩阵、可选量化滤波、Tikhonov 正则化求解岭回归，并
    返回每个降水驱动变量的冲激响应。

    Parameters
    ----------
    p : array-like
        Precipitation series (vector or matrix). Missing values are allowed
        and will be removed together with the affected regression rows.
        降水时间序列，可以是向量或矩阵。缺失值会导致相关的回归行被删除。
    q : array-like
        Discharge series aligned with ``p``. / 与 ``p`` 对齐的流量序列。
    wt : array-like, optional
        Optional weights. Defaults to uniform weights.
        可选的观测权重，默认为全 1。
    m : int, optional
        Maximum lag (inclusive) for the RRD. / 最大时滞。
    nk : int, optional
        Number of knots for a broken-stick summarisation. If <=0 no
        summarisation is produced. / 折线化表示的节点数量，<=0 表示不计算。
    nu : float, optional
        Tikhonov regularisation strength. Values near zero correspond to an
        ordinary least squares fit. / Tikhonov 正则化参数。
    fq : float, optional
        Running quantile (0 <= fq < 1) removed from the discharge series to
        mitigate slow drifts. / 流量序列的滑动分位滤波参数。
    dt : float, optional
        Time step of the input series. / 输入序列时间步长。
    agg : int, optional
        Aggregate input time step by summing precipitation and averaging
        discharge over ``agg`` consecutive samples. / 时间序列聚合系数。
    labels : sequence of str, optional
        Column labels for precipitation inputs. If omitted, columns are
        auto-numbered. / 降水变量列名。
    """

    p_matrix, labels = _prepare_precipitation(p, labels)
    q_vec = _to_numpy(q)

    if agg > 1:
        p_matrix, q_vec, wt = _aggregate(p_matrix, q_vec, wt, agg)

    if wt is None:
        wt_vec = np.ones_like(q_vec, dtype=float)
    else:
        wt_vec = _to_numpy(wt)
        if len(wt_vec) != len(q_vec):
            raise ValueError("weights must have the same length as q")

    q_filtered = _apply_quantile_filter(q_vec, fq, window=4 * m + 1)

    design, response, weights = _build_design_matrix(
        p_matrix, q_filtered, wt_vec, m
    )

    beta, stderr, fitted, residuals = _solve_rrd(
        design, response, weights, nu, m, p_matrix.shape[1]
    )

    lags = np.arange(beta.shape[0]) * dt
    rrd_df = _to_rrd_dataframe(beta, labels)
    stderr_df = _to_rrd_dataframe(stderr, labels)

    knots = None
    if nk and nk > 0:
        knots = _broken_stick(beta.reshape(-1, p_matrix.shape[1]), nk, dt, labels)

    return ERRAResult(
        lags=lags,
        rrd=rrd_df,
        stderr=stderr_df,
        fitted=fitted,
        residuals=residuals,
        weights=weights,
        design_shape=design.shape,
        reg_strength=nu,
        lag_knots=knots,
    )


# ---------------------------------------------------------------------------
# Helper functions / 辅助函数
# ---------------------------------------------------------------------------


def _to_numpy(arr: Sequence[float] | pd.Series | np.ndarray) -> np.ndarray:
    """Convert input to a 1-D NumPy array. / 转换为一维 NumPy 数组。"""

    if isinstance(arr, np.ndarray):
        return np.asarray(arr, dtype=float)
    if isinstance(arr, pd.Series):
        return arr.to_numpy(dtype=float)
    return np.asarray(list(arr), dtype=float)


def _prepare_precipitation(
    p: Iterable[Sequence[float]] | pd.DataFrame | pd.Series | np.ndarray,
    labels: Optional[Sequence[str]],
) -> Tuple[np.ndarray, Sequence[str]]:
    """Ensure precipitation input is a 2-D NumPy array.

    保证降水输入为二维数组，并生成列名称。"""

    if isinstance(p, pd.DataFrame):
        matrix = p.to_numpy(dtype=float)
        cols = list(p.columns)
    elif isinstance(p, pd.Series):
        matrix = p.to_numpy(dtype=float)[:, None]
        cols = [p.name or "p1"]
    else:
        arr = np.asarray(p, dtype=float)
        if arr.ndim == 1:
            matrix = arr[:, None]
        elif arr.ndim == 2:
            matrix = arr
        else:
            raise ValueError("precipitation input must be 1-D or 2-D")
        cols = [f"p{i+1}" for i in range(matrix.shape[1])]

    if labels is not None:
        if len(labels) != matrix.shape[1]:
            raise ValueError("labels must match the number of precipitation columns")
        cols = list(labels)

    return matrix, cols


def _aggregate(
    p: np.ndarray,
    q: np.ndarray,
    wt: Optional[Sequence[float]],
    agg: int,
) -> Tuple[np.ndarray, np.ndarray, Optional[np.ndarray]]:
    """Aggregate time series by ``agg`` samples.

    将时间序列按照 ``agg`` 长度聚合。降水累加，流量和权重取平均。"""

    n = p.shape[0]
    trimmed = n - (n % agg)
    if trimmed == 0:
        raise ValueError("time series shorter than aggregation window")

    p_trim = p[:trimmed]
    q_trim = q[:trimmed]

    p_agg = p_trim.reshape(trimmed // agg, agg, p.shape[1]).sum(axis=1)
    q_agg = q_trim.reshape(trimmed // agg, agg).mean(axis=1)

    if wt is None:
        wt_agg = None
    else:
        wt_arr = _to_numpy(wt)[:trimmed]
        wt_agg = wt_arr.reshape(trimmed // agg, agg).mean(axis=1)

    return p_agg, q_agg, wt_agg


def _apply_quantile_filter(q: np.ndarray, fq: float, window: int) -> np.ndarray:
    """Apply a running quantile filter to discharge.

    对流量序列应用滑动分位滤波，以去除季节性或趋势。"""

    if fq <= 0.0:
        return q
    if fq >= 1.0:
        raise ValueError("fq must be in [0, 1)")

    series = pd.Series(q)
    window = max(3, min(window, len(series)))
    if fq == 0.5:
        trend = series.rolling(window, center=True, min_periods=1).median()
    else:
        trend = series.rolling(window, center=True, min_periods=1).quantile(fq)
    return (series - trend).to_numpy(dtype=float)


def _build_design_matrix(
    p: np.ndarray,
    q: np.ndarray,
    wt: np.ndarray,
    m: int,
) -> Tuple[np.ndarray, np.ndarray, np.ndarray]:
    """Construct the regression matrix with lags up to ``m``.

    根据最大时滞 ``m`` 构建设计矩阵，并删除包含缺测值的行。"""

    n, k = p.shape
    cols = k * (m + 1)
    design = np.zeros((n, cols), dtype=float)

    for lag in range(m + 1):
        shifted = np.roll(p, shift=lag, axis=0)
        shifted[:lag, :] = np.nan
        design[:, lag * k : (lag + 1) * k] = shifted

    valid = (~np.isnan(design).any(axis=1)) & (~np.isnan(q))
    design = design[valid]
    response = q[valid]
    weights = wt[valid]

    # Normalise weights to avoid numerical issues / 归一化权重
    mean_w = weights.mean()
    if not np.isfinite(mean_w) or mean_w <= 0:
        weights = np.ones_like(weights)
    else:
        weights = weights / mean_w

    return design, response, weights


def _tikhonov_matrix(size: int) -> np.ndarray:
    """Construct a 2nd-order difference operator for Tikhonov smoothing."""

    if size <= 2:
        return np.eye(size, dtype=float)
    L = np.zeros((size - 2, size), dtype=float)
    for i in range(size - 2):
        L[i, i : i + 3] = np.array([1.0, -2.0, 1.0])
    return L


def _solve_rrd(
    design: np.ndarray,
    response: np.ndarray,
    weights: np.ndarray,
    nu: float,
    m: int,
    k: int,
) -> Tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray]:
    """Solve the weighted ridge regression for RRD coefficients.

    求解加权岭回归。返回的系数以列主序排列：lag0 var1, lag0 var2, ..."""

    W_sqrt = np.sqrt(weights)[:, None]
    Xw = design * W_sqrt
    yw = response * W_sqrt[:, 0]

    beta_size = (m + 1) * k
    XtX = Xw.T @ Xw
    Xty = Xw.T @ yw

    if nu > 0.0:
        L = _tikhonov_matrix(m + 1)
        L_big = np.kron(np.eye(k), L)
        reg = nu * (L_big.T @ L_big)
        XtX = XtX + reg

    try:
        beta = np.linalg.solve(XtX, Xty)
    except np.linalg.LinAlgError as exc:
        raise np.linalg.LinAlgError(
            "Regression matrix is singular; consider increasing nu or reducing m"
        ) from exc

    fitted = design @ beta
    residuals = response - fitted

    dof = max(len(response) - beta_size, 1)
    sigma2 = float((weights * residuals**2).sum() / dof)
    try:
        inv_xtx = np.linalg.inv(XtX)
    except np.linalg.LinAlgError:
        inv_xtx = np.linalg.pinv(XtX)
    cov = inv_xtx * sigma2
    stderr = np.sqrt(np.maximum(np.diag(cov), 0.0))

    beta = beta.reshape(m + 1, k)
    stderr = stderr.reshape(m + 1, k)

    return beta, stderr, fitted, residuals


def _to_rrd_dataframe(arr: np.ndarray, labels: Sequence[str]) -> pd.DataFrame:
    """Convert coefficient array to DataFrame with bilingual headers."""

    data = {label: arr[:, i] for i, label in enumerate(labels)}
    df = pd.DataFrame(data)
    df.index.name = "lag"
    return df


def _broken_stick(
    beta: np.ndarray,
    nk: int,
    dt: float,
    labels: Sequence[str],
) -> pd.DataFrame:
    """Produce a broken-stick (piecewise-linear) summary of the RRD.

    生成折线化的 RRD 表示，便于与原始 R 代码的 ``nk`` 参数对应。"""

    m_plus_1, k = beta.shape
    if nk < 2:
        raise ValueError("nk must be >= 2 for broken-stick representation")

    lags = np.arange(m_plus_1) * dt
    max_lag = lags[-1]
    if nk > m_plus_1:
        nk = m_plus_1
    span_end = max(max_lag, dt)
    if nk <= 2:
        knot_lags = np.array([0.0, max_lag])
    else:
        inner = np.geomspace(dt, span_end, num=nk - 2, endpoint=True)
        knot_lags = np.concatenate(([0.0], inner, [max_lag]))
    knot_lags = np.unique(np.clip(knot_lags, 0.0, max_lag))

    rows = {}
    for j, label in enumerate(labels):
        rows[label] = np.interp(knot_lags, lags, beta[:, j])

    result = pd.DataFrame(rows, index=knot_lags)
    result.index.name = "lag"
    return result


# ---------------------------------------------------------------------------
# Plotting functions / 绘图函数
# ---------------------------------------------------------------------------


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
    
    Parameters
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
    _plot_rrd_with_error_bars(result, output_dir, filename_prefix, 
                             figsize, dpi, show_plots, save_plots, use_chinese)
    
    # 2. Fitted vs observed (if observed data provided)
    if observed_q is not None:
        _plot_fitted_vs_observed(result, observed_q, output_dir, filename_prefix,
                                figsize, dpi, show_plots, save_plots, use_chinese)
    
    # 3. Residuals analysis / 残差分析
    _plot_residuals_analysis(result, output_dir, filename_prefix,
                            figsize, dpi, show_plots, save_plots, use_chinese)
    
    # 4. Broken-stick representation (if available)
    if result.lag_knots is not None:
        _plot_broken_stick(result, output_dir, filename_prefix,
                          figsize, dpi, show_plots, save_plots, use_chinese)
    
    if save_plots:
        print(f"\n图片已保存到 / Figures saved to: {output_dir}")


def _configure_chinese_fonts():
    """Configure matplotlib to support Chinese fonts."""
    import matplotlib.pyplot as plt
    import matplotlib.font_manager as fm
    
    # Try to find available Chinese fonts
    chinese_fonts = ['SimHei', 'SimSun', 'Microsoft YaHei', 'WenQuanYi Zen Hei',
                     'Noto Sans CJK SC', 'Source Han Sans SC']
    
    available_fonts = [f.name for f in fm.fontManager.ttflist]
    chinese_font = None
    
    for font in chinese_fonts:
        if font in available_fonts:
            chinese_font = font
            break
    
    if chinese_font:
        plt.rcParams['font.sans-serif'] = [chinese_font] + plt.rcParams['font.sans-serif']
        plt.rcParams['axes.unicode_minus'] = False  # Fix negative sign display


def _plot_rrd_with_error_bars(result, output_dir, filename_prefix, 
                             figsize, dpi, show_plots, save_plots, use_chinese):
    """Plot RRD with error bars."""
    import matplotlib.pyplot as plt
    
    plt.figure(figsize=figsize)
    
    for col in result.rrd.columns:
        rrd_values = result.rrd[col].values
        stderr_values = result.stderr[col].values
        
        plt.errorbar(result.lags, rrd_values, yerr=stderr_values, 
                    label=f'RRD {col}', capsize=3, alpha=0.8)
    
    if use_chinese:
        plt.xlabel('Lag (time units) / 时滞 (时间单位)')
        plt.ylabel('RRD Coefficient / RRD系数')
        plt.title('Runoff Response Distribution with Error Bars\n径流响应分布 (带误差棒)')
    else:
        plt.xlabel('Lag (time units)')
        plt.ylabel('RRD Coefficient')
        plt.title('Runoff Response Distribution with Error Bars')
    
    plt.legend()
    plt.grid(True, alpha=0.3)
    plt.tight_layout()
    
    if save_plots:
        fig_path = output_dir / f"{filename_prefix}_rrd_with_errors.png"
        plt.savefig(fig_path, dpi=dpi, bbox_inches='tight')
        if use_chinese:
            print(f"保存RRD图 / Saved RRD plot: {fig_path}")
        else:
            print(f"Saved RRD plot: {fig_path}")
    
    if show_plots:
        plt.show()
    else:
        plt.close()


def _plot_fitted_vs_observed(result, observed_q, output_dir, filename_prefix,
                            figsize, dpi, show_plots, save_plots, use_chinese):
    """Plot fitted vs observed discharge."""
    import matplotlib.pyplot as plt
    
    if isinstance(observed_q, pd.Series):
        observed_q = observed_q.values
    
    plt.figure(figsize=(figsize[0], figsize[1] * 1.5))
    
    # Time series comparison (subset for visibility)
    plot_length = min(1000, len(result.fitted), len(observed_q))
    
    plt.subplot(2, 1, 1)
    time_idx = np.arange(plot_length)
    
    if use_chinese:
        fitted_label = 'Fitted / 拟合值'
        observed_label = 'Observed / 观测值'
        ylabel = 'Discharge / 流量'
        title1 = f'Fitted vs Observed Discharge (First {plot_length} points)\n拟合值与观测值对比 (前{plot_length}个点)'
    else:
        fitted_label = 'Fitted'
        observed_label = 'Observed'
        ylabel = 'Discharge'
        title1 = f'Fitted vs Observed Discharge (First {plot_length} points)'
    
    plt.plot(time_idx, result.fitted[:plot_length], 
             'b-', label=fitted_label, alpha=0.8, linewidth=1)
    plt.plot(time_idx, observed_q[:plot_length], 
             'r-', label=observed_label, alpha=0.6, linewidth=1)
    plt.ylabel(ylabel)
    plt.title(title1)
    plt.legend()
    plt.grid(True, alpha=0.3)
    
    # Scatter plot
    plt.subplot(2, 1, 2)
    min_len = min(len(result.fitted), len(observed_q))
    plt.scatter(observed_q[:min_len], result.fitted[:min_len], 
               alpha=0.5, s=1)
    
    # Add 1:1 line
    min_val = min(np.min(observed_q[:min_len]), np.min(result.fitted[:min_len]))
    max_val = max(np.max(observed_q[:min_len]), np.max(result.fitted[:min_len]))
    plt.plot([min_val, max_val], [min_val, max_val], 'k--', label='1:1 Line')
    
    # Calculate R²
    r2 = np.corrcoef(observed_q[:min_len], result.fitted[:min_len])[0, 1] ** 2
    plt.text(0.05, 0.95, f'R² = {r2:.3f}', transform=plt.gca().transAxes,
             bbox=dict(boxstyle='round', facecolor='white', alpha=0.8))
    
    if use_chinese:
        xlabel = 'Observed Discharge / 观测流量'
        ylabel = 'Fitted Discharge / 拟合流量'
        title2 = 'Scatter Plot: Fitted vs Observed\n散点图：拟合值与观测值'
    else:
        xlabel = 'Observed Discharge'
        ylabel = 'Fitted Discharge'
        title2 = 'Scatter Plot: Fitted vs Observed'
    
    plt.xlabel(xlabel)
    plt.ylabel(ylabel)
    plt.title(title2)
    plt.legend()
    plt.grid(True, alpha=0.3)
    
    plt.tight_layout()
    
    if save_plots:
        fig_path = output_dir / f"{filename_prefix}_fitted_vs_observed.png"
        plt.savefig(fig_path, dpi=dpi, bbox_inches='tight')
        if use_chinese:
            print(f"保存拟合对比图 / Saved fitted vs observed plot: {fig_path}")
        else:
            print(f"Saved fitted vs observed plot: {fig_path}")
    
    if show_plots:
        plt.show()
    else:
        plt.close()


def _plot_residuals_analysis(result, output_dir, filename_prefix,
                            figsize, dpi, show_plots, save_plots, use_chinese):
    """Plot residuals analysis."""
    import matplotlib.pyplot as plt
    from scipy import stats
    
    plt.figure(figsize=(figsize[0] * 1.2, figsize[1] * 1.3))
    
    # Time series of residuals
    plt.subplot(2, 2, 1)
    plt.plot(result.residuals, 'b-', alpha=0.7, linewidth=0.8)
    plt.axhline(y=0, color='r', linestyle='--', alpha=0.7)
    
    if use_chinese:
        plt.xlabel('Time Index / 时间索引')
        plt.ylabel('Residuals / 残差')
        plt.title('Residual Time Series\n残差时间序列')
    else:
        plt.xlabel('Time Index')
        plt.ylabel('Residuals')
        plt.title('Residual Time Series')
    plt.grid(True, alpha=0.3)
    
    # Histogram of residuals
    plt.subplot(2, 2, 2)
    plt.hist(result.residuals, bins=50, alpha=0.7, edgecolor='black', linewidth=0.5)
    
    if use_chinese:
        plt.xlabel('Residuals / 残差')
        plt.ylabel('Frequency / 频次')
        plt.title('Residual Distribution\n残差分布')
    else:
        plt.xlabel('Residuals')
        plt.ylabel('Frequency')
        plt.title('Residual Distribution')
    plt.grid(True, alpha=0.3)
    
    # Q-Q plot
    plt.subplot(2, 2, 3)
    stats.probplot(result.residuals, dist="norm", plot=plt)
    if use_chinese:
        plt.title('Q-Q Plot (Normal)\n正态Q-Q图')
    else:
        plt.title('Q-Q Plot (Normal)')
    plt.grid(True, alpha=0.3)
    
    # Residuals vs fitted
    plt.subplot(2, 2, 4)
    plt.scatter(result.fitted, result.residuals, alpha=0.5, s=1)
    plt.axhline(y=0, color='r', linestyle='--', alpha=0.7)
    
    if use_chinese:
        plt.xlabel('Fitted Values / 拟合值')
        plt.ylabel('Residuals / 残差')
        plt.title('Residuals vs Fitted\n残差与拟合值')
    else:
        plt.xlabel('Fitted Values')
        plt.ylabel('Residuals')
        plt.title('Residuals vs Fitted')
    plt.grid(True, alpha=0.3)
    
    plt.tight_layout()
    
    if save_plots:
        fig_path = output_dir / f"{filename_prefix}_residuals_analysis.png"
        plt.savefig(fig_path, dpi=dpi, bbox_inches='tight')
        if use_chinese:
            print(f"保存残差分析图 / Saved residuals analysis: {fig_path}")
        else:
            print(f"Saved residuals analysis: {fig_path}")
    
    if show_plots:
        plt.show()
    else:
        plt.close()


def _plot_broken_stick(result, output_dir, filename_prefix,
                      figsize, dpi, show_plots, save_plots, use_chinese):
    """Plot broken-stick representation."""
    import matplotlib.pyplot as plt
    
    plt.figure(figsize=figsize)
    
    # Plot original RRD
    for col in result.rrd.columns:
        plt.plot(result.lags, result.rrd[col].values, 
                '-', alpha=0.6, label=f'Original RRD {col}')
    
    # Plot broken-stick representation
    for col in result.lag_knots.columns:
        plt.plot(result.lag_knots.index, result.lag_knots[col].values,
                'o-', linewidth=2, markersize=4, 
                label=f'Broken-stick {col}')
    
    if use_chinese:
        plt.xlabel('Lag (time units) / 时滞 (时间单位)')
        plt.ylabel('RRD Coefficient / RRD系数')
        plt.title('Broken-stick Representation of RRD\nRRD的折线表示')
    else:
        plt.xlabel('Lag (time units)')
        plt.ylabel('RRD Coefficient')
        plt.title('Broken-stick Representation of RRD')
    
    plt.legend()
    plt.grid(True, alpha=0.3)
    plt.tight_layout()
    
    if save_plots:
        fig_path = output_dir / f"{filename_prefix}_broken_stick.png"
        plt.savefig(fig_path, dpi=dpi, bbox_inches='tight')
        if use_chinese:
            print(f"保存折线图 / Saved broken-stick plot: {fig_path}")
        else:
            print(f"Saved broken-stick plot: {fig_path}")
    
    if show_plots:
        plt.show()
    else:
        plt.close()
