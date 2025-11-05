"""Core ERRA implementation with advanced features.

核心 ERRA 实现，包含高级功能

This module provides the main erra() function and supporting classes for
rainfall-runoff analysis, including support for:
- Nonlinear response analysis (xknots)
- Non-stationary analysis (split_params)
- Robust estimation (IRLS)
- Broken-stick lag representation (nk)
- Tikhonov regularization (nu)

本模块提供主要的 erra() 函数和支持类，用于降雨-径流分析，包括：
- 非线性响应分析 (xknots)
- 非平稳分析 (split_params)
- 鲁棒估计 (IRLS)
- 断棍时滞表示 (nk)
- Tikhonov 正则化 (nu)
"""

from __future__ import annotations

from dataclasses import dataclass
from typing import Dict, Iterable, Literal, Optional, Sequence, Tuple

import numpy as np
import pandas as pd

from .nonlin import betaprime_to_nrf, create_nrf_labels, create_xprime_matrix
from .splitting import make_split_sets, validate_split_params


@dataclass
class ERRAResult:
    """Container for ERRA outputs / ERRA 结果容器。

    Attributes / 属性
    -----------------
    lags : np.ndarray
        Lag times corresponding to RRD coefficients (time units follow input data)
        RRD 系数对应的时滞（单位与输入时间步一致）
    rrd : pd.DataFrame
        Estimated runoff response distributions for each precipitation driver
        每个降水驱动变量的径流响应分布估计值
    stderr : pd.DataFrame
        Standard errors of the RRD coefficients (same layout as rrd)
        RRD 系数的标准误差，与 rrd 具有相同结构
    fitted : np.ndarray
        Fitted discharge series / 拟合的流量序列
    residuals : np.ndarray
        Residuals of the regression / 回归残差
    weights : np.ndarray
        Observation weights used in regression / 回归中使用的观测权重
    design_shape : tuple
        Shape of the design matrix / 设计矩阵的形状
    reg_strength : float
        Tikhonov regularization strength / Tikhonov 正则化强度
    lag_knots : pd.DataFrame, optional
        Broken-stick representation when nk > 0 / 当 nk>0 时的折线表示
    nrf : pd.DataFrame, optional
        Nonlinear Response Functions at each xknot (when xknots is used)
        每个 xknot 的非线性响应函数（使用 xknots 时）
    nrf_stderr : pd.DataFrame, optional
        Standard errors of NRF coefficients
        NRF 系数的标准误差
    xknot_values : np.ndarray, optional
        Actual xknot values used in analysis
        分析中使用的实际 xknot 值
    split_labels : list of str, optional
        Labels for split subsets (when split_params is used)
        分割子集的标签（使用 split_params 时）
    split_criteria : pd.DataFrame, optional
        Documentation of splitting criteria
        分割标准的文档
    """

    lags: np.ndarray
    rrd: pd.DataFrame
    stderr: pd.DataFrame
    fitted: np.ndarray
    residuals: np.ndarray
    weights: np.ndarray
    design_shape: Tuple[int, int]
    reg_strength: float
    lag_knots: Optional[pd.DataFrame] = None
    nrf: Optional[pd.DataFrame] = None
    nrf_stderr: Optional[pd.DataFrame] = None
    xknot_values: Optional[np.ndarray] = None
    split_labels: Optional[list] = None
    split_criteria: Optional[pd.DataFrame] = None

    def to_dataframe(self) -> pd.DataFrame:
        """Return a tidy DataFrame with lag-indexed RRD values.

        返回一个以时滞为索引的整理数据表。
        """
        data = {"lag": self.lags}
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
    xknots: Optional[np.ndarray | list] = None,
    xknot_type: Literal[
        "values", "percentiles", "cumsum", "sqsum", "even"
    ] = "percentiles",
    show_top_xknot: bool = False,
    split_params: Optional[Dict] = None,
    robust: bool = False,
    robust_maxiter: int = 10,
    robust_tolerance: float = 1e-4,
) -> ERRAResult:
    """Estimate runoff response distributions using ERRA methodology.

    利用 ERRA 方法估算径流响应分布。

    This function implements the complete ERRA framework including:
    - Linear RRD estimation via ridge regression
    - Nonlinear response analysis via intensity-based splitting (xknots)
    - Non-stationary analysis via covariate-based splitting (split_params)
    - Robust estimation via Iteratively Reweighted Least Squares (IRLS)
    - Broken-stick lag representation for computational efficiency

    此函数实现完整的 ERRA 框架，包括：
    - 通过岭回归进行线性 RRD 估计
    - 通过基于强度的分割进行非线性响应分析 (xknots)
    - 通过基于协变量的分割进行非平稳分析 (split_params)
    - 通过迭代重加权最小二乘法 (IRLS) 进行鲁棒估计
    - 用于计算效率的断棍时滞表示

    Parameters / 参数
    ----------
    p : array-like
        Precipitation series (vector or matrix). Missing values are allowed
        and will be removed together with the affected regression rows.
        降水时间序列，可以是向量或矩阵。缺失值会导致相关的回归行被删除。

    q : array-like
        Discharge series aligned with p. / 与 p 对齐的流量序列。

    wt : array-like, optional
        Optional weights. Defaults to uniform weights.
        可选的观测权重，默认为全 1。

    m : int, default=60
        Maximum lag (inclusive) for the RRD. Typical values: 30-120 depending
        on time resolution and catchment response time.
        最大时滞。典型值：30-120，取决于时间分辨率和流域响应时间。

        Recommended values / 推荐值:
        - Hourly data: m=60-120 (captures 2.5-5 day response)
        - Daily data: m=30-60 (captures 1-2 month response)

    nk : int, default=0
        Number of knots for a broken-stick summarization of the RRD over lag times.
        If <=0, no summarization is produced. If >0, nk knots are placed at lags
        0, 1, and a geometric progression between 1 and m. This allows efficient
        representation of RRDs that decay exponentially with lag.
        RRD 在时滞上的断棍总结的节点数量。
        如果 <=0，则不进行总结。如果 >0，则在时滞 0、1 和 1 到 m 之间的几何级数处放置 nk 个节点。
        这允许有效表示随时滞呈指数衰减的 RRD。

        Recommended values / 推荐值:
        - nk=0: Full resolution (default) / 完全分辨率（默认）
        - nk=10-20: Good compromise for most applications / 大多数应用的良好折衷
        - nk=5-10: For rapid initial screening / 用于快速初步筛选

    nu : float, default=0.0
        Tikhonov regularization strength (0 ≤ nu < 1). Higher values impose
        stronger smoothness constraints on the RRD. Use when data is noisy
        or when the design matrix is near-singular.
        Tikhonov 正则化参数 (0 ≤ nu < 1)。较高的值对 RRD 施加更强的平滑约束。
        当数据噪声较大或设计矩阵接近奇异时使用。

        Recommended values / 推荐值:
        - nu=0: No regularization (default for clean data) / 无正则化（清洁数据的默认值）
        - nu=0.01-0.1: Light smoothing for noisy data / 噪声数据的轻度平滑
        - nu=0.1-0.5: Moderate smoothing for very noisy data / 非常噪声数据的中度平滑

    fq : float, default=0.0
        Running quantile (0 ≤ fq < 1) removed from discharge series to
        mitigate slow drifts or seasonal patterns. A running quantile filter
        with window of 4m+1 timesteps is applied.
        从流量序列中移除的滑动分位数 (0 ≤ fq < 1)，以减轻缓慢漂移或季节模式。
        应用窗口为 4m+1 时间步的滑动分位数滤波器。

        Recommended values / 推荐值:
        - fq=0: No filtering (default) / 无滤波（默认）
        - fq=0.1-0.3: Remove low-flow baseline / 移除低流量基线
        - fq=0.5: Remove median (strong detrending) / 移除中位数（强去趋势）

    dt : float, default=1.0
        Time step of the input series in desired units (hours, days, etc.).
        输入序列的时间步长，以所需单位（小时、天等）表示。

    agg : int, default=1
        Aggregate input time step by summing precipitation and averaging
        discharge over agg consecutive samples. Useful for reducing computational
        cost or for addressing temporal autocorrelation.
        通过对连续 agg 个样本求和降水和平均流量来聚合输入时间步。
        用于减少计算成本或解决时间自相关。

    labels : sequence of str, optional
        Column labels for precipitation inputs. If omitted, columns are
        auto-numbered as "p1", "p2", etc.
        降水变量列名。如果省略，列将自动编号为"p1"、"p2"等。

    xknots : array-like, optional
        Knot values for nonlinear response analysis. These define intensity
        thresholds that split precipitation into segments, allowing estimation
        of how response varies with precipitation intensity.
        用于非线性响应分析的节点值。这些定义将降水分割成段的强度阈值，
        允许估计响应如何随降水强度变化。

        For xknot_type="percentiles" (default), typical values:
        对于 xknot_type="percentiles"（默认），典型值：
        - xknots=[50, 80, 95]: Split at median, 80th, and 95th percentiles
        - xknots=[33, 67]: Tertile split (low, medium, high)
        - xknots=[20, 40, 60, 80]: Quintile split

        If None, linear analysis is performed (single RRD per driver).
        如果为 None，则执行线性分析（每个驱动一个 RRD）。

    xknot_type : str, default="percentiles"
        How to interpret xknots:
        如何解释 xknots：
        - "percentiles": Percentiles of precipitation distribution
        - "cumsum": Percentiles of cumulative sum (splits total precipitation)
        - "sqsum": Percentiles of cumulative sum of squares (splits leverage)
        - "values": Direct precipitation values
        - "even": Evenly spaced knots

    show_top_xknot : bool, default=False
        Whether to include the highest xknot (at max precipitation) in output.
        Generally set to False as estimates at the highest knot tend to be
        unreliable due to sparse data.
        是否在输出中包含最高的 xknot（在最大降水处）。
        通常设置为 False，因为最高节点的估计由于数据稀疏往往不可靠。

    split_params : dict, optional
        Parameters for splitting precipitation based on external covariates
        (e.g., antecedent wetness, temperature) to analyze non-stationary behavior.
        基于外部协变量（例如前期湿度、温度）分割降水以分析非平稳行为的参数。

        Dictionary must contain / 字典必须包含:
        - crit: list of criterion arrays / 标准数组列表
        - crit_label: list of criterion labels / 标准标签列表
        - crit_lag: list of lag values / 滞后值列表
        - pct_breakpts: list of booleans / 布尔值列表
        - breakpts: list of breakpoint lists / 断点列表列表
        - thresh: list of threshold values / 阈值列表
        - by_bin: list of booleans / 布尔值列表

        Example / 示例:
        split_params = {
            'crit': [antecedent_q],
            'crit_label': ['Wetness'],
            'crit_lag': [1],
            'pct_breakpts': [True],
            'breakpts': [[50, 90]],  # Split at 50th and 90th percentiles
            'thresh': [0],
            'by_bin': [True]
        }

    robust : bool, default=False
        Use robust estimation (Iteratively Reweighted Least Squares, IRLS)
        to reduce influence of outliers. Useful when discharge contains
        measurement errors or extreme events.
        使用鲁棒估计（迭代重加权最小二乘法，IRLS）以减少异常值的影响。
        当流量包含测量误差或极端事件时很有用。

    robust_maxiter : int, default=10
        Maximum iterations for IRLS algorithm.
        IRLS 算法的最大迭代次数。

    robust_tolerance : float, default=1e-4
        Convergence tolerance for IRLS (relative change in coefficients).
        IRLS 的收敛容差（系数的相对变化）。

    Returns / 返回
    -------
    ERRAResult
        Object containing all analysis results including RRD, standard errors,
        fitted values, residuals, and (if applicable) NRF and split results.
        包含所有分析结果的对象，包括 RRD、标准误差、拟合值、残差，
        以及（如适用）NRF 和分割结果。

    Notes / 注意事项
    -----
    **Tikhonov Regularization (nu parameter) / Tikhonov 正则化 (nu 参数)**

    The nu parameter controls the strength of Tikhonov-Phillips regularization,
    which imposes a smoothness penalty on the second derivative of the RRD.
    This is implemented via a second-order difference operator L such that
    the regularized objective is:

    nu 参数控制 Tikhonov-Phillips 正则化的强度，它对 RRD 的二阶导数施加平滑惩罚。
    这通过二阶差分算子 L 实现，使得正则化目标为：

        minimize: ||y - Xβ||² + nu × ||Lβ||²

    where L is the second-order difference matrix. This prevents overfitting
    and produces smoother RRDs, especially useful when m is large relative
    to the effective sample size.

    其中 L 是二阶差分矩阵。这可以防止过拟合并产生更平滑的 RRD，
    特别是当 m 相对于有效样本大小较大时很有用。

    **ARMA Noise Correction / ARMA 噪声校正**

    The R implementation includes automatic selection of autoregressive (AR)
    order to correct for serially correlated noise. The current Python
    implementation does not include this feature, but it can be partially
    addressed by:
    - Increasing agg to aggregate timesteps
    - Using fq to remove slow-varying trends
    - Interpreting residual autocorrelation cautiously

    R 实现包括自动选择自回归 (AR) 阶数以校正序列相关噪声。
    当前的 Python 实现不包括此功能，但可以通过以下方式部分解决：
    - 增加 agg 以聚合时间步
    - 使用 fq 移除缓慢变化的趋势
    - 谨慎解释残差自相关

    References / 参考文献
    ----------
    Kirchner, J.W. (2024). Characterizing nonlinear, nonstationary, and
    heterogeneous hydrologic behavior using Ensemble Rainfall-Runoff Analysis
    (ERRA): proof of concept. Hydrology and Earth System Sciences, 28,
    4427-4454. https://doi.org/10.5194/hess-28-4427-2024

    Examples / 示例
    --------
    **Basic linear analysis / 基本线性分析:**

    >>> result = erra(p=precipitation, q=discharge, m=60, nu=0.1)

    **Nonlinear analysis / 非线性分析:**

    >>> result = erra(
    ...     p=precipitation,
    ...     q=discharge,
    ...     m=60,
    ...     xknots=[50, 80, 95],
    ...     xknot_type='percentiles'
    ... )

    **Non-stationary analysis / 非平稳分析:**

    >>> split_params = {
    ...     'crit': [antecedent_discharge],
    ...     'crit_label': ['Wetness'],
    ...     'crit_lag': [1],
    ...     'pct_breakpts': [True],
    ...     'breakpts': [[50, 90]],
    ...     'thresh': [0],
    ...     'by_bin': [True]
    ... }
    >>> result = erra(
    ...     p=precipitation,
    ...     q=discharge,
    ...     m=60,
    ...     split_params=split_params
    ... )
    """
    # Prepare inputs
    p_matrix, col_labels = _prepare_precipitation(p, labels)
    q_vec = _convert_to_numpy_array(q)

    # Apply aggregation if requested
    if agg > 1:
        p_matrix, q_vec, wt = _aggregate(p_matrix, q_vec, wt, agg)

    # Handle weights
    if wt is None:
        wt_vec = np.ones_like(q_vec, dtype=float)
    else:
        wt_vec = _convert_to_numpy_array(wt)
        if len(wt_vec) != len(q_vec):
            raise ValueError("weights must have the same length as q")

    # Apply splitting if requested
    split_labels = None
    split_criteria = None
    if split_params is not None:
        validate_split_params(split_params)
        p_matrix, q_vec, wt_vec, split_labels, split_criteria = make_split_sets(
            p_matrix, q_vec, split_params, wt_vec
        )
        # Update labels to reflect splitting
        col_labels = split_labels

    # Apply nonlinear transformation if requested
    xknot_values_full = None
    seg_wtd_meanx = None
    nrf = None
    nrf_stderr = None
    is_nonlinear = xknots is not None

    if is_nonlinear:
        xknots_arr = np.asarray(xknots, dtype=float)
        p_xprime, xknot_values_full, seg_wtd_meanx = create_xprime_matrix(
            p_matrix, xknots_arr, xknot_type
        )
        # Update p_matrix to x-prime matrix
        p_matrix_original = p_matrix.copy()
        p_matrix = p_xprime
        n_drivers_original = p_matrix_original.shape[1]
    else:
        n_drivers_original = p_matrix.shape[1]

    # Apply quantile filter to discharge
    q_filtered = _apply_quantile_filter(q_vec, fq, window=4 * m + 1)

    # Build design matrix
    design, response, weights = _build_design_matrix(p_matrix, q_filtered, wt_vec, m)

    # Solve regression (with optional robust estimation)
    if robust:
        beta, stderr, fitted, residuals = _solve_rrd_robust(
            design,
            response,
            weights,
            nu,
            m,
            p_matrix.shape[1],
            maxiter=robust_maxiter,
            tolerance=robust_tolerance,
        )
    else:
        beta, stderr, fitted, residuals = _solve_rrd(
            design, response, weights, nu, m, p_matrix.shape[1]
        )

    # Convert results to appropriate format
    lags = np.arange(beta.shape[0]) * dt

    # If nonlinear, convert beta-prime to NRF
    if is_nonlinear:
        # beta and stderr are currently for x-prime
        # Convert to NRF at knots
        nrf_array, rrd_array = betaprime_to_nrf(
            beta, xknot_values_full, seg_wtd_meanx, n_drivers_original
        )

        # Create labels for NRF columns
        nrf_labels = create_nrf_labels(
            col_labels[:n_drivers_original], xknot_values_full, show_top_xknot
        )

        # Create DataFrames
        nrf = _to_rrd_dataframe(nrf_array, nrf_labels)
        rrd_df = _to_rrd_dataframe(rrd_array, col_labels[:n_drivers_original])

        # For now, use NRF stderr as approximation
        # (Proper error propagation would require covariance matrices)
        nrf_stderr = _to_rrd_dataframe(stderr, nrf_labels)
        stderr_df = _to_rrd_dataframe(
            np.sqrt(np.mean(stderr**2, axis=1, keepdims=True))
            * np.ones((stderr.shape[0], n_drivers_original)),
            col_labels[:n_drivers_original],
        )

    else:
        rrd_df = _to_rrd_dataframe(beta, col_labels)
        stderr_df = _to_rrd_dataframe(stderr, col_labels)

    # Compute broken-stick representation if requested
    knots = None
    if nk and nk > 0:
        if is_nonlinear:
            # Use RRD (averaged over intensities) for broken-stick
            knots = _broken_stick(rrd_array, nk, dt, col_labels[:n_drivers_original])
        else:
            knots = _broken_stick(beta, nk, dt, col_labels)

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
        nrf=nrf,
        nrf_stderr=nrf_stderr,
        xknot_values=xknot_values_full,
        split_labels=split_labels,
        split_criteria=split_criteria,
    )


# ---------------------------------------------------------------------------
# Helper functions / 辅助函数
# ---------------------------------------------------------------------------

# Magic number constants for better readability
_MAD_TO_STD_FACTOR = 1.4826  # Makes MAD consistent with std for normal distribution
_HUBER_TUNING_CONSTANT = 1.345  # Standard choice for 95% efficiency


def _convert_to_numpy_array(
    arr: Sequence[float] | pd.Series | np.ndarray,
) -> np.ndarray:
    """Convert input to a 1-D NumPy array. / 转换为一维 NumPy 数组。"""
    if isinstance(arr, np.ndarray):
        return np.asarray(arr, dtype=float).ravel()
    if isinstance(arr, pd.Series):
        return arr.to_numpy(dtype=float)
    return np.asarray(list(arr), dtype=float)


def _prepare_precipitation(
    p: Iterable[Sequence[float]] | pd.DataFrame | pd.Series | np.ndarray,
    labels: Optional[Sequence[str]],
) -> Tuple[np.ndarray, Sequence[str]]:
    """Ensure precipitation input is a 2-D NumPy array.

    保证降水输入为二维数组，并生成列名称。
    """
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
    """Aggregate time series by agg samples.

    将时间序列按照 agg 长度聚合。降水累加，流量和权重取平均。
    """
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
        wt_arr = _convert_to_numpy_array(wt)[:trimmed]
        wt_agg = wt_arr.reshape(trimmed // agg, agg).mean(axis=1)

    return p_agg, q_agg, wt_agg


def _apply_quantile_filter(q: np.ndarray, fq: float, window: int) -> np.ndarray:
    """Apply a running quantile filter to discharge.

    对流量序列应用滑动分位滤波，以去除季节性或趋势。
    """
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
    """Construct the regression matrix with lags up to m.

    根据最大时滞 m 构建设计矩阵，并删除包含缺测值的行。
    """
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

    # Normalize weights to avoid numerical issues / 归一化权重
    mean_w = weights.mean()
    if not np.isfinite(mean_w) or mean_w <= 0:
        weights = np.ones_like(weights)
    else:
        weights = weights / mean_w

    return design, response, weights


def _create_tikhonov_regularization_matrix(lag_count: int) -> np.ndarray:
    """Construct a 2nd-order difference operator for Tikhonov smoothing.

    构建用于 Tikhonov 平滑的二阶差分算子。
    """
    if lag_count <= 2:
        return np.eye(lag_count, dtype=float)
    L = np.zeros((lag_count - 2, lag_count), dtype=float)
    for i in range(lag_count - 2):
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

    求解加权岭回归。返回的系数以列主序排列：lag0 var1, lag0 var2, ...
    """
    W_sqrt = np.sqrt(weights)[:, None]
    Xw = design * W_sqrt
    yw = response * W_sqrt[:, 0]

    beta_size = (m + 1) * k
    XtX = Xw.T @ Xw
    Xty = Xw.T @ yw

    if nu > 0.0:
        L = _create_tikhonov_regularization_matrix(m + 1)
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


def _solve_rrd_robust(
    design: np.ndarray,
    response: np.ndarray,
    weights: np.ndarray,
    nu: float,
    m: int,
    k: int,
    maxiter: int = 10,
    tolerance: float = 1e-4,
) -> Tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray]:
    """Solve using Iteratively Reweighted Least Squares (IRLS).

    使用迭代重加权最小二乘法 (IRLS) 求解。

    This implements robust regression by iteratively downweighting outliers
    based on their residuals using Huber weights.

    通过基于残差使用 Huber 权重迭代降低异常值的权重来实现鲁棒回归。
    """
    # Start with OLS solution
    beta, stderr, fitted, residuals = _solve_rrd(design, response, weights, nu, m, k)

    rob_weights = np.ones_like(weights)

    for iteration in range(maxiter):
        # Calculate robust weights using Huber function
        # Scale by median absolute deviation (MAD)
        mad = np.median(np.abs(residuals - np.median(residuals)))
        if mad < 1e-10:
            break  # Converged or degenerate

        # Scale residuals by MAD (using factor to make MAD consistent with std)
        scaled_resid = residuals / (_MAD_TO_STD_FACTOR * mad)

        # Huber weights: 1 for |r| <= k, k/|r| for |r| > k
        # Using standard tuning constant for 95% efficiency
        rob_weights = np.where(
            np.abs(scaled_resid) <= _HUBER_TUNING_CONSTANT,
            np.ones_like(scaled_resid),
            _HUBER_TUNING_CONSTANT / np.abs(scaled_resid),
        )

        # Combine original weights with robust weights
        combined_weights = weights * rob_weights

        # Re-solve with new weights
        beta_new, stderr, fitted, residuals = _solve_rrd(
            design, response, combined_weights, nu, m, k
        )

        # Check convergence
        beta_change = np.max(np.abs(beta_new - beta) / (np.abs(beta) + 1e-10))
        beta = beta_new

        if beta_change < tolerance:
            break

    return beta, stderr, fitted, residuals


def _to_rrd_dataframe(arr: np.ndarray, labels: Sequence[str]) -> pd.DataFrame:
    """Convert coefficient array to DataFrame with bilingual headers.

    将系数数组转换为带有双语标题的 DataFrame。
    """
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

    生成折线化的 RRD 表示，便于与原始 R 代码的 nk 参数对应。
    """
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
