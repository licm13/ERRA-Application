"""Nonlinear response function analysis.

非线性响应函数分析

This module implements nonlinear impulse response function (NRF) estimation,
allowing ERRA to capture how runoff response varies with precipitation intensity.

本模块实现非线性脉冲响应函数（NRF）估计，
使 ERRA 能够捕捉径流响应如何随降水强度变化。

The key concept is to split precipitation into intensity-based segments using
"xknots" (intensity knots), estimate separate responses for each segment, and
then reconstruct the full nonlinear response function.

关键概念是使用"xknots"（强度节点）将降水分割为基于强度的片段，
为每个片段估计单独的响应，然后重建完整的非线性响应函数。
"""

from __future__ import annotations

from typing import List, Literal, Tuple

import numpy as np

# Constant for minimum precipitation threshold
_MIN_PRECIPITATION_VALUE = 0  # Exclude zero precipitation in knot calculations
_EPSILON_WEIGHT = 1e-10  # Small value to prevent division by zero


def create_xprime_matrix(
    p: np.ndarray,
    xknots: np.ndarray,
    xknot_type: Literal[
        "values", "percentiles", "cumsum", "sqsum", "even"
    ] = "percentiles",
) -> Tuple[np.ndarray, np.ndarray, np.ndarray]:
    """Create x' (x-prime) matrix for nonlinear analysis.

    为非线性分析创建 x' (x-prime) 矩阵。

    This function implements the core transformation from precipitation (p) to
    incremental precipitation (x') based on intensity knots. This allows the
    regression to estimate how response varies across different precipitation
    intensity ranges.

    此函数实现从降水 (p) 到基于强度节点的增量降水 (x') 的核心转换。
    这允许回归估计响应如何在不同降水强度范围内变化。

    The transformation follows Eq. 43 from Kirchner (2022):
    x'_i(t) = max(0, min(p(t) - k_i, k_{i+1} - k_i))

    where k_i are the knot points and x'_i represents the increment of
    precipitation between knots k_i and k_{i+1}.

    变换遵循 Kirchner (2022) 的方程 43：
    x'_i(t) = max(0, min(p(t) - k_i, k_{i+1} - k_i))

    其中 k_i 是节点，x'_i 表示节点 k_i 和 k_{i+1} 之间的降水增量。

    Parameters / 参数
    ----------
    p : np.ndarray
        Precipitation matrix (n_timesteps × n_drivers)
        降水矩阵（时间步数 × 驱动变量数）
    xknots : np.ndarray
        Knot values or percentiles, shape (n_knots,) or (n_knots, n_drivers)
        节点值或百分位数，形状 (n_knots,) 或 (n_knots, n_drivers)
    xknot_type : str
        How to interpret xknots:
        如何解释 xknots：
        - "values": Direct precipitation values (直接降水值)
        - "percentiles": Percentiles of p distribution (p 分布的百分位数)
        - "cumsum": Percentiles of cumulative sum of p (p 累积和的百分位数)
        - "sqsum": Percentiles of cumulative sum of p² (p² 累积和的百分位数)
        - "even": Evenly spaced knots (均匀间隔节点)

    Returns / 返回
    -------
    xprime : np.ndarray
        Transformed precipitation matrix (n_timesteps × n_drivers × n_knots)
        转换后的降水矩阵（时间步数 × 驱动变量数 × 节点数）
    knot_values : np.ndarray
        Actual knot values used (including min and max)
        实际使用的节点值（包括最小值和最大值）
    seg_wtd_meanx : np.ndarray
        Segment-weighted mean precipitation in each interval
        每个区间中段加权平均降水
    """
    n_timesteps, n_drivers = p.shape

    # Convert xknots to array if needed
    xknots = np.asarray(xknots, dtype=float)

    # Ensure xknots is 2D
    if xknots.ndim == 1:
        xknots = np.tile(xknots[:, np.newaxis], (1, n_drivers))
    elif xknots.shape[1] != n_drivers:
        raise ValueError(
            f"xknots has {xknots.shape[1]} columns but p has {n_drivers} drivers"
        )

    n_xknots = xknots.shape[0]

    # Calculate actual knot values based on type
    knot_values = np.zeros((n_xknots + 2, n_drivers))  # +2 for min and max

    for i in range(n_drivers):
        p_col = p[:, i]
        p_nonzero = p_col[p_col > _MIN_PRECIPITATION_VALUE]  # Exclude zeros

        if len(p_nonzero) == 0:
            raise ValueError(f"Driver {i} has no positive values")

        if xknot_type == "values":
            # Direct values
            kpts = xknots[:, i]
        elif xknot_type == "percentiles":
            # Percentiles of p
            kpts = np.percentile(p_nonzero, xknots[:, i])
        elif xknot_type == "cumsum":
            # Percentiles of cumulative sum
            p_sorted = np.sort(p_nonzero)
            cumsum = np.cumsum(p_sorted)
            kpts = []
            for pct in xknots[:, i]:
                target = cumsum[-1] * pct / 100
                idx = np.argmin(np.abs(cumsum - target))
                kpts.append(p_sorted[idx])
            kpts = np.array(kpts)
        elif xknot_type == "sqsum":
            # Percentiles of cumulative sum of squares
            p_sorted = np.sort(p_nonzero)
            cumsum_sq = np.cumsum(p_sorted**2)
            kpts = []
            for pct in xknots[:, i]:
                target = cumsum_sq[-1] * pct / 100
                idx = np.argmin(np.abs(cumsum_sq - target))
                kpts.append(p_sorted[idx])
            kpts = np.array(kpts)
        elif xknot_type == "even":
            # Evenly spaced (xknots should specify [n_knots, min_points])
            # For simplicity, use regular percentiles
            # Full implementation would ensure min_points in each interval
            n_knots = int(xknots[0, i])
            percentiles = np.linspace(0, 100, n_knots + 2)[1:-1]
            kpts = np.percentile(p_nonzero, percentiles)
        else:
            raise ValueError(f"Unknown xknot_type: {xknot_type}")

        # Add min (0) and max
        knot_values[:, i] = np.concatenate([[0], kpts, [np.max(p_col)]])

    # Create x-prime matrix
    # For each driver and each knot interval, create an x-prime column
    n_segments = n_xknots + 1  # Number of intervals between knots
    xprime = np.zeros((n_timesteps, n_drivers * n_segments))

    # Calculate segment-weighted means
    seg_wtd_meanx = np.zeros((n_segments, n_drivers))

    for driver_idx in range(n_drivers):
        p_col = p[:, driver_idx]
        kpts = knot_values[:, driver_idx]

        for seg_idx in range(n_segments):
            k_lower = kpts[seg_idx]
            k_upper = kpts[seg_idx + 1]
            interval_width = k_upper - k_lower

            # Calculate x-prime for this segment
            # x'[seg] = max(0, min(p - k_lower, interval_width))
            xp_col = np.maximum(0, np.minimum(p_col - k_lower, interval_width))

            col_idx = driver_idx * n_segments + seg_idx
            xprime[:, col_idx] = xp_col

            # Calculate weighted mean precipitation in this segment
            in_segment = (p_col > k_lower) & (p_col <= k_upper)
            if np.any(in_segment):
                # Weight by precipitation intensity
                p_in_seg = p_col[in_segment]
                seg_wtd_meanx[seg_idx, driver_idx] = np.sum(p_in_seg**2) / np.sum(
                    p_in_seg
                )
            else:
                seg_wtd_meanx[seg_idx, driver_idx] = (k_lower + k_upper) / 2

    return xprime, knot_values, seg_wtd_meanx


def betaprime_to_nrf(
    betaprime: np.ndarray,
    knot_values: np.ndarray,
    seg_wtd_meanx: np.ndarray,
    n_drivers: int,
) -> Tuple[np.ndarray, np.ndarray]:
    """Convert β' (beta-prime) coefficients to NRF (Nonlinear Response Functions).

    将 β'（beta-prime）系数转换为 NRF（非线性响应函数）。

    After regressing q on x', we obtain β' coefficients that represent
    the response per unit increment of x'. To get the actual response
    as a function of precipitation intensity, we need to convert these
    to NRF values at each knot.

    在对 x' 回归 q 后，我们获得 β' 系数，表示每单位 x' 增量的响应。
    为了获得作为降水强度函数的实际响应，我们需要将这些转换为每个节点的 NRF 值。

    The conversion follows:
    NRF(k_i) = Σ_{j=1}^{i} β'_j × (k_j - k_{j-1})

    转换遵循：
    NRF(k_i) = Σ_{j=1}^{i} β'_j × (k_j - k_{j-1})

    Parameters / 参数
    ----------
    betaprime : np.ndarray
        Beta-prime coefficients from regression, shape (m+1, n_drivers × n_segments)
        回归得到的 beta-prime 系数，形状 (m+1, n_drivers × n_segments)
    knot_values : np.ndarray
        Knot values including min and max, shape (n_knots+2, n_drivers)
        包括最小值和最大值的节点值，形状 (n_knots+2, n_drivers)
    seg_wtd_meanx : np.ndarray
        Segment-weighted mean precipitation, shape (n_segments, n_drivers)
        段加权平均降水，形状 (n_segments, n_drivers)
    n_drivers : int
        Number of precipitation drivers
        降水驱动变量数量

    Returns / 返回
    -------
    nrf : np.ndarray
        Nonlinear Response Functions at each knot, shape (m+1, n_drivers × n_segments)
        每个节点的非线性响应函数，形状 (m+1, n_drivers × n_segments)
    rrd : np.ndarray
        Runoff Response Distribution (average NRF), shape (m+1, n_drivers)
        径流响应分布（平均 NRF），形状 (m+1, n_drivers)
    """
    m_plus_1 = betaprime.shape[0]
    n_segments = knot_values.shape[0] - 1

    # Initialize NRF matrix
    nrf = np.zeros_like(betaprime)

    # Convert β' to NRF for each driver
    for driver_idx in range(n_drivers):
        kpts = knot_values[:, driver_idx]

        for lag_idx in range(m_plus_1):
            # For this lag, get all β' values for this driver
            bp_start = driver_idx * n_segments
            bp_end = (driver_idx + 1) * n_segments
            bp_lag = betaprime[lag_idx, bp_start:bp_end]

            # Calculate cumulative NRF at each knot
            # NRF[i] = Σ_{j=0}^{i} β'[j] × (k_{j+1} - k_j)
            for seg_idx in range(n_segments):
                cumulative_response = np.sum(
                    bp_lag[: seg_idx + 1] * np.diff(kpts[: seg_idx + 2])
                )
                nrf[lag_idx, bp_start + seg_idx] = cumulative_response

    # Calculate average RRD (weighted by precipitation)
    rrd = np.zeros((m_plus_1, n_drivers))

    for driver_idx in range(n_drivers):
        bp_start = driver_idx * n_segments
        bp_end = (driver_idx + 1) * n_segments

        # Weight by segment-weighted mean precipitation
        weights = seg_wtd_meanx[:, driver_idx]
        weights = weights / (np.sum(weights) + _EPSILON_WEIGHT)  # Normalize

        for lag_idx in range(m_plus_1):
            rrd[lag_idx, driver_idx] = np.sum(nrf[lag_idx, bp_start:bp_end] * weights)

    return nrf, rrd


def create_nrf_labels(
    base_labels: List[str],
    knot_values: np.ndarray,
    show_top_knot: bool = False,
) -> List[str]:
    """Create labels for NRF columns.

    为 NRF 列创建标签。

    Parameters / 参数
    ----------
    base_labels : list of str
        Original precipitation driver labels
        原始降水驱动标签
    knot_values : np.ndarray
        Knot values, shape (n_knots+2, n_drivers)
        节点值，形状 (n_knots+2, n_drivers)
    show_top_knot : bool
        Whether to include the top knot in labels
        是否在标签中包含顶部节点

    Returns / 返回
    -------
    labels : list of str
        Labels for each NRF column
        每个 NRF 列的标签
    """
    n_segments = knot_values.shape[0] - 1

    labels = []
    for driver_idx, base_label in enumerate(base_labels):
        kpts = knot_values[:, driver_idx]

        n_segs_to_show = n_segments if show_top_knot else n_segments - 1

        for seg_idx in range(n_segs_to_show):
            k_lower = kpts[seg_idx]
            k_upper = kpts[seg_idx + 1]
            label = f"{base_label}_k{seg_idx}({k_lower:.2f}-{k_upper:.2f})"
            labels.append(label)

    return labels
