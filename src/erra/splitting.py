"""Data splitting utilities for non-stationary analysis.

非平稳性分析的数据分割工具

This module provides functions to split precipitation time series based on
external covariates (e.g., antecedent wetness, temperature) for analyzing
non-stationary hydrologic behavior.

本模块提供基于外部协变量（如前期湿度、温度）分割降水时间序列的功能，
用于分析非平稳水文行为。
"""

from __future__ import annotations

from typing import Dict, List, Optional, Tuple

import numpy as np
import pandas as pd

import logging


logger = logging.getLogger(__name__)


def make_split_sets(
    p: np.ndarray,
    q: np.ndarray,
    split_params: Dict,
    wt: Optional[np.ndarray] = None,
    *,
    min_bin_size: int = 15,
) -> Tuple[np.ndarray, np.ndarray, np.ndarray, List[str], pd.DataFrame]:
    """Split precipitation time series based on external criteria.

    根据外部标准分割降水时间序列。

    This function implements the splitting logic from the R ERRA script,
    allowing users to analyze how runoff response varies under different
    conditions (e.g., wet vs. dry antecedent conditions, different temperatures).

    此函数实现 R ERRA 脚本中的分割逻辑，允许用户分析不同条件下径流响应
    如何变化（例如，湿润与干燥的前期条件、不同温度）。

    Parameters / 参数
    ----------
    p : np.ndarray
        Precipitation matrix (n_timesteps × n_drivers)
        降水矩阵（时间步数 × 驱动变量数）
    q : np.ndarray
        Discharge vector (n_timesteps,)
        流量向量（时间步数，）
    split_params : dict
        Dictionary containing splitting parameters:
        包含分割参数的字典：
        - crit: list of criterion arrays (协变量数组列表)
        - crit_label: list of criterion labels (协变量标签列表)
        - crit_lag: list of lag values for each criterion (每个协变量的滞后值)
        - pct_breakpts: list of booleans (是否使用百分位断点)
        - breakpts: list of breakpoint lists (断点列表)
        - thresh: list of threshold values (阈值列表)
        - by_bin: list of booleans (是否在每个 bin 内设置断点)
    wt : np.ndarray, optional
        Observation weights (观测权重)
    min_bin_size : int, optional
        Minimum number of observations recommended for each split subset.
        每个分组建议的最小样本数量。当数据不足时将记录日志警告。

    Returns / 返回
    -------
    p_split : np.ndarray
        Split precipitation matrix with interleaved subsets
        带交错子集的分割降水矩阵
    q_out : np.ndarray
        Discharge vector (unchanged)
        流量向量（未更改）
    wt_out : np.ndarray
        Updated weights
        更新后的权重
    set_labels : list of str
        Labels for each split subset
        每个分割子集的标签
    criteria_df : pd.DataFrame
        DataFrame documenting the splitting criteria for each subset
        记录每个子集分割标准的 DataFrame

    Examples / 示例
    --------
    >>> # Split precipitation by antecedent discharge (wetness proxy)
    >>> # 根据前期流量（湿润度代理）分割降水
    >>> split_params = {
    ...     'crit': [q],
    ...     'crit_label': ['Antecedent Q'],
    ...     'crit_lag': [1],  # Use lagged discharge
    ...     'pct_breakpts': [True],
    ...     'breakpts': [[50, 90]],  # Split at 50th and 90th percentiles
    ...     'thresh': [0],  # Ignore zero discharge
    ...     'by_bin': [True]
    ... }
    >>> p_split, q_out, wt_out, labels, criteria = make_split_sets(
    ...     p, q, split_params
    ... )
    """
    if not isinstance(p, np.ndarray):
        raise TypeError("p must be a numpy.ndarray")
    if not isinstance(q, np.ndarray):
        raise TypeError("q must be a numpy.ndarray")
    if wt is not None and not isinstance(wt, np.ndarray):
        raise TypeError("wt must be a numpy.ndarray when provided")

    if p.ndim == 1:
        p = p.reshape(-1, 1)
    if q.ndim != 1:
        raise ValueError("q must be one-dimensional")

    n = len(q)
    if p.shape[0] != n:
        raise ValueError(
            "p and q must have the same number of observations "
            f"(got p rows={p.shape[0]}, q length={n})"
        )

    if min_bin_size < 1:
        raise ValueError("min_bin_size must be >= 1")

    n_drivers = p.shape[1]

    # Extract split parameters
    crit_list = split_params["crit"]
    crit_label = split_params["crit_label"]
    crit_lag = split_params["crit_lag"]
    pct_breakpts = split_params["pct_breakpts"]
    breakpts = split_params["breakpts"]
    thresh = split_params["thresh"]
    by_bin = split_params["by_bin"]

    n_crit = len(crit_list)

    # Validate inputs
    if not (
        len(crit_label) == n_crit
        and len(crit_lag) == n_crit
        and len(pct_breakpts) == n_crit
        and len(breakpts) == n_crit
        and len(thresh) == n_crit
        and len(by_bin) == n_crit
    ):
        raise ValueError("All split_params lists must have the same length")

    # Apply lags to criteria
    lagged_crit = []
    for i, (crit, lag) in enumerate(zip(crit_list, crit_lag)):
        crit_arr = np.asarray(crit, dtype=float)
        if len(crit_arr) != n:
            raise ValueError(f"Criterion {i} must have same length as q")

        if lag >= 1:
            crit_lagged = np.concatenate([np.full(lag, np.nan), crit_arr[:-lag]])
        else:
            crit_lagged = crit_arr.copy()

        lagged_crit.append(crit_lagged)

    # Calculate breakpoints and create bin assignments
    bin_assignments = np.zeros((n, n_crit), dtype=int)
    breakpt_values: List[Dict[Optional[Tuple[int, ...]], np.ndarray]] = []

    logger.debug("Starting split computation for %d criteria", n_crit)

    for i in range(n_crit):
        crit_values = lagged_crit[i]
        bp_list = breakpts[i]

        if len(bp_list) == 0:
            global_bp_vals = np.array([], dtype=float)
            valid_mask = (~np.isnan(crit_values)) & (crit_values > thresh[i])
        elif pct_breakpts[i]:
            valid_mask = (~np.isnan(crit_values)) & (crit_values > thresh[i])
            valid_values = crit_values[valid_mask]

            if valid_values.size == 0:
                logger.warning(
                    "Criterion '%s' has no values above threshold %.3f; using NaNs.",
                    crit_label[i],
                    thresh[i],
                )
                global_bp_vals = np.full(len(bp_list), np.nan)
            elif np.unique(valid_values).size == 1:
                single_val = float(np.unique(valid_values)[0])
                logger.info(
                    "Criterion '%s' valid data collapses to %.3f; using repeated breakpoints.",
                    crit_label[i],
                    single_val,
                )
                global_bp_vals = np.full(len(bp_list), single_val)
            else:
                global_bp_vals = np.percentile(valid_values, bp_list)
        else:
            valid_mask = ~np.isnan(crit_values)
            global_bp_vals = np.array(bp_list, dtype=float)

        if not by_bin[i] or i == 0:
            bp_map: Dict[Optional[Tuple[int, ...]], np.ndarray] = {(): global_bp_vals}
            bins = np.digitize(crit_values, global_bp_vals)
        else:
            bp_map = {}
            prev_assignments = bin_assignments[:, :i]
            combos, inverse_idx = np.unique(
                prev_assignments, axis=0, return_inverse=True
            )
            bins = np.zeros(n, dtype=int)

            for combo_idx, combo in enumerate(combos):
                combo_key = tuple(int(x) for x in combo)
                combo_mask = inverse_idx == combo_idx
                combo_valid_mask = combo_mask & valid_mask
                combo_count = int(combo_valid_mask.sum())

                if pct_breakpts[i] and combo_count > 1:
                    combo_values = crit_values[combo_valid_mask]
                    if combo_count < min_bin_size:
                        logger.warning(
                            "Criterion '%s' bin %s has only %d valid observations; consider merging bins or adjusting thresholds.",
                            crit_label[i],
                            combo_key,
                            combo_count,
                        )

                    if np.unique(combo_values).size == 1:
                        bp_vals = np.full(len(bp_list), float(combo_values[0]))
                    else:
                        bp_vals = np.percentile(combo_values, bp_list)
                else:
                    if combo_count <= 1:
                        logger.warning(
                            "Criterion '%s' bin %s lacks sufficient data (count=%d); falling back to global breakpoints.",
                            crit_label[i],
                            combo_key,
                            combo_count,
                        )
                    bp_vals = global_bp_vals

                bp_map[combo_key] = bp_vals
                bins[combo_mask] = np.digitize(crit_values[combo_mask], bp_vals)

            if () not in bp_map:
                bp_map[()] = global_bp_vals

        breakpt_values.append(bp_map)
        bin_assignments[:, i] = bins

    # Create unique bin combinations and labels
    # Each row gets a unique combination of bin indices
    n_bins_per_crit = [len(bp) + 1 for bp in breakpts]
    total_bins = np.prod(n_bins_per_crit)

    # Create split precipitation matrix
    # Each original p column is split into total_bins subsets
    p_split = np.zeros((n, n_drivers * total_bins))

    # Create set labels and criteria dataframe
    set_labels = []
    criteria_records = []

    bin_idx = 0
    subset_counts: Dict[str, int] = {}
    for bin_combo in _generate_bin_combinations(n_bins_per_crit):
        # Create mask for this bin combination
        mask = np.ones(n, dtype=bool)
        label_parts = []
        criteria_row = {}

        for i, bin_num in enumerate(bin_combo):
            mask &= bin_assignments[:, i] == bin_num
            label_parts.append(f"{crit_label[i]}_{bin_num}")

            # Record criteria bounds for this bin
            bp_map = breakpt_values[i]
            prev_combo = tuple(bin_combo[:i])
            bp_vals = _resolve_breakpoints(bp_map, prev_combo)

            if bin_num == 0:
                lower = -np.inf
                upper = bp_vals[0] if len(bp_vals) > 0 else np.inf
            elif bin_num == len(bp_vals):
                lower = bp_vals[-1]
                upper = np.inf
            else:
                lower = bp_vals[bin_num - 1]
                upper = bp_vals[bin_num]

            criteria_row[f"{crit_label[i]}_lower"] = lower
            criteria_row[f"{crit_label[i]}_upper"] = upper

            # Calculate mean value in this bin
            bin_mask = bin_assignments[:, i] == bin_num
            if np.any(bin_mask):
                criteria_row[f"{crit_label[i]}_mean"] = np.nanmean(
                    lagged_crit[i][bin_mask]
                )
            else:
                criteria_row[f"{crit_label[i]}_mean"] = np.nan

        set_label = "_".join(label_parts)
        set_labels.append(set_label)
        criteria_records.append(criteria_row)

        observation_count = int(mask.sum())
        subset_counts[set_label] = observation_count
        if observation_count == 0:
            logger.info(
                "Split subset %s is empty. Consider revising breakpoints or thresholds.",
                set_label,
            )
        elif observation_count < min_bin_size:
            logger.warning(
                "Split subset %s has only %d observations; statistical estimates may be noisy.",
                set_label,
                observation_count,
            )

        # Fill in precipitation for this bin combination
        for driver_idx in range(n_drivers):
            col_idx = bin_idx * n_drivers + driver_idx
            p_split[mask, col_idx] = p[mask, driver_idx]

        bin_idx += 1

    # Create criteria dataframe
    criteria_df = pd.DataFrame(criteria_records, index=set_labels)
    if subset_counts:
        criteria_df["count"] = pd.Series(subset_counts)
        logger.debug("Subset observation counts: %s", subset_counts)

    # Handle weights
    if wt is None:
        wt_out = np.ones(n)
    else:
        if len(wt) != n:
            raise ValueError(
                "Weights must have the same length as q "
                f"(got {len(wt)} vs {n})"
            )
        wt_out = wt.copy()

    return p_split, q, wt_out, set_labels, criteria_df


def _generate_bin_combinations(n_bins_per_crit: List[int]):
    """Generate all combinations of bin indices.

    生成所有 bin 索引组合。

    Parameters / 参数
    ----------
    n_bins_per_crit : list of int
        Number of bins for each criterion
        每个标准的 bin 数量

    Yields / 生成
    ------
    tuple
        Each combination of bin indices
        每个 bin 索引组合
    """
    if len(n_bins_per_crit) == 0:
        yield ()
        return

    for i in range(n_bins_per_crit[0]):
        if len(n_bins_per_crit) == 1:
            yield (i,)
        else:
            for rest in _generate_bin_combinations(n_bins_per_crit[1:]):
                yield (i,) + rest


def _resolve_breakpoints(
    bp_map: Dict[Optional[Tuple[int, ...]], np.ndarray],
    prev_combo: Tuple[int, ...],
) -> np.ndarray:
    """Resolve the breakpoint array for a specific combination of previous bins."""

    if prev_combo in bp_map:
        return bp_map[prev_combo]

    if () in bp_map:
        return bp_map[()]

    if None in bp_map:
        return bp_map[None]

    if bp_map:
        first_key = next(iter(bp_map))
        return bp_map[first_key]

    return np.array([])


def validate_split_params(split_params: Dict) -> None:
    """Validate split_params dictionary structure.

    验证 split_params 字典结构。

    Parameters / 参数
    ----------
    split_params : dict
        Splitting parameters to validate
        要验证的分割参数

    Raises / 抛出
    ------
    ValueError
        If split_params is invalid
        如果 split_params 无效
    """
    required_keys = [
        "crit",
        "crit_label",
        "crit_lag",
        "pct_breakpts",
        "breakpts",
        "thresh",
        "by_bin",
    ]

    for key in required_keys:
        if key not in split_params:
            raise ValueError(f"split_params must contain '{key}'")

    n = len(split_params["crit"])
    for key in required_keys:
        if len(split_params[key]) != n:
            raise ValueError(
                f"All split_params lists must have same length "
                f"(got {key}={len(split_params[key])}, expected {n})"
            )
