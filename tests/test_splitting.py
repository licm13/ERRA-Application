import logging
from pathlib import Path
import sys

import numpy as np

sys.path.append(str(Path(__file__).resolve().parents[1] / "src"))

from erra.splitting import make_split_sets


def _base_split_params(crit1, crit2):
    return {
        "crit": [crit1, crit2],
        "crit_label": ["C1", "C2"],
        "crit_lag": [0, 0],
        "pct_breakpts": [True, True],
        "breakpts": [[50], [50]],
        "thresh": [0.0, 0.0],
        "by_bin": [False, False],
    }


def test_make_split_sets_global_percentiles_matches_expected():
    p = np.arange(24, dtype=float).reshape(12, 2)
    q = np.linspace(0.1, 1.2, 12)
    crit1 = np.linspace(0.0, 11.0, 12)
    crit2 = np.linspace(10.0, 21.0, 12)

    params = _base_split_params(crit1, crit2)
    params["pct_breakpts"][1] = False
    params["breakpts"][1] = [14.5]

    p_split, q_out, wt_out, labels, criteria = make_split_sets(p, q, params)

    assert np.shares_memory(q_out, q)
    assert np.all(wt_out == 1.0)

    c1_bins = np.digitize(crit1, np.percentile(crit1, [50]))
    c2_bins = np.digitize(crit2, np.array(params["breakpts"][1], dtype=float))
    counts = {(i, j): 0 for i in range(2) for j in range(2)}
    for b1, b2 in zip(c1_bins, c2_bins):
        counts[(b1, b2)] += 1

    for combo, expected in counts.items():
        label = f"C1_{combo[0]}_C2_{combo[1]}"
        assert label in labels
        assert int(criteria.loc[label, "count"]) == expected

    global_breaks_c1 = np.percentile(crit1[crit1 > 0], [50])
    np.testing.assert_allclose(
        criteria.loc["C1_0_C2_0", "C1_upper"], global_breaks_c1[0]
    )


def test_make_split_sets_hierarchical_percentiles():
    p = np.tile(np.arange(6, dtype=float), (2, 1)).T
    q = np.linspace(0.0, 5.0, 6)
    crit1 = np.array([0, 0, 1, 1, 2, 2], dtype=float)
    crit2 = np.array([1, 2, 3, 4, 10, 20], dtype=float)

    params = _base_split_params(crit1, crit2)
    params["breakpts"][0] = [33, 66]
    params["by_bin"][1] = True

    _, _, _, labels, criteria = make_split_sets(p, q, params, min_bin_size=1)

    expected_breaks = {0: 1.5, 1: 3.5, 2: 15.0}
    for idx, break_value in expected_breaks.items():
        label = f"C1_{idx}_C2_1"
        assert label in labels
        np.testing.assert_allclose(
            criteria.loc[label, "C2_lower"], break_value
        )


def test_make_split_sets_fallback_to_global_percentiles(caplog):
    p = np.arange(20, dtype=float).reshape(10, 2)
    q = np.linspace(0.0, 0.9, 10)
    crit1 = np.linspace(0.0, 9.0, 10)
    crit2 = np.array([0.1, 0.2, 0.3, 0.4, np.nan, np.nan, np.nan, np.nan, np.nan, np.nan])

    params = _base_split_params(crit1, crit2)
    params["by_bin"][1] = True

    with caplog.at_level(logging.WARNING):
        _, _, _, labels, criteria = make_split_sets(
            p, q, params, min_bin_size=2
        )

    warning_messages = " ".join(rec.message for rec in caplog.records)
    assert "lacks sufficient data" in warning_messages

    global_breaks = np.percentile(crit2[~np.isnan(crit2)], [50])
    label = "C1_1_C2_1"
    if label in labels:
        np.testing.assert_allclose(
            criteria.loc[label, "C2_lower"], global_breaks[0]
        )
