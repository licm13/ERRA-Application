"""Download and preprocess data for the Gao et al. (2025) case study.

This script retrieves daily discharge from the USGS National Water
Information System and daily precipitation from the NOAA Global
Historical Climatology Network (GHCN). The data are transformed into the
three ERRA precipitation drivers used in the Gao et al. example:
convective bursts, stratiform rain, and a slow groundwater recharge
proxy.

Data sources / 数据来源
---------------------
- USGS NWIS Daily Discharge (parameter 00060, cubic feet per second)
- NOAA GHCN Daily Summaries (PRCP data type, metric units)
"""

from __future__ import annotations

import json
import logging
from pathlib import Path
from typing import Dict
from urllib.error import HTTPError, URLError
from urllib.parse import urlencode
from urllib.request import urlopen

import numpy as np
import pandas as pd

LOGGER = logging.getLogger("gao2025_data_prep")

CONFIG: Dict[str, str] = {
    "usgs_site": "11477000",  # Russian River near Healdsburg, CA
    "noaa_station": "USC00043110",  # Healdsburg cooperative station
    "start_date": "2010-01-01",
    "end_date": "2020-12-31",
}

DATA_DIR = Path(__file__).resolve().parent.parent / "data"
RAW_DIR = DATA_DIR / "raw"
PROCESSED_DIR = DATA_DIR / "processed"


def _download_json(url: str) -> Dict:
    with urlopen(url) as response:
        text = response.read().decode("utf-8")
    return json.loads(text)


def fetch_usgs_streamflow(site: str, start: str, end: str) -> pd.Series:
    """Fetch daily discharge (cfs) and convert to m³/s."""

    params = {
        "format": "json",
        "sites": site,
        "startDT": start,
        "endDT": end,
        "parameterCd": "00060",
        "siteStatus": "all",
    }
    url = "https://waterservices.usgs.gov/nwis/dv/?" + urlencode(params)
    LOGGER.info("Downloading USGS streamflow: %s", url)
    data = _download_json(url)

    series_data = []
    try:
        values = data["value"]["timeSeries"][0]["values"][0]["value"]
    except (KeyError, IndexError) as exc:
        raise RuntimeError("Unexpected USGS JSON structure") from exc

    for entry in values:
        value = entry.get("value")
        if value in (None, ""):
            continue
        date = entry["dateTime"][0:10]
        try:
            discharge_cfs = float(value)
        except ValueError:
            continue
        discharge_cms = discharge_cfs * 0.0283168466
        series_data.append((date, discharge_cms))

    if not series_data:
        raise RuntimeError("No streamflow records downloaded")

    df = pd.DataFrame(series_data, columns=["date", "discharge_cms"])
    df["date"] = pd.to_datetime(df["date"])
    df = df.set_index("date").sort_index()
    return df["discharge_cms"]


def fetch_noaa_precip(station: str, start: str, end: str) -> pd.Series:
    params = {
        "dataset": "daily-summaries",
        "stations": station,
        "startDate": start,
        "endDate": end,
        "dataTypes": "PRCP",
        "units": "metric",
        "format": "json",
    }
    url = "https://www.ncei.noaa.gov/access/services/data/v1?" + urlencode(params)
    LOGGER.info("Downloading NOAA precipitation: %s", url)
    with urlopen(url) as response:
        text = response.read().decode("utf-8")

    records = json.loads(text)
    rows = []
    for item in records:
        date = pd.to_datetime(item["DATE"])
        prcp = item.get("PRCP")
        if prcp in (None, ""):
            continue
        rows.append((date, float(prcp)))

    if not rows:
        raise RuntimeError("No precipitation records downloaded")

    df = pd.DataFrame(rows, columns=["date", "prcp_mm"])
    df = df.set_index("date").sort_index()
    return df["prcp_mm"]


def build_erra_inputs(precip_mm: pd.Series, discharge_cms: pd.Series) -> pd.DataFrame:
    data = pd.concat([precip_mm, discharge_cms], axis=1, join="inner").dropna()
    data.columns = ["precip_mm", "discharge_cms"]

    prcp = data["precip_mm"].clip(lower=0)
    convective_threshold = prcp.quantile(0.85)
    convective = (prcp - convective_threshold).clip(lower=0)
    stratiform = prcp - convective
    recharge = prcp.rolling(window=14, min_periods=1).mean()

    drivers = pd.DataFrame(
        {
            "convective_bursts": convective,
            "stratiform_rain": stratiform,
            "recharge_proxy": recharge,
            "discharge": data["discharge_cms"],
            "weights": np.ones(len(data), dtype=float),
        },
        index=data.index,
    )
    return drivers


def main() -> None:
    logging.basicConfig(level=logging.INFO, format="%(levelname)s %(message)s")
    RAW_DIR.mkdir(parents=True, exist_ok=True)
    PROCESSED_DIR.mkdir(parents=True, exist_ok=True)

    try:
        discharge = fetch_usgs_streamflow(
            CONFIG["usgs_site"], CONFIG["start_date"], CONFIG["end_date"]
        )
        precip = fetch_noaa_precip(
            CONFIG["noaa_station"], CONFIG["start_date"], CONFIG["end_date"]
        )
    except (HTTPError, URLError, RuntimeError) as exc:
        LOGGER.error("Failed to download data: %s", exc)
        raise SystemExit(1) from exc

    raw_stream_path = RAW_DIR / "gao2025_usgs_streamflow.csv"
    raw_precip_path = RAW_DIR / "gao2025_noaa_precip.csv"
    discharge.to_frame(name="discharge_cms").to_csv(raw_stream_path)
    precip.to_frame(name="precip_mm").to_csv(raw_precip_path)
    LOGGER.info("Saved raw data to %s and %s", raw_stream_path, raw_precip_path)

    processed = build_erra_inputs(precip, discharge)
    processed_path = PROCESSED_DIR / "gao2025_inputs.csv"
    processed.to_csv(processed_path, index_label="date")
    LOGGER.info("Processed ERRA inputs saved to %s", processed_path)


if __name__ == "__main__":
    main()
