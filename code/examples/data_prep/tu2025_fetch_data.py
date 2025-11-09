"""Download data for the Tu et al. (2025) permafrost transition scenario."""

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

LOGGER = logging.getLogger("tu2025_data_prep")

CONFIG: Dict[str, str] = {
    "usgs_site": "15515500",  # Beaver Creek near Fairbanks, AK
    "noaa_station": "USW00026411",  # Fairbanks International Airport
    "start_date": "2010-01-01",
    "end_date": "2020-12-31",
}

DATA_DIR = Path(__file__).resolve().parent.parent / "data"
RAW_DIR = DATA_DIR / "raw"
PROCESSED_DIR = DATA_DIR / "processed"


def _download_json(url: str):
    with urlopen(url) as response:
        return json.loads(response.read().decode("utf-8"))


def fetch_usgs_streamflow(site: str, start: str, end: str) -> pd.Series:
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
    values = data["value"]["timeSeries"][0]["values"][0]["value"]

    rows = []
    for entry in values:
        value = entry.get("value")
        if value in (None, ""):
            continue
        try:
            discharge_cfs = float(value)
        except ValueError:
            continue
        rows.append((entry["dateTime"][0:10], discharge_cfs * 0.0283168466))

    if not rows:
        raise RuntimeError("No streamflow records downloaded")

    df = pd.DataFrame(rows, columns=["date", "discharge_cms"])
    df["date"] = pd.to_datetime(df["date"])
    df = df.set_index("date").sort_index()
    return df["discharge_cms"]


def fetch_noaa_climate(station: str, start: str, end: str) -> pd.DataFrame:
    params = {
        "dataset": "daily-summaries",
        "stations": station,
        "startDate": start,
        "endDate": end,
        "dataTypes": "PRCP,TAVG,SNWD",
        "units": "metric",
        "format": "json",
    }
    url = "https://www.ncei.noaa.gov/access/services/data/v1?" + urlencode(params)
    LOGGER.info("Downloading NOAA climate: %s", url)
    with urlopen(url) as response:
        records = json.loads(response.read().decode("utf-8"))

    rows = []
    for item in records:
        date = pd.to_datetime(item["DATE"])
        prcp = item.get("PRCP")
        tavg = item.get("TAVG")
        snwd = item.get("SNWD")
        rows.append((date, float(prcp) if prcp not in (None, "") else np.nan,
                     float(tavg) if tavg not in (None, "") else np.nan,
                     float(snwd) if snwd not in (None, "") else np.nan))

    df = pd.DataFrame(rows, columns=["date", "prcp_mm", "tavg_c", "snwd_mm"])
    df = df.set_index("date").sort_index()
    return df


def build_erra_inputs(climate: pd.DataFrame, discharge_cms: pd.Series) -> pd.DataFrame:
    data = climate.join(discharge_cms, how="inner").dropna(subset=["discharge_cms"])
    data["prcp_mm"] = data["prcp_mm"].clip(lower=0)

    rainfall = data.loc[data["tavg_c"] > 1.0, "prcp_mm"].fillna(0.0)
    rainfall = rainfall.reindex(data.index, fill_value=0.0)

    snow_depth = data["snwd_mm"].fillna(method="ffill").fillna(0.0)
    snowmelt = (-snow_depth.diff()).clip(lower=0.0)
    snowmelt = snowmelt.fillna(0.0)

    positive_degree_days = (data["tavg_c"].clip(lower=0) / 10.0).rolling(
        window=30, min_periods=1
    ).sum()
    active_layer = positive_degree_days / positive_degree_days.max()

    drivers = pd.DataFrame(
        {
            "rainfall": rainfall,
            "snowmelt": snowmelt,
            "active_layer": active_layer,
            "discharge": data["discharge_cms"],
            "weights": 1.0 + 0.1 * active_layer.fillna(0.0),
        },
        index=data.index,
    ).fillna(0.0)
    return drivers


def main() -> None:
    logging.basicConfig(level=logging.INFO, format="%(levelname)s %(message)s")
    RAW_DIR.mkdir(parents=True, exist_ok=True)
    PROCESSED_DIR.mkdir(parents=True, exist_ok=True)

    try:
        discharge = fetch_usgs_streamflow(
            CONFIG["usgs_site"], CONFIG["start_date"], CONFIG["end_date"]
        )
        climate = fetch_noaa_climate(
            CONFIG["noaa_station"], CONFIG["start_date"], CONFIG["end_date"]
        )
    except (HTTPError, URLError, RuntimeError) as exc:
        LOGGER.error("Failed to download data: %s", exc)
        raise SystemExit(1) from exc

    discharge.to_frame(name="discharge_cms").to_csv(
        RAW_DIR / "tu2025_usgs_streamflow.csv"
    )
    climate.to_csv(RAW_DIR / "tu2025_noaa_climate.csv")

    processed = build_erra_inputs(climate, discharge)
    processed_path = PROCESSED_DIR / "tu2025_inputs.csv"
    processed.to_csv(processed_path, index_label="date")
    LOGGER.info("Processed ERRA inputs saved to %s", processed_path)


if __name__ == "__main__":
    main()
