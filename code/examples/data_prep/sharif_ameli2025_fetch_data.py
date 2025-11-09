"""Download data for the Sharif & Ameli (2025) functional simplicity case.

This script mirrors the Gao downloader but targets a semiarid basin in
Arizona to emulate the wet/dry regime transitions discussed by Sharif &
Ameli (2025). The precipitation drivers are derived from NOAA GHCN
records and the discharge driver from USGS NWIS.
"""

from __future__ import annotations

import json
import logging
from pathlib import Path
from typing import Dict
from urllib.error import HTTPError, URLError
from urllib.parse import urlencode
from urllib.request import urlopen

import pandas as pd

LOGGER = logging.getLogger("sharif_ameli_data_prep")

CONFIG: Dict[str, str] = {
    "usgs_site": "09506000",  # Sabino Creek near Tucson, AZ
    "noaa_station": "USW00023160",  # Tucson International Airport
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
        records = json.loads(response.read().decode("utf-8"))

    rows = []
    for item in records:
        prcp = item.get("PRCP")
        if prcp in (None, ""):
            continue
        rows.append((pd.to_datetime(item["DATE"]), float(prcp)))

    if not rows:
        raise RuntimeError("No precipitation records downloaded")

    df = pd.DataFrame(rows, columns=["date", "prcp_mm"])
    df = df.set_index("date").sort_index()
    return df["prcp_mm"]


def build_erra_inputs(precip_mm: pd.Series, discharge_cms: pd.Series) -> pd.DataFrame:
    data = pd.concat([precip_mm, discharge_cms], axis=1, join="inner").dropna()
    data.columns = ["precip_mm", "discharge_cms"]

    prcp = data["precip_mm"].clip(lower=0)
    high_threshold = prcp.quantile(0.9)
    burst = (prcp - high_threshold).clip(lower=0)

    antecedent = prcp.ewm(alpha=0.25, adjust=False).mean()
    fractured = 0.6 * burst + 0.4 * data["discharge_cms"].rolling(
        window=7, min_periods=1
    ).mean()
    weights = 1.0 + 0.05 * (antecedent / antecedent.max()).fillna(0.0)

    drivers = pd.DataFrame(
        {
            "burst_rain": burst,
            "antecedent_moisture": antecedent,
            "fractured_bedrock": fractured,
            "discharge": data["discharge_cms"],
            "weights": weights,
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

    discharge.to_frame(name="discharge_cms").to_csv(
        RAW_DIR / "sharif_ameli_usgs_streamflow.csv"
    )
    precip.to_frame(name="precip_mm").to_csv(
        RAW_DIR / "sharif_ameli_noaa_precip.csv"
    )

    processed = build_erra_inputs(precip, discharge)
    processed_path = PROCESSED_DIR / "sharif_ameli2025_inputs.csv"
    processed.to_csv(processed_path, index_label="date")
    LOGGER.info("Processed ERRA inputs saved to %s", processed_path)


if __name__ == "__main__":
    main()
