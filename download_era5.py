#!/usr/bin/env python3
"""Download the necessary ERA5 data to run the atmospheric RTM.

The ERA5 dataset: https://doi.org/10.24381/cds.bd0915c6

The API reference: https://cds.climate.copernicus.eu/api-how-to

Note that an API key is required to use CDS. After registering for an account
with CDS, run this script with the following two environment variables set:

- CDS_UID
- CDS_API_KEY
"""

import os
import sys
from datetime import date
from pathlib import Path

import cdsapi

CDS_URL = "https://cds.climate.copernicus.eu/api/v2"
try:
    cds_uid = os.environ["CDS_UID"]
    cds_api_key = os.environ["CDS_API_KEY"]
except KeyError:
    print(
        "For CDS authentication, both 'CDS_UID' and "
        "'CDS_API_KEY' environment variables must be set"
    )
    sys.exit(1)

cds_key = ":".join([cds_uid, cds_api_key])

c = cdsapi.Client(url=CDS_URL, key=cds_key, verify=True)

d = date(2020, 1, 1)
out_file = Path(f"era5_{d.isoformat()}.nc")

if not out_file.exists():
    c.retrieve(
        "reanalysis-era5-pressure-levels",
        {
            "product_type": "reanalysis",
            "format": "netcdf",
            "variable": ["temperature", "specific_humidity", "relative_humidity"],
            "pressure_level": ["700", "1000"],
            "year": d.strftime("%Y"),
            "month": d.strftime("%m"),
            "day": d.strftime("%d"),
            "time": ["00:00", "12:00"],
        },
        str(out_file),
    )
    print(f"Data downloaded: {out_file}")
else:
    print(f"File already exists: {out_file}")
