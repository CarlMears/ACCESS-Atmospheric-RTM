#!/usr/bin/env python3

"""Download the necessary ERA5 data to run the atmospheric RTM.

For reference: https://cds.climate.copernicus.eu/api-how-to

Note that an API key is required to use CDS. After registering, run this script
with the following two environment variables set:

- CDS_UID
- CDS_API_KEY
"""

import os
import sys

import cdsapi

CDS_URL = "https://cds.climate.copernicus.eu/api/v2"
try:
    cds_uid = os.environ["CDS_UID"]
    cds_api_key = os.environ["CDS_API_KEY"]
except KeyError:
    print(
        "For CDS authentication, both 'CDS_UID' and 'CDS_API_KEY' environment variables must be set"
    )
    sys.exit(1)

cds_key = ":".join([cds_uid, cds_api_key])

c = cdsapi.Client(url=CDS_URL, key=cds_key, verify=True)

c.retrieve(
    "reanalysis-era5-pressure-levels",
    {
        "variable": "temperature",
        "pressure_level": "1000",
        "product_type": "reanalysis",
        "year": "2020",
        "month": "01",
        "day": "01",
        "time": "12:00",
        "format": "netcdf",
    },
    "download.nc",
)
