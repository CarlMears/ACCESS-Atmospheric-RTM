#!/usr/bin/env python3
"""Download the necessary ERA5 data to run the atmospheric RTM.

The ERA5 datasets:

- https://doi.org/10.24381/cds.bd0915c6
- https://doi.org/10.24381/cds.adbb2d47

The API reference: https://cds.climate.copernicus.eu/api-how-to

Note that an API key is required to use CDS. After registering for an account
with CDS, run this script with the following two environment variables set:

- CDS_UID
- CDS_API_KEY
"""

import argparse
import os
from datetime import date, timedelta
from pathlib import Path

import cdsapi

# CDS API
CDS_URL = "https://cds.climate.copernicus.eu/api/v2"

# These pressure levels (in hPa) are the same as the NCEP GDAS data but they are
# a subset of the available ERA5 levels. (In other words, more could be added if
# desired.)
LEVELS = [
    "1000",
    "975",
    "950",
    "925",
    "900",
    "850",
    "800",
    "750",
    "700",
    "650",
    "600",
    "550",
    "500",
    "450",
    "400",
    "350",
    "300",
    "250",
    "200",
    "150",
    "100",
    "70",
    "50",
    "30",
    "20",
    "10",
]

# The ERA5 variables for the pressure-level analysis product
ATM_VARIABLES = [
    "temperature",
    "geopotential",
    "relative_humidity",
    "specific_cloud_liquid_water_content",
]


# The ERA5 variables for the surface-level analysis product
SURFACE_VARIABLES = [
    "surface_pressure",
    "geopotential",
    "total_column_cloud_liquid_water",
    "total_column_water_vapour",
    "2m_temperature",
    "2m_dewpoint_temperature",
]

# The hours per day to download
HOURS = ["00:00", "12:00"]


def main() -> None:
    """Parse arguments and download the data."""
    parser = argparse.ArgumentParser(description="Download relevant ERA5 data")
    parser.add_argument("start_date", type=date.fromisoformat, help="starting day")
    parser.add_argument("end_date", type=date.fromisoformat, help="ending day")
    parser.add_argument(
        "--out-dir",
        type=Path,
        default=Path.cwd(),
        metavar="DIR",
        help="directory to write files",
    )
    args = parser.parse_args()

    if args.end_date < args.start_date:
        parser.error("Ending date must be after starting date")

    try:
        cds_uid = os.environ["CDS_UID"]
        cds_api_key = os.environ["CDS_API_KEY"]
    except KeyError:
        parser.exit(
            1,
            "For CDS authentication, both 'CDS_UID' and "
            "'CDS_API_KEY' environment variables must be set",
        )

    cds_key = ":".join([cds_uid, cds_api_key])
    c = cdsapi.Client(url=CDS_URL, key=cds_key, verify=True)

    cur_day: date = args.start_date
    end_day: date = args.end_date
    while cur_day <= end_day:
        download_day(c, cur_day, args.out_dir)
        cur_day += timedelta(days=1)


def download_day(client: cdsapi.Client, cur_day: date, out_dir: Path) -> None:
    """Download the ERA5 datasets of interest for a single day."""
    out_surface = out_dir / Path(f"era5_surface_{cur_day.isoformat()}.nc")
    out_levels = out_dir / Path(f"era5_levels_{cur_day.isoformat()}.nc")

    if not out_surface.exists():
        client.retrieve(
            "reanalysis-era5-single-levels",
            {
                "product_type": "reanalysis",
                "format": "netcdf",
                "variable": SURFACE_VARIABLES,
                "time": HOURS,
                "year": cur_day.strftime("%Y"),
                "month": cur_day.strftime("%m"),
                "day": cur_day.strftime("%d"),
            },
            str(out_surface),
        )
        print(f"Data downloaded: {out_surface}")
    else:
        print(f"File already exists: {out_surface}")

    if not out_levels.exists():
        client.retrieve(
            "reanalysis-era5-pressure-levels",
            {
                "product_type": "reanalysis",
                "format": "netcdf",
                "variable": ATM_VARIABLES,
                "pressure_level": LEVELS,
                "time": HOURS,
                "year": cur_day.strftime("%Y"),
                "month": cur_day.strftime("%m"),
                "day": cur_day.strftime("%d"),
            },
            str(out_levels),
        )
        print(f"Data downloaded: {out_levels}")
    else:
        print(f"File already exists: {out_levels}")


if __name__ == "__main__":
    main()
