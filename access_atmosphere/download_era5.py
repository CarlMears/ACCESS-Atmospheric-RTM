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


class Era5Downloader:
    """Download ERA5 data from CDS."""

    def __init__(
        self, uid: str, api_key: str, out_dir: Path, *, url: str = CDS_URL
    ) -> None:
        """Initialize the downloader with the CDS credentials and output directory."""
        self.out_dir = out_dir
        cds_key = ":".join([uid, api_key])
        self.client = cdsapi.Client(url=url, key=cds_key, verify=True)

    def download_day(self, day: date, verbose: bool = False) -> None:
        """Download the ERA5 datasets of interest for a single day.

        For a given day, two netCDF files are downloaded, one for the surface
        and one for profile levels. If the files already exist in the output
        directory, they are not re-downloaded.
        """
        out_surface = self.out_dir / Path(f"era5_surface_{day.isoformat()}.nc")
        out_levels = self.out_dir / Path(f"era5_levels_{day.isoformat()}.nc")

        if not out_surface.exists():
            self.client.retrieve(
                "reanalysis-era5-single-levels",
                {
                    "product_type": "reanalysis",
                    "format": "netcdf",
                    "variable": SURFACE_VARIABLES,
                    "time": HOURS,
                    "year": day.strftime("%Y"),
                    "month": day.strftime("%m"),
                    "day": day.strftime("%d"),
                },
                str(out_surface),
            )
            if verbose:
                print(f"Data downloaded: {out_surface}")
        else:
            if verbose:
                print(f"File already exists: {out_surface}")

        if not out_levels.exists():
            self.client.retrieve(
                "reanalysis-era5-pressure-levels",
                {
                    "product_type": "reanalysis",
                    "format": "netcdf",
                    "variable": ATM_VARIABLES,
                    "pressure_level": LEVELS,
                    "time": HOURS,
                    "year": day.strftime("%Y"),
                    "month": day.strftime("%m"),
                    "day": day.strftime("%d"),
                },
                str(out_levels),
            )
            if verbose:
                print(f"Data downloaded: {out_levels}")
        else:
            if verbose:
                print(f"File already exists: {out_levels}")


if __name__ == "__main__":
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

    downloader = Era5Downloader(cds_uid, cds_api_key, args.out_dir)
    cur_day: date = args.start_date
    end_day: date = args.end_date
    while cur_day <= end_day:
        downloader.download_day(cur_day)
        cur_day += timedelta(days=1)
