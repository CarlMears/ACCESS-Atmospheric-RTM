"""Process an entire daily file.

Example usage:

python -m access_atmosphere.process \
    era5_surface_2020-01-01.nc \
    era5_levels_2020-01-01.nc \
    access_era5_2020-01-01.nc
"""

import argparse
import sys
from datetime import datetime, timezone
from pathlib import Path
from time import perf_counter_ns
from typing import NamedTuple, Optional, Sequence

import numpy as np
from netCDF4 import Dataset, getlibversion, num2date
from numpy.typing import NDArray

from . import era5, rtm

# Reference frequencies (in GHz) to use
REF_FREQ: NDArray[np.float32] = np.array(
    [1.41, 6.8, 10.7, 18.7, 23.8, 37.0, 89.0], np.float32
)

# Reference Earth incidence angle to use for each reference frequency (in degrees)
REF_EIA: NDArray[np.float32] = np.array(
    [40.0, 53.0, 53.0, 53.0, 53.0, 53.0, 53.0], np.float32
)


class RtmDailyData(NamedTuple):
    """Output values after computing the atmospheric RTM.

    This is for a full day of data, so it is more "full" than should normally be
    computed.
    """

    # Latitude in degrees North, dimensioned as (num_lats, ). They should be in
    # ascending order (e.g., -90 to 90).
    lats: NDArray[np.float32]

    # Longitude in degrees East, dimensioned as (num_lons, ). They should be in
    # ascending order (e.g., 0 to 360).
    lons: NDArray[np.float32]

    # Hours since 1900-01-01, dimensioned as (num_time, ).
    time: NDArray[np.int32]

    # Microwave frequency in GHz, dimensioned as (num_freq, ).
    frequency: NDArray[np.float32]

    # Earth incidence angle in degrees, dimensioned as (num_freq, ).
    incidence: NDArray[np.float32]

    # Total column water vapor in kg/m^2, dimensioned as (time, lats, lons)
    columnar_water_vapor: NDArray[np.float32]

    # Total column cloud liquid water in kg/m^2, dimensioned as (time, lats, lons)
    columnar_cloud_liquid: NDArray[np.float32]

    # Atmospheric transmissivity, unitless, dimensioned as (time, lats, lons,
    # freq).
    transmissivity: NDArray[np.float32]

    # Atmospheric upwelling brightness temperature in K, dimensioned as (time,
    # lats, lons, freq).
    tb_up: NDArray[np.float32]

    # Atmospheric downwelling brightness temperature in K, dimensioned as (time,
    # lats, lons, freq).
    tb_down: NDArray[np.float32]


def write_rtm_data(rtm_data: RtmDailyData, rtm_output: Path) -> None:
    """Write the daily RTM data to a netCDF output file."""
    timestamp = datetime.now(timezone.utc).isoformat(" ", "seconds")
    nc_version = getlibversion().partition(" ")[0]

    time_units = "hours since 1900-01-01 00:00:00Z"
    data_times = num2date(rtm_data.time, time_units)
    time_start = min(data_times).isoformat(" ", "seconds")
    time_end = max(data_times).isoformat(" ", "seconds")

    with Dataset(rtm_output, "w") as f:
        # ----------
        # Global attributes
        f.setncattr_string("Conventions", "CF-1.9,ACDD-1.3")
        f.setncattr_string("title", "ACCESS RTM output")
        f.setncattr_string("institution", "REMSS")
        f.setncattr_string("history", f"{timestamp} created: {' '.join(sys.argv[1:])}")
        f.setncattr_string("netcdf_version_id", nc_version)
        f.setncattr_string("date_created", timestamp)
        f.setncattr_string("creator_name", "Remote Sensing Systems")
        f.setncattr_string("creator_email", "support@remss.com")
        f.setncattr_string("creator_url", "http://www.remss.com")
        f.setncattr("geospatial_lat_min", np.float32(-90.0))
        f.setncattr("geospatial_lat_max", np.float32(90.0))
        f.setncattr("geospatial_lon_min", np.float32(0.0))
        f.setncattr("geospatial_lon_max", np.float32(360.0))
        f.setncattr_string("time_coverage_start", time_start)
        f.setncattr_string("time_coverage_end", time_end)
        f.setncattr_string("standard_name_vocabulary", "CF Standard Name Table v78")

        # ----------
        # Dimensions
        f.createDimension("time", len(rtm_data.time))
        f.createDimension("lat", len(rtm_data.lats))
        f.createDimension("lon", len(rtm_data.lons))
        f.createDimension("freq", len(rtm_data.frequency))

        # ----------
        # Coordinate variables
        v = f.createVariable("time", np.int32, ("time",))
        v[:] = rtm_data.time
        v.setncattr_string("standard_name", "time")
        v.setncattr_string("axis", "T")
        v.setncattr_string("units", time_units)

        v = f.createVariable("lat", np.float32, ("lat",))
        v[:] = rtm_data.lats
        v.setncattr_string("standard_name", "latitude")
        v.setncattr_string("axis", "Y")
        v.setncattr_string("units", "degrees_north")

        v = f.createVariable("lon", np.float32, ("lon",))
        v[:] = rtm_data.lons
        v.setncattr_string("standard_name", "longitude")
        v.setncattr_string("axis", "X")
        v.setncattr_string("units", "degrees_east")

        v = f.createVariable("freq", np.float32, ("freq",))
        v[:] = rtm_data.frequency
        v.setncattr_string("standard_name", "sensor_band_central_radiation_frequency")
        v.setncattr_string("long_name", "frequency")
        v.setncattr_string("units", "GHz")

        v = f.createVariable("eia", np.float32, ("freq",))
        v[:] = rtm_data.incidence
        v.setncattr_string("standard_name", "sensor_zenith_angle")
        v.setncattr_string("long_name", "incidence angle")
        v.setncattr_string("units", "degree")

        # ----------
        # Variables
        v = f.createVariable("col_vapor", np.float32, ("time", "lat", "lon"), zlib=True)
        v[...] = rtm_data.columnar_water_vapor
        v.setncattr_string("standard_name", "atmosphere_mass_content_of_water_vapor")
        v.setncattr_string("long_name", "columnar water vapor")
        v.setncattr_string("units", "kg m-2")
        v.setncattr_string("coordinates", "lat lon")

        v = f.createVariable("col_water", np.float32, ("time", "lat", "lon"), zlib=True)
        v[...] = rtm_data.columnar_cloud_liquid
        v.setncattr_string(
            "standard_name", "atmosphere_mass_content_of_cloud_liquid_water"
        )
        v.setncattr_string("long_name", "columnar liquid cloud content")
        v.setncattr_string("units", "kg m-2")
        v.setncattr_string("coordinates", "lat lon")

        v = f.createVariable(
            "tran", np.float32, ("time", "lat", "lon", "freq"), zlib=True
        )
        v[...] = rtm_data.transmissivity
        v.setncattr_string("long_name", "atmospheric transmissivity")
        v.setncattr_string("coordinates", "lat lon")

        v = f.createVariable(
            "tb_up", np.float32, ("time", "lat", "lon", "freq"), zlib=True
        )
        v[...] = rtm_data.tb_up
        v.setncattr_string("long_name", "upwelling brightness temperature")
        v.setncattr_string("units", "kelvin")
        v.setncattr_string("coordinates", "lat lon")

        v = f.createVariable(
            "tb_down", np.float32, ("time", "lat", "lon", "freq"), zlib=True
        )
        v[...] = rtm_data.tb_down
        v.setncattr_string("long_name", "downwelling brightness temperature")
        v.setncattr_string("units", "kelvin")
        v.setncattr_string("coordinates", "lat lon")


def combine_rtm(hourly: Sequence[RtmDailyData]) -> RtmDailyData:
    """Combine subsets of daily RTM data together."""
    if len(hourly) == 0:
        raise Exception("At least one input is needed to combine")

    # Check that all coordinate variables match except for the "time" variable
    if not all(np.allclose(data.lats, hourly[0].lats) for data in hourly):
        raise Exception("Mismatched latitudes")
    if not all(np.allclose(data.lons, hourly[0].lons) for data in hourly):
        raise Exception("Mismatched longitudes")
    if not all(np.allclose(data.frequency, hourly[0].frequency) for data in hourly):
        raise Exception("Mismatched frequencies")
    if not all(np.allclose(data.incidence, hourly[0].incidence) for data in hourly):
        raise Exception("Mismatched incidence angles")

    full_data = RtmDailyData(
        hourly[0].lats,
        hourly[0].lons,
        np.concatenate([data.time for data in hourly]),
        hourly[0].frequency,
        hourly[0].incidence,
        np.concatenate([data.columnar_water_vapor for data in hourly], axis=0),
        np.concatenate([data.columnar_cloud_liquid for data in hourly], axis=0),
        np.concatenate([data.transmissivity for data in hourly], axis=0),
        np.concatenate([data.tb_up for data in hourly], axis=0),
        np.concatenate([data.tb_down for data in hourly], axis=0),
    )

    return full_data


def run_rtm(era5_data: era5.Era5DailyData, verbose: bool = False) -> RtmDailyData:
    """Run the RTM on ERA5 daily data."""
    if verbose:
        print("Running RTM over all data")

    tick = perf_counter_ns()

    # The ERA5 data is organized by lat/lon, but we need to vectorize that down
    # for the RTM and then reshape the output when finished
    shape_4d = era5_data.temperature.shape
    num_time, num_lat, num_lon, num_levels = shape_4d[0:4]
    num_points = num_time * num_lat * num_lon
    atmo_results = rtm.compute(
        era5_data.levels,
        np.reshape(era5_data.temperature, (num_points, num_levels)),
        np.reshape(era5_data.height, (num_points, num_levels)),
        np.reshape(era5_data.specific_humidity, (num_points, num_levels)),
        np.reshape(era5_data.liquid_content, (num_points, num_levels)),
        np.ravel(era5_data.surface_temperature),
        np.ravel(era5_data.surface_height),
        np.ravel(era5_data.surface_dewpoint),
        np.ravel(era5_data.surface_pressure),
        REF_EIA,
        REF_FREQ,
    )

    tock = perf_counter_ns()
    duration_seconds = (tock - tick) * 1e-9
    if verbose:
        print(f"Finished RTM in {duration_seconds:0.2f} s")

    # Now the output values need to be un-vectorized
    num_freq = len(REF_FREQ)
    tran = np.reshape(atmo_results.tran, (num_time, num_lat, num_lon, num_freq))
    tb_up = np.reshape(atmo_results.tb_up, (num_time, num_lat, num_lon, num_freq))
    tb_down = np.reshape(atmo_results.tb_down, (num_time, num_lat, num_lon, num_freq))

    return RtmDailyData(
        era5_data.lats,
        era5_data.lons,
        era5_data.time,
        REF_FREQ,
        REF_EIA,
        era5_data.columnar_water_vapor,
        era5_data.columnar_cloud_liquid,
        tran,
        tb_up,
        tb_down,
    )


def convert_all(
    era5_surface_input: Path,
    era5_levels_input: Path,
    rtm_output: Path,
    time_indices: Optional[Sequence[int]],
    one_pass: bool,
    verbose: bool = False,
) -> None:
    """Read the ERA5 profile/surface files and run the RTM and write its output."""
    all_time_indices: Sequence[int] | list[int]
    if time_indices is None:
        all_time_indices = era5.read_time_indices(era5_surface_input, era5_levels_input)
    else:
        all_time_indices = time_indices

    if one_pass:
        era5_data = era5.read_era5_data(
            era5_surface_input, era5_levels_input, all_time_indices, verbose=verbose
        )
        rtm_data = run_rtm(era5_data, verbose)
    else:
        # Read the ERA5 data one hour at a time and process just that much
        hourly_rtm_data = []
        for time_index in all_time_indices:
            era5_data = era5.read_era5_data(
                era5_surface_input, era5_levels_input, [time_index], verbose=verbose
            )
            hourly_rtm_data.append(run_rtm(era5_data, verbose))

        # Accumulate all the hourly data together
        rtm_data = combine_rtm(hourly_rtm_data)

    if verbose:
        print(f"Writing output data: {rtm_output}")
    write_rtm_data(rtm_data, rtm_output)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Run the RTM on one day of ERA5 data")
    parser.add_argument("era5_surface", type=Path, help="ERA5 surface data for one day")
    parser.add_argument("era5_levels", type=Path, help="ERA5 levels data for one day")
    parser.add_argument(
        "rtm_out",
        type=Path,
        help="RTM daily output file",
    )
    parser.add_argument(
        "--time",
        action="append",
        type=int,
        help="Time index to use in the ERA5 data, can be specified multiple times",
    )
    parser.add_argument(
        "--one-pass",
        action="store_true",
        help=(
            "Read all ERA5 time indices at once instead of iterating over them. "
            "This reduces the time to complete at the cost of "
            "increased memory consumption."
        ),
    )
    args = parser.parse_args()

    print(f"ERA5 surface file: {args.era5_surface}")
    print(f"ERA5 levels file: {args.era5_levels}")
    print(f"RTM output file: {args.rtm_out}")
    if args.one_pass:
        print("One-pass mode enabled")
    convert_all(
        args.era5_surface,
        args.era5_levels,
        args.rtm_out,
        args.time,
        args.one_pass,
        verbose=True,
    )
