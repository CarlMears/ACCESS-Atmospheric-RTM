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
from typing import NamedTuple, cast

import numpy as np
from netCDF4 import Dataset, getlibversion, num2date
from numpy.typing import NDArray

from . import rtm

# Reference frequencies (in GHz) to use
REF_FREQ = np.array([1.41, 6.8, 10.7, 18.7, 23.8, 37.0, 89.0], np.float32)

# Reference Earth incidence angle to use for each reference frequency (in degrees)
REF_EIA = np.array([40.0, 53.0, 53.0, 53.0, 53.0, 53.0, 53.0], np.float32)


def buck_vap(temperature: NDArray[np.float32]) -> NDArray[np.float32]:
    """Buck equation.

    Use the Buck equation to convert temperature in Kelvin into water vapor
    saturation pressure in hPa. The equation is from [1], which cites Buck 1996.

    To convert to water vapor partial pressure, multiply the result by the
    relative humidity.

    [1] https://en.wikipedia.org/wiki/Arden_Buck_equation
    """
    # Temperature in degrees Celsius
    temp_c: NDArray[np.float32] = temperature - 273.15
    return cast(
        NDArray[np.float32],
        6.1121 * np.exp((18.678 - temp_c / 234.5) * (temp_c / (257.14 + temp_c))),
    )


class Era5DailyData(NamedTuple):
    """The daily data for ERA5 surface and levels data."""

    # Pressure levels in hPa, with shape (num_levels, ). They should be in
    # descending order (e.g., 1000 to 10).
    levels: NDArray[np.float32]

    # Latitude in degrees North, dimensioned as (num_lats, ). They should be in
    # ascending order (e.g., -90 to 90).
    lats: NDArray[np.float32]

    # Longitude in degrees East, dimensioned as (num_lons, ). They should be in
    # ascending order (e.g., -180 to 180).
    lons: NDArray[np.float32]

    # Hours since 1900-01-01, dimensioned as (num_time, ).
    time: NDArray[np.int32]

    # Profile air temperature in Kelvin, dimensioned as (time, lats, lons, levels)
    temperature: NDArray[np.float32]

    # Profile relative humidity in percentage (0 to 100), dimensioned as (time,
    # lats, lons, levels)
    relative_humidity: NDArray[np.float32]

    # Geopotential height profile in meters, dimensioned as (time, lats, lons, levels)
    height: NDArray[np.float32]

    # Profile specific liquid water content (from clouds) in kg/kg, dimensioned
    # as (time, lats, lons, levels)
    liquid_content: NDArray[np.float32]

    # Surface pressure in hPa, dimensioned as (time, lats, lons)
    surface_pressure: NDArray[np.float32]

    # 2-meter air temperature in kelvin, dimensioned as (time, lats, lons)
    surface_temperature: NDArray[np.float32]

    # 2-meter relative humidity in percentage, dimensioned as (time, lats, lons)
    surface_relative_humidity: NDArray[np.float32]

    # Geopotential height at the surface in meters, dimensioned as (time, lats, lons)
    surface_height: NDArray[np.float32]

    # Total column water vapor in kg/m^2, dimensioned as (time, lats, lons)
    columnar_water_vapor: NDArray[np.float32]

    # Total column cloud liquid water in kg/m^2, dimensioned as (time, lats, lons)
    columnar_cloud_liquid: NDArray[np.float32]


def read_era5_data(
    surface_file: Path, levels_file: Path, verbose: bool = False
) -> Era5DailyData:
    """Read the pair of ERA5 surface/levels files."""
    if verbose:
        print(f"Reading surface data: {surface_file}")
    with Dataset(surface_file, "r") as f:
        lats = f["latitude"][:]
        lons = f["longitude"][:]
        time = f["time"][:]
        surface_pressure = f["sp"][:]
        surface_temperature = f["t2m"][:]
        surface_dewpoint = f["d2m"][:]
        surface_height = f["z"][:]
        columnar_water_vapor = f["tcwv"][:]
        columnar_cloud_liquid = f["tclw"][:]

    if verbose:
        print(f"Reading profiles data: {levels_file}")
    with Dataset(levels_file, "r") as f:
        # TODO: Probably the lats/lons/time coordinate variables should be
        # checked to ensure the levels file matches the surface file...but for
        # now we'll just be really trusting
        levels = f["level"][:]
        temperature = f["t"][:]
        relative_humidity = f["r"][:]
        height = f["z"][:]
        liquid_content = f["clwc"][:]

    if verbose:
        print("Post-processing ERA5 data")

    # By default, netCDF4 returns masked arrays for all the variables above.
    # However, there shouldn't be any values that are actually masked in the
    # ERA5 data. Check that assumption and then convert everything to "vanilla"
    # ndarrays and also cast the float64 values to float32 values.
    if any(
        np.ma.count_masked(a) > 0
        for a in (
            lats,
            lons,
            time,
            surface_pressure,
            surface_temperature,
            surface_dewpoint,
            surface_height,
            columnar_water_vapor,
            columnar_cloud_liquid,
            levels,
            temperature,
            relative_humidity,
            height,
            liquid_content,
        )
    ):
        raise Exception("Masked input values detected")
    lats = np.ma.getdata(lats)
    lons = np.ma.getdata(lons)
    time = np.ma.getdata(time)
    surface_pressure = np.ma.getdata(surface_pressure).astype(np.float32)
    surface_temperature = np.ma.getdata(surface_temperature).astype(np.float32)
    surface_dewpoint = np.ma.getdata(surface_dewpoint).astype(np.float32)
    surface_height = np.ma.getdata(surface_height).astype(np.float32)
    columnar_water_vapor = np.ma.getdata(columnar_water_vapor).astype(np.float32)
    columnar_cloud_liquid = np.ma.getdata(columnar_cloud_liquid).astype(np.float32)
    levels = np.ma.getdata(levels).astype(np.float32)
    temperature = np.ma.getdata(temperature).astype(np.float32)
    relative_humidity = np.ma.getdata(relative_humidity).astype(np.float32)
    height = np.ma.getdata(height).astype(np.float32)
    liquid_content = np.ma.getdata(liquid_content).astype(np.float32)

    # The reciprocal of the standard gravity, in units of s^2 / m
    # https://en.wikipedia.org/wiki/Standard_gravity
    INV_STANDARD_GRAVITY = 1 / 9.80665

    # Convert geopotential to geopotential height
    # (https://apps.ecmwf.int/codes/grib/param-db?id=129)
    height *= INV_STANDARD_GRAVITY
    surface_height *= INV_STANDARD_GRAVITY

    # Convert surface pressure from Pa to hPa
    surface_pressure *= 1e-2

    # Convert dewpoint at 2 m to relative humidity at 2 m
    # http://bmcnoldy.rsmas.miami.edu/Humidity.html
    # However, rather than using the "Magnus approximation" I use the Buck equation.
    surface_relative_humidity = cast(
        NDArray[np.float32],
        100.0 * buck_vap(surface_dewpoint) / buck_vap(surface_temperature),
    )

    # The 4d arrays should be reordered from (time, levels, lat, lon) to (time,
    # lat, lon, levels)
    temperature = np.moveaxis(temperature, 1, -1)
    relative_humidity = np.moveaxis(relative_humidity, 1, -1)
    height = np.moveaxis(height, 1, -1)
    liquid_content = np.moveaxis(liquid_content, 1, -1)

    # The latitudes/longitudes need to be adjusted. In the ERA5 files, the
    # latitudes are in *descending* order from 90 to -90, and while the longitudes
    # are in ascending order, they go from 0 to 360. The desired output is that
    # the latitudes go from -90 to 90 and longitudes from -180 to 180.
    lats = -lats
    lons = lons - 180.0
    half_lon_len = len(lons) // 2

    def flip_and_roll_4d(a: NDArray[np.float32]) -> NDArray[np.float32]:
        return np.roll(a[:, ::-1, :, :], half_lon_len, axis=2)

    def flip_and_roll_3d(a: NDArray[np.float32]) -> NDArray[np.float32]:
        return np.roll(a[:, ::-1, :], half_lon_len, axis=2)

    temperature = flip_and_roll_4d(temperature)
    relative_humidity = flip_and_roll_4d(relative_humidity)
    height = flip_and_roll_4d(height)
    liquid_content = flip_and_roll_4d(liquid_content)
    surface_pressure = flip_and_roll_3d(surface_pressure)
    surface_temperature = flip_and_roll_3d(surface_temperature)
    surface_relative_humidity = flip_and_roll_3d(surface_relative_humidity)
    surface_height = flip_and_roll_3d(surface_height)
    columnar_water_vapor = flip_and_roll_3d(columnar_water_vapor)
    columnar_cloud_liquid = flip_and_roll_3d(columnar_cloud_liquid)

    return Era5DailyData(
        levels,
        lats,
        lons,
        time,
        temperature,
        relative_humidity,
        height,
        liquid_content,
        surface_pressure,
        surface_temperature,
        surface_relative_humidity,
        surface_height,
        columnar_water_vapor,
        columnar_cloud_liquid,
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
    # ascending order (e.g., -180 to 180).
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

    time_units = "hours since 1990-01-01 00:00:00Z"
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
        f.setncattr("geospatial_lon_min", np.float32(-180.0))
        f.setncattr("geospatial_lon_max", np.float32(180.0))
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


def convert(
    era5_surface_input: Path,
    era5_levels_input: Path,
    rtm_output: Path,
    verbose: bool = False,
) -> None:
    """Read ERA5 profile/surface files and run the RTM and write its output."""
    era5_data = read_era5_data(era5_surface_input, era5_levels_input, verbose)

    if verbose:
        print("Running RTM over all data")
    # The ERA5 data is organized by lat/lon, but we need to vectorize that down
    # for the RTM and then reshape the output when finished
    shape_4d = era5_data.temperature.shape
    num_time, num_lat, num_lon, num_levels = shape_4d[0:4]
    num_points = num_time * num_lat * num_lon
    atmo_results = rtm.compute(
        era5_data.levels,
        np.reshape(era5_data.temperature, (num_points, num_levels)),
        np.reshape(era5_data.height, (num_points, num_levels)),
        np.reshape(era5_data.relative_humidity, (num_points, num_levels)),
        np.reshape(era5_data.liquid_content, (num_points, num_levels)),
        np.ravel(era5_data.surface_temperature),
        np.ravel(era5_data.surface_height),
        np.ravel(era5_data.surface_relative_humidity),
        np.ravel(era5_data.surface_pressure),
        REF_EIA,
        REF_FREQ,
    )

    # Now the output values need to be un-vectorized
    num_freq = len(REF_FREQ)
    tran = np.reshape(atmo_results.tran, (num_time, num_lat, num_lon, num_freq))
    tb_up = np.reshape(atmo_results.tb_up, (num_time, num_lat, num_lon, num_freq))
    tb_down = np.reshape(atmo_results.tb_down, (num_time, num_lat, num_lon, num_freq))

    rtm_data = RtmDailyData(
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
    args = parser.parse_args()

    print(f"ERA5 surface file: {args.era5_surface}")
    print(f"ERA5 levels file: {args.era5_levels}")
    print(f"RTM output file: {args.rtm_out}")
    convert(args.era5_surface, args.era5_levels, args.rtm_out, verbose=True)
