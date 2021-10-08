"""Read ERA5 daily surface and profile data."""

from pathlib import Path
from typing import NamedTuple, SupportsIndex, Union, cast

import numpy as np
from numpy.typing import NDArray
from netCDF4 import Dataset


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

    # Profile air temperature in kelvin, dimensioned as (time, lats, lons, levels)
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


def buck_vap(temperature: NDArray[np.float32]) -> NDArray[np.float32]:
    """Buck equation.

    Use the Buck equation to convert temperature in kelvin into water vapor
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


def read_era5_data(
    surface_file: Path,
    levels_file: Path,
    time_subset: Union[slice, SupportsIndex] = slice(None),
    verbose: bool = False,
) -> Era5DailyData:
    """Read the pair of ERA5 surface/levels files.

    Optionally, a subset of the time values can be read.
    """
    if verbose:
        print(f"Reading surface data: {surface_file}")
    with Dataset(surface_file, "r") as f:
        lats = f["latitude"][:]
        lons = f["longitude"][:]
        time = f["time"][time_subset]
        surface_pressure = f["sp"][time_subset, :, :]
        surface_temperature = f["t2m"][time_subset, :, :]
        surface_dewpoint = f["d2m"][time_subset, :, :]
        surface_height = f["z"][time_subset, :, :]
        columnar_water_vapor = f["tcwv"][time_subset, :, :]
        columnar_cloud_liquid = f["tclw"][time_subset, :, :]

    if verbose:
        print(f"Reading profiles data: {levels_file}")
    with Dataset(levels_file, "r") as f:
        # TODO: Probably the lats/lons/time coordinate variables should be
        # checked to ensure the levels file matches the surface file...but for
        # now we'll just be really trusting
        levels = f["level"][:]
        temperature = f["t"][time_subset, :, :, :]
        relative_humidity = f["r"][time_subset, :, :, :]
        height = f["z"][time_subset, :, :, :]
        liquid_content = f["clwc"][time_subset, :, :, :]

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
