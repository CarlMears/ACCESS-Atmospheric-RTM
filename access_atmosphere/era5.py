"""Read ERA5 daily surface and profile data."""

from pathlib import Path
from typing import NamedTuple, Optional, Sequence, Union, cast

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

    # Profile specific humidity in kg/kg, dimensioned as (time, lats, lons,
    # levels)
    specific_humidity: NDArray[np.float32]

    # Geopotential height profile in meters, dimensioned as (time, lats, lons, levels)
    height: NDArray[np.float32]

    # Profile specific liquid water content (from clouds) in kg/kg, dimensioned
    # as (time, lats, lons, levels)
    liquid_content: NDArray[np.float32]

    # Surface pressure in hPa, dimensioned as (time, lats, lons)
    surface_pressure: NDArray[np.float32]

    # 2-meter air temperature in kelvin, dimensioned as (time, lats, lons)
    surface_temperature: NDArray[np.float32]

    # 2-meter dewpoint in kelvin, dimensioned as (time, lats, lons)
    surface_dewpoint: NDArray[np.float32]

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
    time_subset: Optional[Sequence[int]] = None,
    verbose: bool = False,
) -> Era5DailyData:
    """Read the pair of ERA5 surface/levels files.

    Optionally, a subset of the time values can be read.
    """
    if verbose:
        print(f"Reading surface data: {surface_file}")
        if time_subset is not None:
            print(f"Subsetting hours to: {time_subset}")

    times: Union[slice, Sequence[int]]
    if time_subset is None:
        times = slice(None)
    else:
        times = time_subset
    with Dataset(surface_file, "r") as f:
        lats = f["latitude"][:]
        lons = f["longitude"][:]
        time = f["time"][times]
        surface_pressure = f["sp"][times, :, :]
        surface_temperature = f["t2m"][times, :, :]
        surface_dewpoint = f["d2m"][times, :, :]
        surface_height = f["z"][times, :, :]
        columnar_water_vapor = f["tcwv"][times, :, :]
        columnar_cloud_liquid = f["tclw"][times, :, :]

    if verbose:
        print(f"Reading profiles data: {levels_file}")
    with Dataset(levels_file, "r") as f:
        # TODO: Probably the lats/lons/time coordinate variables should be
        # checked to ensure the levels file matches the surface file...but for
        # now we'll just be really trusting
        levels = f["level"][:]
        temperature = f["t"][times, :, :, :]
        specific_humidity = f["q"][times, :, :, :]
        height = f["z"][times, :, :, :]
        liquid_content = f["clwc"][times, :, :, :]

    if verbose:
        print(f"Post-processing ERA5 data ({len(levels)} pressure levels)")

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
            specific_humidity,
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
    specific_humidity = np.ma.getdata(specific_humidity).astype(np.float32)
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

    # The 4d arrays should be reordered from (time, levels, lat, lon) to (time,
    # lat, lon, levels)
    temperature = np.moveaxis(temperature, 1, -1)
    specific_humidity = np.moveaxis(specific_humidity, 1, -1)
    height = np.moveaxis(height, 1, -1)
    liquid_content = np.moveaxis(liquid_content, 1, -1)

    # The latitudes need to be adjusted. In the ERA5 files, the latitudes are in
    # *descending* order from 90 to -90, and the longitudes are in ascending
    # order from 0 to 360. Leave the longitudes alone, but flip the latitudes so
    # they go from -90 to 90.
    lats = -lats

    temperature = np.flip(temperature, 1)
    specific_humidity = np.flip(specific_humidity, 1)
    height = np.flip(height, 1)
    liquid_content = np.flip(liquid_content, 1)
    surface_pressure = np.flip(surface_pressure, 1)
    surface_temperature = np.flip(surface_temperature, 1)
    surface_dewpoint = np.flip(surface_dewpoint, 1)
    surface_height = np.flip(surface_height, 1)
    columnar_water_vapor = np.flip(columnar_water_vapor, 1)
    columnar_cloud_liquid = np.flip(columnar_cloud_liquid, 1)

    return Era5DailyData(
        levels,
        lats,
        lons,
        time,
        temperature,
        specific_humidity,
        height,
        liquid_content,
        surface_pressure,
        surface_temperature,
        surface_dewpoint,
        surface_height,
        columnar_water_vapor,
        columnar_cloud_liquid,
    )
