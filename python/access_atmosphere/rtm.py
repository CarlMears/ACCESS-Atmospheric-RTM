"""Perform the radiative transfer model calculations.

This is implemented as a Python extension, written in Rust.
"""

from typing import Optional

import numpy as np
from numpy.typing import NDArray

from .access_atmosphere import AtmoParameters, compute_rtm


def compute(
    pressure: NDArray[np.float32],
    temperature: NDArray[np.float32],
    height: NDArray[np.float32],
    specific_humidity: NDArray[np.float32],
    liquid_content: NDArray[np.float32],
    surface_temperature: NDArray[np.float32],
    surface_height: NDArray[np.float32],
    surface_dewpoint: NDArray[np.float32],
    surface_pressure: NDArray[np.float32],
    incidence_angle: NDArray[np.float32],
    frequency: NDArray[np.float32],
    num_threads: Optional[int] = None,
) -> AtmoParameters:
    """Compute the radiative transfer model for the atmosphere.

    This is a wrapper function to ensure all the inputs have consistent and
    expected shapes. All inputs are numpy arrays and are either 1d or 2d. The
    `pressure` parameter is the pressure levels in hPa and has shape
    (`num_levels`, ). It is treated as a constant (i.e., not a function of
    `num_points`).

    `pressure`: pressure levels, in hPa

    The following are input profiles and have shape (`num_points`,
    `num_levels`):

    `temperature`: physical temperature in K

    `height`: geometric height above the geoid in m

    `specific_humidity`: specific humidity in kg/kg

    `liquid_content`: liquid water content (from clouds) in kg/kg

    The following are surface parameters and have shape (`num_points`, ):

    `surface_temperature`: 2 meter air temperature in K

    `surface_height`: geopotential height at the surface in m

    `surface_dewpoint`: 2 meter dewpoint in K

    `surface_pressure`: surface pressure in hPa

    The following are RTM parameters and have shape (`num_freq`, ):

    `incidence_angle`: Earth incidence angle in degrees

    `frequency`: microwave frequency in GHz

    The returned atmospheric parameters are each dimensioned as (`num_points`,
    `num_freq`).

    The number of worker threads is controlled by `num_threads`. It must be a
    positive integer, or `None` to automatically choose the number of threads.
    """
    # Check all inputs
    if pressure.ndim != 1 or temperature.ndim != 2 or incidence_angle.ndim != 1:
        raise Exception("Unexpected input shapes")

    num_levels = len(pressure)
    num_points = temperature.shape[0]
    num_freq = len(incidence_angle)

    if (
        not all(
            a.shape == (num_points, num_levels)
            for a in (temperature, height, specific_humidity, liquid_content)
        )
        or not all(
            a.shape == (num_points,)
            for a in (
                surface_temperature,
                surface_height,
                surface_dewpoint,
                surface_pressure,
            )
        )
        or not all(
            a.shape == (num_freq,)
            for a in (
                incidence_angle,
                frequency,
            )
        )
    ):
        raise Exception("Unexpected input shapes")

    if num_threads is not None and num_threads < 0:
        raise Exception("Number of threads must be non-negative")

    return compute_rtm(
        num_freq,
        num_levels,
        num_points,
        pressure,
        temperature,
        height,
        specific_humidity,
        liquid_content,
        surface_temperature,
        surface_height,
        surface_dewpoint,
        surface_pressure,
        incidence_angle,
        frequency,
        num_threads,
    )
