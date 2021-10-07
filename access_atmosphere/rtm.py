"""Perform the radiative transfer model calculations.

This is implemented as a Python extension, written in Fortran.
"""

from typing import NamedTuple

import numpy as np
from numpy.typing import NDArray

from .rtm_f90 import access_rtm


class AtmoParameters(NamedTuple):
    """Atmospheric radiative parameters.

    This is the output of the RTM.
    """

    # Transmissivity, from 0 to 1
    tran: NDArray[np.float32]
    # Upwelling TB, in K
    tb_up: NDArray[np.float32]
    # Downwelling TB, in K
    tb_down: NDArray[np.float32]


def compute(
    pressure: NDArray[np.float32],
    temperature: NDArray[np.float32],
    height: NDArray[np.float32],
    relative_humidity: NDArray[np.float32],
    liquid_content: NDArray[np.float32],
    surface_temperature: NDArray[np.float32],
    surface_height: NDArray[np.float32],
    surface_relative_humidity: NDArray[np.float32],
    surface_pressure: NDArray[np.float32],
    incidence_angle: NDArray[np.float32],
    frequency: NDArray[np.float32],
) -> AtmoParameters:
    """Compute the radiative transfer model for the atmosphere.

    All inputs are numpy arrays and are 1d or 2d. The "pressure" parameter is
    the pressure levels in hPa and has shape (num_levels, ). It is treated as a
    constant (i.e., not a function of num_points).

    p: pressure levels, in hPa

    The following are input profiles and have shape (num_points, num_levels):

    temperature: physical temperature in K

    height: geometric height above the geoid in m

    relative_humidity: relative humidity in %

    liquid_content: liquid water content (from clouds) in kg/kg

    The following are surface parameters and have shape (num_points, ):

    surface_temperature: 2 meter air temperature in K

    surface_height: geopotential height at the surface in m

    surface_relative_humidity: 2 meter relative humidity in %

    surface_pressure: surface pressure in hPa

    The following are RTM parameters and have shape (num_freq, ):

    incidence_angle: Earth incidence angle in degrees

    frequency: microwave frequency in GHz

    The returned atmospheric parameters are dimensioned as (num_points,
    num_freq).
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
            for a in (temperature, height, relative_humidity, liquid_content)
        )
        or not all(
            a.shape == (num_points,)
            for a in (
                surface_temperature,
                surface_height,
                surface_relative_humidity,
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

    # Why are the 2d arrays transposed? Because f2py doesn't transpose the
    # arrays. It should be the case that the un-transposed values are
    # C-contiguous but the transposed values are F-contiguous and therefore the
    # best option for the Fortran code.
    tran, tb_up, tb_down = access_rtm.compute_rtm(
        pressure,
        temperature.T,
        height.T,
        relative_humidity.T,
        liquid_content.T,
        surface_temperature,
        surface_height,
        surface_relative_humidity,
        surface_pressure,
        incidence_angle,
        frequency,
        num_points,
        num_freq,
    )

    # Similarly, the output arrays need to be transposed to get them back into C
    # order
    return AtmoParameters(tran.T, tb_up.T, tb_down.T)
