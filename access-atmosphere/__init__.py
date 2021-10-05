"""Atmospheric RTM for the ACCESS project."""


from typing import NamedTuple


class AtmoParameters(NamedTuple):
    """Atmospheric radiative parameters.

    This is the output of the RTM.
    """

    # Transmissivity, from 0 to 1
    tran: float
    # Upwelling TB, in K
    tb_up: float
    # Downwelling TB, in K
    tb_down: float


class AtmoRtm:
    """The atmospheric RTM."""

    @classmethod
    def compute(p, t, pv, rhol, z, ibegin, eia, freq) -> AtmoParameters:
        """Compute the radiative transfer model for the atmosphere.

        The following are input profiles and are arrays of the same length:

        p: pressure levels, in hPa

        t: temperature, in K

        pv: water vapor pressure, in hPa

        rhol: liquid cloud density, in g/m^3

        z: geometric height above the geoid, in m

        The following are scalar inputs:

        ibegin: the index in the profiles to begin at

        eia: Earth incidence angle, in degrees

        freq: Microwave frequency, in GHz

        """
        # TODO
        return AtmoParameters(0.5, 200.0, 100.0)
