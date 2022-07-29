from typing import final

import numpy as np
from numpy.typing import NDArray

@final
class AtmoParameters:
    """Atmospheric radiative parameters.

    This is the output of the RTM.
    """

    @property
    def tran(self) -> NDArray[np.float32]:
        """Transmissivity, from 0 to 1.

        Dimensioned as (`num_points`, `num_freq`).
        """
        ...
    @property
    def tb_up(self) -> NDArray[np.float32]:
        """Upwelling TB, in K.

        Dimensioned as (`num_points`, `num_freq`).
        """
        ...
    @property
    def tb_down(self) -> NDArray[np.float32]:
        """Downwelling TB, in K.

        Dimensioned as (`num_points`, `num_freq`).
        """
        ...

def compute_rtm(
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
    num_levels: int,
    num_points: int,
    num_freq: int,
) -> AtmoParameters: ...
