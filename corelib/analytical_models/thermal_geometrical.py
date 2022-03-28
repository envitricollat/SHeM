import numpy as np

from corelib.analytical_models.isentropic_nozzle import isentropic_intensity


def TG_intensity(pressure, temperature, nozzle_radius, detection_rad, xs, a):
    # compute the thermal geometrical intensity as defined in
    # Neutral Helium Microscopy (SHeM): A Review
    # https://arxiv.org/abs/2111.12582
    ITG = (
        np.pi
        * isentropic_intensity(pressure, temperature, nozzle_radius)
        * (detection_rad / (xs + a)) ** 2
    )
    return ITG
