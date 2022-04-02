import numpy as np

from corelib.analytical_models.thermal_geometrical import TG_intensity as TGI


def sikora_intensity(
    pressure, temperature, nozzle_radius, detection_rad, xs, a, s, rs, rf
):
    # compute sikora intensity as defined in
    # Neutral Helium Microscopy (SHeM): A Review
    # https://arxiv.org/abs/2111.12582
    # s is the speed ratio
    # rs is the skimmer radius
    # rf is the quitting surface position
    ITG = TGI(pressure, temperature, nozzle_radius, detection_rad, xs, a)
    intensity = ITG * (
        1 - np.exp(-(s ** 2) * (rs * (rf + a) / (rf * (rf - xs + a))) ** 2)
    )
    return intensity
