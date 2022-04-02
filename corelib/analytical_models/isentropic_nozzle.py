import numpy as np
from scipy.constants import k as k_b

# [kg] He mass
m = 4 * 1.66053886e-27

gamma = 5 / 3
F_g = np.sqrt(gamma / (gamma + 1)) * (2 / (gamma + 1)) ** (1 / (gamma - 1))


def isentropic_intensity(pressure, temperature, nozzle_radius):
    # compute the isentropic intensity as defined in
    # Neutral Helium Microscopy (SHeM): A Review
    # https://arxiv.org/abs/2111.12582
    Intensity = (
        (pressure / (k_b * temperature))
        * np.sqrt(2 * k_b * temperature / m)
        * np.pi
        * (nozzle_radius ** 2)
        * F_g
    )
    return Intensity
