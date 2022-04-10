# a simple pyhton script that calculates the intensity
# in a microscope using an analytical intensity equation according
# to certain parameters
# pressure in pascals
import itertools

import numpy as np
import pandas as pd

from corelib.analytical_models.sikora import sikora_intensity
from corelib.analytical_models.thermal_geometrical import TG_intensity

pressure = 100e5
# temperature in kelvins
temperature = 300
# nozzle radius in meters
nozzle_radius = 5e-6
# detection radius in meters
det_rad = 100e-6
# distance to the skimmer
xs = 492e-3
# distance to the detector
a = 1
# skimmer radius
rs = 0.5e-3
# detector efficiency
efficiency = 2.1e-6

# intensity in counts/s
Intensity = TG_intensity(pressure, temperature, nozzle_radius, det_rad, xs, a)
# multiply by detector efficiency
Intensity = Intensity * efficiency
# given a list of speed ratios and pressures, get the sikora matrix
S_r = [125, 300]
T = 295
P = np.arange(1, 100, 1)
intensities = pd.DataFrame(index=P, columns=[str(element) for element in S_r])
for el in itertools.product(S_r, P):
    intensities.loc[el[1], str(el[0])] = sikora_intensity(
        el[1], temperature, nozzle_radius, det_rad, xs, a, el[0], rs, xs
    )
print(Intensity)
print(intensities)
