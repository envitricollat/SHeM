import pytest

from corelib.analytical_models.isentropic_nozzle import isentropic_intensity
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
# speed ratio
S_r = 300


@pytest.mark.unit
def testTG_smaller_IO():
    # test that the thermal geometrical intensity is smaller than IO
    ITG = TG_intensity(pressure, temperature, nozzle_radius, det_rad, xs, a)
    IIO = isentropic_intensity(pressure, temperature, nozzle_radius)
    assert ITG < IIO


@pytest.mark.unit
def testSikoraSmaller():
    # test that the sikora intensity is smaller than the TG
    ITG = TG_intensity(pressure, temperature, nozzle_radius, det_rad, xs, a)
    ISikora = sikora_intensity(
        pressure, temperature, nozzle_radius, det_rad, xs, a, S_r, rs, xs
    )
    assert ISikora < ITG
