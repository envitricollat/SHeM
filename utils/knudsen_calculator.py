# this is a simple python script that calculates
# the denominator of the following relation
# Knudsen_number = distance_to_source^2 / denominator
# using all simple assumptions [pauly's eqn 2.45 from
# "Atom Molecule and Cluster Beams I]
# all units in ISU
# skimmer radius
import numpy as np
from scipy.constants import k as k_b

r_s = 0.5e-3
# nozzle radius
r_nz = 50e-6
# kinetic cross section
sigma = 260e-12
# gamma (for helium)
gamma = 5 / 3
# pressure (Pa)
P_0 = 100e3
# temperature (K)
T_0 = 300
# effective point source (distance from nozzle)
ps = 2.5
# approximate speed ratio
S = 90
# potential (in this case lennard jones)
n_p = 13
#
F_gamma = (1 + (gamma - 1) / 2) ** (-1 / (gamma - 1))
denominator = (
    sigma * (P_0 / (k_b * T_0)) * F_gamma * (ps * r_nz) ** 2 * r_s * np.sqrt(2)
)
#
d_K_1 = np.sqrt(denominator)
print("denominator: ", denominator)
# print some helpers
print("distance at K=1 (m): ", d_K_1)
# modified knudsen number calculations
correction = ((2 / 5) * S ** 2) ** (-2 / (n_p - 1))
d_K_1_mod = np.sqrt(denominator / correction)
print("correction", correction)
print("distance at K=1 mod (m): ", d_K_1_mod)
