# a simple python script that calculates the intensity
# in a microscope using an analytical intensity equation according
# to certain parameters
# pressure in pascals
# this replicates the calculations published in [TO BE DISCLOSED]

import numpy as np
import pandas as pd
from scipy.constants import k as k_b

from corelib.analytical_models.sikora import sikora_intensity

# physical constant
m = 4 * 1.66053886e-27
# results dir
rdir = "./results/"
# load expansion parameters at the skimmer
expansion_highP = np.loadtxt(
    "./numerical_data/papercals/nz200_dist_nz_hp.dat", unpack=True
)
expansion_lowP = np.loadtxt(
    "./numerical_data/papercals/nz200_dist_nz_lp.dat", unpack=True
)
expansion_all = np.vstack([expansion_highP.T, expansion_lowP.T])
del expansion_highP, expansion_lowP
# save into a dataframe
expansion_df = pd.DataFrame(
    expansion_all,
    columns=[
        "pressure",
        "pressurexnozzle",
        "temperature",
        "S_par",
        "Qs",
        "v_He",
        "Tpar",
        "Tper",
    ],
)
# add column with perpendicular speed ratio
expansion_df["S_per"] = np.sqrt(
    m * expansion_df["v_He"] ** 2 / (2 * k_b * expansion_df["Tper"])
)
# nozzle radius in meters
nozzle_radius = 100e-6
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

# given a list of speed ratios and pressures, get the sikora matrix
# we are just going to add to the current dataframe columns
# we are assuming Qs=Skimmer
expansion_df["intensity_S_par"] = expansion_df.apply(
    lambda x: sikora_intensity(
        x["pressure"],
        x["temperature"],
        nozzle_radius,
        det_rad,
        xs,
        a,
        x["S_par"],
        rs,
        xs,
    ),
    axis=1,
)
expansion_df["intensity_S_per"] = expansion_df.apply(
    lambda x: sikora_intensity(
        x["pressure"],
        x["temperature"],
        nozzle_radius,
        det_rad,
        xs,
        a,
        x["S_per"],
        rs,
        xs,
    ),
    axis=1,
)
# saving this into the results directory
expansion_df.to_csv(rdir + "expansion_df.csv")
