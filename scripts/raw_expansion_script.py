import numpy as np
import pandas as pd
from scipy.constants import k as k_b

# this is a raw script that ports the original fortran code
# by @Gianangelo Bracco to python
from corelib.translated_fortran_scripts import (
    diff_tempwradius,
    interpolate_c,
    load_data_benchmark,
    rho_real,
    test_solutions_equal,
)

msuk = 4.814053e-4
r_l = 2.5
gamma = 5 / 3
h_0 = 1e-3
pas = 1.001
acc = 1e-2
# a different value of k_B was used in the fortran numerical calculations,
# dummy print to keep the import
print(k_b)
# we over-write for debug (and we will for testing)
# k_b = 1.380662e-23
# define the values of p_0dv and T_0v (source properties)

p0d_v = [
    1,
    2,
    4,
    6,
    8,
    12,
    16,
    20,
    25,
    30,
    35,
    40,
    50,
    60,
    70,
    80,
    90,
    100,
    110,
    120,
    130,
    150,
    160,
    180,
    200,
]
t0_v = [311, 131]
n_press = len(p0d_v)
n_temp = len(t0_v)

# read the values of Omega(T) (the collision integral)
potential_type = "LJ_re"  # lennard-jones potential with a real gas correction
if potential_type == "LJ_re":
    path = "../numerical_data/" + "omega_" + "LJ" + ".dat"
else:
    path = "../numerical_data/" + "omega_" + potential_type + ".dat"

omega = pd.read_table(path, sep=r"\s+", header=None, names=["dat_T", "dat_0"])
[dat_T, dat_0] = [omega["dat_T"].values, omega["dat_0"].values]
# load numerical benchmark and save in dictionary
loaded_dict = load_data_benchmark(p0d_v, t0_v, potential_type)

c = interpolate_c(dat_T, dat_0)
# second routine: obtain the names of the output files
# within the loop to solve for each combination
global_df = pd.DataFrame(
    index=np.arange(n_temp * n_press),
    columns=["pressure", "temperature", "speed_ratio"],
)

global_dictionary_result = {}
count = 0
sr = 0

column_names = ["r", "n_r", "u_r", "t_l1", "t_l2"]
for temperature_index in range(n_temp):
    for k in range(n_press):
        t0 = t0_v[temperature_index]
        p0d = p0d_v[k]
        benchmark_df = loaded_dict["t" + str(t0) + "p" + str(p0d)]
        experiment_df = pd.DataFrame(None, columns=column_names)
        # correction for a real gas - only valid for helium
        rho_r = rho_real(t0, p0d)
        # L*bar in Joule/molecola
        dE_real = 2 * p0d / rho_r / 1e-2 / 6.02214179e23 / k_b
        # internal energy 3 kT plus PV
        T_Ent_E = 3.0 * t0 + dE_real
        name = "t" + str(t0) + "p" + str(p0d) + "LJ_re" + ".dat"
        # third routine: set initial conditions for an spherical approximation
        M = (r_l ** (gamma - 1.0)) * (
            3.232 - 0.7563 / r_l + 0.3937 / (r_l ** 2) - 0.0729 / (r_l ** 3)
        )
        # initial conditions checked and ok
        t_l0 = 1 / (1 + (gamma - 1.0) / 2 * M ** 2)
        t_l = t0 * t_l0
        u_l = M * np.sqrt(gamma * t_l / msuk)
        n_l = p0d / (k_b * t0) * t_l0 ** (1.0 / (gamma - 1.0))
        # fourth routine: solve the system of differential partial equations
        fi = n_l * u_l * r_l ** 2
        r = r_l
        u_r = u_l
        n_r = n_l
        t_r = [t_l, t_l]
        j = 0
        experiment_df.loc[j, :] = [r_l, n_l, u_l, t_l, t_l]

        flag = True
        while flag:
            h = h_0 * pas ** j
            tp_r = diff_tempwradius(r, t_r, c, dat_T, n_r, u_r)
            t_rn = t_r + h * tp_r / 6.0
            r = r + h / 2.0
            for _i in [0, 1]:
                t_r1 = t_r + h / 2.0 * tp_r
                tp_r = diff_tempwradius(r, t_r1, c, dat_T, n_r, u_r)
                t_rn = t_rn + h * tp_r / 3.0
            r = r + h / 2.0
            t_r1 = t_r + h * tp_r
            tp_r = diff_tempwradius(r, t_r, c, dat_T, n_r, u_r)
            t_r = t_rn + h * tp_r / 6.0
            u_r = np.sqrt((T_Ent_E - 3.0 * t_r[1] - 2.0 * t_r[0]) / msuk)
            n_r = fi / (u_r * r ** 2)
            rapp = t_r[0] / t_r[1]
            if rapp < acc:
                flag = False
                # speed ratio
                sr = np.sqrt(msuk * u_r ** 2.0 / (2.0 * t_r[1]))
            j = j + 1
            experiment_df.loc[j, :] = [r, n_r, u_r, t_r[0], t_r[1]]

        if benchmark_df is not None:
            test_solutions_equal(experiment_df, benchmark_df)
        else:
            print("There is no benchmark dataframe!")
        print("ALl tests have passed!")
        print("System of differential eqns finished. Parameters:")
        print("Pressure: " + str(p0d))
        print("Temperature: " + str(t0))
        print("Speed ratio:" + str(sr))
        # todo: include tests for this variables
        global_df.loc[count, "pressure"] = p0d
        global_df.loc[count, "temperature"] = t0
        global_df.loc[count, "speed_ratio"] = sr
        count = count + 1
print("loop completed!")
global_df.to_csv("../experimental_data/overall_results.csv")
