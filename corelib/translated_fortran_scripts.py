import os
from itertools import product

import numpy as np
import pandas as pd
from pandas._testing import assert_series_equal

ROOT_DIR = os.path.abspath(os.curdir)


def calculate_third_degree_coefficients(x, y):
    c = [0, 0, 0]
    x10 = x[1] - x[0]
    xq01 = x[0] ** 2 - x[1] ** 2
    x20 = x[2] - x[0]
    xq02 = x[0] ** 2 - x[2] ** 2

    y10 = y[1] - y[0]
    y20 = y[2] - y[0]
    p = y20 / x20 - y10 / x10
    q = xq01 / x10 - xq02 / x20
    c[2] = p / q
    c[1] = (c[2] * xq02 + y20) / x20
    c[0] = y[0] - c[2] * x[0] ** 2 - c[1] * x[0]
    return c


def load_data_benchmark(pressure_list, temperature_list, potential_type="LJ"):
    dict_results = {}
    col_list = ["r", "n_r", "u_r", "t_l1", "t_l2"]
    for temperature, pressure in product(temperature_list, pressure_list):
        key = (temperature, pressure)
        try:
            path = (
                ROOT_DIR
                + "/numerical_data/benchmark/t"
                + str(temperature)
                + "p"
                + str(pressure)
                + potential_type
                + ".dat"
            )
            dict_results[key] = pd.read_table(
                path, sep=r"\s+", header=None, names=col_list
            )
        except FileNotFoundError:
            dict_results[key] = None
    return dict_results


def rho_real(T, P):
    # this function calculates the real rho for helium according to
    # McCarty and Arp (1990) Advances in cryogenic engineering.
    # C R. W. Fast. New York, Plenum Press. 35: 1465-1475.
    # C **  Echelle provisoire de temperature de 1976. Bureau
    # C international des poids et mesures, F-92310 Sevres, France
    # C *** Hands, B. A. and V. D. Arp (1981) A correlation of thermal
    # C conductivity data for helium. Cryogenics 21: 697-703.
    # todo: proper function docstring
    # P is in Bar, T in Kelvin
    rho = 0.1 * P
    dr = rho * 1e-2
    ind = 0
    indd = 0

    while ind == 0:
        Pc = EOS(T, rho)
        if Pc < P:
            rho = rho + dr
            indd = 0
        else:
            if indd == 0:
                dr = 0.25 * dr
                indd = 1
            rho = rho - dr
            if dr < 1e-10:
                ind = 1
    return rho


def rho_real_python(rho, T, P):
    # this function calculates the real rho for helium according to
    # McCarty and Arp (1990) Advances in cryogenic engineering.
    # C R. W. Fast. New York, Plenum Press. 35: 1465-1475.
    # C **  Echelle provisoire de temperature de 1976. Bureau
    # C international des poids et mesures, F-92310 Sevres, France
    # C *** Hands, B. A. and V. D. Arp (1981) A correlation of thermal
    # C conductivity data for helium. Cryogenics 21: 697-703.
    # todo: proper function docstring
    # P is in Bar, T in Kelvin
    assert 0.1 * P == rho
    rho = 0.1 * P
    dr = rho * 1e-2
    ind = 0
    indd = 0

    while ind == 0:
        Pc = EOS(T, rho)
        if Pc < P:
            rho = rho + dr
            indd = 0
        else:
            if indd == 0:
                dr = 0.25 * dr
                indd = 1
            rho = rho - dr
            if dr < 1e-10:
                ind = 1
    return rho


def EOS(temp, rho):
    # temperature in K, density in kg/m3, pressure in bar, rho in mol/L
    # The validity range of Hands and McCarty fits (used here)
    # !extends from 2 K up to 1500 K, for pressure up to 2 GPa
    # !Limited accuracy in the critical region (~ 5.19 K, ~ 69 kg/m3).
    # mol / L-> kg / m ** 3
    dens = rho * 4.0026
    tp2 = temp ** 2
    tp3 = temp ** 3
    tp4 = temp ** 4

    e_ = 7.636186157005e-18
    f = 6.097223119177e-19
    e_1 = 3.848665703556e-18
    tempo2 = (2.367205158795324e-7) * (dens ** 11) * (
        6.872567403738e-15 / tp3 - 4.595138561035e-15 / tp2
    ) + 1.477581743701978e-8 * dens ** 13 * (e_1 / tp4 - e_ / tp3 - f / tp2)
    e_2 = 0.00006075816692925613
    e_3 = 3.850153114958e-8
    e_4 = 2.067693644675999e-8
    e_5 = 3.792453641033501e-6
    e_6 = 1.888462892389e-12
    e_7 = 1.399040627e-11
    tempo3 = (
        e_2 * dens ** 7 * (e_3 / tp3 - e_4 / tp2)
        + (e_5) * dens ** 9 * (-((e_6) / tp4) - e_7 / tp2)
        + tempo2
    )
    e_8 = 0.01559457081650663
    e_9 = 0.009166109232806
    e_10 = 0.006885401367690001
    e_11 = 0.000973394851465435
    e_12 = 0.00003315398880031
    e_13 = 6.544314242937e-6
    tempo1 = (
        e_8 * dens ** 3 * (e_9 / tp3 - e_10 / tp2)
        + e_11 * dens ** 5 * (-(e_12 / tp4) - e_13 / tp2)
        + tempo3
    )
    e_14 = 0.003896110232475551
    e_15 = 2.6997269279e-6
    f1 = 0.00003954146691114
    e_16 = 5.093547838380999e-9
    tempo3 = e_14 * dens ** 4 * (e_15 - f1 / temp - e_16 * temp)
    e_17 = 0.00004708238429298
    e_18 = 0.002410763742104
    e_19 = 0.001132915232587
    e_20 = 1.454229259623e-6
    tempo3 = tempo3 + e_8 * dens ** 3 * (
        -e_17 + e_18 / (tp2) + e_19 / temp + e_20 * temp
    )
    e_21 = 0.06241882915014949
    e_22 = 0.007139657549318
    e_23 = 0.01589302471561998
    e_24 = 0.009728903861441
    e_25 = 0.001260692007853
    e_26 = 0.00004558980227431
    tempo3 = tempo3 + e_21 * dens ** 2 * (
        -e_22 - e_23 / tp2 + e_24 / temp + e_25 * temp ** 0.5 + e_26 * temp
    )
    e_27 = 1.348439408105817e-18
    e_28 = 6.304713842604079e-15
    e_29 = 0.002077227302253535
    tempo3 = (
        tempo3
        - (e_27 * dens ** 9) / tp2
        - (e_28 * dens ** 7) / temp
        + e_29 * dens * temp
    )
    e_30 = 1.510671273545713e-12
    e_31 = 0.0002061897349793637
    e_32 = 0.00001517967494360068
    e_33 = 3.298960057070999e-11
    e_34 = 6.446881346447997e-13
    e_35 = 0.0002431906389510405
    e_36 = 5.501158366750001e-8
    e_37 = 1.050712335784999e-8
    eos = 10 * (
        e_30 * dens ** 5
        + tempo1 / np.exp(e_31 * dens ** 2)
        + e_32 * dens ** 8 * (e_33 / tp2 + e_34 / temp)
        + e_35 * dens ** 6 * (-(e_36 / tp2) + e_37 / temp)
        + tempo3
    )
    return eos


def diff_tempwradius(r, t_r, c, dat_T, n_r, u_r):
    tp_r = [0, 0]
    msuk = 4.814053e-4
    ris = integrate_custom(t_r, c, dat_T)
    ris = 2 * (n_r / u_r) * ris / (t_r[0] * np.sqrt(t_r[1]))
    t1 = u_r ** 2.0 * msuk
    tp_r[0] = -2.0 * t_r[0] / r + ris
    tp_r[1] = (-2.0 * ris * t1 + 2.0 * t_r[1] * (-2.0 * t_r[0] / r + ris)) / (
        t1 - 3.0 * t_r[1]
    )
    return np.array(tp_r)


def integrate_custom(t_r, c, dat_T):
    # this integral is by far the slower part of the implementation
    # so the first that I will refactor and optimise
    A_O = 2.48004e-20
    ris = 0.0
    dchi = 1.0e-2
    chi = dchi
    flag1 = True
    i = 0
    while flag1:
        a = (t_r[1] - t_r[0]) / t_r[1]
        t_eff = t_r[0] / (1 - a * chi ** 2.0)
        try:
            i = np.where(dat_T > t_eff)[0][0] + 1
        except IndexError:
            print("error in the integral calculation")
        # flag2 = True
        # i = 0
        # while flag2:
        #     if dat_T[i] > t_eff:
        #         flag2 = False
        #     i = i + 1
        if i == 1:
            omega = 0
        # note that we do not integrate over the last row
        else:
            omega = A_O * (
                c[i - 1, 2] * t_eff ** 2.0 + c[i - 1, 1] * t_eff + c[i - 1, 0]
            )
        fun = t_eff ** (5.0 / 2.0) * omega * (3.0 * chi ** 2.0 - 1.0)
        ris = ris + fun * dchi
        chi = chi + dchi
        if chi > 1:
            flag1 = False
    return ris


def check_solutions_equal(
    experiment_df: pd.DataFrame,
    benchmark_df: pd.DataFrame,
    abs_tolerance_speed=1,
    abs_tolerance_temp=0.01,
    rel_tolerance=0.002,
):
    """ """
    last_index = min([benchmark_df.index.max(), experiment_df.index.max()])
    experiment_df = experiment_df[:last_index]
    benchmark_df = benchmark_df[:last_index]
    experiment_df = experiment_df.astype("float64")
    # add relative tolerance indensities
    # due to large swings in absolute value
    # depending on source pressure
    assert_series_equal(
        experiment_df["n_r"],
        benchmark_df["n_r"],
        check_exact=False,
        rtol=rel_tolerance,
        check_dtype=False,
    )
    # the power function of fortran ** at large exponents
    # does not coincide with python's and therefore a
    # discrepancy accumulates at large j's
    # - therefore the more leniant test
    # based on relative tolerance
    assert_series_equal(
        experiment_df["r"],
        benchmark_df["r"],
        check_exact=False,
        rtol=rel_tolerance,
    )
    # allow for 1 m/s discrepancy in the speed
    # (way below experimental error)
    assert_series_equal(
        experiment_df["u_r"],
        benchmark_df["u_r"],
        check_exact=False,
        atol=abs_tolerance_speed,
    )
    # allow for 0.01 K discrepancy in the temperatures
    # (numeric error of the fortran routine)
    assert_series_equal(
        experiment_df["t_l1"],
        benchmark_df["t_l1"],
        check_exact=False,
        atol=abs_tolerance_temp,
    )
    assert_series_equal(
        experiment_df["t_l2"],
        benchmark_df["t_l2"],
        check_exact=False,
        atol=abs_tolerance_temp,
    )


def interpolate_c(dat_T, dat_0):
    # initialise null variables as done in the original fortran code
    cd_old = [None, None, None]
    cd = cd_old
    x = [None, None, None]
    y = [None, None, None]
    c = np.zeros([len(dat_T) - 1, 3])
    # first routine: interpolate omega with T using a second order polinomial
    # [todo: vectorise]
    # todo: call fortran subroutines from python and test them with unit tests
    # todo: simplify using python functions and vectorisation
    for kk in np.arange(1, len(dat_T) - 1, 1):
        for dimension in np.arange(0, 3, 1):
            cd_old[dimension] = cd[dimension]
            x[dimension] = dat_T[kk - 1 + dimension]
            y[dimension] = dat_0[kk - 1 + dimension]
        cd = calculate_third_degree_coefficients(x, y)
        if kk == 1:
            for dimension in np.arange(0, 3, 1):
                c[kk - 1, dimension] = cd[dimension]
        else:
            for dimension in np.arange(0, 3, 1):
                c[kk - 1, dimension] = (cd[dimension] + cd_old[dimension]) / 2
                if kk == len(dat_T) - 1:
                    c[kk, dimension] = cd[dimension]
    return c
