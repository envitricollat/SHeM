import numpy as np
import pandas as pd
from scipy.constants import k as k_b

from corelib.translated_fortran_scripts import diff_tempwradius, interpolate_c

msuk = 4.814053e-4
h_0 = 1e-3
pas = 1.001


class BoltzmannSpherical:
    def __init__(
        self,
        interaction_potential: str = "LJ_re",
        velocity_distribution: str = "Maxwellian",
        gamma=5 / 3,
        r_l=2.5,
        rho_correction=lambda rho, x, y: rho,
    ):
        self.interaction_potential = interaction_potential
        self.velocity_distribution = velocity_distribution
        self.gamma = gamma
        self.rho_correction = rho_correction
        self.r_l = r_l

    def solve_expansions(self, temperatures: list, pressures: list):
        self.get_collision_integral()
        for temperature, pressure in zip(temperatures, pressures):
            self.solve_expansion(temperature, pressure)
        self.save_results()

    def solve_expansion(
        self,
        temperature: [int, float],
        pressure: [int, float],
        qs_condition=1e-2,
        stopstep=None,
    ):
        self.initialise_expansion(temperature, pressure)
        while self.condition > qs_condition:
            self.expansion_step()

    def initialise_expansion(self, temp, press):
        rho = 0.1 * press
        rho_r = self.rho_correction(rho, temp, press)
        # L*bar in Joule/molecola
        dE_real = 2 * press / rho_r / 1e-2 / 6.02214179e23 / k_b
        # internal energy 3 kT plus PV
        self.T_Ent_E = 3.0 * temp + dE_real
        # third routine: set initial conditions for an spherical approximation
        M = (self.r_l ** (self.gamma - 1)) * (
            3.232
            - 0.7563 / self.r_l
            + 0.3937 / (self.r_l ** 2)
            - 0.0729 / (self.r_l ** 3)
        )
        t_l0 = 1 / (1 + (self.gamma - 1) / 2 * M ** 2)
        # save initial conditions
        self.t_l = temp * t_l0
        self.t_r = [self.t_l, self.t_l]
        self.u_l = M * np.sqrt(self.gamma * self.t_l / msuk)
        self.n_l = press / (k_b * temp) * t_l0 ** (1 / (self.gamma - 1))
        self.fi = self.n_l * self.u_l * self.r_l ** 2
        # assign initial conditions to variables
        self.r = self.r_l
        self.u_r = self.u_l
        self.n_r = self.n_l
        self.condition = 1
        self.j = 0
        # check that the rho makes sense
        self.assert_unchanged_rho_physical(rho, rho_r, temp)

    def assert_unchanged_rho_physical(self, rho, rho_r, temperature):
        if rho_r == rho:
            assert self.T_Ent_E == 5 * temperature
        else:
            pass

    def save_results(self):
        return NotImplementedError

    def get_collision_integral(self):
        # read the values of Omega(T) (the collision integral)
        potential_type = (
            self.interaction_potential
        )  # lennard-jones potential with a real gas correction
        if potential_type == "LJ_re":
            path = "../numerical_data/" + "omega_" + "LJ" + ".dat"
        else:
            path = "../numerical_data/" + "omega_" + potential_type + ".dat"
        cols = ["dat_T", "dat_0"]
        omega = pd.read_table(path, sep=r"\s+", header=None, names=cols)
        self.dat_T = omega["dat_T"].values
        self.dat_0 = omega["dat_0"].values
        self.c = interpolate_c(self.dat_T, self.dat_0)

    def expansion_step(self):
        # todo: change to use a memory class
        h = h_0 * pas ** self.j
        tp_r = diff_tempwradius(
            self.r, self.t_r, self.c, self.dat_T, self.n_r, self.u_r
        )
        t_rn = self.t_r + h * tp_r / 6.0
        self.r = self.r + h / 2.0
        for _i in [0, 1]:
            self.t_r1 = self.t_r + h / 2.0 * tp_r
            tp_r = diff_tempwradius(
                self.r, self.t_r1, self.c, self.dat_T, self.n_r, self.u_r
            )
            t_rn = t_rn + h * tp_r / 3.0
        self.r = self.r + h / 2.0
        self.t_r1 = self.t_r + h * tp_r
        tp_r = diff_tempwradius(
            self.r, self.t_r, self.c, self.dat_T, self.n_r, self.u_r
        )
        t_r = t_rn + h * tp_r / 6.0
        self.u_r = np.sqrt((self.T_Ent_E - 3.0 * t_r[1] - 2.0 * t_r[0]) / msuk)
        self.n_r = self.fi / (self.u_r * self.r ** 2)
        self.condition = t_r[0] / t_r[1]
