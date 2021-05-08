import numpy as np
from scipy.constants import k as k_b

msuk = 4.814053e-4
h_0 = 1e-3
pas = 1.001


class BoltzmannSpherical:
    def __init__(
        self,
        interaction_potential: str = "LJ",
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
        for temperature, pressure in zip(temperatures, pressures):
            self.solve_expansion(temperature, pressure)
        self.save_results()

    def solve_expansion(
        self, temperature: [int, float], pressure: [int, float], accuracy=1e-2
    ):
        self.initialise_expansion(temperature, pressure)

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
        # initial conditions checked and ok
        t_l0 = 1 / (1 + (self.gamma - 1) / 2 * M ** 2)
        self.t_l = temp * t_l0
        self.u_l = M * np.sqrt(self.gamma * self.t_l / msuk)
        self.n_l = press / (k_b * temp) * t_l0 ** (1 / (self.gamma - 1))
        # check that the rho makes sense
        self.assert_unchanged_rho_physical(rho, rho_r, temp)

    def assert_unchanged_rho_physical(self, rho, rho_r, temperature):
        if rho_r == rho:
            assert self.T_Ent_E == 5 * temperature
        else:
            pass

    def save_results(self):
        return NotImplementedError

    def expansion_step(self):
        return NotImplementedError
