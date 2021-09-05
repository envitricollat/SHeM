import pytest

from corelib.solvers.boltzmann_spherical import BoltzmannSpherical


@pytest.mark.parametrize("pressure", [1, 6])
@pytest.mark.unit
def compare_benchmark_results(pressure):
    spherical_solver = BoltzmannSpherical()
    spherical_solver.get_collision_integral()
    spherical_solver.solve_expansion(311, pressure)
    assert True
