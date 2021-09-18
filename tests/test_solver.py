import pytest

from corelib.solvers.boltzmann_spherical import BoltzmannSpherical
from corelib.translated_fortran_scripts import (
    check_solutions_equal,
    load_data_benchmark,
    rho_real_python,
)


@pytest.mark.parametrize("pressure", [1, 6])
@pytest.mark.unit
def test_load_data(pressure):
    potential_type = "LJ_re"
    loaded_dict = load_data_benchmark(
        [
            pressure,
        ],
        [
            311,
        ],
        potential_type,
    )
    for kk in loaded_dict.keys:
        assert loaded_dict[kk] is not None


@pytest.mark.parametrize("pressure", [1, 6])
@pytest.mark.integration
def test_benchmark_expansion(pressure):
    # load benchmark data
    potential_type = "LJ_re"
    key = (311, pressure)
    loaded_dict = load_data_benchmark(
        [
            pressure,
        ],
        [
            311,
        ],
        potential_type,
    )
    spherical_solver = BoltzmannSpherical(rho_correction=rho_real_python)
    spherical_solver.get_collision_integral()
    spherical_solver.solve_expansion(311, pressure)
    check_solutions_equal(spherical_solver.dict_results[key], loaded_dict[key])


# @pytest.mark.integration
# def test_benchmark_expansions():
#     pressure_list = [1,6]
#     ## load benchmark data
#     potential_type = "LJ_re"
#     loaded_dict = load_data_benchmark(pressure_list, [311], potential_type)
#     spherical_solver = BoltzmannSpherical(rho_correction=rho_real_python)
#     spherical_solver.solve_expansions([311], pressure_list)
#     for pressure in pressure_list:
#         tup = (311,pressure)
#         assert check_solutions_equal(spherical_solver.dict_results[tup],
#         loaded_dict[tup])
