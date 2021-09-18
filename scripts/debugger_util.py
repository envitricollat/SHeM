from corelib.solvers.boltzmann_spherical import BoltzmannSpherical
from corelib.translated_fortran_scripts import (
    check_solutions_equal,
    load_data_benchmark,
    rho_real_python,
)

pressure = 2
potential_type = "LJ_re"
key = (311, pressure)
loaded_dict = load_data_benchmark([pressure], [311], potential_type)
spherical_solver = BoltzmannSpherical(rho_correction=rho_real_python)
spherical_solver.get_collision_integral()
spherical_solver.solve_expansion(311, pressure)
print("we are here")
check_solutions_equal(spherical_solver.dict_results[key], loaded_dict[key])
