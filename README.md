# SHeM
![Python Versions Supported](https://img.shields.io/badge/python-3.8+-blue.svg)
[![codecov](https://codecov.io/gh/envitricollat/SHeM/branch/main/graph/badge.svg?token=NQ0YI56RT5)](https://codecov.io/gh/envitricollat/SHeM)
![Continuous Integration](https://github.com/envitricollat/SHeM/actions/workflows/python-app.yml/badge.svg)

This repo contains utilities for Scanning Helium Microscopy, Molecular beams, and resolving the supersonic expansion of Helium Gas into vacuum.

This repo aims to combine physics with software engineering by covering most of its code with unit and integration tests. 
Any development proposed here needs to fulfill these software tests in order to be merged.

## Usage
To solve a supersonic expansion of gas into vacuum under a spherical assumption replicate the following steps:
1. Install the project dependencies with `pip install -r requirements.txt`
2. Define the pressure range for which you want to solve the expansion `pressure_list = [1, 6]`
3. Define the temperature range (in Kelvin) for which you want to solve the expansion `temperature_list = [200, 314]`
4. Initialise the solver class `spherical_solver = BoltzmannSpherical()`
5. Call the class to solve the expansions `spherical_solver.solve_expansions(temperature_list, pressure_list)`
