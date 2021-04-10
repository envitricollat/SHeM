from scipy.constants import k as k_b
import numpy as np
import pandas as pd
# this is a raw script that ports the original fortran code
# by @Gianangelo Bracco to python

msuk = 4.814053e-4
r_l = 2.5
gamma = 5/3
h_0 = 1e-3
pas = 1.001
acc = 1e-2

# define the values of p_0dv and T_0v (source properties)

p0d_v = [1,2,4,6,8,12,16,20,25,30,35,40,50,60,70,80,90,100,110,120,130,150,160,180,200]
t0_v = [300,131,132]
n_press = len(p0d_v)
n_temp = len(t0_v)

# read the values of Omega(T) (the collision integral)
potential_type = 'LJ' # lennard-jones potential
path = 'experimental_data/'+'omega_'+potential_type
omega = pd.read_table(path)
