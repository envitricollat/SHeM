from scipy.constants import k as k_b
import numpy as np
import pandas as pd
# this is a raw script that ports the original fortran code
# by @Gianangelo Bracco to python
def calculate_third_degree_coefficients(x,y):
    c = [0,0,0]
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
path = '../numerical_data/'+'omega_'+potential_type+'.dat'
omega = pd.read_table(path, sep="\s+",header=None, names=['dat_T', 'dat_0'])
[dat_T,dat_0] = [omega['dat_T'].values,omega['dat_0'].values]

# initialise null variables as done in the original fortran code
cd_old = [None,None,None]
cd = cd_old
x=[None,None,None]
y=[None,None,None]
c = np.zeros([len(dat_T)-1,3])

# first routine: interpolate omega with T using a second order polinomial [todo: vectorise]
# todo: call fortran subroutines from python and test them with unit tests
# todo: simplify using python functions and vectorisation
for kk in np.arange(1,len(dat_T)-1,1):
    for dimension in np.arange(0,3,1):
        cd_old[dimension]=cd[dimension]
        x[dimension] = dat_T[kk - 1 + dimension]
        y[dimension] = dat_0[kk - 1 + dimension]
    cd = calculate_third_degree_coefficients(x,y)
    if kk==1:
        for dimension in np.arange(0, 3, 1):
            c[kk-1,dimension] = cd[dimension]
    else:
        for dimension in np.arange(0, 3, 1):
            c[kk-1,dimension] = (cd[dimension]+cd_old[dimension])/2
            if kk==len(dat_T)-2:
                c[kk,dimension] = cd[dimension]
print(c)
# second routine: obtain the names of the output files
# third routine: set initial conditions for an spherical approximation
# fourth routine: solve the system of differential partial equations




