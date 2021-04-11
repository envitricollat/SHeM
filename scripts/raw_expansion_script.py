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

def rho_real(T,P):
    # this function calculates the real rho for helium according to
    # McCarty and Arp (1990) Advances in cryogenic engineering.
    # C R. W. Fast. New York, Plenum Press. 35: 1465-1475.
    # C **  Echelle provisoire de temp�rature de 1976. Bureau
    # C international des poids et mesures, F-92310 Sevres, France
    # C *** Hands, B. A. and V. D. Arp (1981) A correlation of thermal
    # C conductivity data for helium. Cryogenics 21: 697-703.
    # todo: function docstring
    # P is in Bar, T in Kelvin
    rho = 0.1e0*P
    dr = rho*1e-2
    ind = 0
    indd = 0

    while ind==0:
        Pc = EOS(T,rho)
        if Pc<P:
            rho = rho+dr
            indd = 0
        else:
            if indd==0:
                dr = 0.25*dr
                indd = 1
            rho = rho-dr
            if dr<1e-10:
                ind = 1
    return rho

def EOS(temp,rho):
    # temperature in K, density in kg/m3, pressure in bar, rho in mol/L
    # The validity range of Hands and McCarty fits (used here)
    # !extends from 2 K up to 1500 K, for pressure up to 2 GPa
    # !Limited accuracy in the critical region (~ 5.19 K, ~ 69 kg/m3).
    # mol / L-> kg / m ** 3
    dens = rho*4.0026
    tp2 = temp**2
    tp3 = temp**3
    tp4 = temp**4
    tempo2 = (2.367205158795324e-7)*(dens**11)*(6.872567403738e-15 / tp3- 4.595138561035e-15/tp2)+1.477581743701978e-8 * dens ** 13 *(3.848665703556e-18/tp4 - 7.636186157005e-18/tp3-6.097223119177e-19/tp2)
    tempo3 = 0.00006075816692925613e0 * dens ** 7 *(3.850153114958e-8 / tp3 - 2.067693644675999e-8 / tp2)+(3.792453641033501e-6)*dens ** 9 * (-((1.888462892389e-12) / tp4)  - 1.399040627e-11 / tp2) + tempo2
    tempo1 = (0.01559457081650663e0 * dens ** 3* (0.009166109232806e0 / tp3- 0.006885401367690001e0 / tp2) + 0.000973394851465435e0 *dens ** 5 * (-(0.00003315398880031e0 / tp4) - 6.544314242937e-6 / tp2)+ tempo3)
    tempo3 = 0.003896110232475551e0*dens**4*(2.6997269279e-6  - 0.00003954146691114/temp - 5.093547838380999e-9*temp)
    tempo3 = tempo3 + 0.01559457081650663e0 * dens ** 3 *(-0.00004708238429298e0 + 0.002410763742104e0 / (tp2)+ 0.001132915232587e0 / temp + 1.454229259623e-6 * temp)
    tempo3 = tempo3 + 0.06241882915014949e0 * dens ** 2 *(-0.007139657549318e0 - 0.01589302471561998e0 / tp2+ 0.009728903861441e0 / temp+ 0.001260692007853e0 * temp ** 0.5e0+ 0.00004558980227431e0 * temp)
    tempo3 = tempo3 - (1.348439408105817e-18 * dens ** 9) / tp2 - (6.304713842604079e-15 * dens ** 7) / temp+ 0.002077227302253535e0 * dens * temp
    EOS = 10 * (1.510671273545713e-12 * dens ** 5 +tempo1 / np.exp(0.0002061897349793637e0 * dens ** 2)+ 0.00001517967494360068e0 * dens ** 8 *(3.298960057070999e-11 / tp2 + 6.446881346447997e-13 / temp)+ 0.0002431906389510405e0 * dens ** 6 * (-(5.501158366750001e-8 / tp2)+ 1.050712335784999e-8 / temp) + tempo3)
    return EOS

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
# second routine: obtain the names of the output files - saves them in a matrix
for l in range(n_temp):
    for k in range(n_press):
        t0 = t0_v[l]
        p0d = p0d_v[k]
        # correction for a real gas - only valid for helium
        rho_r = rho_real(t0,p0d)

# third routine: set initial conditions for an spherical approximation
# fourth routine: solve the system of differential partial equations




