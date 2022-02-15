""" 
when changing parameters, need to update all tables! 
""" 
from numpy import pi

#-- cosmology params --#
om0 = 0.27
ode0 = 0.73
obh2 = 0.0223
w = -1
h0 = 0.7
omb = 0.0469
sigma8 = 0.82
n_s = 0.96

#-- physical constants --#
Ggrav_cgs = 6.67259e-8   # cm^3/g/s^2
mp = 1.6726231e-24  # in g
kb = 1.380658e-16  # erg/K
sigma_T = 6.6524574e-25  # cm^2
mec2 = 0.51099891 * 1.60217657e-6  # erg  
c_cgs = 2.9979e10
hbar = 1.0546e-27 


#--additional cosmology params --#
Tcmb = 2.726
rho0_CMB = (kb * Tcmb)**4 / (hbar * c_cgs)**3 * pi**2 / c_cgs**2 / 15
rho0_Gamma = 1.68 * rho0_CMB  # for high z where neutrinos are relativistic, not for z=0!!!


#-- additional physical params --#
mu_e = 1.14   # for fully ionized premordial gas
mu = 0.59

#-- units conversion --#
Msun_g = 1.9891e33
Mpc_cm = 3.08567758e24
rho_MsunMpc_cgs =  6.770e-41  # from h^2 Msun/Mpc^3 to h^2 g/cm^3
eV_erg = 1.60217657e-12
hubbleunit2cgs = 1e5 / Mpc_cm # km/s/Mpc to s^-1
K_eV = 8.621738e-5

#
Ggrav_MsunMpckms = Ggrav_cgs / 1e10 / Mpc_cm*Msun_g  # in unit of (km/s)^2 * Mpc/ Msun

