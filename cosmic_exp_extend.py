"""
----- functions for computing cosmic time, distance, volume, Hubble parameter -----
with radiation, curvature...  e.g. for early universe
"""

import os
import scipy.integrate
from scipy.interpolate import splrep, splev
import numpy as np
from scipy import vectorize

from param import om0, h0, Ggrav_MsunMpckms, Msun_g, Mpc_cm, rho0_Gamma, hubbleunit2cgs

rho_crit0 = 3*100**2/8./np.pi/Ggrav_MsunMpckms   # in unit of h^2 Msun/Mpc^3 

# **
def rho_crit0_cgs(h0=h0):
    """
    critical density at z=0 in cgs unit
    """
    return rho_crit0 * h0**2 * Msun_g / Mpc_cm**3

def omGamma(h0=h0):
    return rho0_Gamma / rho_crit0_cgs(h0)

def Ez(z, om0=om0, omk=0., omGamma=0.):
    omL = 1. - om0 - omGamma - omk 
    return (om0 * (1.+z)**3 + omGamma * (1.+z)**4 + omL + omk * (1.+z)**2)**0.5

# **
def Hz(z, om0=om0, omk=0, omGamma=0., h0=h0):  
    """ 
    Hubble parameter as function of redshift in unit of km/s/Mpc 
    """
    Ezz = Ez(z, om0, omk, omGamma)
    H0 = h0 * 100.
    return Ezz * H0

# **
def Hz_in_cgs(z, om0=om0, omk=0, omGamma=0., h0=h0):  
    """ 
    in unit of s-1
    """ 
    return Hz(z,om0,omk,omGamma,h0) * hubbleunit2cgs

def _int4t(a, om0=om0, omk=0, omGamma=0., h0=h0):   
    """ 
    integrand for time: dt = da/H/a
    """ 
    z = 1. / a - 1.
    return 1. / Hz_in_cgs(z, om0, omk, omGamma, h0) / a      #1/H/a

# ** 
def propertime(z, om0=om0, omk=0, omGamma=0., h0=h0):  
    """ 
    proper time in unit of seconds. age of univ: ~ 4.3e17 s
    """ 
    a = 1. / (1.+z)
    return scipy.integrate.quad(_int4t, 0, a, (om0,omk,omGamma,h0))[0]

propertimes = vectorize(propertime)

def tfromz_fast(z, tcal):
    return splev(z, tcal)

# **
def tfromz(z):
    """
    proper time from redshift
    using presaved table for fast computation for many redshifts
    """
    tcal = get_tcal()
    return splev(z, tcal)

def zfromt_fast(t, zcal):
    if t < zcal[0].min():
        print('only for z<'+str(zcal[1].max())+'!')
    elif t > zcal[0].max():
        print('only for z>'+str(zcal[1].min())+'!')
    else:
        return splev(t,zcal)

def get_tcal():  
    """ 
    z -> t
    """ 
    if os.path.isfile('tszs.tab'):
        data = np.loadtxt("tszs.tab")
    else:
        save_tszs()
        data = np.loadtxt("tszs.tab")
    return splrep(data[:,1], data[:,0])  # z,t

def zfromt_logzp1cal(t, logzp1cal):
    tt = np.array(t)
    if np.log10(tt.min()) < logzp1cal[0].min():
        print('only for z<'+str(10**(logzp1cal[1].max())-1.)+'!')
    elif np.log10(tt.max()) > logzp1cal[0].max():
        print('only for z>'+str(10**(logzp1cal[1].min())-1.)+'!')
    else:
        return 10**(splev(np.log10(t), logzp1cal)) - 1.

def get_logzp1cal_zi_zf(zi=20, zf=-0.5, nz=100, om0=om0, omk=0, omGamma=0., h0=h0, filename=None):
    """ 
    t -> z
    """
    cos = str(om0)+'_'+str(omk)+'_'+str(omGamma)+'_'+str(h0)
    if filename==None: filename = 'logtslogzp1s_'+str(zi)+'_'+str(zf)+'_'+str(nz)+'_'+cos+'.tab'
    if os.path.isfile(filename): data = np.loadtxt(filename)
    else: data = save_tszs_zi_zf(zi,zf,nz,om0,omk,omGamma,h0)
    return splrep(data[:,0], data[:,1])  # logt, logzp1 

def save_tszs_zi_zf(zi, zf, nz, om0=om0, omk=0, omGamma=0., h0=h0, filename=None):
    cos = str(om0)+'_'+str(omk)+'_'+str(omGamma)+'_'+str(h0)
    logzp1s = np.linspace(np.log10(zi+1.),np.log10(zf+1.),nz)
    zs = 10**logzp1s - 1. 
    logts = np.log10(np.array([propertime(z, om0, omk, omGamma, h0) for z in zs]))
    data = np.array(zip(logts, logzp1s))
    if filename==None:
        np.savetxt('logtslogzp1s_'+str(zi)+'_'+str(zf)+'_'+str(nz)+'_'+cos+'.tab',data)  #shape(60,2)
    else:
        np.savetxt(filename,data)
    return data

def _int4chi(a, om0=om0, omk=0, omGamma=0.):
    """ 
    integrand for chi: dchi = c*da/H/a^2, in unit of h^-1 Mpc
    """
    z = 1./a - 1.
    return 2998. / Ez(z,om0,omk,omGamma) / a**2


# **
def chi(z, om0=om0, omk=0, omGamma=0.):
    """
    comoving distance, h^-1 Mpc
    """
    a = 1. / (1.+z)
    return scipy.integrate.quad(_int4chi, a, 1., (om0,omk,omGamma))[0]

# **
def Da(z, om0=om0, omk=0, omGamma=0.):
    """
    angular diameter distance, h^-1 Mpc, for flat universe (fK(chi)=chi)
    """
    return chi(z,om0,omk,omGamma) / (1.+z)

# **
def drprop_dz(z, om0=om0, omk=0, omGamma=0.):
    """ 
    proper separation along l.o.s. dr_prop = a dchi = c*a/H dz, in unit of h^-1 Mpc
    """
    return c_cgs / 1.e7 / Ez(z, om0, omk, omGamma) / (1. + z)

# **
def dVprop_dz(z, fsky=1, om0=om0, omk=0, omGamma=0.):
    """ 
    proper volume element dV_prop = Da^2 * omega * dr_prop, in unit of (h^-1 Mpc)^3
    """
    return Da(z, om0, omk, omGamma)**2 * drprop_dz(z, om0, omk, omGamma) * 4 * np.pi * fsky

# **
def Vprop(z, fsky=1, om0=om0, omk=0, omGamma=0.):
    """ 
    proper volume V_prop from z=0 to z. = \int Da^2 * omega * dr_prop, in unit of (h^-1 Mpc)^3
    """
    return scipy.integrate.quad(dVprop_dz, 0., z, (fsky,om0,omk,omGamma))[0]
    

# **
def rho_crit(z, om0=om0, omk=0, omGamma=0.):
    """ 
    (physical) critical density of the universe, in unit of h^2 Msun/Mpc^3
    """
    return Ez(z,om0,omk,omGamma)**2 * rho_crit0

# **
def rho_crit_cgs(z, om0=om0, omk=0, omGamma=0., h0=h0):
    return rho_crit(z, om0, omk, omGamma) * h0**2 * Msun_g / Mpc_cm**3








