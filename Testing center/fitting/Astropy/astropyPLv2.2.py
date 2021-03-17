#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Feb  8 11:15:29 2021

@author: asif
"""

import numpy as np
import matplotlib.pyplot as plt
from astropy.modeling import models,fitting

def nu_break_evln(t,nu_0,m=2,t0=1):
    return nu_0*np.power(t/t0,-m)

def F_0_evln(t,F_0,m=2,t0=1):
    return F_0*np.power(t/t0,-m)

def calc_flux(nu,F_0t,nu_bt,alpha_1t,alpha_2t):
    """Flux as flux[t_i,nu_i] for given arrays of ts and nus. Give other paramters as arrays 
    i.e F_0t[ti]= amplitude for time t_i"""
    flux=np.zeros((n,len(nu)))           
    f=models.SmoothlyBrokenPowerLaw1D()                         
    for i in range(n):
        f.amplitude=F_0t[i]
        f.x_break=nu_bt[i]
        f.alpha_1=alpha_1t[i]
        f.alpha_2=alpha_2t[i]
        f.delta=delta_t[i]
        flux[i]=f(nu)
    
    return flux

def f_nu(t,F_0,m_f,nu_b0,m_b,alpha_1t,alpha_2t,delta_t):
    F_0t=F_0_evln(t,F_0,m_f)
    nu_bt=nu_break_evln(t,nu_b0,m_b)
    flux=calc_flux(nu,F_0t,nu_bt,alpha_1t,alpha_2t)          
    return flux

def IdxEvln(F,freqs):
    "Fits a power law to F[i1,i2] for all i1 with i2 being freqs"
    
    s_idx=[]
    for ff in F:
        model_x=np.log(freqs)
        model_y=np.log(ff)
        temp = np.polyfit(model_x,model_y, 1)
        s_idx.append(temp)
    return np.array(s_idx)

# def calc_phIdx_evln(t,nu,F_0t,nu_bt,alpha_1t,alpha_2t):
    # flux=f_nu(t,F_0,m_f,nu_b0,m_b,alpha_1t,alpha_2t,delta_t)

# t = np.logspace(0.7, 6.3, 500)
t=np.logspace(1,8,num=7)
n=len(t)
F_0=1e4
m_f=2
F_0t=F_0_evln(t,F_0,m_f)
nu_b0=1e28
m_b=2
nu_bt=nu_break_evln(t,nu_b0,m_b)
delta_t=0.4*np.ones(n)
alpha_1t=1*np.ones(n)
alpha_2t=3*np.ones(n)


nu=2.4e17*np.logspace(-0.5,1.,100)


# f = models.SmoothlyBrokenPowerLaw1D(amplitude=1, x_break=2000,
                                    # alpha_1=1, alpha_2=4)
# flux=calc_flux(nu,F_0t,nu_bt,alpha_1t,alpha_2t)
flux=f_nu(t, F_0, m_f, nu_b0, m_b, alpha_1t, alpha_2t, delta_t)

# ph=[]
# for f in flux:
#     s=IdxEvln(f,nu)
#     ph.append(s)

s=IdxEvln(flux,nu)
ph=-s[:,0]+1


plt.figure()
plt.title("amplitude=1, x_break=20, alpha_1=-2, alpha_2=2")

# f.delta = 0.5
plt.loglog(t, flux, '--', label='delta=0.5')
