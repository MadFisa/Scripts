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
    """break frequency evolution accordin to a powerlaw"""
    return nu_0*np.power(t/t0,-m)

def F_0_evln(t,F_0,m=2,t0=1):
    """Flux evolution accordin to a powerlaw"""
    return F_0*np.power(t/t0,-m)

def calc_flux(nu,F_0t,nu_bt,alpha_1t,alpha_2t,delta_t):
    """Flux as flux[t_i,nu_i] for given arrays of ts and nus. Give other paramters as arrays 
    i.e F_0t[ti]= amplitude for time t_i"""
    n=len(F_0t)
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
    # if not isinstance(t, np.ndarray):
    #     t=np.array([t])
    F_0t=F_0_evln(t,F_0,m_f)
    nu_bt=nu_break_evln(t,nu_b0,m_b)
    flux=calc_flux(nu,F_0t,nu_bt,alpha_1t,alpha_2t,delta_t)          
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

def calc_phIdx_evln(t,nu,F_0,nu_b,alpha_1t,alpha_2t,delta_t):
    flux=f_nu(t,F_0,m_f,nu_b0,m_b,alpha_1t,alpha_2t,delta_t)
    s=IdxEvln(flux,nu)
    ph=-s[:,0]+1
    return ph

##################################Read light curve data########################

GRB="111209"
ph_file_name="ph_idx_"+"GRB "+ GRB+".txt"
ph_file=np.loadtxt(ph_file_name)
t_data=ph_file[:,0]
ph_data=ph_file[:,3]
ph_err=2*ph_file[:,4] 
nu=2.4e17*np.logspace(-0.5,1.,100)
##############################################################################


# def model_f(t,nu,beta1,beta2,p_F,p_b,F_0,nu_b0,s=2.2):
#     """Calculates model flux for given set of parameters"""
#     F_0t=F_0_evln(t,F_0,p_F)
#     nu_bt=nu_break(t,nu_b0,p_b)
#     nu=10*2.4e17 #10 Kev To Hz
#     F=F_arb(nu,nu_bt,F_0t,beta1,beta2,s)
#     return F

def model_ph(t,F_0,nu_b0,alpha_1,alpha_2):
    """Calculates photon index as observed in window of freqs nu as per model model_f"""
    if not isinstance(t, np.ndarray):
        t=np.array([t])
    t=np.array([t])
    nu=2.4e17*np.logspace(-0.5,1.,100)#XRT range
    n=len(t)
    delta_t=0.4*np.ones(n)
    alpha_1t=alpha_1*np.ones(n)
    alpha_2t=alpha_2*np.ones(n)
    ph_idx_t=calc_phIdx_evln(t,nu,F_0,nu_b0,alpha_1t,alpha_2t,delta_t)
    
    return ph_idx_t
t=1e6
# t = np.array([1e7])
# t=np.logspace(1,8,num=7)
# t = np.logspace(0.7, 6.3, 500)
# n=len(t)
F_0=1e4
m_f=2
# F_0t=F_0_evln(t,F_0,m_f)
nu_b0=1e28
m_b=2
# nu_bt=nu_break_evln(t,nu_b0,m_b)
delta=0.4
alpha_1=1
alpha_2=3

ph=model_ph(t,F_0,nu_b0,alpha_1,alpha_2)


plt.figure()

plt.scatter(t,ph)
plt.xscale("log")
# plt.yscale("log")
