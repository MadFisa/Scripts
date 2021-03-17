#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Feb  8 11:15:29 2021

@author: asif
"""

import numpy as np
import matplotlib.pyplot as plt
from astropy.modeling import models,fitting
from scipy.optimize import curve_fit,least_squares

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

def f_nu(t,nu,F_0,m_f,nu_b0,m_b,alpha_1,alpha_2,delta_t):
    # if not isinstance(t, np.ndarray):
    #     t=np.array([t])
    n=len(t)
    F_0t=F_0_evln(t,F_0,m_f)
    nu_bt=nu_break_evln(t,nu_b0,m_b)
    alpha_1t=alpha_1*np.ones(n)
    alpha_2t=alpha_2*np.ones(n)
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

def calc_phIdx_evln(t,nu,F_0,m_f,nu_b0,m_b,alpha_1t,alpha_2t,delta_t):
    flux=f_nu(t,nu,F_0,m_f,nu_b0,m_b,alpha_1t,alpha_2t,delta_t)
    s=IdxEvln(flux,nu)
    ph=-s[:,0]+1
    return ph

##################################Read light curve data#######################
GRB="190114"
lc_file_name="flux_"+"GRB "+ GRB+".txt"
lc_file=np.loadtxt(lc_file_name)
t_data=lc_file[:,0]
lc_data=lc_file[:,3]
lc_err=2*lc_file[:,4] 
##############################################################################


def model_f(t,F_0,m_f,nu_b0,m_b,alpha_1,alpha_2):
    """Calculates model flux for given set of parameters"""
    if not isinstance(t, np.ndarray):
        t=np.array([t])
    nu=np.array([2.4e17*10])#10 kev
    n=len(t)
    delta_t=0.4*np.ones(n)
    F= f_nu(t,nu,F_0,m_f,nu_b0,m_b,alpha_1,alpha_2,delta_t)
    return F[:,0]

def model_ph(t,F_0,m_f,nu_b0,m_b,alpha_1,alpha_2):
    """Calculates photon index as observed in window of freqs nu as per model model_f"""
    if not isinstance(t, np.ndarray):
        t=np.array([t])
    nu=2.4e17*np.logspace(-0.5,1.,100)#XRT range
    n=len(t)
    delta_t=0.4*np.ones(n)
    ph_idx_t=calc_phIdx_evln(t,nu,F_0,m_f,nu_b0,m_b,alpha_1,alpha_2,delta_t)
    
    return ph_idx_t

parameters = ["F_0","beta1","beta2","p_F","p_b","nu_b0"]
n_params = len(parameters)
# t=1e6
# F_0=1e4
# m_f=2
# nu_b0=1e28
# m_b=2
# delta=0.4
# alpha_1=1
# alpha_2=3

# ph=model_ph(t,F_0,m_f,nu_b0,m_b,alpha_1,alpha_2)
LB=[0,0,0,0,0,0] #Lower bound for parameters
UB=[np.inf,10,np.inf,10,10,10] #Upperbound for parameters


# popt, pcov = curve_fit(model_ph, t_data, ph_data)
popt, pcov = curve_fit(model_f, t_data, lc_data,bounds=(LB,UB),method='trf')

nu=2.4e17*np.logspace(-0.5,1.,100)
t_plot=np.logspace(4,7)

n=len(t_plot)
lc_fit=model_f(t_plot,*popt)
F_0,m_f,nu_b0,m_b,alpha_1,alpha_2=popt
delta_t=0.4*np.ones(n)
alpha_1t=alpha_1*np.ones(n)
alpha_2t=alpha_2*np.ones(n)
# f_n=f_nu(t_plot,nu,F_0,m_f,nu_b0,m_b,alpha_1t,alpha_2t,delta_t)

plt.figure()

plt.scatter(t_data,lc_data,label="data")
plt.plot(t_plot,lc_fit,label="fit")
plt.legend()
plt.xscale("log")
plt.yscale("log")
