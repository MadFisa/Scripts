#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Feb 11 12:28:26 2021

@author: asif
"""

import numpy as np
import matplotlib.pyplot as plt
from astropy.modeling import models,fitting
from scipy.optimize import curve_fit,least_squares
plt.style.use('seaborn-whitegrid')

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
GRB="130907"
lc_file_name="flux_"+"GRB "+ GRB+".txt"
lc_file=np.loadtxt(lc_file_name)
t_data=lc_file[:,0]
lc_data=lc_file[:,3]
lc_err=lc_file[:,4] 
log_lc_err=(lc_err/lc_data)
##############################################################################

def model_f(log_t,log_F_0,m_f,log_nu_b0,m_b,alpha_1,alpha_2):
    """Calculates model flux for given set of parameters"""
    if not isinstance(log_t, np.ndarray):
        log_t=np.array([log_t])
    nu=np.array([2.4e17*10])#10 kev
    n=len(log_t)
    delta_t=0.4*np.ones(n)
    F_0=np.power(10.,log_F_0)
    nu_b0=np.power(10.,log_nu_b0)
    F= f_nu(np.power(10.,log_t),nu,F_0,m_f,nu_b0,m_b,alpha_1,alpha_2,delta_t)
    # print (F)
    return np.log10(F[:,0])
def model_f1(log_t,log_F_0,log_nu_b0,m_b,alpha_1,alpha_2):
    """Calculates model flux for given set of parameters"""
    if not isinstance(log_t, np.ndarray):
        log_t=np.array([log_t])
    nu=np.array([2.4e17*10])#10 kev
    n=len(log_t)
    delta_t=0.4*np.ones(n)
    F_0=np.power(10.,log_F_0)
    nu_b0=np.power(10.,log_nu_b0)
    F= f_nu(np.power(10.,log_t),nu,F_0,0,nu_b0,m_b,alpha_1,alpha_2,delta_t)
    # print (F)
    return np.log10(F[:,0])

# LB=[-20,0,0,0,-10,-10] #Lower bound for parameters
# UB=[20,10,np.inf,10,10,10] #Upperbound for parameters
LB=[-20,0,0,-10,-10] #Lower bound for parameters
UB=[20,np.inf,10,10,10] #Upperbound for parameters
nu=np.array([2.4e17*10])
nu_b0=2.4e17*10*1e5

lg_F_0=-8.
lg_nu_b0=np.log10(nu_b0)
m_f=0
m_b=1.
alpha_1=4.
alpha_2=2.
t_b=np.power(nu_b0/nu, 1/m_b)
pars=[lg_F_0,lg_nu_b0,m_b,alpha_1,alpha_2]
# pp=[ 6.2975156 ,  2.88708825, 30.61078379,  2.61705714,  1.38926202,
#         -0.39287496]
pp=[ 6.2975156 , 30.61078379,  2.61705714,  1.38926202,
        -0.39287496]
# popt2, pcov2 = curve_fit(model_f, np.log10(t_data), np.log10(lc_data),p0=pars,bounds=(LB,UB))
# popt, pcov = curve_fit(model_f, np.log10(t_data), np.log10(lc_data),p0=(-5,2,15,2,3,4),bounds=(LB,UB))
# popt, pcov = curve_fit(model_f, np.log10(t_data), np.log10(lc_data),bounds=(LB,UB))
# popt, pcov = curve_fit(model_f, np.log10(t_data), np.log10(lc_data),p0=[lg_F_0,0.,lg_nu_b0,m_b,alpha_1,alpha_2])
popt, pcov = curve_fit(model_f1, np.log10(t_data), np.log10(lc_data),p0=pp,bounds=(LB,UB))

# popt, pcov = curve_fit(model_f,np.log10(t_data), np.log10(lc_data),sigma=log_lc_err,p0=(-5,2,15,2,3,4),bounds=(LB,UB))
popt2, pcov2 = curve_fit(model_f1,np.log10(t_data), np.log10(lc_data),sigma=log_lc_err,p0=(-5,15,2,3,4),bounds=(LB,UB))
# popt, pcov = curve_fit(model_f,np.log10(t_data), np.log10(lc_data),sigma=np.log10(lc_err))
# nu=2.4e17*np.logspace(-0.5,1.,100)

t_plot=np.linspace(4,6,200)
# n=len(t_plot)


# popt2, pcov2 = curve_fit(model_f1, np.log10(t_data), np.log10(lc_data))
# popt2, pcov2 = curve_fit(model_f1, np.log10(t_data), np.log10(lc_data),p0=pars)
lc_fit2=model_f1(t_plot,*popt2)
# lc_fit1=model_f1(t_plot,*pars)
lc_fit=model_f1(t_plot,*popt)




plt.figure(figsize=(20,10))

# plt.title("[lg_F_0,lg_nu_b0,m_b,alpha_1,alpha_2]="+ str(pars))
plt.errorbar(np.log10(t_data),np.log10(lc_data),fmt='.k',yerr=log_lc_err,ecolor='lightgray', elinewidth=3, capsize=0,label="data")
plt.plot(t_plot,lc_fit,label="fit")
# plt.plot(t_plot,lc_fit1,label="manual_fit")
plt.plot(t_plot,lc_fit2,label="p_fit")
plt.legend()
plt.show()
plt.savefig("manual fit GRB"+GRB+".png")

