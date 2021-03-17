#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Feb  2 22:11:29 2021

@author: asif
"""

import numpy
import numpy as np
# import scipy.stats, scipy
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit,least_squares

def F_arb(nu,nu_b,F_0,beta1,beta2,s):
    """ retruns F[t_i,nu_i] broken power law. Give nu_b and F_0 as lists even 
        if single value"""
    F=[]
    for Fi,nu_bi in zip(F_0,nu_b):        
        temp=np.power(nu/nu_bi,-beta1*s)+np.power(nu/nu_bi,-beta2*s)
        F.append(Fi*np.power(temp,-1/s))
    return np.array(F)

def nu_break(t,nu_0,m=-2,t0=1):
    return nu_0*np.power(t/t0,m)

def F_0_evln(t,F_0,m=-2,t0=1):
    return F_0*np.power(t/t0,m)

def IdxEvln(F,freqs):
    "Fits a power law to F[i1,i2] for all i1 with i2 being freqs"
    
    s_idx=[]
    for ff in F:
        model_x=np.log(freqs)
        model_y=np.log(ff)
        temp = np.polyfit(model_x,model_y, 1)
        s_idx.append(temp)
    return np.array(s_idx)

##################################Read light curve data########################

GRB="190114"
lc_file_name="flux_"+"GRB "+ GRB+".txt"
lc_file=np.loadtxt(lc_file_name)
t_data=lc_file[:,0]
lc_data=lc_file[:,3]
lc_err=2*lc_file[:,4] 
nu=2.4e17*np.logspace(-0.5,1.,100)

################################Multinest Stuff################################

# parameters = ["beta1","beta2","p_F","p_b","F_0",r"$nu_{b_0}$"]
# n_params = len(parameters)

# def model_f(t,nu,beta1,beta2,p_F,p_b,F_0,nu_b0,s=2.2):
#     """Calculates model flux for given set of parameters"""
#     F_0t=F_0_evln(t,F_0,p_F)
#     nu_bt=nu_break(t,nu_b0,p_b)
#     nu=10*2.4e17 #10 Kev To Hz
#     F=F_arb(nu,nu_bt,F_0t,beta1,beta2,s)
#     return F

def calc_brokenPL(x,A,x_b,beta_1,beta_2):
    """A broken power law"""
    F=np.zeros(len(x))
    mask1= x<=x_b
    mask2= x>x_b
    if (x[mask1].size>0):
        F[mask1]=A*np.power((x/x_b)[mask1],-beta_1)
    if (x[mask2].size>0):
        F[mask2]=A*np.power((x/x_b)[mask2],-beta_2)
    return F
def PL(x,A,beta_1):
    """A power law"""
    F=A*np.power((x),-beta_1)

    return F

# def calc_brokenPL(x,A,x_b,beta_1,beta_2):
#     """A broken power law"""
#     if x<=x_b:
#         F=A*(x/x_b)**(-beta_1)
#     if x>x_b:
#         F=A*(x/x_b)**(-beta_1)
#     return F

# popt, pcov = curve_fit(calc_brokenPL, t_data, lc_data,bounds=[(0,0,0,0),(np.inf,np.inf,10,10)],maxfev=100000)
popt, pcov = curve_fit(calc_brokenPL, t_data, lc_data)
# popt, pcov = curve_fit(PL, t_data, lc_data,maxfev=10000)#,bounds=[(0,0,0),(np.inf,np.inf,10)])
# popt, pcov = curve_fit(PL, t_data, lc_data,maxfev=100000,bounds=[(0,0),(np.inf,10)])

t_plot=np.logspace(4,7)
lc_fit=calc_brokenPL(t_plot,*popt)
# lc_fit=PL(t_plot,*popt)



plt.figure()
plt.scatter(t_data,lc_data)
plt.plot(t_plot,lc_fit)
plt.xscale("log")
plt.yscale("log")
