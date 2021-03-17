#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sat Feb 13 22:01:20 2021

@author: asif
"""

import numpy as np
import matplotlib.pyplot as plt
from astropy.modeling import models
import os
import pymultinest
import json

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
lc_err=2*lc_file[:,4] 
log_tdata=np.log10(t_data)
log_lc_data=np.log10(lc_data)
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


# def prior(cube):
#     """Priors for parameters are beta1,beta2,p_F,p_b,F_0,nu_b0 and possibly s"""
#     cube[0] = cube[0]*15 - 5   #F_0
#     cube[1] = cube[1]*-10  + 5     #m_f
#     cube[2] = cube[2]*40   # nu_b
#     cube[3] = cube[3]*-20 + 10         # m_b 
#     cube[4] =  cube[4]*-20 + 10 # alpha_1
#     cube[5] =  cube[5]*-20 + 10# alpha_2
#     # if ndim < 7:
#     #     return
#     #    	cube[7] =cube[7]*4# s=uniform prior between 0 and 4
#     return cube
def prior(cube):
    """Priors for parameters are beta1,beta2,p_F,p_b,F_0,nu_b0 and possibly s"""
    cube[0] = -cube[0]*10 + 10  #F_0
    # cube[1] = cube[1]*-10  + 5     #m_f
    cube[1] = cube[1]*20 + 20   # nu_b
    cube[2] = cube[2]*5         # m_b 
    cube[3] =  cube[3]*5 # alpha_1
    cube[4] =  cube[4]*5# alpha_2
    # if ndim < 7:
    #     return
    #    	cube[7] =cube[7]*4# s=uniform prior between 0 and 4
    return cube



# def loglike(cube):
#     parms=np.ones(n_params)
#     for i in range(n_params):
#       parms[i] = cube[i]
#     log_F_0,m_f,log_nu_b0,m_b,alpha_1,alpha_2=parms
#     lc_model=model_f(log_tdata,log_F_0,m_f,log_nu_b0,m_b,alpha_1,alpha_2)
#     # likelihood = ((-0.5*(lc_model- lc_data)/lc_err)**2).sum()
#     likelihood = ((-0.5*(lc_model- log_lc_data))**2).sum()
#     # likelihood = exp(-0.5 * ((pos - positions) / width)**2) / (2*pi*width**2)**0.5
#     return likelihood
def loglike(cube):
    # parms=np.ones(n_params)
    # for i in range(n_params):
    #   parms[i] = cube[i]
    # log_F_0,log_nu_b0,m_b,alpha_1,alpha_2=parms
    lc_model=model_f1(log_tdata,*cube)
    likelihood = -0.5*(((lc_model- lc_data)/lc_err)**2).sum()
    # likelihood = ((-0.5*(lc_model- log_lc_data))**2).sum()
    # likelihood = exp(-0.5 * ((pos - positions) / width)**2) / (2*pi*width**2)**0.5
    return likelihood


# parameters = ["F_0","m_f","nu_b","m_b","alpha_1","alpha_2"]
parameters = ["F_0","nu_b","m_b","alpha_1","alpha_2"]
n_params = len(parameters)
# t=1e6
# F_0=1e4
# m_f=2
# nu_b0=1e28
# m_b=2
# delta=0.4
# alpha_1=1
# alpha_2=3

# ph=model_f(t,F_0,m_f,nu_b0,m_b,alpha_1,alpha_2)
LB=[0,0,0,0,0,0] #Lower bound for parameters
UB=[np.inf,10,np.inf,10,10,10] #Upperbound for parameters

if not os.path.exists('out'):
    os.makedirs('out')

# pymultinest.run(loglike, prior, n_params, outputfiles_basename='out/'+GRB + '_fit_',
# 	resume = False, verbose = True)
result=pymultinest.solve(LogLikelihood=loglike,Prior=prior,n_dims= n_params, outputfiles_basename='out/'+GRB + '_fit_',
	resume=True,verbose = True)
json.dump(parameters, open('out/'+GRB + '_fit_'+'params.json', 'w')) 

############################Plotting and stuff###################################################
# plot the distribution of a posteriori possible models
# plt.figure() 
# plt.scatter(t_data,lc_data)
# # plt.plot(x, ydata, '+ ', color='red', label='data')
# a = pymultinest.Analyzer(outputfiles_basename='out/'+GRB + '_fit_', n_params = n_params)
# for (beta1,beta2,p_F,p_b,F_0,nu_b0) in a.get_equal_weighted_posterior()[::100,:-1]:
# 	plt.plot(t_data, model_f(t_data,nu,beta1,beta2,p_F,p_b,F_0,nu_b0), '-', color='blue', alpha=0.3, label='data')

# plt.xscale('log')
# plt.savefig('out/'+GRB + '_fit_posterior.png')

# a_lnZ = a.get_stats()
##############################################################################################
print()
print('evidence: %(logZ).1f +- %(logZerr).1f' % result)
print()
print('parameter values:')

par_fit=[]
for name, col in zip(parameters, result['samples'].transpose()):
    par_fit.append([col.mean(), col.std()])
    print('%15s : %.3f +- %.3f' % (name, col.mean(), col.std()))
with open('GRB %sparams.json' % GRB, 'w') as f:
	json.dump(parameters, f, indent=2)

par_fit=np.array(par_fit)
# log_t,log_F_0,m_f,log_nu_b0,m_b,alpha_1,alpha_2=par_fit[:,0]
# t_plot=np.logspace(4,5)
lc_fit=model_f1(log_tdata,*par_fit[:,0])
    
plt.figure()
plt.scatter(t_data,lc_data)
plt.scatter(t_data,lc_fit)
# plt.scatter(t_plot,lc_fit)
plt.xscale("log")
# plt.yscale("log")