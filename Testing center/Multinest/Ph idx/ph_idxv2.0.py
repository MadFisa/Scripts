#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Jan 19 16:12:09 2021

@author: asif
"""

import json
import numpy
import numpy as np
import scipy.stats, scipy
import pymultinest
import matplotlib.pyplot as plt

def F_arb(nu,nu_b,F_0,beta1=3.,beta2=2.,s=2.):
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


###############################Reading data from ph_idx files##############################

###############################Start of multinest stuff##########################################


def model_f(t,nu,beta1,beta2,p_F,p_b,F_0,nu_b0,s=2.2):
    """Calculates model flux for given set of parameters"""
    F_0t=F_0_evln(t,F_0,p_F)
    nu_bt=nu_break(t,nu_b0,p_b)
    F=F_arb(nu,nu_bt,F_0t,beta1,beta2)
    return F

def model_ph(t,nu,beta1,beta2,p_F,p_b,F_0,nu_b0,s=2.2):
    F=model_f(t,nu,beta1,beta2,p_F,p_b,F_0,nu_b0)
    sp_idx_m=-IdxEvln(F,nu)[:,0]
    ph_idx_m=sp_idx_m+1
    return ph_idx_m

def prior(cube, ndim, nparams):
    """parameters are beta1,beta2,p_F,p_b,F_0,nu_b0 and possibly s"""
    cube[0] = cube[0]*-10            #beta1 =uniform prior between 0:-10
    cube[1] = cube[2]*-10            #beta2 = uniform prior between 0:-10
    cube[2] = cube[3]*-10            #P_F  =uniform prior between 0:-10
    cube[3] = cube[4]*-10            #P_b = uniform prior between 0:-10
    cube[4] = 10**(cube[5]*10 - 5) # F_0 = log-uniform prior between 10^-5 and 10^5
    cube[5] = 10**(cube[2]*10 + 20) # nu_b0=log-uniform prior between 10^10 and 10^30
    if ndim < 7:
        return
       	cube[7] =cube[3]*4# s=uniform prior between 0 and 4

def loglike(cube, ndim, nparams):
	# get the current parameter (is between 0:2 now)
    t,nu,beta1,beta2,p_F,p_b,F_0,nu_b0=cube
    ph_model=model_f(t_idx,nu,beta1,beta2,p_F,p_b,F_0,nu_b0)
    likelihood = (((ph_model- ph_data)/ph_err)**2).sum()
    # likelihood = exp(-0.5 * ((pos - positions) / width)**2) / (2*pi*width**2)**0.5
    return np.log(likelihood)

parameters = ["beta1","beta2","p_F","p_b","F_0",r"$nu_{b_0}$"]
n_params = len(parameters)

pymultinest.run(loglike, prior, n_params, outputfiles_basename='out/',
	resume = False, verbose = True)
json.dump(parameters, open('out/params.json', 'w')) 

t=np.logspace(2,6,1000)



bu=1e6
beta1=-1.5
beta2=-2.5
p_F=-1.5
p_b=-0.5
F_0=1
nu_b0=1e23
ph=model_ph(t,nu,beta1,beta2,p_F,p_b,F_0,nu_b0)
ff=model_f(t,nu,beta1,beta2,p_F,p_b,F_0,nu_b0)

plt.figure()
plt.loglog(t,ff)
plt.figure()
plt.loglog(t,ph)
# plt.xscale('log')
    
    