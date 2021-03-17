#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Jan 20 16:09:59 2021

@author: asif
"""

#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Jan 19 16:12:09 2021

@author: asif
"""

import json
import numpy
import numpy as np
# import scipy.stats, scipy
import pymultinest
import matplotlib.pyplot as plt
import os


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

##################################Read light curve data########################

GRB="080207"
lc_file_name="flux_"+"GRB "+ GRB+".txt"
lc_file=np.loadtxt(lc_file_name)
t_data=lc_file[:,0]
lc_data=lc_file[:,3]
lc_err=lc_file[:,4] +lc_file[:,5]
nu=2.4e17*np.logspace(-0.5,1.,100)

################################Multinest Stuff################################

parameters = ["beta1","beta2","p_F","p_b","F_0",r"$nu_{b_0}$"]
n_params = len(parameters)

def model_f(t,nu,beta1,beta2,p_F,p_b,F_0,nu_b0,s=2.2):
    """Calculates model flux for given set of parameters"""
    F_0t=F_0_evln(t,F_0,p_F)
    nu_bt=nu_break(t,nu_b0,p_b)
    nu=10*2.4e17 #10 Kev To Hz
    F=F_arb(nu,nu_bt,F_0t,beta1,beta2)
    return F

# def model_ph(t,nu,beta1,beta2,p_F,p_b,F_0,nu_b0,s=2.2):
#     """Calculates photon index as observed in window of freqs nu as per model model_f"""
#     F=model_f(t,nu,beta1,beta2,p_F,p_b,F_0,nu_b0)
#     sp_idx_m=-IdxEvln(F,nu)[:,0]
#     ph_idx_m=sp_idx_m+1
#     return ph_idx_m

def prior(cube):
    """Priors for parameters are beta1,beta2,p_F,p_b,F_0,nu_b0 and possibly s"""
    cube[0] = cube[0]*-15          #beta1 =uniform prior between 0:-10
    cube[1] = cube[1]*-15           #beta2 = uniform prior between 0:-10
    cube[2] = cube[2]*-10             #P_F  =uniform prior between 0:-10
    cube[3] = cube[3]*-10           #P_b = uniform prior between 0:-10
    cube[4] = 10**(cube[4]*15 - 5) # F_0 = log-uniform prior between 10^-5 and 10^5
    cube[5] = 10**(cube[5]*10 + 20) # nu_b0=log-uniform prior between 10^10 and 10^30
    # if ndim < 7:
    #     return
    #    	cube[7] =cube[7]*4# s=uniform prior between 0 and 4
    return cube

def loglike(cube):
    parms=np.ones(n_params)
    for i in range(n_params):
      parms[i] = cube[i]
    beta1,beta2,p_F,p_b,F_0,nu_b0=parms
    lc_model=model_f(t_data,nu,beta1,beta2,p_F,p_b,F_0,nu_b0)
    # ph_model=model_ph(t_data,nu,beta1,beta2,p_F,p_b,F_0,nu_b0)
    # likelihood = ((-0.5*(ph_model- ph_data)/ph_err)**2).sum()
    likelihood = ((-0.5*(lc_model- lc_data))**2).sum()
    # likelihood = exp(-0.5 * ((pos - positions) / width)**2) / (2*pi*width**2)**0.5
    return np.log(likelihood)

if not os.path.exists('out'):
    os.makedirs('out')

# pymultinest.run(loglike, prior, n_params, outputfiles_basename='out/'+GRB + '_fit_',
# 	resume = False, verbose = True)
result=pymultinest.solve(LogLikelihood=loglike,Prior=prior,n_dims= n_params, outputfiles_basename='out/'+GRB + '_fit_',
	resume=False,verbose = True)
json.dump(parameters, open('out/'+GRB + '_fit_'+'params.json', 'w')) 

###################Viewing the results########################################
print()
print('evidence: %(logZ).1f +- %(logZerr).1f' % result)
print()
print('parameter values:')
for name, col in zip(parameters, result['samples'].transpose()):
	print('%15s : %.3f +- %.3f' % (name, col.mean(), col.std()))
with open('GRB %sparams.json' % GRB, 'w') as f:
	json.dump(parameters, f, indent=2)

