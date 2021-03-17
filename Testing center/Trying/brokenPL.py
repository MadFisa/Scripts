#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sat Jan 23 21:26:38 2021

@author: asif
"""

import numpy as np
import matplotlib.pyplot as plt
import os
import json
import pymultinest

def calc_brokenPL(x,A,x_b,beta_1,beta_2):
    """A broken power law"""
    F=np.zeros(len(x))
    mask1= x<=x_b
    mask2= x>x_b
    F[mask1]=A*np.power((x/x_b)[mask1],beta_1)
    F[mask2]=A*np.power((x/x_b)[mask2],beta_2)
    return F

def nu_break(t,nu_0,m=-2,t0=1):
    return nu_0*np.power(t/t0,m)

def F_0_evln(t,F_0,m=-2,t0=1):
    return F_0*np.power(t/t0,m)




##################################Read light curve data########################

GRB="111209"
lc_file_name="flux_"+"GRB "+ GRB+".txt"
lc_file=np.loadtxt(lc_file_name)
t_data=lc_file[:,0]
lc_data=lc_file[:,3]
lc_err=2*lc_file[:,4] 
nu=2.4e17*np.logspace(-0.5,1.,100)

# nu=np.logspace(2,7,200)
# A=3.
# nu_b=1e4
# F=calc_brokenPL(nu,A,nu_b,-2,-4)
# plt.loglog(nu,F)
################################Multinest Stuff################################

parameters = ["F_0","t_b","beta_1","beta_2"]
n_params = len(parameters)



def prior(cube):
    """Priors for parameters are F_0,t_b,beta_1,beta_2 """
    cube[0] =10**(cube[0]*4 - 2)      
    cube[1] = 10**(cube[1]*5 +2)     
    cube[2]=-cube[2]*10 
    cube[3]=-cube[3]*10 
    return cube

def loglike(cube):
    global p
    parms=np.ones(n_params)
    for i in range(n_params):
      parms[i] = cube[i]
      
    F_0,t_b,beta_1,beta_2=parms
    lc_model=calc_brokenPL(t_data,F_0,t_b,beta_1,beta_2)
    # ph_model=model_ph(t_data,nu,beta1,beta2,p_F,p_b,F_0,nu_b0)
    # likelihood = ((-0.5*(ph_model- ph_data)/ph_err)**2).sum()
    loglikelihood = (-0.5*( (lc_model- lc_data))**2).sum()
    loglikelihood=loglikelihood - 0.5*sum(np.log(2.*np.pi*lc_err**2))
    # likelihood = exp(-0.5 * ((pos - positions) / width)**2) / (2*pi*width**2)**0.5
    return loglikelihood

dire="PL"
if not os.path.exists(dire):
    os.makedirs(dire)

# pymultinest.run(loglike, prior, n_params, outputfiles_basename='out/'+GRB + '_fit_',
# 	resume = False, verbose = True)
n_live=500
tol=0.3
result=pymultinest.solve(LogLikelihood=loglike,Prior=prior,n_dims= n_params, outputfiles_basename=dire+'/'+GRB + '_fit_',
	resume=False,verbose = True,n_live_points=n_live,sampling_efficiency=0.3)
json.dump(parameters, open(dire+'/'+GRB + '_fit_'+'params.json', 'w')) 

print()
print('evidence: %(logZ).1f +- %(logZerr).1f' % result)
print()
print('parameter values:')
par_fit=[]
for name, col in zip(parameters, result['samples'].transpose()):
    par_fit.append([col.mean(), col.std()])
    print('%15s : %.3e +- %.3e' % (name, col.mean(), col.std()))
with open('GRB %sparams.json' % GRB, 'w') as f:
	json.dump(parameters, f, indent=2)

par_fit=np.array(par_fit)
F_0,t_b,beta_1,beta_2=par_fit[:,0]
t_plot=np.logspace(4,7)
lc_fit=calc_brokenPL(t_plot,F_0,t_b,beta_1,beta_2)



plt.figure()
plt.scatter(t_data,lc_data)
plt.plot(t_plot,lc_fit)
plt.xscale("log")
plt.yscale("log")

