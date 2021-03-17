#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Jan 21 20:07:25 2021

@author: asif
"""

import numpy as np
import pymultinest
import matplotlib.pyplot as plt
import os
import json

#The idea is to simulate data of a wave with some noise thrown in and run
#MCMC using pymultinest to see what comes out
####################Generating some simulated data of a wave###################

np.random.seed(5)
noise=0.5
A=3
omega=8
phi=np.pi/3
t=np.linspace(0,10)
y_model=A*np.cos(omega*t+phi)
y_data=np.random.normal(y_model,noise)
plt.figure()
plt.plot(t,y_model)
plt.scatter(t,y_data,color='r')

###############################MCMC############################################

params=['A','omega','phi']
n_params=len(params)

def prior(cube):
    """Priors for parameters are A,omega and phi"""
    cube[0] = cube[0]*60       #A =uniform prior between 0:10
    cube[1] = cube[1]*10     #omega = uniform prior between 0:10
    cube[2] = cube[2]*2*np.pi # phi=uniform prior between 0 and 2pi
    return cube

def model(t,A,omega,phi):
    return A*np.cos(omega*t+phi)

def loglike(cube):
    parms=np.ones(n_params)
    for i in range(n_params):
      parms[i] = cube[i]
    A,omega,phi=parms
    y_model1=model(t,A,omega,phi)
    likelihood=(-0.5*((y_model1- y_data)/noise)**2).sum()
    return likelihood

if not os.path.exists('out'):
    os.makedirs('out')

GRB='Wave'
result=pymultinest.solve(LogLikelihood=loglike,Prior=prior,n_dims= n_params, outputfiles_basename='out/'+GRB + '_fit_',
	resume=False,verbose = True)
# json.dump(params, open('out/'+GRB + '_fit_'+'params.json', 'w')) 
#################################plot stuff####################################
print()
print('evidence: %(logZ).1f +- %(logZerr).1f' % result)
print()
print('parameter values:')

par_fit=[]
for name, col in zip(params, result['samples'].transpose()):
    par_fit.append([col.mean(), col.std()])
    print('%15s : %.3f +- %.3f' % (name, col.mean(), col.std()))
with open('GRB %sparams.json' % GRB, 'w') as f:
	json.dump(params, f, indent=2)

t_fit=np.linspace(0,10,200)
par_fit=np.array(par_fit)
A,omega,phi=par_fit[:,0]
# t_plot=np.logspace(4,5)
y_fit=model(t_fit,A,omega,phi)
    

# plt.scatter(t_data,y_data)
plt.plot(t_fit,y_fit,color='g')
# plt.scatter(t_plot,ph_fit)
# plt.xscale("log")
# plt.yscale("log")
    
