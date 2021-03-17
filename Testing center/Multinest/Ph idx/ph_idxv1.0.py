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
    """ retruns F[t_i,nu_i] broken power law"""
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


def model_f(t,nu,beta1,beta2,p_F,p_b,F_0,nu_b0,s=2.2):
    F_0t=F_0_evln(t,F_0,p_F)
    nu_bt=nu_break(t,nu_b0,p_b)
    F=F_arb(nu,nu_bt,F_0t,beta1,beta2)
    return F

t=np.logspace(2,6)
nu=1e15


bu=1e6
beta1=-1.5
beta2=-2.5
p_F=-1.5
p_b=-0.5
F_0=1e4
nu_b0=1e23
ff=model_f(t,nu,beta1,beta2,p_F,p_b,F_0,nu_b0)

plt.figure()
plt.loglog(t,ff)
# plt.xscale('log')
    
    