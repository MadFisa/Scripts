#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Dec 11 09:51:07 2020

@author: asif
"""
import numpy as np
import matplotlib.pyplot as plt

# def F_v(nu,nu_b,F_0=1.,beta1=3.,beta2=2.,s=2.):
#     F=np.power(nu/nu_b,(beta1-beta2)*s)
#     return F_0*np.power(1+F,-1/s)

def F_v(nu,nu_b,F_0=1.,beta1=3.,beta2=2.,s=2.):
    """F[ti,nui]"""
    F=[]
    for nu_bi in nu_b:        
        temp=np.power(nu/nu_bi,-beta1*s)+np.power(nu/nu_bi,-beta2*s)
        F.append(F_0*np.power(temp,-1/s))
    return np.array(F)

def nu_break(t,nu_0,t0=1,m=-2):
    return nu_0*np.power(t/t0,m)

nu=np.logspace(1,30,1000)
t=np.logspace(0,7,100)
nu_b=nu_break(t, 1e15)
beta1=-1/3
beta2=-1/2
s=4

F=F_v(nu,nu_b,1,beta1,beta2,s)  
plt.loglog(nu,F[10,:])
plt.loglog(nu,F[20,:])
plt.loglog(nu,F[30,:])

plt.figure()
plt.loglog(t,F[:,278])
