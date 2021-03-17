#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sun Feb  7 22:58:25 2021

@author: asif
"""

import numpy as np
import matplotlib.pyplot as plt
from astropy.modeling import models,fitting
def nu_break_evln(t,nu_0,m=2,t0=1):
    return nu_0*np.power(t/t0,-m)

def F_0_evln(t,F_0,m=2,t0=1):
    return F_0*np.power(t/t0,-m)

t = np.logspace(0.7, 6.3, 500)
n=len(t)
F_0t=F_0_evln(t,1e4)
nu_bt=nu_break_evln(t,1e10,m=2)
delta_t=0.4*np.ones(n)
alpha_1t=1*np.ones(n)
alpha_2t=3*np.ones(n)

nu=2.4e17*np.logspace(-0.5,1.,100)

# f = models.SmoothlyBrokenPowerLaw1D(amplitude=1, x_break=2000,
                                    # alpha_1=1, alpha_2=4)
flux=np.zeros((n,len(nu)))                                    
for i in range(n):
    f=models.SmoothlyBrokenPowerLaw1D(amplitude=F_0t[i], x_break=nu_bt[i],
                                    alpha_1=alpha_1t[i], alpha_2=alpha_2t[i])
    flux[i]=f(nu)



plt.figure()
plt.title("amplitude=1, x_break=20, alpha_1=-2, alpha_2=2")

# f.delta = 0.5
plt.loglog(t, flux, '--', label='delta=0.5')

# f.delta = 0.3
# plt.loglog(x, f(x), '-.', label='delta=0.3')

# f.delta = 0.1
# plt.loglog(x, f(x), label='delta=0.1')

# plt.axis([x.min(), x.max(), 0.1, 1.1])
# plt.legend(loc='lower center')
# plt.grid(True)
# plt.show()
