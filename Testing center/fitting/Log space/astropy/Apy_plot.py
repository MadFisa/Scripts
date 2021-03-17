#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Feb 12 16:36:10 2021

@author: asif
"""

import numpy as np
import matplotlib.pyplot as plt
from astropy.modeling import models
# from scipy.optimize import curve_fit,least_squares
from mpl_toolkits.mplot3d import axes3d
from matplotlib import cm

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

def model_f(log_t,log_nu,log_F_0,m_f,log_nu_b0,m_b,alpha_1,alpha_2):
    """Calculates model flux for given set of parameters"""
    if not isinstance(log_t, np.ndarray):
        log_t=np.array([log_t])
    # nu=np.array([2.4e17*10])#10 kev
    n=len(log_t)
    delta_t=0.4*np.ones(n)
    F_0=np.power(10.,log_F_0)
    nu_b0=np.power(10.,log_nu_b0)
    F= f_nu(np.power(10.,log_t),np.power(10,log_nu),F_0,m_f,nu_b0,m_b,alpha_1,alpha_2,delta_t)
    # print (F)
    return np.log10(F)

# def model_f1(log_t,log_F_0,log_nu_b0,m_b,alpha_1,alpha_2):
#     """Calculates model flux for given set of parameters"""
#     if not isinstance(log_t, np.ndarray):
#         log_t=np.array([log_t])
#     nu=np.array([2.4e17*10])#10 kev
#     n=len(log_t)
#     delta_t=0.4*np.ones(n)
#     F_0=np.power(10.,log_F_0)
#     nu_b0=np.power(10.,log_nu_b0)
#     F= f_nu(np.power(10.,log_t),nu,F_0,0,nu_b0,m_b,alpha_1,alpha_2,delta_t)
#     # print (F)
#     return np.log10(F[:,0])

nu_b0=2.4e17*10*1e5
lg_F_0=-8.
lg_nu_b0=np.log10(nu_b0)
m_f=0
m_b=1.
alpha_1=4.
alpha_2=2.
# t_b=np.power(nu_b0/nu, 1/m_b)
pars=[lg_F_0,m_f,lg_nu_b0,m_b,alpha_1,alpha_2]

t_plot=np.linspace(4,7,200)
nu_plot=np.linspace(np.log10(2.4e17),np.log10(2.4e17*10),200)
T, Nu = np.meshgrid(t_plot, nu_plot)

F=model_f(t_plot,nu_plot,*pars)

fig = plt.figure()
ax = fig.add_subplot(111, projection='3d')
# ax.plot_wireframe(T, Nu, F.T, rstride=10, cstride=10,color='green')
# ax.plot_surface(T, Nu, F.T, rstride=10, cstride=10,cmap='winter', edgecolor='none')
# ax.plot_surface(T, Nu, F.T, rstride=10, cstride=10,color='b')
# ax.plot_surface(T, Nu, F.T, cmap=cm.coolwarm,linewidth=0, antialiased=False)
ax.plot_surface(T, Nu, F.T, rstride=8, cstride=8, alpha=0.3)
cset = ax.contour(T, Nu, F.T, zdir='z', offset=-100, cmap=cm.coolwarm)
cset = ax.contour(T, Nu, F.T, zdir='x', offset=t_plot[0], cmap=cm.coolwarm)
cset = ax.contour(T, Nu, F.T, zdir='y', offset=nu_plot[0], cmap=cm.coolwarm)
# cset = ax.contour(T, Nu, F.T, zdir='z', offset=-100, cmap=cm.coolwarm)
# cset = ax.contour(T, Nu, F.T, zdir='x', offset=-40, cmap=cm.coolwarm)
# cset = ax.contour(T, Nu, F.T, zdir='y', offset=40, cmap=cm.coolwarm)

ax.set_xlabel('T')
# ax.set_xlim(-40, 40)
ax.set_ylabel('Nu')
# ax.set_ylim(-40, 40)
ax.set_zlabel('F')
# ax.set_zlim(-100, 100)

# plt.figure(figsize=(20,10))
# plt.plot(t_plot,lc,label="fit")
# plt.legend()
# plt.show()

