# -*- coding: utf-8 -*-
"""
Created on Fri Aug 21 02:14:53 2020

@author: Asif
"""

import numpy as np
import matplotlib.pyplot as plt

e=1.60217662e-19 #Charge of an electron
m_e=9.10938356e-31#Mass of an electron
m_p=1.6726219e-27#mass of a proton
c=299792458.#velocity of light
epsilon_0= 8.854187817e-12 #Pemitivit of vaccum

sig_T=6.6524587158e-29

E_N=1#Energy of supernova explsion
n_N=1
n_p=1#
n_e=1
epsilon_B=1
n_1=1
Y=1
p=3


def calc_GAMMA(t,E,n):
    GAMMA= 260*np.power(E,1/8)*np.power(n,-1/8)*np.power(t,-3/8)
    return GAMMA

def calc_R(t,E,n):
    Rt=3.2e6*np.power(E,1/4)*np.power(n,-1/4)*np.power(t,1/4)
    return Rt

def calc_B(GAMMA):
    B=(8*np.pi)*epsilon_B*4*np.power(GAMMA,2)*n_1*m_p*(c**2)
    return np.sqrt(B)

def calc_freq(GAMMA,gamma,B):
    B_GAMMA=calc_B(GAMMA)
    freq=(3/np.pi)*GAMMA*np.power(gamma,2)*e*B_GAMMA/(m_e*c)
    return freq

def calc_gamma_m(GAMMA,p):
    # if p>2:
    #     g= (p-2)/(p-1)
    # else:
    #     g= (p-2)/(p-1)
    g= (p-2)/(p-1)
    
    gamma_m=g*(GAMMA-1)*((m_p/m_e))*((n_p/n_e))
    return gamma_m

def calc_gamma_c(t,GAMMA):
    B=calc_B(GAMMA)
    gamma_c=6*np.pi*m_e*c/(sig_T*GAMMA*((B)**2)*t*(1+Y))
    return gamma_c
    
    
def f1(freq,freq_c,freq_m,F_0):
        return F_0*np.power(freq_c/freq_m,-(p-1)/2)*np.power(freq/freq_c,-(p)/2)
def f2(freq,freq_c,freq_m,F_0):
    return F_0*np.power(freq/freq_m,-(p-1)/2)        

def calc_flux(t,freqs,F_0):
    GAMMA=calc_GAMMA(t,E_N,n_N)
    gamma_m=calc_gamma_m(GAMMA,p)
    gamma_c=calc_gamma_c(t,GAMMA)
    B=calc_B(GAMMA)
    freq_m=calc_freq(GAMMA,gamma_m,B)
    freq_c=calc_freq(GAMMA,gamma_c,B)

    F=[]
    for f_c,f_m in zip(freq_c,freq_m):
        mask1=freqs>f_c
        mask2= (f_m<freqs) & (freqs<f_c)
        mask3= freqs<f_m
        temp1=f1(freqs[mask1],f_c,f_m,F_0)
        temp2=f2(freqs[mask2],f_c,f_m,F_0)
        temp3=np.zeros(len(freqs[mask3]))
        F.append(np.concatenate((temp3,temp2,temp1)) )
        
    
    return(gamma_m,gamma_c,GAMMA,B,freq_m,freq_c,np.array(F))
    # freqs=np.logspace(1,20)
    # F=[]
    # f_0=10
    
    

t=np.linspace(1e-3,1e-1)
freqs=np.logspace(1,35,100)
(gamma_m,gamma_c,GAMMA,B,freq_m,freq_c,F)= calc_flux(t,freqs,10)

plt.loglog(freqs,F[40,:])
    
    
    