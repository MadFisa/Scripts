# -*- coding: utf-8 -*-
"""
Created on Fri Aug 21 19:28:41 2020

@author: Asif
"""
import numpy as np
import matplotlib.pyplot as plt



e=1. #Charge of an electron
m_e=9.10938356e-28#Mass of an electron
m_p=1.6726219e-24#mass of a proton
c=299792458e2#velocity of light
# epsilon_0= 8.854187817e-12 #Pemitivity of vaccum

sig_T=6.6524587158e-25

E_N=1e43#Energy of supernova explsion
n_N=1
n_p=1#
n_e=1
epsilon_B=1
n_1=1
Y=1
p=3

E_N=1e43#ergs Energy of GRB(GAMMA)
n_N=1.#Number density of nova?Shockwave?(GAMMA) (n_p)
n_p=1.#proton number density(gamm_m)
n_e=1.#(gamma_m)
F_0=1.#Flux Normalisation(Flux calculation)
epsilon_B=0.05#Magnetic Field energy density(Magnetic field)
epsilon_e=0.15#Electric Field energy density(gamma_m)
n_1=1.#What's this??(In magnetic filed calculation) (n_p)
Y=1.#IC correction factor??(gamma_c calculation)
p=3#Power index

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
    gamma_c=6*np.pi*m_e*c/(sig_T*((B)**2)*t*(1+Y))
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
    
    

t=np.logspace(-3,10)
freqs=np.logspace(1,35,100)
(gamma_m,gamma_c,GAMMA,B,freq_m,freq_c,F)= calc_flux(t,freqs,10)

# plt.loglog(freqs,F[40,:])
    
    
    