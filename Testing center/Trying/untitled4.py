#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sat Jan 23 13:56:22 2021

@author: asif
"""
import numpy as np
import matplotlib.pyplot as plt

e=4.80325e-10 #Charge of an electron
m_e=9.10938356e-28#Mass of an electron
m_p=1.6726219e-24#mass of a proton
c=299792458e2#velocity of light
# epsilon_0= 8.854187817e-12 #Pemitivity of vaccum

sig_T=6.6524587158e-25#Thompson cross section

##Model parameters####

# E_52=1e50/1e52#ergs Energy of GRB(GAMMA)
n_N=1.#Number density of nova?Shockwave?(GAMMA) (n_p)
n_p=1.#proton number density(gamm_m)
n_e=1.#(gamma_m)
# F_0=1.#Flux Normalisation(Flux calculation)
epsilon_B=0.05#Magnetic Field energy density(Magnetic field)
epsilon_e=0.015#Electric Field energy density(gamma_m)
n_1=1.#What's this??(In magnetic filed calculation) (n_p)
Y=1.#IC correction factor??(gamma_c calculation)
p=3#Power index
h=6.6261e-27# Plank's constant

"""Try to input arguments as np arrays always"""

def calc_GAMMA(t,E,n):
    """Calculates Bulk  Lorentz factor at time t from Energy E(IN E_52ergs), and number density n"""
    GAMMA= 260*np.power(E,1/8)*np.power(n,-1/8)*np.power(t,-3/8)
    # GAMMA[GAMMA<1]=np.NAN
    return GAMMA

def calc_R(t,E,n):
    """Calculates Radius? at time t from Energy E(IN E_52ergs), and number density n"""
    Rt=3.2e6*np.power(E,1/4)*np.power(n,-1/4)*np.power(t,1/4)
    return Rt

def calc_B(GAMMA):
    """Calculates magnetic field for given bulk lorentz factor GAMMA"""
    B=(8*np.pi)*epsilon_B*4*np.power(GAMMA,2)*n_1*m_p*(c**2)
    return np.sqrt(B)

def calc_freq(GAMMA,gamma):
    """Gyration frequncy for given bulg lorenz factor GAMMA and enrgy of electron GAMMA"""
    B_GAMMA=calc_B(GAMMA)
    freq=(3/4*np.pi)*GAMMA*np.power(gamma,2)*e*B_GAMMA/(m_e*c)
    return freq

def calc_gamma_m(GAMMA,p):
    """Lower limit gamma for powerlaw distribution with index p. Currently assumes p>2"""
    # if p>2:
    #     g= (p-2)/(p-1)
    # else:
    #     g= (p-2)/(p-1)
    g= (p-2)/(p-1)
    
    gamma_m=g*epsilon_e*(GAMMA-1)*(m_p/m_e)*(n_p/n_e)
    return gamma_m

def calc_gamma_c(t,GAMMA):
    """cooling gamma for given bulk lorentz factor"""
    
    B=calc_B(GAMMA)
    # print(B)
    gamma_c=(6*np.pi*m_e*c)/(sig_T*GAMMA*t*np.power(B,2)*(1+Y))
    return gamma_c
    
    
def f1(freq,freq_c,freq_m,F_0):
    """Asymptotic function for f>>f_c"""    
    return F_0*np.power(freq_c/freq_m,-(p-1)/2)*np.power(freq/freq_c,-(p)/2)

def f4(freq,freq_c,freq_m,F_0):
    """Asymptotic function for f>>f_c"""    
    return F_0*np.power(freq_c/freq_m,(1)/2)*np.power(freq/freq_m,-(p)/2)

def f2(freq,freq_m,F_0):
    """Asymptotic function for f_m<f<f_c"""
    return F_0*np.power(freq/freq_m,-(p-1)/2) 

def f5(freq,freq_c,F_0):
    """Asymptotic function for f_m<f<f_c"""    
    return F_0*np.power(freq/freq_c,-1/2)        

def f3(freqs,freq_m,F_0):   
    """Asymptotic function for f_m>>f"""    
    return F_0*np.power(freqs/freq_m,1/3)

def f6(freqs,freq_c,F_0):
    """Asymptotic function for f_m>>f"""    
    return F_0*np.power(freqs/freq_c,1/3)

def calc_flux_sharpo(t,freqs,F_0,E_52):
    
    """Calculate flux for times t(Pass a np.array even if single element) and frequencies with normalised flux F_0
    as piecewise asymptotic functions. Indexing for Flux is F[time,frequency]"""
    
    GAMMA=calc_GAMMA(t,E_52,n_N) #Calculates bulk lorentz factor
    gamma_m=calc_gamma_m(GAMMA,p) #Calculates lower gamma factor
    gamma_c=calc_gamma_c(t,GAMMA) #Calculates cooling gamma factor
    B=calc_B(GAMMA) #Magnetic field for given lorents facote
    freq_m=calc_freq(GAMMA,gamma_m)
    freq_c=calc_freq(GAMMA,gamma_c)
    
    #Checking for different regimes
    mask_slow=freq_m<freq_c
    mask_Fast=freq_m>freq_c
    F=[]
    
    for f_c,f_m in zip(freq_c[mask_Fast],freq_m[mask_Fast]):

        #Fast cooling
        mask4=freqs>f_m 
        mask5= (f_c<freqs) & (freqs<f_m)
        mask6= freqs<f_c
        
        #Caclculating piecewise functions
        temp4=f4(freqs[mask4],f_c,f_m,F_0)
        temp5=f5(freqs[mask5],f_c,F_0)
        temp6=f6(freqs[mask6],f_c,F_0)
        #Appending the piecewise values(Assumes frequencies are ordered in increasing order?)
        F.append(np.concatenate((temp6,temp5,temp4)) )
    
    for f_c,f_m in zip(freq_c[mask_slow],freq_m[mask_slow]):
        #Slow Cooling
        mask1=freqs>f_c 
        mask2= (f_m<freqs) & (freqs<f_c)
        mask3= freqs<f_m
        
        #Caclculating piecewise functions
        temp1=f1(freqs[mask1],f_c,f_m,F_0)
        temp2=f2(freqs[mask2],f_m,F_0)
        temp3=f3(freqs[mask3],f_m,F_0)
        #Appending the piecewise values(Assumes frequencies are ordered in increasing order?)
        F.append(np.concatenate((temp3,temp2,temp1)) )
            


    
    return(gamma_m,gamma_c,GAMMA,B,freq_m,freq_c,np.array(F))
##################################Read light curve data########################

GRB="111209"
lc_file_name="flux_"+"GRB "+ GRB+".txt"
lc_file=np.loadtxt(lc_file_name)
t_data=lc_file[:,0]
lc_data=lc_file[:,3]
lc_err=2*lc_file[:,4] 
nu=np.array([2.4e17*10])
t=np.logspace(2,6,300)
nu=np.array([2.4e17*10])

# gamma_m,gamma_c,GAMMA,B,freq_m,freq_c,F=calc_flux_sharpo(t,nu,1,1e1)

F_0,E_52=(1,1e1)
gamma_m,gamma_c,GAMMA,B,freq_m,freq_c,lc_model=calc_flux_sharpo(t_data,nu,F_0,E_52)
# ph_model=model_ph(t_data,nu,beta1,beta2,p_F,p_b,F_0,nu_b0)
# likelihood = ((-0.5*(ph_model- ph_data)/ph_err)**2).sum()
loglikelihood = (-0.5*( (lc_model- lc_data))**2).sum()
loglikelihood1=loglikelihood + 0.5*sum(np.log(2.*np.pi*lc_err**2))