# -*- coding: utf-8 -*-
"""
Created on Fri Aug 21 19:28:41 2020

@author: Asif
"""
import numpy as np
import matplotlib.pyplot as plt



e=4.80325e-10 #Charge of an electron
m_e=9.10938356e-28#Mass of an electron
m_p=1.6726219e-24#mass of a proton
c=299792458e2#velocity of light
# epsilon_0= 8.854187817e-12 #Pemitivity of vaccum

sig_T=6.6524587158e-25


E_52=1e43/1e52#ergs Energy of GRB(GAMMA)
n_N=1.#Number density of nova?Shockwave?(GAMMA) (n_p)
n_p=1.#proton number density(gamm_m)
n_e=1.#(gamma_m)
F_0=100.#Flux Normalisation(Flux calculation)
epsilon_B=0.05#Magnetic Field energy density(Magnetic field)
epsilon_e=0.15#Electric Field energy density(gamma_m)
n_1=1.#What's this??(In magnetic filed calculation) (n_p)
Y=1.#IC correction factor??(gamma_c calculation)
p=3#Power index

"""Try to input arguments as np arrays always"""

def calc_GAMMA(t,E,n):
    """Calculates Bulk  Lorentz factor at time t from Energy E(IN E_52ergs), and number density n"""
    GAMMA= 260*np.power(E,1/8)*np.power(n,-1/8)*np.power(t,-3/8)
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
def f2(freq,freq_c,freq_m,F_0):
    """Asymptotic function for f_m<f<f_c"""
    return F_0*np.power(freq/freq_m,-(p-1)/2)        
def f3(freqs,freq_m,F_0):
    """Asymptotic function for f_m>>f"""
    return F_0*np.power(freqs/freq_m,1/3)

def calc_flux_sharpo(t,freqs,F_0):
    """Calculate flux for times t(Pass a np.array even if single element) and frequencies with normalised flux F_0
    as piecewise asymptotic functions. Indexing for Flux is F[time,frequency]"""
    GAMMA=calc_GAMMA(t,E_52,n_N) #Calculates bulk lorentz factor
    gamma_m=calc_gamma_m(GAMMA,p) #Calculates lower gamma factor
    gamma_c=calc_gamma_c(t,GAMMA) #Calculates cooling gamma factor
    B=calc_B(GAMMA) #Magnetic field for given lorents facote
    freq_m=calc_freq(GAMMA,gamma_m)
    freq_c=calc_freq(GAMMA,gamma_c)

    F=[]
    for f_c,f_m in zip(freq_c,freq_m):
        #Checking for different regimes
        
        mask1=freqs>f_c 
        mask2= (f_m<freqs) & (freqs<f_c)
        mask3= freqs<f_m
        
        #Caclculating piecewise functions
        temp1=f1(freqs[mask1],f_c,f_m,F_0)
        temp2=f2(freqs[mask2],f_c,f_m,F_0)
        temp3=f3(freqs[mask3],f_m,F_0)
        #Appending the piecewise values(Assumes frequencies are ordered in increasing order?)
        F.append(np.concatenate((temp3,temp2,temp1)) )
        
    
    return(gamma_m,gamma_c,GAMMA,B,freq_m,freq_c,np.array(F))
    # freqs=np.logspace(1,20)
    # F=[]
    # f_0=10           

def smoothie(freq,freq_c,freq_m,F_0,s1,s2):
    """Calculates a smoothened version of flux with smoothening parameters s1,s2,s3"""
    y1=F_0*np.power(freq/freq_m,1/3)     
    y2=lambda s:(1+np.power(freq/freq_m,(1/3-(-((p-1)/2)))*s))
    y3=lambda s:(1+np.power(freq/freq_c,((1/2))*s) )
    
    temp1=y1
    temp2=np.power(y2(s1),-1/s1)
    temp3=np.power(y3(s2),-1/s2)
    y= temp1*temp2*temp3
    # y=temp1
    return y

def calc_flux_smoothie(t,freqs,F_0):
    """Calculates a smoothened version of flux with smoothening parameters s1,s2,s3 for tmies t and frequncies freqs
    (Pass as numpy arrays even if single elements).Indexing for flux is F[time,frequency]"""
    GAMMA=calc_GAMMA(t,E_52,n_N)
    gamma_m=calc_gamma_m(GAMMA,p)
    gamma_c=calc_gamma_c(t,GAMMA)
    B=calc_B(GAMMA)
    freq_m=calc_freq(GAMMA,gamma_m)
    freq_c=calc_freq(GAMMA,gamma_c)
    s1=1.1
    s2=1.1
    F=[]
    for f_c,f_m in zip(freq_c,freq_m):
        temp=smoothie(freqs,f_c,f_m,F_0,s1,s2)
        # temp=smoothie2(freqs,f_c,f_m,F_0,1.5,1.6)
        F.append(temp)
        # mask1=freqs>f_c
        # mask2= (f_m<freqs) & (freqs<f_c)
        # mask3= freqs<f_m
        # temp1=f1(freqs[mask1],f_c,f_m,F_0)
        # temp2=f2(freqs[mask2],f_c,f_m,F_0)
        # temp3=f3(freqs[mask3],f_m,F_0)
        # F.append(np.concatenate((temp3,temp2,temp1)) )
        
    
    return(gamma_m,gamma_c,GAMMA,B,freq_m,freq_c,np.array(F))
    # freqs=np.logspace(1,20)
    # F=[]
    # f_0=10           

t=np.logspace(-3,2)
freqs=np.logspace(1,35,100)
(gamma_m1,gamma_c1,GAMMA1,B1,freq_m1,freq_c1,F1)= calc_flux_sharpo(t,freqs,F_0)
(gamma_m2,gamma_c2,GAMMA2,B2,freq_m2,freq_c2,F2)= calc_flux_smoothie(t,freqs,F_0)



referance_f=1e9


referance_F1= calc_flux_sharpo(t,np.array([referance_f]),F_0)[-1][:,0]
referance_F2=calc_flux_smoothie(t,referance_f,F_0)[-1]
scale=referance_F1/referance_F2
F3=F2*np.array([scale]).T

ti=40   

plt.figure()
plt.xlabel("Frequency in Hertz")
plt.ylabel("Normalised flux")
plt.loglog(freqs,F1[ti,:],label='sharpo')
plt.loglog(freqs,F2[ti,:],label='smoothie')
plt.loglog(freqs,F3[ti,:],label='scaled smoothie')
plt.legend()
# for ti in range(0,len(t),10):    
#     plt.loglog(freqs,F[ti,:],label='Time = % 10.3E'%(t[ti]))
#     plt.legend()

plt.savefig("F vs Freq.png",dpi=600)

# plt.figure()
# plt.xlabel("Time in s")
# plt.ylabel("Normalised flux")
# for fi in range(0,len(freqs),20):    
#     plt.loglog(t,F[:,fi],label='Freq = % 10.3E'%(freqs[fi]))
#     plt.legend()
# plt.savefig("F vs time.png",dpi=600)

    
    
    
    