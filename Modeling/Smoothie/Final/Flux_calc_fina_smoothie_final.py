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
    freq=(3/4*np.pi)*GAMMA*np.power(gamma,2)*e*B_GAMMA/(m_e*c)
    return freq

def calc_gamma_m(GAMMA,p):
    # if p>2:
    #     g= (p-2)/(p-1)
    # else:
    #     g= (p-2)/(p-1)
    g= (p-2)/(p-1)
    
    gamma_m=g*epsilon_e*(GAMMA-1)*(m_p/m_e)*(n_p/n_e)
    return gamma_m

def calc_gamma_c(t,GAMMA):
    B=calc_B(GAMMA)
    # print(B)
    gamma_c=(6*np.pi*m_e*c)/(sig_T*GAMMA*t*np.power(B,2)*(1+Y))
    return gamma_c
    
    
def f1(freq,freq_c,freq_m,F_0):
        return F_0*np.power(freq_c/freq_m,-(p-1)/2)*np.power(freq/freq_c,-(p)/2)
def f2(freq,freq_c,freq_m,F_0):
    return F_0*np.power(freq/freq_m,-(p-1)/2)        
def f3(freqs,freq_m,F_0):
    return F_0*np.power(freqs/freq_m,1/3)

def smoothie1(freq,freq_c,freq_m,F_0,s):
    y1=F_0*np.power(freq_c/freq_m,((p-1)*s)/2)*np.power(freq/freq_c,(p*s)/2)
    y2=F_0*np.power(freq/freq_m,((p-1)*s)/2)
    # y1=F_0*np.power(freq_c/freq_m,-((p-1)*s)/2)*np.power(freq/freq_c,-(p*s)/2)
    # y2=F_0*np.power(freq/freq_m,-((p-1)*s)/2)
    y= np.power(y1+y2,-1/s)
    return y
def smoothie2(freq,freq_c,freq_m,F_0,s1,s2):
    y1=lambda s:F_0*np.power(freq_c/freq_m,((p-1)*s)/2)*np.power(freq/freq_c,(p*s)/2)
    y2=lambda s:F_0*np.power(freq/freq_m,((p-1)*s)/2)
    y3=lambda s:F_0*np.power(freqs/freq_m,-1/3)
    temp1=np.power(y1(s1)+y2(s1),-1/s1)
    temp2=np.power(y3(s2)+y2(s2),-1/s2)
    y= temp1*temp2
    # y=temp2
    return y
def smoothie3(freq,freq_c,freq_m,F_0,s1,s2,s3):
    y1=lambda s:F_0*np.power(freq/freq_m,-1/3*s)    
    y2=lambda s:(1+np.power(freq/freq_m,(1/3-(-((p-1)/2)))*s))
    y3=lambda s:(1+np.power(freq/freq_c,((1/2))*s) )
    
    temp1=np.power(y1(s1),-1/s1)
    temp2=np.power(y2(s2),-1/s2)
    temp3=np.power(y3(s3),-1/s3)
    y= temp1*temp2*temp3
    # y=temp1
    return y

def calc_flux_sharpo(t,freqs,F_0):
    GAMMA=calc_GAMMA(t,E_52,n_N)
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
        temp3=f3(freqs[mask3],f_m,F_0)
        F.append(np.concatenate((temp3,temp2,temp1)) )
        
    
    return(gamma_m,gamma_c,GAMMA,B,freq_m,freq_c,np.array(F))
    # freqs=np.logspace(1,20)
    # F=[]
    # f_0=10           


def calc_flux_smoothie(t,freqs,F_0):
    GAMMA=calc_GAMMA(t,E_52,n_N)
    gamma_m=calc_gamma_m(GAMMA,p)
    gamma_c=calc_gamma_c(t,GAMMA)
    B=calc_B(GAMMA)
    freq_m=calc_freq(GAMMA,gamma_m,B)
    freq_c=calc_freq(GAMMA,gamma_c,B)
    s1=0.1
    s2=0.1
    s3=0.1
    F=[]
    for f_c,f_m in zip(freq_c,freq_m):
        temp=smoothie3(freqs,f_c,f_m,F_0,1.2,1.1,1.2)
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
(gamma_m1,gamma_c1,GAMMA1,B1,freq_m1,freq_c1,F1)= calc_flux_sharpo(t,freqs,10)
(gamma_m2,gamma_c2,GAMMA1,B2,freq_m2,freq_c2,F2)= calc_flux_smoothie(t,freqs,10)

ti=10
plt.figure()
plt.xlabel("Frequency in Hertz")
plt.ylabel("Normalised flux")
plt.loglog(freqs,F1[ti,:],label='sharpo')
plt.loglog(freqs,F2[ti,:],label='smoothie')
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

    
    
    
    