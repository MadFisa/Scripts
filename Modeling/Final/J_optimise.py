# -*- coding: utf-8 -*-
"""
Created on Sun Sep  6 23:42:36 2020

@author: Asif
"""


import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit


###Constants, Uses CGS units####

e=4.80325e-10 #Charge of an electron
m_e=9.10938356e-28#Mass of an electron
m_p=1.6726219e-24#mass of a proton
c=299792458e2#velocity of light
# epsilon_0= 8.854187817e-12 #Pemitivity of vaccum

sig_T=6.6524587158e-25#Thompson cross section

##Model parameters####

E_52=1e50/1e52#ergs Energy of GRB(GAMMA)
n_N=1.#Number density of nova?Shockwave?(GAMMA) (n_p)
n_p=1.#proton number density(gamm_m)
n_e=1.#(gamma_m)
F_0=1.#Flux Normalisation(Flux calculation)
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

# def SlowCooling(freqs,freq_c,freq_m,F_0):
    
#     F=[]
#     for f_c,f_m in zip(freq_c,freq_m):
#         #Checking for different regimes
        
#         mask1=freqs>f_c 
#         mask2= (f_m<freqs) & (freqs<f_c)
#         mask3= freqs<f_m
        
#         #Caclculating piecewise functions
#         temp1=f1(freqs[mask1],f_c,f_m,F_0)
#         temp2=f2(freqs[mask2],f_c,f_m,F_0)
#         temp3=f3(freqs[mask3],f_m,F_0)
#         #Appending the piecewise values(Assumes frequencies are ordered in increasing order?)
#         F.append(np.concatenate((temp3,temp2,temp1)) )
    

    
#     return F

# def FastCooling(freqs,freq_c,freq_m,F_0):
    
#     F=[]
#     for f_c,f_m in zip(freq_c,freq_m):
#         #Checking for different regimes
        
#         mask4=freqs>f_m 
#         mask5= (f_c<freqs) & (freqs<f_m)
#         mask6= freqs<f_c
        
#         #Caclculating piecewise functions
#         temp4=f4(freqs[mask4],f_c,f_m,F_0)
#         temp5=f5(freqs[mask5],f_c,f_m,F_0)
#         temp6=f6(freqs[mask6],f_c,F_0)
#         #Appending the piecewise values(Assumes frequencies are ordered in increasing order?)
#         F.append(np.concatenate((temp6,temp5,temp4)) )
    

#     return F
    


def calc_flux_sharpo(t,freqs,F_0):
    
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

def Slow_smoothie(freq,freq_c,freq_m,F_0,s1,s2):
    """Calculates a smoothened version of slow cooling flux with smoothening parameters s1,s2,s3"""
    
    y1=F_0*np.power(freq/freq_m,1/3)     
    y2=lambda s:(1+np.power(freq/freq_m,(1/3-(-((p-1)/2)))*s))
    y3=lambda s:(1+np.power(freq/freq_c,((1/2))*s) )
    
    temp1=y1
    temp2=np.power(y2(s1),-1/s1)
    temp3=np.power(y3(s2),-1/s2)
    y= temp1*temp2*temp3
    # y=temp1
    return y

def Fast_smoothie(freq,freq_c,freq_m,F_0,s1,s2):
    """Calculates a smoothened version of Fast cooling flux with smoothening parameters s1,s2,s3"""
    
    y1=F_0*np.power(freq/freq_c,1/3)     
    y2=lambda s:(1+np.power(freq/freq_m,(1/3-(-1/2))*s))
    y3=lambda s:(1+np.power(freq/freq_c,((-1/2)-(-p/2))*s) )
    
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
    s1=5
    s2=5
        #Checking for different regimes
    mask_slow=freq_m<freq_c
    mask_Fast=freq_m>freq_c
    F=[]
    for f_c,f_m in zip(freq_c[mask_Fast],freq_m[mask_Fast]):
        #Fast cooling
            F.append(Fast_smoothie(freqs,f_c,f_m,F_0,s1,s2))
    for f_c,f_m in zip(freq_c[mask_slow],freq_m[mask_slow]):
        #slow cooling
            F.append(Slow_smoothie(freqs,f_c,f_m,F_0,s1,s2))

    
    return(gamma_m,gamma_c,GAMMA,B,freq_m,freq_c,np.array(F)) 


def powerLaw(x,p,a):
    """Defines a power law for fitting functions"""
    return a*np.power(x,p)

def powerFit(x_data,y_data,p=0.5,a=1.):
    """Fits a power law for given data with initial guesse p and a"""
    popt,pcov=curve_fit(powerLaw,x_data,y_data,p0=[p,a],maxfev=2000)
    return popt,pcov

def linearLaw(x,m,b):
    """Defines a linear law for fitting functions"""
    return (m*x)+b

def IdxEvln(F,freqs):
    "Fits a power law to F[i1,i2] for all i1 with i2 being freqs"
    
    s_idx=[]
    for ff in F:
        model_x=np.log(freqs)
        model_y=np.log(ff)
        temp = np.polyfit(model_x,model_y, 1)
        s_idx.append(temp)
    return np.array(s_idx)

# def lineFit(x_data,y_data,p=3,a=1):
#     """Fits a power law for given data with initial guesse p and a"""
#     popt,pcov=curve_fit(powerLaw,x_data,y_data,p0=[p,a],maxfev=2000)
#     return popt,pcov
    

t_temp=np.logspace(-2,8,4000)#Defines time
GAMMA_temp=calc_GAMMA(t_temp,E_52,n_N)

t=t_temp[GAMMA_temp>1]
freqs=np.logspace(1,35,1000)#Defines frequencies
# freqs=2.4e17*np.logspace(-0.5,1.,100)
#Defines frequencies
(gamma_m1,gamma_c1,GAMMA1,B1,freq_m1,freq_c1,F1)= calc_flux_sharpo(t,freqs,F_0)
(gamma_m2,gamma_c2,GAMMA2,B2,freq_m2,freq_c2,F2)= calc_flux_smoothie(t,freqs,F_0)



ti=3000  


fig, axs = plt.subplots(2,sharex=True)
# axs[0].xlabel("Frequency in Hertz")
# axs[0].ylabel("Normalised flux")
axs[0].loglog(freqs,F1[ti,:],label='sharpo')
axs[0].loglog(freqs,F2[ti,:],label='smoothie')
# plt.loglog(freqs,F3[ti,:],label='scaled smoothie')
axs[1].loglog(freqs,np.abs( F2[ti,:]-F1[ti,:]))
axs[0].legend()

# for ti in range(0,len(t),1000):    
#     plt.loglog(freqs,F1[ti,:],label='Time = % 10.3E'%(t[ti]))
#     plt.legend()




# plt.savefig("F vs Freq.png",dpi=600)

# plt.figure()
# plt.xlabel("Time in s")
# plt.ylabel("Normalised flux")
# for fi in range(0,len(freqs),20):    
#     plt.loglog(t,F1[:,fi],label='Freq = % 10.3E'%(freqs[fi]))
#     plt.legend()
# # plt.savefig("F vs time.png",dpi=600)


# window_f1=500 #Window frequency 1 index
# window_f2=700 #Window frequency 2 index
window_f1=0 #Window frequency 1 index
window_f2=-1 #Window frequency 2 index

window_t1=0
window_t2=-1

obs_F=F2[:,window_f1:window_f2]
obs_freqs=freqs[window_f1:window_f2]
# sp_idx=IdxEvln(obs_F, obs_freqs)

E=h*freqs
n_ph=F2/E
obs_n_ph=n_ph[window_t1:window_t2,window_f1:window_f2]
obs_E=E[window_f1:window_f2]
obs_t=t[window_t1:window_t2]
ph_idx=IdxEvln(obs_n_ph,obs_E)

plt.figure()
plt.plot(obs_t,-ph_idx[:,0])
plt.xscale('log')
plt.xlabel("Time")
plt.ylabel("Photon index")
# # plt.savefig("photon_index.png")


fig1, axs1 = plt.subplots(3,sharex=True)
# axs[0].xlabel("Frequency in Hertz")
# axs[0].ylabel("Normalised flux")
axs1[0].plot(obs_t,-ph_idx[:,0],label='photon index')

# plt.loglog(freqs,F3[ti,:],label='scaled smoothie')
axs1[1].loglog(obs_t,gamma_m2[window_t1:window_t2],label='gamma_m')
axs1[1].loglog(obs_t,gamma_c2[window_t1:window_t2],label='gamma_c')
axs1[2].loglog(obs_t,freq_m2[window_t1:window_t2],label='f_m')
axs1[2].loglog(obs_t,freq_c2[window_t1:window_t2],label='f_c')
axs1[0].legend()
axs1[1].legend()
axs1[2].legend()




# plt.figure()
# plt.plot(t,sp_idx[:,0])
# plt.xscale('log')


# model_x=np.log(freqs[window_f1:window_f2])
# model_y=np.log(F2[ti,window_f1:window_f2])

# m, b = np.polyfit(model_x, model_y, 1)
# fit=np.exp(linearLaw(model_x,m,b))

# plt.loglog(freqs,F2[ti,:],'--',label='smoothie',)
# plt.loglog(freqs[window_f1:window_f2],fit,label='fit')
# plt.loglog(freqs[window_f1:window_f2],F2[ti,window_f1:window_f2],'--',label='model')
# plt.legend()





    
    
    
    