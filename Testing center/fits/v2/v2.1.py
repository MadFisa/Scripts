#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sat Dec 12 01:08:19 2020

@author: asif
"""


import numpy as np
import matplotlib.pyplot as plt
plt.style.use('seaborn-whitegrid')

def Read_lc_qdp(qdp_file_name):
    with open(qdp_file_name,'r') as reader:
        lines=reader.readlines()
        modes=[]
        data=[]
        b=iter(lines[9:])
        data=[]
        block=[]
        for line in b:
            if  line[0].isalpha():
                data.append(block)
                block=[]            
                line=next(b)            
                modes.append(line.strip().split(' ')[1])
                line=next(b)
                line=next(b)
            block.append(line.strip().split('\t'))    
            # print(line)
        data.append(block)
        data.pop(0)
        return [modes,data]

def F_v(nu,nu_b,F_0,beta1=3.,beta2=2.,s=2.):
    """F[t_i,nu_i]"""
    F=[]
    for Fi,nu_bi in zip(F_0,nu_b):        
        temp=np.power(nu/nu_bi,-beta1*s)+np.power(nu/nu_bi,-beta2*s)
        F.append(Fi*np.power(temp,-1/s))
    return np.array(F)

def nu_break(t,nu_0,t0=1,m=-2):
    return nu_0*np.power(t/t0,m)

def F_0_evln(t,F_0,t0=1,m=-2):
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


lc_filename="curve.qdp"
GRB="100621"
modes,data=Read_lc_qdp(lc_filename)


##############################################Light curve#########################################
PC_data=np.array(data[-1],dtype=np.float)
WT_data=np.array(data[1],dtype=np.float)
SL_data=np.array(data[0],dtype=np.float)

t=np.concatenate((SL_data[:,0],WT_data[:,0],PC_data[:,0]))
terr=(-np.concatenate((SL_data[:,2],WT_data[:,2],PC_data[:,2])) , np.concatenate((SL_data[:,1],WT_data[:,1],PC_data[:,1])))
t=np.concatenate((SL_data[:,0],WT_data[:,0],PC_data[:,0]))
ferr=(-np.concatenate((SL_data[:,4],WT_data[:,4],PC_data[:,4])) , np.concatenate((SL_data[:,5],WT_data[:,5],PC_data[:,5])))
f=np.concatenate((SL_data[:,3],WT_data[:,3],PC_data[:,3]))



plt.errorbar(t,f,xerr=terr,yerr=ferr,fmt='.k',ecolor='lightgray', elinewidth=3, capsize=0)

#####################################################################################################
nu=np.logspace(14,25,1000)
t_m=np.logspace(2,6,100)
t_b0=10
t_f0=1

m_b=-1.5
m_f=-0.5

nu_b=nu_break(t_m,1e23,t_b0,m_b)
beta1=-1.5
beta2=-2.5
s=3
F_0=F_0_evln(t_m, 1e4,t_f0,m_f)

freqs_m=2.4e17*np.logspace(-0.5,1.,100)
##################################################################################
# F_f=F_v(nu,nu_b,F_0,beta1,beta2,s)
nu_obs=1e21
F_f1=F_v([nu_obs],nu_b,F_0,beta1,beta2,s)  

plt.plot(t_m,F_f1)

plt.xscale('log')
plt.xlabel("Time (s)")
plt.yscale('log')
plt.title("GRB 100621A")
# plt.savefig("lc_fit.png",dpi=800)

##################################################################################
# spec=np.loadtxt("figj1.qdp",skiprows=3)
nu_b_obs=nu_break(1e5,1e28,1e2,m_b)
# F_f2=F_v(nu,[nu_b_obs],0.06,beta1,beta2,s).T  

# plt.figure()
# plt.errorbar(spec[:,0],spec[:,2],xerr=spec[:,1],yerr=spec[:,3],fmt='.k',ecolor='lightgray', elinewidth=3, capsize=0)
# plt.plot(nu,F_f2)
# plt.xscale('log')
# plt.xlabel("Freq (hz)")
# plt.yscale('log')
# plt.title("GRB 100621A")
# plt.savefig("spec_fit.png",dpi=800)
##################################################################################

F=F_v(freqs_m,nu_b,F_0,beta1,beta2,s)
sp_idx_m=-IdxEvln(F,freqs_m)[:,0]
ph_idx_m=sp_idx_m+1


plt.figure()
plt.plot(t_m,ph_idx_m)




ph=np.loadtxt("ph_idx_GRB 100621.txt")

t=ph[:,0]
t_n=ph[:,1]
t_p=ph[:,2]

ph_idx=ph[:,3]
ph_n=ph[:,4]
ph_p=ph[:,5]


plt.errorbar(t,ph_idx,xerr=(t_n,t_p),yerr=(ph_n,ph_p),fmt='.k',ecolor='lightgray', elinewidth=3, capsize=0)
plt.title("GRB "+ GRB)
plt.xscale('log')
plt.xlabel("Time (s)")
plt.ylabel('Photon Index')


