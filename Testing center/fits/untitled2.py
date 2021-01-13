#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sat Dec 12 17:44:15 2020

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

def F_v(nu,nu_b,F_0=1.,beta1=3.,beta2=2.,s=2.):
    """F[t_i,nu_i]"""
    F=[]
    for nu_bi in nu_b:        
        temp=np.power(nu/nu_bi,-beta1*s)+np.power(nu/nu_bi,-beta2*s)
        F.append(F_0*np.power(temp,-1/s))
    return np.array(F)

def nu_break(t,nu_0,t0=1,m=-2):
    return nu_0*np.power(t/t0,m)



lc_filename="curve.qdp"
modes,data=Read_lc_qdp(lc_filename)



PC_data=np.array(data[-1],dtype=np.float)
WT_data=np.array(data[1],dtype=np.float)
SL_data=np.array(data[0],dtype=np.float)

t=np.concatenate((SL_data[:,0],WT_data[:,0],PC_data[:,0]))
terr=(-np.concatenate((SL_data[:,2],WT_data[:,2],PC_data[:,2])) , np.concatenate((SL_data[:,1],WT_data[:,1],PC_data[:,1])))
t=np.concatenate((SL_data[:,0],WT_data[:,0],PC_data[:,0]))
ferr=(-np.concatenate((SL_data[:,4],WT_data[:,4],PC_data[:,4])) , np.concatenate((SL_data[:,5],WT_data[:,5],PC_data[:,5])))
f=np.concatenate((SL_data[:,3],WT_data[:,3],PC_data[:,3]))


plt.errorbar(t,f,xerr=terr,yerr=ferr,fmt='.k',ecolor='lightgray', elinewidth=3, capsize=0)

nu=np.logspace(17,18.5,1000)
t_f=np.logspace(2,6,100)
m=-6
nu_b=nu_break(t_f, 1e30,m)
beta1=-0.5
beta2=-0.8
s=3
F_0=1e-6

F_f=F_v(nu,nu_b,F_0,beta1,beta2,s)
# nu_obs=1e21
# F_f1=F_v([nu_obs],nu_b,0.06,beta1,beta2,s)  
F_f1=np.sum(F_f,axis=1)

plt.plot(t_f,F_f1)

plt.xscale('log')
plt.xlabel("Time (s)")
plt.yscale('log')
plt.title("GRB 100621A")
# plt.savefig("lc_fit.png",dpi=800)


spec=np.loadtxt("figj1.qdp",skiprows=3)
nu_b_obs=nu_break(1e6,1e28,1e2,m)
F_f2=F_v(nu,[nu_b_obs],F_0,beta1,beta2,s).T  

plt.figure()
plt.errorbar(spec[:,0],spec[:,2],xerr=spec[:,1],yerr=spec[:,3],fmt='.k',ecolor='lightgray', elinewidth=3, capsize=0)
plt.plot(nu,F_f2)
plt.xscale('log')
plt.xlabel("Freq (hz)")
plt.yscale('log')
plt.title("GRB 100621A")
# plt.savefig("spec_fit.png",dpi=800)