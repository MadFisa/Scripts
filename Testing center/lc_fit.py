#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Dec 11 11:53:04 2020

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
plt.xscale('log')
plt.xlabel("Time (s)")
plt.yscale('log')
