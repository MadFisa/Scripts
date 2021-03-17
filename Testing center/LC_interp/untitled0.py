#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Jan 12 13:10:16 2021

@author: asif
"""


import numpy as np
import matplotlib.pyplot as plt 
# import linecache
from scipy.interpolate import interp1d
plt.style.use('seaborn-whitegrid')

def Read_t_file(file_name):
    """Reads time from files for swift website spectra time slicing"""
    t=[]
    
    with open(file_name,'r') as reader:
        temp=reader.readline().strip().split()[-1].split('-')
        t.append(temp[0])
        t.append(temp[1])
        for line in reader.readlines():
            t.append(line.strip().split()[-1].split('-')[-1])
    
    return np.array(t,dtype=np.float32)


def Log_reader(base_name,no_of_specs):    
    """Takes in base name and number of spectra as argument. Returns Photon index and normalisation
    The data is of the format [lower lim, upper lim, -ve error, +ve error] """    
    ph_idx=[]
    norm=[]
    
    for i in range(no_of_specs):    
        log_name=base_name+str(i)+".log"
        with open(log_name,'r') as reader:
            log=reader.readlines()
        
        line=log[-11]
        temp=line.strip().split(' ')
        temp=[i for i in temp if i !='']
        temp1=temp[-1][1:-1].split(',')
        ph_idx.append(temp[2:-1]+temp1)
        
        
        line=log[-7]
        temp=line.strip().split(' ')
        temp=[i for i in temp if i !='']
        temp1=temp[-1][1:-1].split(',')
        norm.append(temp[2:-1]+temp1)
    
    return [np.array(ph_idx,dtype=np.float64),np.array(norm,dtype=np.float64)]


# time_file_name="tbin_bb.txt"
# t_bin=Read_t_file(time_file_name)
# no_of_specs=7
# base_name="spec"

# ph,nm=Log_reader(base_name,no_of_specs) #Data directly from log

# ph_idx=ph[:,0] - ph[:,2]
# norm=nm[:,0] - nm[:,2]

# ph_idx=ph[:,0] - ph[:,2]
# norm=nm[:,0] - nm[:,2]

t,t_n,t_p,ph,ph_n,ph_p=np.loadtxt("ph_idx_GRB 080207.txt").T
GRB="080207"
# t=data[:,0]
# t_n=data[:,0]
# t_p=data[:,0]
ph_interp = lambda tt: interp1d(np.log10(t),ph)(np.log10(tt))


t_new=np.logspace(3.8,4.6,50)
ph_new=ph_interp(t_new)

plt.errorbar(t,ph,xerr=(t_n,t_p),yerr=(ph_n,ph_p),fmt='.k',ecolor='lightgray', elinewidth=3, capsize=0)
plt.plot(t_new,ph_new)
plt.scatter(t_new,ph_new)
plt.title("GRB "+ GRB)
plt.xscale('log')
plt.xlabel("Time (s)")
plt.ylabel('Photon Index')
# plt.savefig(GRB+"_phIdx.png",dpi=800)
