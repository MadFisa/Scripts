#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Nov  2 23:27:12 2020

@author: asif
"""
import numpy as np
import matplotlib.pyplot as plt 
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

no_of_specs=7
base_name="spec"    
file_name="tbin_bb.txt"
GRB="080207"

ph,nm=Log_reader(base_name,no_of_specs)

t_bin=Read_t_file(file_name)

ph_idx=ph[:,0] - ph[:,2]
t=[]
t_n=[]
t_p=[]

for i in range(len(t_bin)-1):
    t.append(np.sqrt(t_bin[i]*t_bin[i+1]))
    t_n.append(t[i]-t_bin[i])
    t_p.append(t_bin[i+1]-t[i])
    
    
    

plt.errorbar(t,ph_idx,xerr=(t_n,t_p),yerr=(-ph[:,2],ph[:,3]),fmt='.k',ecolor='lightgray', elinewidth=3, capsize=0)
plt.title("GRB "+ GRB)
plt.xscale('log')
plt.xlabel("Time (s)")
plt.ylabel('Photon Index')
plt.savefig(GRB+"_phIdx.png",dpi=800)

final_data=np.array([t,t_n,t_p,ph_idx,-ph[:,2],ph[:,3]]).T
np.savetxt("ph_idx_"+"GRB "+ GRB+".txt",final_data,header="Time Time_-error Time_+ve_error Photon_index -ve_error +ve_error")