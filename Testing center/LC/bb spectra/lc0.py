#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Dec 30 18:48:42 2020

@author: asif
"""

import numpy as np
import matplotlib.pyplot as plt 
import linecache
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

def read_corr(no_of_specs,base_name="spec"):
    """Reads File count rate and correction factor"""
    count_r=[]
    corr_fact=[]
    for i in range(no_of_specs):
        path="./"+base_name+str(i)+"/"+base_name+str(i)+"pc_fit.fit"
        # path="./spec0/spec0pc_fit.fit"
        
        temp_count_r=linecache.getline(path,9).strip().replace('File count rate: ','')
        count_r.append(temp_count_r)
        temp_corr_fact=linecache.getline(path,10).strip().replace('Corr factor:\t','')
        corr_fact.append(temp_corr_fact)
    
    count_r=np.array(count_r,dtype=np.float32)
    corr_fact=np.array(corr_fact,dtype=np.float32)
    return [count_r,corr_fact]

time_file_name="tbin_bb.txt"
GRB="080207"
text_file="ph_idx_"+"GRB "+ GRB+".txt"
t_bin=Read_t_file(time_file_name)
no_of_specs=7
base_name="spec"

data=np.genfromtxt(text_file,skip_header=1)
t_data=data[:,:3]
ph_data=data[:,3:]

[file_count_r,corr_fact]=read_corr(no_of_specs,base_name)


