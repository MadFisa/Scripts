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

def create_blocks(data,edges):
    """Breaks (sorted)data into blocks according to edges and returns indices for np.split"""
    idx=[]
    for e in edges[1:-1]:
        temp=np.where(data<e)
        idx.append(temp[0][-1]+1)

    return idx
    
# def calc_CF()

lc_filename="curve_mod.qdp"
modes,data=Read_lc_qdp(lc_filename)
PC_data=np.array(data[-1],dtype=np.float32)



time_file_name="tbin_bb.txt"
GRB="080207"
text_file="ph_idx_"+"GRB "+ GRB+".txt"
t_bin=Read_t_file(time_file_name)
no_of_specs=7
base_name="spec"

ph,nm=Log_reader(base_name,no_of_specs) #Data directly from log

ph_idx=ph[:,0] - ph[:,2]
norm=nm[:,0] - nm[:,2]


# data=np.genfromtxt(text_file,skip_header=1) #Data from plotter
# t_ph=data[:,:3]
# ph_data=data[:,3:]

[file_count_r,corr_fact]=read_corr(no_of_specs,base_name)   

t_data=PC_data[:,:3]
time=t_data[:,0]
cnt_data=PC_data[:,3:]

block_idx=create_blocks(time,t_bin)#Gets indices for splitting into blocks
cnt_block=np.split(cnt_data,block_idx)#Splits the initial data into blocks

cf=(norm*667.5)/(file_count_r*corr_fact)

ff=[]
for i in range(no_of_specs):
    conversion=cf[i]*np.power(10,-(ph_idx[i]-1))
    temp=conversion*cnt_block[i]
    ff.append(temp)

flux=np.vstack(ff)

plt.scatter(time,flux[:,0])
# plt.scatter(time,cnt_data[:,0])
plt.xscale("log")
plt.yscale("log")

