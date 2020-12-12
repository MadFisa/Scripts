id#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sun Oct 18 18:09:33 2020

@author: asif
"""
import numpy as np







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



def Calc_t_bin(data,r_bin=10.,min_counts=100):
    """Takes tdata which has time,count rate. Calculates time bins 
    with uniform binning in log with minimum couts"""
    tstart=data[:,0]
    del_t=np.diff(tstart)
    cr=data[:,3]
    # cr_err=data[:,2]
    del_t=np.insert(del_t, 0, del_t[0])
    
    counts=cr*del_t
    
    t_bin=[0.0,r_bin]
    
    
    while tstart[-1]//t_bin[-1]!=0:
        t_bin.append(t_bin[-1]*r_bin)
    
    total_bins=len(t_bin)
    
    c_bin=[]
    
    for i in range(total_bins-1):
        mask=(tstart >= t_bin[i]) & (tstart<=t_bin[i+1])
        print(counts[mask])
        cum=np.sum(counts[mask])
        c_bin.append(cum)
    
    
    
    indices=[]
    for i in range(total_bins-1):
        if c_bin[i]<min_counts:
            c_bin[i+1]=c_bin[i+1] + c_bin[i]
            indices.append(i)
    
    c_bin=np.delete(c_bin, indices)
    t_bin=np.delete(t_bin, np.array(indices)+1)
    # t_bin=np.append(t_bin,tstart[-1]) 
    return t_bin,c_bin


def CreateBinFile(t_bin,out_file_name="tbin.txt"):
    """Takes bin times as input and create a text file used for time splicing in xrt site"""
    n=len(t_bin)
    f=open(out_file_name,'x')    
    for i in range(n-1):
        tstart=t_bin[i]
        tstop=t_bin[i+1]
        line="spec"+str(i)+" "+str(tstart)+"-"+str(tstop)
        f.write(line+"\n")
    f.close()
    

lc_filename="curve_mod.qdp"
bin_textfile="tbin_2.txt"
r_bin1=2
min_counts1=100

modes,data=Read_lc_qdp(lc_filename)
PC_data=np.array(data[-1],dtype=np.float32)

t_bin,c_bin=Calc_t_bin(PC_data,r_bin1,min_counts1)


CreateBinFile(t_bin,bin_textfile)

    
    