#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Nov 30 17:59:59 2020

@author: asif
"""
import numpy as np
from astropy.stats import bayesian_blocks

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
    

lc_filename="curve.qdp"
bin_textfile="tbin_bb.txt"
modes,data=Read_lc_qdp(lc_filename)
PC_data=np.array(data[-1],dtype=np.float32)
t=PC_data[:,0]
x=PC_data[:,3]
del_t=np.diff(t)
s=(PC_data[:,4]-PC_data[:,5])/2

# edges = bayesian_blocks(t, fitness='events')
# edges = bayesian_blocks(t, x,fitness='measures')
# edges1 = bayesian_blocks(t, x,sigma=s,fitness='measures')
edges = bayesian_blocks(t, x,sigma=s,fitness='measures',p0=0.01)


CreateBinFile(edges,bin_textfile)

# plt.scatter(range(len(edges)),edges)
# Out[19]: <matplotlib.collections.PathCollection at 0x7f303bb3c250>

# plt.scatter(range(len(edges1)),edges1)
# Out[20]: <matplotlib.collections.PathCollection at 0x7f30380653d0>

# plt.scatter(range(len(edges2)),edges2)
# Out[21]: <matplotlib.collections.PathCollection at 0x7f303806c760>

# plt.yscale("log")