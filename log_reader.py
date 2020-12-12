#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Nov  2 22:03:16 2020

@author: asif
"""

import numpy as np

no_of_specs=8
base_name="spec"



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
    
ph,nm=Log_reader(base_name,no_of_specs)