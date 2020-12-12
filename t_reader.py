#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Nov  2 23:27:12 2020

@author: asif
"""
import numpy as np

file_name="tbin_2.txt"


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
