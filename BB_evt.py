#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Nov 30 17:59:59 2020

@author: asif
"""
import numpy as np
from astropy.stats import bayesian_blocks
from astropy.io import fits

dir_list="dir"


def read_times(dir_list):
    """Reads xrt data from given directory list and appends all the times"""
    time=[]
    with open(dir_list,'r') as reader:
        for dire in reader.readlines():
            fits_name="sw"+dire[:-2]+"xpcw3po_cl.evt.gz"
            path="./"+dire[:-1]+"xrt/event/"+fits_name
            with fits.open(path) as hdul:
                time.append(hdul['EVENTS'].data["Time"])
    
    return np.concatenate(time)



time=read_times(dir_list)

edges = bayesian_blocks(time, fitness='events',p0=0.01)
edges1 = bayesian_blocks(time, fitness='event')
# edges = bayesian_blocks(t, fitness='events')
# edges = bayesian_blocks(t, x,fitness='measures')
# edges1 = bayesian_blocks(t, x,sigma=s,fitness='measures')
# edges2 = bayesian_blocks(t, x,sigma=s,fitness='measures',p0=0.01)


# plt.scatter(range(len(edges)),edges)
# Out[19]: <matplotlib.collections.PathCollection at 0x7f303bb3c250>

# plt.scatter(range(len(edges1)),edges1)
# Out[20]: <matplotlib.collections.PathCollection at 0x7f30380653d0>

# plt.scatter(range(len(edges2)),edges2)
# Out[21]: <matplotlib.collections.PathCollection at 0x7f303806c760>

# plt.yscale("log")