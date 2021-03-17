#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Dec 11 21:28:57 2020

@author: asif
"""

import numpy as np
import  matplotlib.pyplot as plt
plt.style.use('seaborn-whitegrid')

spec=np.loadtxt("figj1.qdp",skiprows=3)



plt.errorbar(spec[:,0],spec[:,2],xerr=spec[:,1],yerr=spec[:,3],fmt='.k',ecolor='lightgray', elinewidth=3, capsize=0)
plt.xscale('log')
plt.xlabel("Freq (hz)")
plt.yscale('log')
plt.title("GRB 100621A")
plt.savefig("spec_fit.png",dpi=800)