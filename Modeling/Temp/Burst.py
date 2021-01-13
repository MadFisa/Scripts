# -*- coding: utf-8 -*-
"""
Created on Fri Oct  2 02:19:42 2020

@author: Asif
"""
import numpy as np
import matplotlib.pyplot as plt


ph_idx=np.array([2.10,2.46,2.16,2.07,2.55,3.06])
phPerr=[0.24,0.27,0.27,0.33,0.60,1.17]
phNerr=[0.23,0.25,0.26,0.30,0.51,0.86]
phErr=np.array([phPerr,phNerr])

Nh=[8.8,14.8,7.7,5.2,5.1,7.6]#*1e22
NhPerr=[2.5,3.4,2.6,2.6,3.3,6.2]
NhNerr=[2.2,3.1,2.2,2.1,2.5,4.6]

Enweiph=[2.0,2.05,2.1,2.34,2.6,2.9,3.3]
Enweit=[4.2e3,1.3e4,2e4,5e4,1.2e5,2.3e5,5e5]

time=[(10000,16903),(20216,40000),(40001,80000),(80001,196071),(201028,397611),(402287,883371)]

t=np.array(time,dtype=float)[:,0]


plt.figure()
plt.errorbar(t,ph_idx,yerr=[phNerr,phPerr],fmt='o')
plt.scatter(Enweit,Enweiph,label="Wang 2016",color='r')
plt.legend()
plt.xscale("log")
plt.xlabel("time")
plt.ylabel("Photon index")

plt.figure()
plt.errorbar(t,Nh,yerr=[NhNerr,NhPerr],fmt='o')
plt.xscale("log")
plt.xlabel("time")
plt.ylabel("NH")

plt.figure()
plt.errorbar(ph_idx,Nh,xerr=[phPerr,phNerr],yerr=[NhNerr,NhPerr],fmt='o',color='black',markersize=2,
             ecolor='lightgray', elinewidth=3, capsize=0)
plt.xlabel("photone index")
plt.ylabel("NH")

