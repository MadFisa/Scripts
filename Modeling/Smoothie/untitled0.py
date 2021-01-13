# -*- coding: utf-8 -*-
"""
Created on Sat Aug 29 13:54:42 2020

@author: Asif
"""
import numpy as np
import matplotlib.pyplot as plt


def f1(x,s):
    return np.power(x,-3*s)
    

def f2(x,s):
    return np.power(x,-8*s)

def f(x,x0,s):
    y=f1(x/x0,s)+ f2(x/x0,s)
    return np.power(y,1/s)

x=np.logspace(1,10)
x0=1e5

y=f(x,x0,-3)
plt.loglog(x,y)
