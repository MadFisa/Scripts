#!/usr/bin/env plc_datathon3
# -*- coding: utf-8 -*-
"""
Created on Sun Jan 24 12:57:24 2021

@author: asif
"""

import numpy as np
from astropy.modeling import models,fitting
import matplotlib.pyplot as plt
##################################Read light curve data########################

GRB="190114"
lc_file_name="flux_"+"GRB "+ GRB+".txt"
lc_file=np.loadtxt(lc_file_name)
t_data=lc_file[:,0]
lc_data=lc_file[:,3]
lc_err=2*lc_file[:,4] 
nu=2.4e17*np.logspace(-0.5,1.,100)
#########################Fitting###############################################
# BrokenPL=models.BrokenPowerLaw1D(amplitude=1, x_break=1,alpha_1=1,alpha_2=1,bounds={"x_break": (t_data[0], t_data[-1])})
BrokenPL=models.BrokenPowerLaw1D(amplitude=1, x_break=1,alpha_1=1,alpha_2=1)
# BrokenPL=models.BrokenPowerLaw1D()
fit = fitting.LevMarLSQFitter()
# fit= fitting.LinearLSQFitter()
PL = fit(BrokenPL, t_data, lc_data)
# PL.x_break.min=t_data[0]

print(PL)
t_plot=np.logspace(np.log10(t_data[0]),np.log10(t_data[-1]),100)
plt.figure(figsize=(8,5))
plt.plot(t_data, lc_data, 'ko')
plt.plot(t_plot,PL(t_plot), label='Fit')
plt.xscale("log")
plt.yscale("log")
plt.legend()