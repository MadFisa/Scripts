#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Feb  3 10:50:10 2021

@author: asif
"""

import matplotlib.pyplot as plt

from scipy.optimize import curve_fit
import numpy as np



def func(x, a, b, c):

    return a * np.exp(-b * x) + c

#Define the data to be fit with some noise:


xdata = np.linspace(0, 4, 50)

y = func(xdata, 2.5, 1.3, 0.5)

np.random.seed(1729)

y_noise = 0.2 * np.random.normal(size=xdata.size)

ydata = y + y_noise

plt.plot(xdata, ydata, 'b-', label='data')

#Fit for the parameters a, b, c of the function func:


popt, pcov = curve_fit(func, xdata, ydata)



plt.plot(xdata, func(xdata, *popt), 'r-',

         label='fit: a=%5.3f, b=%5.3f, c=%5.3f' % tuple(popt))

#Constrain the optimization to the region of 0 <= a <= 3, 0 <= b <= 1 and 0 <= c <= 0.5:


popt, pcov = curve_fit(func, xdata, ydata, bounds=(0, [3., 1., 0.5]))



plt.plot(xdata, func(xdata, *popt), 'g--',

         label='fit: a=%5.3f, b=%5.3f, c=%5.3f' % tuple(popt))



plt.xlabel('x')

plt.ylabel('y')

plt.legend()

plt.show()
