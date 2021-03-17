#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sun Feb  7 22:58:25 2021

@author: asif
"""

import numpy as np
import matplotlib.pyplot as plt
from astropy.modeling import models,fitting

x = np.logspace(0.7, 2.3, 500)
f = models.SmoothlyBrokenPowerLaw1D(amplitude=1, x_break=20,
                                    alpha_1=-2, alpha_2=2)

plt.figure()
plt.title("amplitude=1, x_break=20, alpha_1=-2, alpha_2=2")

f.delta = 0.5
plt.loglog(x, f(x), '--', label='delta=0.5')

f.delta = 0.3
plt.loglog(x, f(x), '-.', label='delta=0.3')

f.delta = 0.1
plt.loglog(x, f(x), label='delta=0.1')

plt.axis([x.min(), x.max(), 0.1, 1.1])
plt.legend(loc='lower center')
plt.grid(True)
plt.show()
