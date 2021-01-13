# -*- coding: utf-8 -*-
"""
Created on Sun Aug 23 20:47:58 2020

@author: Asif
"""

import numpy as np
import Flux_calc_final as fc
from matplotlib import pyplot as plt
from matplotlib.animation import FuncAnimation
plt.style.use('seaborn-pastel')


fig = plt.figure()
ax = fig.add_subplot(1, 1, 1)
line, = ax.plot([], [], lw=3)
# ax.set_xscale("log")
ax.set_ylim((1e-23,1))
ax.set_yscale("log")
ax.set_xlim((1e4,1e35))
ax.set_xscale("log")

t=np.logspace(-3,2)
freqs=np.logspace(1,35,100)
F_0=1
(gamma_m,gamma_c,GAMMA,B,freq_m,freq_c,F)= fc.calc_flux(t,freqs,F_0)

def init():
    line.set_data(t, F[:,0])
    return line,
def animate(i):
    x = freqs
    y = F[i,:]
    line.set_data(x, y)
    # ax.scatter(freq_c[i],fc.calc_flux(t[i],freq_c,F_0))
    print(i)
    return line,

anim = FuncAnimation(fig, animate, init_func=init,
                                frames=len(t), interval=80, blit=True)

# anim = FuncAnimation(fig, animate, init_func=init,
#                                frames=200, interval=20, blit=True)


anim.save('F vs freq s.gif', writer='imagemagick')