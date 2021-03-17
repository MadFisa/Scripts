#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Mar  3 09:36:48 2021

@author: asif
"""

import numpy as np
import matplotlib.pyplot as plt

evt=np.loadtxt("mltnst_dim6_ev.dat")
live=np.loadtxt("mltnst_dim6_phys_live.points")
eq=np.loadtxt("mltnst_dim6_post_equal_weights.dat")