#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Thu Jan 24 12:23:37 2019

@author: burt
"""
import numpy as np

start = 0
stop = 3
stepsize = 0.001

time = np.arange(start, stop, stepsize)

alpha_1 = 10
alpha_2 = 1

beta_1 = float(alpha_1)
beta_2 = float(alpha_2)

alpha_prolif = 20
beta_prolif = 20

rate_birth = 2.0 
rate_death = 2.0

# define initial conditions based on number of intermediary states
th1_0 = 0
th2_0 = 0

cells_0 = [th1_0, th2_0]

rate_ifn = 1.0
rate_il4 = 1.0

rate_cytos = [rate_ifn, rate_il4]

fb_th1 = [10, -1, 0]
fb_th2 = [-1, 10, 0]

hill_th1 = [1., 1., 1.]
hill_th2 = [1., 1., 1.]

K_th1 = [1., 1., 1.]
K_th2 = [1., 1., 1.]

conc_il12 = 0