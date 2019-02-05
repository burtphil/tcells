#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Thu Jan 24 12:23:37 2019

@author: burt
"""
import numpy as np

start = 0
stop = 20
stepsize = 0.001

time = np.arange(start, stop, stepsize)

alpha_1 = 10
alpha_2 = 1

mean_th1 = 1.
mean_th2 = 1.
beta_1 = alpha_1 / mean_th1
beta_2 = alpha_2 / mean_th2

mean_prolif = 0.25
alpha_prolif = 10
beta_prolif = alpha_prolif / mean_prolif

rate_birth = 1.0 
rate_diff_th0 = 0.5
rate_death_th0 = 0.1
rate_death = 4.0

# define initial conditions based on number of intermediary states
th1_0 = 0
th2_0 = 0

cells_0 = [th1_0, th2_0]

rate_ifn = 10.0
rate_il4 = 10.0

rate_cytos = [rate_ifn, rate_il4]

fb_th1 = [1, 0.5, 1]
fb_th2 = [0.5, 1, 1]

hill_th1 = [1., 1., 1.]
hill_th2 = [1., 1., 1.]

K_th1 = [1., 1., 1.]
K_th2 = [1., 1., 1.]

conc_il12 = 0