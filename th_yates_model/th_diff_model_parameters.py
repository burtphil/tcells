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

beta_1 = 10.
beta_2 = 1.

alpha_prolif = 20
beta_prolif = 20

rate_birth = 2. 
rate_death = 2.0

# define initial conditions based on number of intermediary states
th1_0 = 0
th2_0 = 0

no_of_states = int(alpha_1 + alpha_2 + alpha_prolif + alpha_prolif)
y0 = np.zeros(no_of_states)

y0[0] = th1_0
y0[alpha_1+alpha_prolif] = th2_0

rate_ifn = 10.0
rate_il4 = 10.0

rate_cytos = [rate_ifn, rate_il4]

fb_th1 = [0, 0, 0]
fb_th2 = [0, 0, 0]

hill_th1 = [1., 1., 1.]
hill_th2 = [1., 1., 1.]

K_th1 = [1., 1., 1.]
K_th2 = [1., 1., 1.]

conc_il12 = 0