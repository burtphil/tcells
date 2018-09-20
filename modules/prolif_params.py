#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Wed Sep 12 16:51:52 2018

@author: burt
"""
import numpy as np
### set up params (roughly Hammer et al params)
initial_cells = 1.
death_rate = 1./24
division_time = 10.
non_dividing_cells = 0.05/np.exp(-death_rate*(6*24))

start = 0
stop = 200
stepsize = 0.01
simulation_time = np.arange(start, stop, stepsize)

alpha = 10.
mean = 55.
beta = alpha/ mean

mu = 137.66
sigma = 0.11
c = 9301
div0 = -86.2
div_t = 12.38
d = 0.047

n = 6