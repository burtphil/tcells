#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Fri Jan 11 10:25:20 2019

@author: burt
"""

import numpy as np
import matplotlib.pyplot as plt
from scipy.integrate import odeint

def prolif(state, t, diff_rate, birth_rate, death_rate):

    dt_state = np.zeros_like(state)
    
    for i in range(len(state)):
        if i == 0:
            dt_state[i] = birth_rate +  2 * diff_rate * state[-1] - (diff_rate + death_rate) * state[i]
        else:
            dt_state[i] = diff_rate * (state[i-1]-state[i])
    
    return dt_state

def prolif2(state, t, diff_rate, birth_rate, death_rate):

    dt_state = birth_rate + (2 - death_rate) * state
    
    return dt_state

stage_no = 10
y0 = np.zeros(stage_no)
y0[0] = 1
diff_rate = 10.
birth_rate = 0
death_rate = 0

sim_time = np.arange(0, 5, 0.01)
state = odeint(prolif, y0, sim_time, args = (diff_rate, birth_rate, death_rate))

all_cells = np.sum(state, axis = 1)

fig, ax = plt.subplots()
plt.plot(sim_time,state)

fig, ax = plt.subplots()
plt.plot(sim_time,all_cells)