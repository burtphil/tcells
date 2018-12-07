#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Tue Oct  2 12:43:02 2018

@author: burt
new model implementation
"""

import numpy as np
from scipy.integrate import odeint
import matplotlib.pyplot as plt

def th_cell_diff(state, t, alpha_1, alpha_2, beta_1, beta_2, th0_influx = 0, degradation = 1):
        
    th_0 = state[0]
    p_1 = 0.5
    p_2 = 0.5
    #print th1_0,th2_0
    # assign th1 states and th2 states from initial vector based on chain length alpha
    th1 = state[1:(alpha_1+2)]
    th2 = state[(alpha_1+2):]
      
    dt_th1 = np.ones_like(th1)
    dt_th2 = np.ones_like(th2)
        
    th_states = [th1,th2]
    dt_th_states = [dt_th1,dt_th2]
    rate = [beta_1,beta_2]
    prob = [p_1,p_2]
        
    ### differential equations depending on the number on intermediary states
    # loop over states of th1 and th2 cells
    for i in range(2):
        th_state = th_states[i]
        dt_state = dt_th_states[i]
        r = rate[i]
        p = prob[i]
            
    # calculate derivatives
        for j in range(len(th_state)):
            if j == 0:
                dt_state[j] = p * th_0 - r * th_state[j]
            elif j != (len(th_state)-1):
                dt_state[j] = r * (th_state[j-1] - th_state[j])
            else:
                dt_state[j] = r * th_state[j-1]- degradation * th_state[j]

    # assign new number of naive cells based on the number of present naive cells that were designated th1_0 or th2_0
    #dt_th0 = state[0]
    dt_th0 = 0
    
    #print type(dt_th0),type(dt_th1)
    # return cell states
    dt_state = np.concatenate(([dt_th0], dt_th1, dt_th2))
    
    return dt_state  

alpha_1 = 10
alpha_2 = 1
beta_1 = 10.
beta_2 = 1.

simulation_time = np.arange(0, 5, 0.01)

initial_cells = 10000
no_of_states = int(alpha_1 + alpha_2 + 3)
ini_cond = np.zeros(no_of_states)
ini_cond[0] = initial_cells

state = odeint(th_cell_diff, ini_cond, simulation_time, args =(alpha_1, alpha_2, beta_1, beta_2,))

labels = ["Thn","th1_0","th1_1","th1_2","th1_3","th2_0","th2_1"]
#plt.plot(simulation_time,state)
#plt.legend(labels)

norm = initial_cells / 100
th0_cells = state[:,0] / norm
th1_cells = state[:,alpha_1 + 1] / norm
th2_cells = state[:,-1] / norm

fig, ax = plt.subplots(1, 1, figsize = (5,5))
ax.plot(simulation_time, th0_cells, "k",
    simulation_time, th1_cells, "tab:blue",
    simulation_time, th2_cells, "tab:red")
ax.set_yticks([0, 50, 100])
ax.set_ylim([0, 100])
ax.set_xlim([simulation_time[0], simulation_time[-1]])
ax.set_xlabel("Time [h]")
ax.set_ylabel("% Th cells")
ax.legend(["Thn", "Th1", "Th2"])
ax.set_title(r"$p_1=p_2, \alpha_1=2, \alpha_2=1, mean = const.$")
