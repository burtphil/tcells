#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Thu Sep 27 17:47:00 2018

@author: burt
"""
import numpy as np
from scipy.integrate import odeint
import matplotlib.pyplot as plt

def th_cell_diff(state,t, alpha_1, alpha_2, beta_1, beta_2, th0_influx = 2.0, degradation = 1.0):
        
    th1 = state[:(alpha_1+1)]
    th2 = state[(alpha_1+1):]
      
    dt_th1 = np.ones_like(th1)
    dt_th2 = np.ones_like(th2)
        
    th_states = [th1, th2]
    dt_th_states = [dt_th1, dt_th2]
    rate = [beta_1, beta_2]
    
    p_1 = 0.5
    p_2 = 0.5
    
    cell_flux = [p_1 * th0_influx, p_2 * th0_influx]
        
    ### differential equations depending on the number on intermediary states
    # loop over states of th1 and th2 cells
    for i in range(2):
        th_state = th_states[i]
        dt_state = dt_th_states[i]
        r = rate[i]
        flux = cell_flux[i]
            
    # calculate derivatives
        for j in range(len(th_state)):
            if j == 0:
                dt_state[j] = flux - r * th_state[j]
                
            elif j != (len(th_state) - 1):
                dt_state[j] = r * (th_state[j-1] - th_state[j])
                
            else:
                dt_state[j] = r * th_state[j-1] - degradation * th_state[j]

    # assign new number of naive cells based on the number of present naive cells that were designated th1_0 or th2_0
    #dt_th0 = state[0]    
    # return cell states
    dt_state = np.concatenate([dt_th1, dt_th2])
    
    return dt_state  

def run_model(parameters):
    """ 
    run ode model with parameters from params file
    needs shape and rate params as input as well as simulation time
    """
    alpha_1, alpha_2, beta_1, beta_2, simulation_time = parameters

    # define initial conditions based on number of intermediary states
    no_of_states = int(alpha_1 + alpha_2 + 2)
    ini_cond = np.zeros(no_of_states)
    print no_of_states
    print alpha_1
    print alpha_2
    
    initial_cells = 1
    ini_cond[0] = initial_cells
    ini_cond[alpha_1+1] = initial_cells

    
    state = odeint(th_cell_diff, ini_cond, simulation_time, 
                   args =(alpha_1, alpha_2, beta_1, beta_2,))
    
    return state

#==============================================================================
# import parameters
#==============================================================================
simulation_time = np.arange(0, 5, 0.01)
alpha_1 = 10
alpha_2 = 1
beta_1 = float(alpha_1)
beta_2 = float(alpha_2)

conceptual_params = [alpha_1, alpha_2, beta_1, beta_2, simulation_time]
initial_cells= 1

#==============================================================================
# run time course simulation
#==============================================================================
th1_th2_model = run_model(parameters = conceptual_params)

state = th1_th2_model

# plot time course
#alpha_params = list(conceptual_params)
alpha_params = list(conceptual_params)
norm = initial_cells / 100

fig, ax = plt.subplots(1, 1, figsize = (5,4))

th1_cells = state[:, alpha_1]
th2_cells = state[:, -1]
#ax.plot(simulation_time, th1_cells)
ax.plot(simulation_time, th1_cells)
ax.plot(simulation_time, th2_cells)

plt.tight_layout()

