#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Thu Sep 27 17:47:00 2018

@author: burt
include yates proliferation assumption in det model
"""
import numpy as np
from scipy.integrate import odeint
import matplotlib.pyplot as plt
from scipy.special import gamma
from scipy.interpolate import interp1d
import seaborn as sns


#==============================================================================
# choose plotting style and saving path
#==============================================================================
sns.set(context = "talk", style = "ticks")
save_path = "/home/burt/Documents/tcell_project/figures/"
#==============================================================================
# import parameters
#==============================================================================
# rates
alpha_1 = 10
alpha_2 = 5
beta_1 = float(alpha_1)
beta_2 = float(alpha_2)
rate_birth = 0
rate_death = 1.0

# other params
simulation_time = np.arange(0, 5, 0.01)

alpha_prolif = 10
beta_prolif = 10.

prolif_params = [simulation_time, 
                     alpha_1, 
                     alpha_2, 
                     beta_1, 
                     beta_2, 
                     alpha_prolif,
                     beta_prolif,
                     rate_birth, 
                     rate_death,
                ]

colors = ['tab:blue', 
          'tab:orange', 
          'tab:green', 
          'tab:red', 
          'tab:purple', 
          'tab:brown', 
          'tab:pink', 
          'tab:gray', 
          'tab:olive', 
          'tab:cyan']
#==============================================================================
# functions
#==============================================================================
def gamma_dist(t, alpha, beta):
    """
    regular gamma distrubtion
    """
    if t >= 0:
        return np.exp(-beta*t)*(beta**alpha)*(t**(alpha-1))/gamma(alpha)
    else:
        return 0
    
y = [gamma_dist(t, 30, 30) for t in simulation_time]
fig, ax = plt.subplots()
ax.plot(simulation_time, y)

def th_cell_diff(state, t, alpha_1, alpha_2, beta_1, beta_2, alpha_prolif, beta_prolif, rate_birth, rate_death):
        
    th1 = state[:(alpha_1+alpha_prolif)]
    th2 = state[(alpha_1+alpha_prolif):]
      
    dt_th1 = np.zeros_like(th1)
    dt_th2 = np.zeros_like(th2)
        
    th_states = [th1, th2]
    dt_th_states = [dt_th1, dt_th2]
    rate = [beta_1, beta_2]
    
    p_1 = 1.0
    p_2 = 0
    
    cell_flux = [p_1 * rate_birth, p_2 * rate_birth]
    
    alphas = [alpha_1, alpha_2]
    rate = [beta_1, beta_2]
    ### differential equations depending on the number on intermediary states
    # loop over states of th1 and th2 cells
    for i in range(2):
        th_state = th_states[i]
        dt_state = dt_th_states[i]
        r = rate[i]
        alpha = alphas[i]
        flux = cell_flux[i]
        #print flux
            
    # calculate derivatives
        for j in range(len(th_state)):
            if j == 0:
                dt_state[j] = flux - r * th_state[j]
                
            elif j < alpha:
                dt_state[j] = r * (th_state[j-1] - th_state[j])
                
            elif j == alpha:
                dt_state[j] = 2 * (r * th_state[j-1] + beta_prolif * th_state[-1]) - (rate_death + beta_prolif) * th_state[j]
            
            else:
                assert j > alpha
                dt_state[j] = beta_prolif * (th_state[j-1] - th_state[j])

    # assign new number of naive cells based on the number of present naive cells that were designated th1_0 or th2_0
    #dt_th0 = state[0]    
    # return cell states
    dt_state = np.concatenate([dt_th1, dt_th2])
    
    return dt_state  

def run_model(simulation_time, 
              alpha_1, 
              alpha_2, 
              beta_1, 
              beta_2, 
              alpha_prolif,
              beta_prolif,
              rate_birth, 
              rate_death,
              ):
    """ 
    run ode model with parameters from params file
    needs shape and rate params as input as well as simulation time
    """

    # define initial conditions based on number of intermediary states
    no_of_states = int(alpha_1 + alpha_2 + alpha_prolif + alpha_prolif)
    ini_cond = np.zeros(no_of_states)
    
    initial_cells = 1
    ini_cond[0] = initial_cells
    ini_cond[alpha_1+alpha_prolif] = 0

    
    state = odeint(th_cell_diff, 
                   ini_cond, 
                   simulation_time, 
                   args =(alpha_1, 
                          alpha_2, 
                          beta_1, 
                          beta_2, 
                          alpha_prolif,
                          beta_prolif,
                          rate_birth, 
                          rate_death)
                   )
    
    return state

#==============================================================================
# run time course simulation (proliferation condtitions)
#==============================================================================
state = run_model(*prolif_params)

# get cells from the first generation
th1_cells = state[:, alpha_1:(alpha_1+alpha_prolif)]
th2_cells = state[:, (alpha_1+alpha_prolif+alpha_2):]


th1_all_cells = np.sum(th1_cells, axis = 1)
th2_all_cells = np.sum(th2_cells, axis = 1)

fig, ax = plt.subplots()
ax.plot(simulation_time, th1_all_cells)
ax.plot(simulation_time, th2_all_cells)