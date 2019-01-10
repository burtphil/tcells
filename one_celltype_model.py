#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Tue Jan  8 11:17:57 2019

@author: burt
new model with influx of naive cells that is modulated by probabilities
double check because branching probabilities might be set to a fixed value within function
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
simulation_time = np.arange(0, 3, 0.01)
alpha_1 = 20
alpha_2 = 1
beta_1 = float(alpha_1)
beta_2 = float(alpha_2)

initial_cells = 1
t_div = 0.1
th0_influx = 0.1
degradation = 1.0
prolif = True
no_prolif = False

#==============================================================================
# summarize params for model input
#==============================================================================
prolif_params = [simulation_time, 
                     alpha_1, 
                     alpha_2, 
                     beta_1, 
                     beta_2, 
                     t_div, 
                     th0_influx, 
                     degradation,
                     prolif,]

no_prolif_params = [simulation_time, 
                     alpha_1, 
                     alpha_2, 
                     beta_1, 
                     beta_2, 
                     t_div, 
                     th0_influx, 
                     degradation,
                     no_prolif,]

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
    
def th_cell_diff(state, t, alpha_1, alpha_2, beta_1, beta_2, t_div, th0_influx, degradation, prolif):
        
    th1 = state[:(alpha_1+1)]
    th2 = state[(alpha_1+1):]
      
    dt_th1 = np.ones_like(th1)
    dt_th2 = np.ones_like(th2)
        
    th_states = [th1, th2]
    dt_th_states = [dt_th1, dt_th2]
    rate = [beta_1, beta_2]
    
    # note that cell_flux also determines model output
    p_1 = 1.0
    p_2 = 0
    
    cell_flux = [p_1 * th0_influx, p_2 * th0_influx]
    
    alphas = [alpha_1, alpha_2]
    betas = [beta_1, beta_2]
    ### differential equations depending on the number on intermediary states
    # loop over states of th1 and th2 cells
    for i in range(2):
        th_state = th_states[i]
        dt_state = dt_th_states[i]
        r = rate[i]
        flux = cell_flux[i]
        #print flux
            
    # calculate derivatives
        for j in range(len(th_state)):
            if j == 0:
                dt_state[j] = flux - r * th_state[j]
                
            elif j != (len(th_state) - 1):
                dt_state[j] = r * (th_state[j-1] - th_state[j])
                
            elif prolif == True:
                dt_state[j] = (2 * r * th_state[j-1]
                - (gamma_dist(t - t_div, alphas[i], betas[i]) * np.exp(-degradation * t_div))
                - degradation * th_state[j])
                
            else:
                dt_state[j] = r * th_state[j-1] - degradation * th_state[j]

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
              t_div, 
              th0_influx, 
              degradation, 
              prolif):
    """ 
    run ode model with parameters from params file
    needs shape and rate params as input as well as simulation time
    """

    # define initial conditions based on number of intermediary states
    no_of_states = int(alpha_1 + alpha_2 + 2)
    ini_cond = np.zeros(no_of_states)
    
    initial_cells = 1
    ini_cond[0] = initial_cells
    ini_cond[alpha_1+1] = initial_cells

    
    state = odeint(th_cell_diff, 
                   ini_cond, 
                   simulation_time, 
                   args =(alpha_1, 
                          alpha_2, 
                          beta_1, 
                          beta_2, 
                          t_div, 
                          th0_influx, 
                          degradation, 
                          prolif))
    
    return state

#==============================================================================
# analysis
#==============================================================================
no_prolif_state = run_model(*no_prolif_params)
th1_no_prolif = no_prolif_state[:, alpha_1]
th2_no_prolif = no_prolif_state[:, -1] 


xlabel = "time"
ylabel = "% Th cells"

# plot only first generation
fig, ax = plt.subplots(1, 1, figsize = (5,4))
ax.plot(simulation_time, th1_no_prolif, c = "tab:blue", label = "Th1", linestyle = "-")
#ax.plot(simulation_time, th2_no_prolif, c = "tab:red", label = "Th2")
#ax.set_title("no proliferation \n" + r"$\alpha_{Th1}=1$, $\alpha_{Th2}=10$")
ax.set_title("no proliferation")
ax.set_xlabel(xlabel)
ax.set_ylabel(ylabel)
ax.legend()
plt.tight_layout()