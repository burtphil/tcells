#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Thu Sep 27 17:47:00 2018

@author: burt
this script uses a multistep process for a cell to diff
from type Thn to Th1 or Th2 with branching
I found that de Boer et al can only be reproduced when no additional birth rate
is assumed
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
alpha_1 = 20
alpha_2 = 1
beta_1 = float(alpha_1)
beta_2 = float(alpha_2)
t_div = 0.2
rate_birth = 1.0
rate_death = 2.0

# other params
simulation_time = np.arange(0, 3, 0.01)
last_gen = 5
prolif = True


prolif_params = [simulation_time, 
                     alpha_1, 
                     alpha_2, 
                     beta_1, 
                     beta_2, 
                     t_div, 
                     rate_birth, 
                     rate_death,
                     prolif,]

no_prolif_params = [simulation_time, 
                     alpha_1, 
                     alpha_2, 
                     beta_1, 
                     beta_2, 
                     t_div, 
                     rate_birth, 
                     rate_death,
                     False,]

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
    
def th_cell_diff(state, t, alpha_1, alpha_2, beta_1, beta_2, t_div, rate_birth, rate_death, prolif):
        
    th1 = state[:(alpha_1+1)]
    th2 = state[(alpha_1+1):]
      
    dt_th1 = np.ones_like(th1)
    dt_th2 = np.ones_like(th2)
        
    th_states = [th1, th2]
    dt_th_states = [dt_th1, dt_th2]
    rate = [beta_1, beta_2]
    
    p_1 = 1.0
    p_2 = 0
    
    cell_flux = [p_1 * rate_birth, p_2 * rate_birth]
    
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
                dt_state[j] = (r * th_state[j-1]
                - (gamma_dist(t - t_div, alphas[i], betas[i]) * np.exp(-rate_death * t_div))
                - rate_death * th_state[j])
                
            else:
                dt_state[j] = r * th_state[j-1] - rate_death * th_state[j]

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
              rate_birth, 
              rate_death, 
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
                          rate_birth, 
                          rate_death, 
                          prolif))
    
    return state

def N_i(t, i, d, t_div, n_1):
    """
    analytical function N_i(t) for cell numbers that have undergone i divisions
    depends on number of cells that have undergone 1 division (n_1)
    """
    scale_off = 1.
    scale_on = (2 * np.exp(-d * t_div))**(i-1)
    return scale_on * n_1(t - (i-1) * t_div)

def next_gens(simulation_time, rate_death, t_div, first_gen_arr, no_of_next_gens):
    """
    calculate next generations based on array of first generation cells 
    and return as list of arrays
    """
    dummy_list = [N_i(simulation_time,
                     i, 
                     d = rate_death, 
                     t_div = t_div, 
                     n_1 = first_gen_arr) for i in range(2, no_of_next_gens)]
    return dummy_list
#==============================================================================
# run time course simulation (proliferation condtitions)
#==============================================================================
state = run_model(*prolif_params)

# get cells from the first generation
th1_cells = state[:, alpha_1]
th2_cells = state[:, -1]

# interpolate the th1 and th2 cells from first generation
th2_n1 = interp1d(simulation_time, th2_cells, axis = 0, kind='cubic', fill_value = 0, bounds_error = False)
th1_n1 = interp1d(simulation_time, th1_cells, axis = 0, kind='cubic', fill_value = 0, bounds_error = False)

# calculate subsequent generations based on interpolated gen1 cells
th1_gens = next_gens(simulation_time, rate_death, t_div, th1_n1, last_gen)
th2_gens = next_gens(simulation_time, rate_death, t_div, th2_n1, last_gen)
#th2_gens = [N_i(simulation_time, i, d = degradation, t_div = t_div, n_1 = th2_n1) for i in range(2, last_gen)]    
#th1_gens = [N_i(simulation_time, i, d = degradation, t_div = t_div, n_1 = th1_n1) for i in range(2, last_gen)]    

# get sum of all generations
th1_all_gens = list(th1_gens)
th1_all_gens.insert(0, th1_cells)
sum_th1 = [sum(x) for x in zip(*th1_all_gens)]

th2_all_gens = list(th2_gens)
th2_all_gens.insert(0, th2_cells)
sum_th2 = [sum(x) for x in zip(*th2_all_gens)]

#==============================================================================
# time course (no proliferation conditions)
#==============================================================================
no_prolif_state = run_model(*no_prolif_params)
th1_no_prolif = no_prolif_state[:, alpha_1]
th2_no_prolif = no_prolif_state[:, -1] 


xlabel = "time"
ylabel = "% Th cells"
"""
# plot only first generation
fig, ax = plt.subplots(1, 1, figsize = (5,4))
ax.plot(simulation_time, th1_no_prolif, c = "tab:blue", label = "Th1", linestyle = "--")
ax.plot(simulation_time, th2_no_prolif, c = "tab:red", label = "Th2")
ax.set_title("no proliferation \n" + r"$\alpha_{Th1}=1$, $\alpha_{Th2}=10$")
ax.set_xlabel(xlabel)
ax.set_ylabel(ylabel)
ax.legend()
plt.tight_layout()
#fig.savefig(save_path + "th1_th2_no_proliferation.pdf", bbox_inches="tight")

#==============================================================================
# plot cells (proliferation conditions)
#==============================================================================
xlabel = "time"
ylabel = "#cells (normalized)"


fig, ax = plt.subplots(1, 3, figsize = (15,4))
# plot only first generation
ax[0].plot(simulation_time, th1_cells, c = "tab:blue", label = "Th1", linestyle = "--")
ax[0].plot(simulation_time, th2_cells, c = "tab:red", label = "Th2")
ax[0].set_title('cells after first division')
ax[0].set_xlabel(xlabel)
ax[0].set_ylabel(ylabel)
ax[0].legend()

# plot all generations individually
#fig, ax = plt.subplots(1, 1, figsize = (5,4))
palette = ['k', 
          'tab:gray',
          'tab:green',
          'tab:orange', 
          ]

# note that due to the append, first generations are last in the list
generations = ['1st div',
               '2nd div',
               '3rd div',
               '4th div',
               ]
for i in range(len(th1_all_gens)):
    col = palette[i]
    th1_gen = th1_all_gens[i]
    th2_gen = th2_all_gens[i]
    label = generations[i]
    
    ax[1].plot(simulation_time, th1_gen, c = col, linestyle = "--")
    ax[1].plot(simulation_time, th2_gen, c = col, label = label)

ax[1].set_title("-- Th1, - Th2")
ax[1].set_xlabel(xlabel)
ax[1].set_ylabel(ylabel)
ax[1].legend()

#plt.tight_layout()

# plot total cell numbers
#fig, ax = plt.subplots(1,1, figsize = (5,4))
ax[2].plot(simulation_time, sum_th1, c = "tab:blue", label = "Th1", linestyle = "--")
ax[2].plot(simulation_time, sum_th2, c = "tab:red", label = "Th2")
ax[2].set_title(r"total cell numbers")
ax[2].legend()
ax[2].set_xlabel(xlabel)
ax[2].set_ylabel(ylabel)
fig.suptitle(r"th1 th2 model with proliferation, $\alpha_{Th1}=1$, $\alpha_{Th2}=10$")
plt.tight_layout()
plt.subplots_adjust(top=0.8)

#fig.savefig(save_path + "th1_th2_proliferation.pdf", bbox_inches="tight")
# plot individual generations together with total cell numbers
fig, ax = plt.subplots(1, 1, figsize = (5,4))
for th1_gen, th2_gen in zip(th1_all_gens, th2_all_gens):
    ax.plot(simulation_time, th1_gen, c = "tab:blue", alpha = 0.5, linestyle = "--")
    ax.plot(simulation_time, th2_gen, c = "tab:red", alpha = 0.5)

ax.plot(simulation_time, sum_th1, c = "tab:blue", linestyle = "--", label = "Th1 (total)")
ax.plot(simulation_time, sum_th2, c = "tab:red", label = "Th2 (total)")
ax.legend()
ax.set_xlabel(xlabel)
ax.set_ylabel(ylabel)

plt.tight_layout()
"""
fig, ax = plt.subplots(1, 1, figsize = (5,4))
for th1_gen, th2_gen in zip(th1_all_gens, th2_all_gens):
    ax.plot(simulation_time, th1_gen)
    #ax.plot(simulation_time, th2_gen, c = "tab:red", alpha = 0.5)

#ax.plot(simulation_time, sum_th1, c = "tab:blue", linestyle = "--", label = "Th1 (total)")
#ax.plot(simulation_time, sum_th2, c = "tab:red", label = "Th2 (total)")
ax.legend()
ax.set_xlabel(xlabel)
ax.set_ylabel(ylabel)
plt.tight_layout()