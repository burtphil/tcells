#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Fri Dec 14 10:39:00 2018

@author: burt
new model: cells only increase due to scaling reaction
plot only first generation
"""
import numpy as np
import matplotlib.pyplot as plt
from scipy.special import gammainc
import seaborn as sns
import pandas as pd
from scipy.integrate import odeint
from scipy.special import gamma
from scipy.interpolate import interp1d

sns.set(context = "talk", style = "ticks")

#==============================================================================
# params
#==============================================================================
start = 0
stop = 2.5
nsteps = 1000
simulation_time = np.linspace(start, stop, nsteps)
last_gen = 2

# stoc params
nstates = 8
nsim = 300
n_dummies = 40
t_div_stoc = 0.2

# rates
alpha_1 = 20
alpha_2 = 1
beta_1 = float(alpha_1)
beta_2 = float(alpha_2)
rate_birth = 0
rate_death = 2.0

# other params
ncells = 1
t_div = t_div_stoc

no_prolif = False
prolif = True

initial_cells = 1

colors = ['tab:blue', 
          'tab:orange', 
          'tab:green', 
          'tab:red', 
          ]

prolif_params = [simulation_time, 
                     alpha_1, 
                     alpha_2, 
                     beta_1, 
                     beta_2, 
                     t_div, 
                     rate_birth, 
                     rate_death,
                     prolif,]

stoc_params = [start, 
               stop, 
               nsteps, 
               ncells, 
               nstates, 
               rate_birth, 
               alpha_1, 
               beta_1, 
               rate_death, 
               n_dummies,
               t_div_stoc,]
#==============================================================================
# functions
#==============================================================================
def N_i(t, i, d, t_div, n_1):
    """
    analytical function N_i(t) for cell numbers that have undergone i divisions
    depends on number of cells that have undergone 1 division (n_1)
    """
    scale_off = 1.
    scale_on = (2 * np.exp(-d * t_div))**(i-1)
    return scale_on * n_1(t - (i-1) * t_div)

def next_gens(simulation_time, rate_deatj, t_div, first_gen_arr, no_of_next_gens):
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

def gamma_cdf(t, alpha, beta):
    dummy = t*beta
    return gammainc(alpha,dummy)      


def make_cells(ncells, nsteps, nstates, rate_death):
    """
    returns an array of cells for given number of
    desired cells, steps and cell states
    initializes random numbers probability
    of cell state transition
    """
    cells = np.zeros((ncells, nsteps, nstates))
    for i in range(ncells):        
        # random number for th0 to th_eff transition
        cells[i, :, 2] = np.random.rand()  
        cells[i, :, 3] = np.random.exponential(1. / rate_death)
    return cells

def stoc_simulation(start, stop, nsteps, ncells, nstates, rate_birth, alpha_diff, beta_diff, rate_death, n_dummies, t_div):
    """
    run a simulation of the th1 th2 models   
    this function uses the a cumulative gamma fct integral to calculate transition probability
    """
    # define state indices for the nstates array
    cell_state_idx = 0
    prob_diff_idx = 1
    rnd_diff_idx = 2
    death_time = 3
    diff_time = 4
    cell_gen = 5
    
    # define cell_state_idx state values
    th0_cell = 0
    dead_cell = 2

    # initialize some dummy dead cells t
    cells = make_cells(n_dummies+ncells, nsteps, nstates, rate_death)
    time = np.linspace(start, stop, nsteps)
    cells[:n_dummies,:,cell_state_idx] = dead_cell

    # loop over each time point
    for i in range(len(time)-1):
        t = time[i]
        t_new = time[i+1]

        # get all cells for current time step
        cell_j = cells[:,i,:]

        # cell number should be number of alive cells not total number              
        cell_number = cell_j.shape[0]
     
        # loop over each cell for current time point
        for j in range(cell_number):
                       
            cell = cell_j[j,:] 

            #check if cell divides for the first time
            if cell[cell_gen] == 0 and cell[cell_state_idx] == th0_cell:
                if cell[rnd_diff_idx] > cell[prob_diff_idx]:
                    cell[prob_diff_idx] = (cell[prob_diff_idx]+
                        (gamma_cdf(t_new, alpha_diff, beta_diff)-
                         gamma_cdf(t, alpha_diff, beta_diff)))
                else:
                    cell[cell_state_idx] = th0_cell
                    cell[diff_time] = t
                    cell[cell_gen] = 1
                                                            
            # check if division time is reached for differentiated cells
            if cell[cell_gen] > 0 and cell[cell_state_idx] == th0_cell:
                if t - cell[diff_time] > t_div:
                    cell[cell_state_idx] = dead_cell
            
                if cell[death_time] < t - cell[diff_time]:
                    cell[cell_state_idx] = dead_cell
            
            cells[j,i+1,:] = cell
        
    return [cells, time]     


def get_cells(cells, time):
    """
    input: cells and time, which is output from the function: run_stochastic_simulation
    use this for the model with precursor state times included
    """
    all_cells = cells[:,:, 0]
          
    th0_cells = np.sum(all_cells == 0, axis = 0)
    th1_cells = np.sum(all_cells == 1, axis = 0)
    dead_cells = np.sum(all_cells == 2, axis = 0)
   
    return th0_cells, th1_cells, dead_cells

def get_gens(cells, time, gen_no):
    """
    input: cells and time, which is output from the function: run_stochastic_simulation
    use this for the model with precursor state times included
    """
    cell_states = cells[:,:, 0]
    cell_gens = cells[:,:,5]
    generations = [np.sum((cell_states == 0) & (cell_gens == i), axis = 0) for i in range(1, gen_no)]      
    return generations

#==============================================================================
# det model functions
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
    
    initial_cells = 1.
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

# =============================================================================
# run stochastic simulation
# =============================================================================
simulation = [stoc_simulation(*stoc_params) for i in range(nsim)]

#==============================================================================
# plot generations in stoc model
#==============================================================================
th_gens = [get_gens(*simu, gen_no = last_gen) for simu in simulation]
# make dummy labels for pandas df
label = [str(i)+"th_gen" for i in range(1, last_gen)]
labels = [label for _ in range(nsim)]
flat_labels = [item for sublist in labels for item in sublist]

# make a df for single step stoc simulation (adjusted time)
flat_cells = [item for sublist in th_gens for item in sublist]
df = pd.DataFrame.from_records(flat_cells).transpose()
df.columns = flat_labels

df["time"] = simulation_time
df_long = df.melt(id_vars = ["time"])

fig, ax = plt.subplots(1, 1, figsize =(5,4))
stoc_plot = sns.lineplot(x = "time", 
                         y = "value", 
                         data = df_long,
                         hue = "variable",
                         ci = None,
                         ax = ax,
                         legend = False)
ax.set_ylabel("$n_{cells}$")
ax.set_xlim(0, simulation_time[-1])
plt.tight_layout()


#==============================================================================
# run time course simulation (proliferation condtitions)
#==============================================================================
state = run_model(*prolif_params)

# get cells from the first generation
th1_cells = state[:, alpha_1]

# interpolate the th1 and th2 cells from first generation
th1_n1 = interp1d(simulation_time, th1_cells, axis = 0, kind='cubic', fill_value = 0, bounds_error = False)

# calculate subsequent generations based on interpolated gen1 cells
th1_gens = next_gens(simulation_time, rate_death, t_div, th1_n1, last_gen)

# get sum of all generations
th1_all_gens = list(th1_gens)
th1_all_gens.insert(0, th1_cells)
sum_th1 = [sum(x) for x in zip(*th1_all_gens)]

xlabel = "time"
ylabel = "% Th cells"

fig, ax = plt.subplots(1, 1, figsize = (5,4))
for th1_gen in th1_all_gens:
    ax.plot(simulation_time, th1_gen)

ax.legend()
ax.set_xlabel(xlabel)
ax.set_ylabel(ylabel)
plt.tight_layout()

#==============================================================================
# stoc and det model together
#==============================================================================
fig, ax = plt.subplots(1, 1, figsize =(5,4))
stoc_plot = sns.lineplot(x = "time", 
                         y = "value", 
                         data = df_long,
                         hue = "variable",
                         ci = "sd",
                         ax = ax,
                         legend = False)

for i, th1_gen in enumerate(th1_all_gens):
    ax.plot(simulation_time, th1_gen, c = colors[i], linestyle = "--")
ax.set_ylabel("$n_{cells}$")
ax.set_xlim(0, simulation_time[-1])
plt.tight_layout()