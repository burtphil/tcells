#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Fri Dec 14 10:39:00 2018

@author: burt
--> thn --> th1 --> death
here, I try to model only first generation and only calculate the
rest through deterministic scaling equation
"""
import numpy as np
import matplotlib.pyplot as plt
from scipy.special import gammainc
import seaborn as sns
import pandas as pd
from det_th1_th2_model import run_model
from scipy.interpolate import interp1d

sns.set(context = "talk", style = "ticks")

#==============================================================================
# params
#==============================================================================
start = 0
stop = 3
nsteps = 10000
simulation_time = np.linspace(start, stop, nsteps)
last_gen = 2

# stoc params
nstates = 6
nsim = 50
n_dummies = 50
t_div = 0.2

# rates
alpha_1 = 20
alpha_2 = 1
beta_1 = float(alpha_1)
beta_2 = float(alpha_2)
rate_birth = 1.0
rate_death = 2.0

# other params
ncells = 1

no_prolif = False
prolif = True

no_prolif_params = [simulation_time, 
                     alpha_1, 
                     alpha_2, 
                     beta_1, 
                     beta_2, 
                     t_div, 
                     rate_birth, 
                     rate_death,
                     no_prolif,]

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
               t_div,]
#==============================================================================
# functions
#==============================================================================
def N_i(t, i, d, t_div, n_1):
    """
    analytical function N_i(t) for cell numbers that have undergone i divisions
    depends on number of cells that have undergone 1 division (n_1)
    """
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

def exp_cdf(t, rate):
    return 1 - np.exp(-t * rate)


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
    revival_time = 5
    
    # define cell_state_idx state values
    th0_cell = 0
    th1_cell = 1
    dead_cell = 2

    # initialize some dummy dead cells t
    cells = make_cells(n_dummies+ncells, nsteps, nstates, rate_death)
    time = np.linspace(start, stop, nsteps)
    cells[:n_dummies,:,cell_state_idx] = dead_cell

    
    # initialize random number and probability for birth process
    rnd_birth = np.random.rand()
    p_birth = 0
    t0 = 0
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

            #check if cell differentiates
            if cell[cell_state_idx] == th0_cell:
                if cell[rnd_diff_idx] > cell[prob_diff_idx]:
                    cell[prob_diff_idx] = (cell[prob_diff_idx]+
                        (gamma_cdf(t_new-cell[revival_time], alpha_diff, beta_diff)-
                         gamma_cdf(t-cell[revival_time], alpha_diff, beta_diff)))
                else:
                    cell[cell_state_idx] = th1_cell
                    cell[diff_time] = t

            
            if cell[cell_state_idx] == th1_cell:
                
                # check if division time is reached for differentiated cells
                if  t_div < t - cell[diff_time]:
                    cell[cell_state_idx] = dead_cell
                
                # check if cell dies
                if cell[death_time] < t - cell[diff_time]:
                    cell[cell_state_idx] = dead_cell
                    
            # update cells                        
            cells[j,i+1,:] = cell
                    
        # is there a birth? if not, cumulatively add p_birth
        if rnd_birth > p_birth:           
            p_birth = p_birth + (exp_cdf(t_new-t0, rate_birth)-
                                 exp_cdf(t-t0, rate_birth))
        
        # if there is a birth, draw new rnd number 
        # and make a new th0 cell out of a dead cell
        else:        
            p_birth = 0
            t0 = t
            rnd_birth = np.random.rand()            
            #search for a dead cell to transform             
            c = cells[:, i, cell_state_idx] == dead_cell
            assert c.any(), "cannot find any dead cells"
            # note that np.where returns a tuple
            cell_idx = np.where(cells[:, i, cell_state_idx] == dead_cell)[0][0]
            #print cell_idx
            cells[cell_idx, i+1, :] = make_cells(1,1, nstates, rate_death)
            cells[cell_idx, i+1, revival_time] = t
            
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
    cell_gens = cells[:,:,7]
    generations = [np.sum((cell_states == 1) & (cell_gens == i), axis = 0) for i in range(1, gen_no)]      
    return generations
# =============================================================================
# run stochastic simulation
# =============================================================================
simulation = [stoc_simulation(*stoc_params) for i in range(nsim)]

thn_th1_cells = [get_cells(*simu)[:2] for simu in simulation]

#==============================================================================
# run deterministic simulation
#==============================================================================
prolif_state = run_model(*prolif_params)
th1_prolif = prolif_state[:, alpha_1]

th1_n1 = interp1d(simulation_time, th1_prolif, axis = 0, kind='cubic', fill_value = 0, bounds_error = False)

# calculate subsequent generations based on interpolated gen1 cells
th1_gens = next_gens(simulation_time, rate_death, t_div, th1_n1, last_gen)

# get sum of all generations (note that this is only for those specified in last gen variable)
th1_all_gens = list(th1_gens)
th1_all_gens.insert(0, th1_prolif)
sum_th1 = [sum(x) for x in zip(*th1_all_gens)]
# =============================================================================
# plot stochastic simulation
# =============================================================================

# make dummy labels for pandas df
label = ["Thn","Th1"]
labels = [label for _ in range(nsim)]
flat_labels = [item for sublist in labels for item in sublist]

# make a df for single step stoc simulation (adjusted time)
flat_cells = [item for sublist in thn_th1_cells for item in sublist]
df = pd.DataFrame.from_records(flat_cells).transpose()
df.columns = flat_labels

df["time"] = simulation_time
df_long = df.melt(id_vars = ["time"])

fig, ax = plt.subplots(1, 1, figsize =(5,4))
stoc_plot = sns.lineplot(x = "time", 
                         y = "value", 
                         data = df_long,
                         hue = "variable",
                         ci = "sd",
                         ax = ax,
                         legend = False)
ax.set_ylabel("$n_{cells}$")
ax.set_xlim(0, simulation_time[-1])
ax.plot(simulation_time, th1_prolif, c = "tab:blue", linestyle = "--")
ax.legend(["Thn stoc", "Th1 stoc", "Th1 det"])
plt.tight_layout()

#==============================================================================
# plot generations in stoc model
#==============================================================================
#colors = ["tab:blue", "tab:orange"]

# plot deterministic generations
#for i, th1_gen in enumerate(th1_all_gens):
#    ax.plot(simulation_time, th1_gen, c = colors[i], linestyle = "--")
#ax.legend(label)
#ax.legend(["1st div stoc", "2nd div stoc", "1st div det", "2nd div det"])
#ax.set_ylim(0,3)
#plt.tight_layout()

