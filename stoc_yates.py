#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Fri Dec 14 10:39:00 2018

@author: burt
--> thn --> th1 --> death
also, th1 cells try to divide but this needs to be fixed
"""
import numpy as np
import matplotlib.pyplot as plt
from scipy.special import gammainc
import seaborn as sns
import pandas as pd
sns.set(context = "talk", style = "ticks")

#==============================================================================
# params
#==============================================================================
start = 0
stop = 3
nsteps = 6000
simulation_time = np.linspace(start, stop, nsteps)
last_gen = 3

# stoc params
nstates = 5
nsim = 50
n_dummies = 100

# diff rate
alpha_1 = 20
alpha_2 = 1
beta_1 = float(alpha_1)
beta_2 = float(alpha_2)

# prolif rate
alpha_prolif = 30
beta_prolif = float(alpha_prolif)

# other rates
rate_birth = 0
rate_death = 2.0


# other params
ncells = 1

stoc_params = [start, 
               stop, 
               nsteps, 
               ncells, 
               nstates, 
               rate_birth, 
               alpha_1, 
               beta_1, 
               alpha_prolif,
               beta_prolif,
               rate_death, 
               n_dummies,
               ]
#==============================================================================
# functions
#==============================================================================
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
        cells[i, :, 4] = np.random.exponential(1. / rate_death)
    return cells

def exp_cdf(t, rate):
    return 1 - np.exp(-t * rate)


def stoc_simulation(start, stop, nsteps, ncells, nstates, rate_birth, alpha_diff, beta_diff, alpha_prolif, beta_prolif, rate_death, n_dummies):
    """
    run a simulation of the th1 th2 models   
    this function uses the a cumulative gamma fct integral to calculate transition probability
    """
    # define state indices for the nstates array
    cell_state = 0
    prob_state_change = 1
    rnd_change = 2
    t_last_change = 3
    death_time = 4
    
    # define cell_state_idx state values
    th0_cell = 0
    th1_cell = 1
    dead_cell = 2

    # initialize some dummy dead cells t
    cells = make_cells(n_dummies+ncells, nsteps, nstates, rate_death)
    time = np.linspace(start, stop, nsteps)
    cells[:n_dummies,:,cell_state] = dead_cell

    
    # initialize random number and probability for birth process
    #rnd_birth = np.random.rand()
    #p_birth = 0
    #t0 = 0
    # loop over each time point
    for i in range(len(time)-1):
        t = time[i]
        t_new = time[i+1]

        # get all cells for current time step
        cell_j = cells[:,i,:]

        # cell number should be number of alive cells not total number              
        cell_number = cell_j.shape[0]
        
        counter = 1      
        # loop over each cell for current time point
        for j in range(cell_number):
            
            # define a counter to index dead cells so that no dead cell is revived
            # for two different born cells
            # index 0 is reserved for birth process, higher indices for diff cells
            
            cell = cell_j[j,:] 

            #check if cell differentiates
            if cell[cell_state] == th0_cell:
                if cell[rnd_change] > cell[prob_state_change]:
                    cell[prob_state_change] = (cell[prob_state_change]+
                        (gamma_cdf(t_new-cell[t_last_change], alpha_diff, beta_diff)-
                         gamma_cdf(t-cell[t_last_change], alpha_diff, beta_diff)))
                else:
                    cell[cell_state] = th1_cell
                    cell[t_last_change] = t
                    cell[prob_state_change] = 0
                    cell[rnd_change] = np.random.rand()
                    cell[death_time] = np.random.exponential(1. / rate_death)

                    #search for a dead cell to transform             
                    c = cells[:, i, cell_state] == dead_cell
                    assert c.any(), "cannot find any dead cells"
                    # note that np.where returns a tuple
                    idx = np.where(cells[:, i, cell_state] == dead_cell)[0][counter]
                    #print cell_idx
                    cells[idx, i+1, :] = make_cells(1,1, nstates, rate_death)
                    cells[idx, i+1, t_last_change] = cell[t_last_change]
                    cells[idx, i+1, cell_state] = cell[cell_state]

                    # increase counter so that same dead cell is not used in this time step
                    counter += 1
                    
                    
            # check if division time is reached for differentiated cells
            if cell[cell_state] == th1_cell:
                
                # check if division time is reached for th1 cells
                if cell[rnd_change] > cell[prob_state_change]:
                    #print cell[t_last_change], t
                    cell[prob_state_change] = (cell[prob_state_change]+
                        (gamma_cdf(t_new-cell[t_last_change], alpha_prolif, beta_prolif)-
                         gamma_cdf(t-cell[t_last_change], alpha_prolif, beta_prolif)))

                    
                else:
                    
                    cell[cell_state] = th1_cell
                    cell[t_last_change] = t
                    cell[prob_state_change] = 0
                    cell[rnd_change] = np.random.rand()
                    cell[death_time] = np.random.exponential(1. / rate_death)

                    #search for a dead cell to transform             
                    c = cells[:, i, cell_state] == dead_cell
                    assert c.any(), "cannot find any dead cells"
                    # note that np.where returns a tuple
                    idx = np.where(cells[:, i, cell_state] == dead_cell)[0][counter]
                    #print cell_idx
                    cells[idx, i+1, :] = make_cells(1,1, nstates, rate_death)
                    cells[idx, i+1, t_last_change] = cell[t_last_change]
                    cells[idx, i+1, cell_state] = cell[cell_state]

                    # increase counter so that same dead cell is not used in this time step
                    counter += 1
                               
                if cell[death_time] < t - cell[t_last_change]:
                    cell[cell_state] = dead_cell

            
            cells[j,i+1,:] = cell
        
            
        # is there a birth? if not, cumulatively add p_birth
        #if rnd_birth > p_birth:           
        #    p_birth = p_birth + (exp_cdf(t_new-t0, rate_birth)-
        #                         exp_cdf(t-t0, rate_birth))
        
        # if there is a birth, draw new rnd number 
        # and make a new th0 cell out of a dead cell
        #else:        
        #    p_birth = 0
        #    t0 = t
        #    rnd_birth = np.random.rand()            
        #    #search for a dead cell to transform             
        #    c = cells[:, i, cell_state_idx] == dead_cell
        #    assert c.any(), "cannot find any dead cells"
        #    # note that np.where returns a tuple
        #    cell_idx = np.where(cells[:, i, cell_state_idx] == dead_cell)[0][0]
        #    #print cell_idx
        #    cells[cell_idx, i+1, :] = make_cells(1,1, nstates)
        #    cells[cell_idx, i+1, revival_time] = t
            
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

# =============================================================================
# run stochastic simulation
# =============================================================================
#th1_cells = np.zeros_like(np.linspace(start, stop, nsteps))

simulation = [stoc_simulation(*stoc_params) for i in range(nsim)]

thn_th1_cells = [get_cells(*simu)[:2] for simu in simulation]

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
ax.legend(label)
#ax.set_title(r"$\rightarrow$ Thn $\rightarrow$ Th diff $\rightarrow$ $\emptyset$")
plt.tight_layout()

