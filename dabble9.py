#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Fri Dec 14 10:39:00 2018

@author: burt
same as dabble 6 but now I try to include prolif in stoc process
"""
import numpy as np
import matplotlib.pyplot as plt
from scipy.special import gammainc
import seaborn as sns
import pandas as pd
import one_celltype_model as det_model
from scipy.interpolate import interp1d

sns.set(context = "talk", style = "ticks")

#==============================================================================
# params
#==============================================================================
start = 0
stop = 3
nsteps = 12000
simulation_time = np.linspace(start, stop, nsteps)
last_gen = 4

# stoc params
nstates = 8
nsim = 100
n_dummies = 50
t_div_stoc = 0.2

# rates
alpha_1 = 20
alpha_2 = 1
beta_1 = float(alpha_1)
beta_2 = float(alpha_2)
rate_birth = 1.0
rate_death = 2.0

# other params
ncells = 1
t_div = t_div_stoc

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
    return scale_off * n_1(t - (i-1) * t_div)

def next_gens(simulation_time, degradation, t_div, first_gen_arr, no_of_next_gens):
    """
    calculate next generations based on array of first generation cells 
    and return as list of arrays
    """
    dummy_list = [N_i(simulation_time,
                     i, 
                     d = degradation, 
                     t_div = t_div, 
                     n_1 = first_gen_arr) for i in range(2, no_of_next_gens)]
    return dummy_list

def gamma_cdf(t, alpha, beta):
    dummy = t*beta
    return gammainc(alpha,dummy)      


def make_cells(ncells, nsteps, nstates):
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
        cells[i, :, 4] = np.random.rand() 
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
    prob_death_idx = 3
    rnd_death_idx = 4
    diff_time = 5
    revival_time = 6
    cell_gen = 7
    
    # define cell_state_idx state values
    th0_cell = 0
    th1_cell = 1
    dead_cell = 2

    # initialize some dummy dead cells t
    cells = make_cells(n_dummies+ncells, nsteps, nstates)
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
            
            # define a counter to index dead cells so that no dead cell is revived
            # for two different born cells
            # index 0 is reserved for birth process, higher indices for diff cells
            counter = 1
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
                    cell[cell_gen] = 1
                    #search for a dead cell to transform             
                    c = cells[:, i, cell_state_idx] == dead_cell
                    assert c.any(), "cannot find any dead cells"
                    # note that np.where returns a tuple
                    cell_idx = np.where(cells[:, i, cell_state_idx] == dead_cell)[0][counter]
                    #print cell_idx
                    cells[cell_idx, i+1, :] = make_cells(1,1, nstates)
                    cells[cell_idx, i+1, diff_time] = cell[diff_time]
                    cells[cell_idx, i+1, cell_state_idx] = cell[cell_state_idx]
                    cells[cell_idx, i+1, cell_gen] = cell[cell_gen]
                    # increase counter so that same dead cell is not used in this time step
                    counter += 1
                    
            # check if division time is reached for differentiated cells
            if cell[cell_state_idx] == th1_cell:
                # check if division time is reached
                assert cell[cell_gen] > 0, "if cell is th1, generation should be greater than 0"

                if t - cell[diff_time] > cell[cell_gen] * t_div:
                    # for current cell, increase generation number and remember div time
                    cell[cell_gen] += 1
                    cell[prob_death_idx] = 0
                    # check if there are dead cells at current time step
                    c = cells[:, i, cell_state_idx] == dead_cell
                    assert c.any(), "cannot find any dead cells"
                    # note that np.where returns a tuple
                    # get index of dead cell that is not already being used at this time step
                    cell_idx = np.where(cells[:, i, cell_state_idx] == dead_cell)[0][counter]
                    cells[cell_idx, i+1, :] = make_cells(1,1, nstates)
                    
                    # inherit diff_time, state_index and generation from mother cell
                    # note that cell generation is already updated in mother cell
                    cells[cell_idx, i+1, diff_time] = cell[diff_time]
                    cells[cell_idx, i+1, cell_state_idx] = cell[cell_state_idx]
                    cells[cell_idx, i+1, cell_gen] = cell[cell_gen]
                    counter += 1

                
            #check if cell dies
            if cell[cell_state_idx] == th1_cell:
                # do I need to assert that diff_time < t? otherwise newl
                if cell[rnd_death_idx] > cell[prob_death_idx]:
                    cell[prob_death_idx] = (cell[prob_death_idx]+
                        (exp_cdf(t_new-(cell[diff_time]+(cell[cell_gen]-1)*t_div), rate_death)-
                         exp_cdf(t-(cell[diff_time]+(cell[cell_gen]-1)*t_div), rate_death)))
                else:
                    cell[cell_state_idx] = dead_cell
                    
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
            cells[cell_idx, i+1, :] = make_cells(1,1, nstates)
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

def get_gens(cells, time):
    """
    input: cells and time, which is output from the function: run_stochastic_simulation
    use this for the model with precursor state times included
    """
    cell_states = cells[:,:, 0]
    cell_gens = cells[:,:,7]      
    first_gen = np.sum((cell_states == 1) & (cell_gens == 1), axis = 0)
    second_gen = np.sum((cell_states == 1) & (cell_gens == 2), axis = 0)
    third_gen = np.sum((cell_states == 1) & (cell_gens == 3), axis = 0)
   
    return first_gen, second_gen, third_gen

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
                         palette = ["k", "tab:blue"],
                         legend = False)
ax.set_ylabel("$n_{cells}$")
ax.set_xlim(0, simulation_time[-1])
ax.legend(["Thn","Th diff"])
ax.set_title(r"$\rightarrow$ Thn $\rightarrow$ Th diff $\rightarrow$ $\emptyset$")
plt.tight_layout()
#fig.savefig("stoc_simulation_new2.pdf", bbox_inches = "tight")

#==============================================================================
# compare stoc and deterministic model
#==============================================================================
#no_prolif_state = det_model.run_model(*no_prolif_params)
#th1_no_prolif = no_prolif_state[:, alpha_1]
#th2_no_prolif = no_prolif_state[:, -1] 

prolif_state = det_model.run_model(*prolif_params)
th1_prolif = prolif_state[:, alpha_1]

th1_n1 = interp1d(simulation_time, th1_prolif, axis = 0, kind='cubic', fill_value = 0, bounds_error = False)

# calculate subsequent generations based on interpolated gen1 cells
th1_gens = next_gens(simulation_time, rate_death, t_div, th1_n1, last_gen)

# get sum of all generations (note that this is only for those specified in last gen variable)
th1_all_gens = list(th1_gens)
th1_all_gens.insert(0, th1_prolif)
sum_th1 = [sum(x) for x in zip(*th1_all_gens)]

#==============================================================================
# plot stoc and det sim together
#==============================================================================
df_th1 = df_long[df_long["variable"] == "Th1"]
fig, ax = plt.subplots(1, 1, figsize =(5,4))
stoc_plot = sns.lineplot(x = "time", 
                         y = "value", 
                         data = df_th1,
                         ci = "sd",
                         ax = ax,
                         palette = ["tab:blue"],
                         legend = False)

#ax.plot(simulation_time, th1_no_prolif, c = "tab:red", label = "det sim", linestyle = "--")
ax.set_ylabel("$n_{cells}$")
ax.set_xlim(0, simulation_time[-1])
ax.legend(["stoc sim", "det sim"])
#ax.set_title(r"$\rightarrow$ Thn $\rightarrow$ Th diff $\rightarrow$ $\emptyset$")
plt.tight_layout()

#==============================================================================
# plot generations in stoc model
#==============================================================================
th_gens = [get_gens(*simu) for simu in simulation]
# make dummy labels for pandas df
label = ["1st gen","2nd gen", "3rd gen"]
labels = [label for _ in range(nsim)]
flat_labels = [item for sublist in labels for item in sublist]

# make a df for single step stoc simulation (adjusted time)
flat_cells = [item for sublist in th_gens for item in sublist]
df = pd.DataFrame.from_records(flat_cells).transpose()
df.columns = flat_labels

df["time"] = simulation_time
df_long = df.melt(id_vars = ["time"])

colors = ["k", "tab:blue", "tab:orange"]

fig, ax = plt.subplots(1, 1, figsize =(5,4))
stoc_plot = sns.lineplot(x = "time", 
                         y = "value", 
                         data = df_long,
                         hue = "variable",
                         ci = "sd",
                         ax = ax,
                         palette = colors,
                         legend = False)
ax.set_ylabel("$n_{cells}$")
ax.set_xlim(0, simulation_time[-1])

# plot deterministic generations
for th1_gen, col in zip(th1_all_gens, colors):
    ax.plot(simulation_time, th1_gen, c = col)
ax.legend(label)
plt.tight_layout()
