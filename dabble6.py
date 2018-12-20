#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Fri Dec 14 10:39:00 2018

@author: burt
"""
#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Tue Dec 11 16:04:49 2018

@author: burt
"""
import numpy as np
import matplotlib.pyplot as plt
from scipy.special import gammainc
import seaborn as sns
import pandas as pd

sns.set(context = "talk", style = "ticks")
#==============================================================================
# functions
#==============================================================================
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

def stoc_simulation(start, stop, nsteps, ncells, nstates):
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
    
    # define cell_state_idx state values
    th0_cell = 0
    th1_cell = 1
    dead_cell = 2

    # initialize some dummy dead cells t
    n_dummies = 100
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
        
        # get number of current th0 and th1 cells
        #n_th0 = sum(cell_j[:, cell_state_idx] == th0_cell)
        #n_th1 = sum(cell_j[:, cell_state_idx] == th1_cell)
        
        # define rates
        rate_birth = 20.
        rate_diff = 1.
        rate_death = 1.
        # loop over each cell for current time point
        for j in range(cell_number):
            cell = cell_j[j,:] 

            #check if cell differentiates
            if cell[cell_state_idx] == th0_cell:
                if cell[rnd_diff_idx] > cell[prob_diff_idx]:
                    cell[prob_diff_idx] = (cell[prob_diff_idx]+
                        (gamma_cdf(t_new-cell[revival_time], 20, 20)-
                         gamma_cdf(t-cell[revival_time], 20, 20)))
                else:
                    cell[cell_state_idx] = th1_cell
                    cell[diff_time] = t
                       
            #check if cell dies
            if cell[cell_state_idx] == th1_cell:
                # do I need to assert that diff_time < t? otherwise newl
                if cell[rnd_death_idx] > cell[prob_death_idx]:
                    cell[prob_death_idx] = (cell[prob_death_idx]+
                        (exp_cdf(t_new-cell[diff_time], rate_death)-
                         exp_cdf(t-cell[diff_time], rate_death)))
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
        
    # make some dummy lists
    th0_cells = []
    th1_cells = []
    dead_cells = []
    
    # for each time step, check how many of the cells belong to each subtype
    for t in range(len(time)):
        x = all_cells[:,t]
        th0_cells.append(len(x[x==0]))
        th1_cells.append(len(x[x==1]))
        dead_cells.append(len(x[x==2]))
   
    return th0_cells, th1_cells, dead_cells

#==============================================================================
# set up params
#==============================================================================
start = 0
stop = 3
nsteps = 30000
ncells = 20
nstates = 7
nsim = 50
# =============================================================================
# run simulation
# =============================================================================
th1_cells = np.zeros_like(np.linspace(start, stop, nsteps))

simulation = [stoc_simulation(start, stop, nsteps, ncells, nstates) for i in range(nsim)]
time = np.linspace(start, stop, nsteps)
thn_th1_cells = [get_cells(*simu)[:2] for simu in simulation]

# =============================================================================
# plot figure
# =============================================================================

# make dummy labels for pandas df
label = ["Thn","Th1"]
labels = [label for _ in range(nsim)]
flat_labels = [item for sublist in labels for item in sublist]

# make a df for single step stoc simulation (adjusted time)
flat_cells = [item for sublist in thn_th1_cells for item in sublist]
df = pd.DataFrame.from_records(flat_cells).transpose()
df.columns = flat_labels

df["time"] = time
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
ax.set_xlim(0,time[-1])
plt.tight_layout()
fig.savefig("stoc_simulation_new2.pdf", bbox_inches = "tight")
