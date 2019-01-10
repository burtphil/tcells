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
    for i in range(cells.shape[0]):        
        # random number for th0 to th_eff transition
        cells[i, :, 2] = np.random.rand()  
    return cells

def exp_cdf(t, rate):
    return 1 - np.exp(-t * rate)

def stoc_simulation(start, stop, nsteps, ncells, nstates):
    """
    run a simulation of the th1 th2 models   
    this function uses the a cumulative gamma fct integral to calculate transition probability
    """
    # initialize some dummy dead cells t
    n_dummies = 50
    
    cells = make_cells(n_dummies+ncells, nsteps, nstates)
    time = np.linspace(start, stop, nsteps)
    
    # make some dummy dead cells that can be used for birth process
    th0_idx = 0
    dead_cell_idx = 1
    cells[:n_dummies,:,0] = dead_cell_idx
    
    # initialize random number for birth process
    rnd_birth = np.random.rand()
    
    # initialize p_birth and birth time, will be increased over time
    p_birth = 0
    t0 = 0

    # rates
    rate_birth = 1.
    rate_death = 0.1        

    # loop over each time point
    for i in range(len(time)-1):
        t = time[i]
        t_new = time[i+1]

        # get all cells and corresponding states for current time step
        cell_j = cells[:,i,:]

        # get number of cells for current time step              
        cell_number = cell_j.shape[0]
        
        # loop over each cell for current time point
        for j in range(cell_number):
            cell = cell_j[j,:] 
                       
            #check if cell dies
            if cell[0] == th0_idx:
                if cell[2] > cell[1]:
                    cell[1] = cell[1]+ (exp_cdf(t_new-cell[3], rate_death)-
                        exp_cdf(t-cell[3], rate_death))
                else:
                    cell[0] = dead_cell_idx
                    
            cells[j,i+1,:] = cell
            
        # is there a birth? if not, cumulatively add p_birth
        if rnd_birth > p_birth:           
            p_birth = p_birth + (exp_cdf(t_new-t0, rate_birth)-exp_cdf(t-t0, rate_birth))
        
        # if there is a birth, draw new rnd number for the next birth time
        # and make a new th0 cell out of a dead cell
        else:          
            p_birth = 0
            t0 = t
            rnd_birth = np.random.rand()            
            #search for a dead cell to transform             
            c = cells[:,i,0] == dead_cell_idx
            assert c.any()
            # note that np.where returns a tuple
            cell_idx = np.where(cells[:,i,0] == dead_cell_idx)[0][0]
            #print cell_idx
            cells[cell_idx, i+1,:] = make_cells(1,1, nstates)
            
            # add survival time
            cells[cell_idx, i+1,3] = t

            
    return [cells, time]     


def get_cells(cells, time):
    """
    input: cells and time, which is output from the function: run_stochastic_simulation
    use this for the model with precursor state times included
    """
    all_cells = cells[:,:, 0]
          
    th1_cells = np.sum(all_cells == 0, axis = 0)
    dead_cells = np.sum(all_cells == 1, axis = 0)
   
    return th1_cells, dead_cells
#==============================================================================
# set up params
#==============================================================================
start = 0
stop = 3
nsteps = 3000
ncells = 10
nstates = 4
nsim = 100

# =============================================================================
# run simulation
# =============================================================================
th1_cells = np.zeros_like(np.linspace(start, stop, nsteps))

simulation = [stoc_simulation(start, stop, nsteps, ncells, nstates) for i in range(nsim)]
cell_numbers = [get_cells(*simu)[0] for simu in simulation]

time = np.linspace(start, stop, nsteps)

# make a df for single step stoc simulation (adjusted time)
flat_cells = [item for sublist in cell_numbers for item in sublist]
df = pd.DataFrame.from_records(cell_numbers).transpose()
labels = ["Th1" for _ in range(nsim)]
df.colums = ["Th1" for _ in range(nsim)]

df["time"] = time
df_long = df.melt(id_vars = ["time"])

# =============================================================================
# plot figure
# =============================================================================
fig, ax = plt.subplots(1, 1, figsize =(5,4))
stoc_plot = sns.lineplot(x = "time", 
                         y = "value", 
                         data = df_long,
                         ci = "sd",
                         ax = ax,
                         legend = False)
ax.set_ylabel("$n_{cells}$")
ax.set_title(str(nsteps)+" steps, "+str(nsim)+" simulations")
ax.set_xlim(0,time[-1])
plt.tight_layout()
#fig.savefig("stoc_simulation.pdf", bbox_inches = "tight")