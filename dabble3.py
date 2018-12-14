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
    cells = np.zeros((ncells, nsteps, nstates))
    for i in range(cells.shape[0]):        
        # random number for th0 to th_eff transition
        cells[i, :, 2] = np.random.rand()  
    return cells

def exp_cdf(t, mean):
    return 1 - np.exp(-t / mean)

def add_cell(cells):
    """
    insert a new naive cell into cell vector if there are dead cells available
    otherwise append a new naive cell to cell vector
    """
    ncells, nsteps, nstates = cells.shape

    # initialize new cell
    new_cell = make_cells(1, nsteps, nstates)    
    dead_cell_idx = 1
    # check if there are any dead cells
    c = cells[:,:,0] == dead_cell_idx   
    if c.any():
        # if yes, insert new cell into dead cell vector
        d = np.where(cells[:,:,0] == dead_cell_idx)[0][0]
        #print new_cell
        cells[d,:,:] = new_cell        
    # else append new cell
    else:
        cells = np.concatenate((cells, new_cell))
    
    return cells

def stoc_simulation(start, stop, nsteps, ncells, nstates):
    """
    run a simulation of the th1 th2 models   
    this function uses the a cumulative gamma fct integral to calculate transition probability
    """
    cells = make_cells(ncells, nsteps, nstates)
    time = np.linspace(start, stop, nsteps)
    
    rnd_birth = np.random.rand()
    p_birth = 0
    t0 = 0
    for i in range(len(time)-1):
        t = time[i]
        t_new = time[i+1]

        # I think if I add a new cell, it should not be totally new for all timesteps
        if rnd_birth > p_birth:           
            p_birth = p_birth + (exp_cdf(t_new-t0, 1.)-exp_cdf(t-t0, 1.))
        else:
            if t == 0:
                print "hi"
            p_birth = 0
            t0 = t
            rnd_birth = np.random.rand()
            cells = add_cell(cells)              
        # cell for each time step
        cell_j = cells[:,i,:]
                
        cell_number = cell_j.shape[0]
           
        for j in range(cell_number):
            cell = cell_j[j,:]                        
            #check if cell dies
            if cell[0] == 0:
                if cell[2] > cell[1]:
                    cell[1] = cell[1]+ (exp_cdf(t_new, cell_number)-exp_cdf(t, cell_number))
                else:
                    cell[0] = 1
              
                    
            cells[j,i+1,:] = cell
        
         
        
        #check if a new cell is added
            
    return [cells, time]     

def get_cells(cells, time):
    """
    input: cells and time, which is output from the function: run_stochastic_simulation
    use this for the model with precursor state times included
    """
    all_cells = cells[:,:, 0]
        
    # make some dummy lists
    th1_cells = []
    dead_cells = []
    
    # for each time step, check how many of the cells belong to each subtype
    for t in range(len(time)):
        x = all_cells[:,t]
        th1_cells.append(len(x[x==0]))
        dead_cells.append(len(x[x==1]))
   
    return th1_cells, dead_cells

#==============================================================================
# set up params and run
#==============================================================================
start = 0
stop = 5
nsteps = 2000
ncells = 1000
nstates = 3
nsim = 20

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

fig, ax = plt.subplots(1, 1, figsize =(5,4))
stoc_plot = sns.lineplot(x = "time", 
                         y = "value", 
                         data = df_long,
                         ci = "sd",
                         ax = ax,
                         legend = False)
ax.set_ylabel("$n_{cells}$")
ax.set_title(str(nsteps)+" steps, "+str(nsim)+ "simulations")
plt.tight_layout()
fig.savefig("stoc_simulation.pdf", bbox_inches = "tight")