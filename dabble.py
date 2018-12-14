#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Tue Dec 11 16:04:49 2018

@author: burt
"""
import numpy as np
import matplotlib.pyplot as plt
from scipy.special import gammainc
#==============================================================================
# parameter settings
#==============================================================================

# simulation for one cell
### assign state indices 
# 0 represents Thn
# 1 represents Th1_0
# 2 represents Th2_0
# 3 represents Th1
# 4 represents Th2

#==============================================================================
# functions
#==============================================================================
def gamma_cdf(t, alpha, beta):
    dummy = t*beta
    return gammainc(alpha,dummy)      


def make_cells(ncells, nsteps, nstates):
    cells = np.zeros((ncells, nsteps, nstates))
    for i in range(cells.shape[0]):
        # rnd_nr for influx
        cells[i, :, 4] = np.random.exponential(1.)
        # random number for th0 to th_eff transition
        cells[i, :, 5] = np.random.rand()
        # exp rnd for death time of th_eff cell
        cells[i, :, 6] = np.random.exponential(1.)
    
    return cells

def add_cell(cells):
    """
    insert a new naive cell into cell vector if there are dead cells available
    otherwise append a new naive cell to cell vector
    """
    ncells, nsteps, nstates = cells.shape

    # initialize new cell
    new_cell = make_cells(1, nsteps, nstates)    
    dead_cell_idx = 2
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

#==============================================================================
# simulations
#==============================================================================
def semi_markov_simulation(start, stop, nsteps, ncells, alpha, beta, nstates):
    """
    run a simulation of the th1 th2 models   
    this function uses the a cumulative gamma fct integral to calculate transition probability
    """
    cells = make_cells(ncells, nsteps, nstates)
    time = np.linspace(start, stop, nsteps)
    
    for i in range(len(time)-1):
        t = time[i]
        t_new = time[i+1]
           
        # cell for each time step
        cell_j = cells[:,i,:]
        
        for j in range(cell_j.shape[0]):
            #print cell_j.shape[0]
            cell = cell_j[j,:]
                        
            # if cell is a naive cell, check if probability to transition is greater than rnd nr
            if cell[0] == 0:
                if cell[5] > cell[1]:
                    cell[1] = cell[1]+(gamma_cdf(t_new, alpha, beta)-gamma_cdf(t, alpha, beta))
                else:
                    cell[0] = 1
                    cell[2] = t
            
            if cell[0] == 1 and cell[6] < (t - cell[2]):
                cell[0] = 2       
                
                    
            cells[j,i+1,:] = cell
        
        # check if a new cell is added
            
    return [cells, time]     

def get_cells(cells, time):
    """
    input: cells and time, which is output from the function: run_stochastic_simulation
    use this for the model with precursor state times included
    """
    all_cells = cells[:,:, 0]
        
    # make some dummy lists
    naive_cells = []
    th1_cells = []
    dead_cells = []
    
    # for each time step, check how many of the cells belong to each subtype
    for t in range(len(time)):
        x = all_cells[:,t]
        naive_cells.append(len(x[x==0]))
        th1_cells.append(len(x[x==1]))
        dead_cells.append(len(x[x==2]))
    
    return [naive_cells,th1_cells, dead_cells]

start = 0
stop = 5
nsteps = 100
ncells = 100
alpha = 10
beta = 10.
nstates = 7

for _ in range(20):
    simulation = semi_markov_simulation(start, stop, nsteps, ncells, alpha, beta, nstates)
    cell_numbers = get_cells(*simulation)

    time = np.linspace(start, stop, nsteps)
    #plt.plot(time, cell_numbers[0])
    plt.plot(time, cell_numbers[1])
