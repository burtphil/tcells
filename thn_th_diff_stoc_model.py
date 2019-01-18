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
from rising_star_model import run_model
from scipy.interpolate import interp1d

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
        # random number for death time
        cells[i, :, 4] = np.random.rand() 
    return cells

def exp_cdf(t, rate):
    return 1 - np.exp(-t * rate)

def stoc_simulation(start,
                    stop,
                    nsteps,
                    ncells,
                    nstates,
                    rate_birth,
                    alpha_diff,
                    beta_diff,
                    rate_death,
                    n_dummies,
                    t_div):
    """
    run a simulation of the th1 th2 models   
    this function uses the a cumulative gamma fct integral to calculate transition probability
    """
    # each cell contains an 8 dimensional array that is accesible through these indices
    cell_state_idx = 0
    prob_diff_idx = 1
    rnd_diff_idx = 2
    prob_death_idx = 3
    rnd_death_idx = 4
    diff_time = 5
    revival_time = 6
    cell_gen = 7
    
    # cell_state_idx can take following values
    th0_cell = 0
    th1_cell = 1
    dead_cell = 2
    
    # make th naive cells, each cell is initialized with a random number to check
    # when it differentiates and a random number to check when it dies
    # these numbers are compared to dynamically updated probabilities

    cells = make_cells(n_dummies+ncells, nsteps, nstates)
    time = np.linspace(start, stop, nsteps)
    # initialize some dummy dead cells to revive later
    cells[:n_dummies,:,cell_state_idx] = dead_cell
   
    # initialize random number and probability for birth process
    rnd_birth = np.random.rand()
    p_birth = 0
    t0 = 0
    
    ######### Here come the actual iterations #################################
    # loop over each time point
    for i in range(len(time)-1):
        t = time[i]
        t_new = time[i+1]

        # get all cells and their number for current time step
        cell_j = cells[:,i,:]
        cell_number = cell_j.shape[0]
              
        # loop over each cell for current time point
        for j in range(cell_number):
            
            # define a counter to index dead cells so that no dead cell is revived
            # for two different cells at the same time
            # index 0 is reserved for birth process, higher indices for diff cells
            counter = 1
            cell = cell_j[j,:] 

            #check if cell differentiates (and thereby undergoes first division)
            if cell[cell_state_idx] == th0_cell:
                
                # if probability to differentiate is too small, increase probability
                if cell[rnd_diff_idx] > cell[prob_diff_idx]:
                    cell[prob_diff_idx] = (cell[prob_diff_idx]+
                        (gamma_cdf(t_new-cell[revival_time], alpha_diff, beta_diff)-
                         gamma_cdf(t-cell[revival_time], alpha_diff, beta_diff)))
                else:
                    # update the state index and remember the time and generation
                    cell[cell_state_idx] = th1_cell
                    cell[diff_time] = t
                    cell[cell_gen] = 1
                    
                    #search for a dead cell to transform             
                    c = cells[:, i, cell_state_idx] == dead_cell
                    assert c.any(), "cannot find any dead cells"
                    # note that np.where returns a tuple
                    cell_idx = np.where(cells[:, i, cell_state_idx] == dead_cell)[0][counter]
                    
                    # make a daughter th1 cell out of a dead cell
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
                    # and set death probability back to 0
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


#==============================================================================
# params
#==============================================================================
start = 0
stop = 3
nsteps = 12000
simulation_time = np.linspace(start, stop, nsteps)
last_gen = 3

# stoc params
nstates = 8
nsim = 50
n_dummies = 20
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
