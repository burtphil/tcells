#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Wed Oct 17 17:10:28 2018

@author: burt
"""

import numpy as np
import matplotlib.pyplot as plt
from scipy.special import gammainc
import os
os.chdir("/home/burt/Documents/code/th_cell_differentiation")
import modules.th1_th2_conceptual_parameters as cparams
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

def survival_fct(t, alpha, beta):
    return 1-gamma_cdf(t,alpha,beta)

def draw_fate(probs):     
    """
    sample a precursor fate
    """
    if np.random.rand() < probs[0]:
        return 1
    else:
        return 2

def draw_time(mean = 1.):
    """
    time of a poisson process
    """
    time = np.random.exponential(mean)
    return time


def cytokine_prod(n_th, rate, base_rate = 0):
    """
    calculate cytokine concentration based on number of cells, the production rate and
    a basal cytokine concentration
    note that if neg and pos feedback should be symmetric, base_rate should be set to 
    the half saturation konstant
    """
    return n_th*rate+base_rate

def count_cells(cell_t, cell_idx):
    """
    count thn, th1, th2, cells for a given time point
    """
    n_cells = len(cell_t[cell_t[:,0]==cell_idx])
    return n_cells

def p_norm(p1,p2):
    """
    return normalized probabilities
    """
    p1 = p1/(p1+p2)
    p2 = 1-p1
    return [p1,p2]

def p_gamma(conc_ifn, conc_il4, conc_il12, hill, K, fb_strength):

    assert conc_ifn >= 0, "ifn conc is "+str(conc_ifn)+", concentrations must be non-negative."
    assert conc_il4 >= 0, "il4 conc is "+str(conc_il4)+", concentrations must be non-negative."
    assert conc_il12 >= 0, "il12conc is "+str(conc_il12)+", concentrations must be non-negative."
    
    # watch out, this is tricky because if concentrations are zero and negative feedback (hill coeff < 0)
    # then you divide by zero, best avoid zero concentrations through base cytokine rate
    ifn_prob = (fb_strength[0] * (conc_ifn**hill[0]) + K[0]) / ((conc_ifn**hill[0]) + K[0])
    il4_prob = (fb_strength[1] * (conc_il4**hill[1]) + K[1]) / ((conc_il4**hill[1]) + K[1])
    il12_prob = (fb_strength[2] * (conc_il12**hill[2]) + K[2]) / ((conc_il12**hill[2]) + K[2])    
    
    prob_th_diff = ifn_prob*il4_prob*il12_prob

    assert prob_th_diff >= 0, "probability is "+str(prob_th_diff)+" must be non-negative."
    return prob_th_diff

#==============================================================================
# simulations
#==============================================================================

def run_stochastic_simulation(start, stop, nsteps, ncells, feedback_strength,
                              alpha_th1 = cparams.alpha_th1,
                              alpha_th2 = cparams.alpha_th2,
                              beta_th1 = cparams.beta_th1,
                              beta_th2 = cparams.beta_th2,
                              hill_1 = cparams.stoc_hill_1,
                              hill_2 = cparams.stoc_hill_2,
                              K = cparams.K,
                              nstates = 3,
                              ):
    """
    run a simulation of the th1 th2 models   
    this function uses the a cumulative gamma fct integral to calculate transition probability
    """
    
    # 3d cell vector that contains number of cells, number of time steps and number of states of each cell
    cells = np.zeros((ncells,nsteps,nstates))    
    time = np.linspace(start,stop,nsteps)
    
    # draw two random numbers for each cell, one for the time when cell fate is decided (from exp distribution)
    # and one, to compare with response time to decide time of state transition
    numbers_rnd_1 = [np.random.rand() for i in range(ncells)]
    fate_times = [np.random.exponential(cparams.precursor_rate) for i in range(ncells)]
    
    for i in range(len(time)-1):
        # for every time point loop over every cell
        t = time[i]
        t_new = time[i+1]
        cell_j = cells[:,i,:]
        
        # get numbers of differentiated cells for each time point
        n_th1 = count_cells(cell_j, cell_idx = cparams.th1_idx)
        n_th2 = count_cells(cell_j, cell_idx = cparams.th2_idx)
        
        # get cytokine concentrations based on differentiated cell numbers
        ifn = cytokine_prod(n_th1, rate = cparams.rate_ifn, base_rate = cparams.kd_ifn)
        il4 = cytokine_prod(n_th2, rate = cparams.rate_il4, base_rate = cparams.kd_il4)
        
        # calculate branching probabilities based on cytokine concentrations
        p_1 = p_gamma(ifn, il4, cparams.conc_il12, hill = hill_1, K = K, fb_strength = feedback_strength[0])
        p_2 = p_gamma(ifn, il4, cparams.conc_il12, hill = hill_2, K = K, fb_strength = feedback_strength[1])
        
        #p_1 = 0.7
        #p_2 = 0.3
        
        probs = p_norm(p_1,p_2)

        for idx, cell in enumerate(cell_j):
            
            # get random numbers calculated above for each cell
            n_rnd_1 = numbers_rnd_1[idx]
            fate_time = fate_times[idx]
            
            # check if there is a fate decision, by checking if cell is a naive cell 
            # and check if time value is greater than random nr drawn from exp distr.
            if cell[0] == cparams.thn_idx and fate_time < t:
                # calculate probabilities which fate to choose and
                # assign new cell fate
                cell[0] = draw_fate(probs)
                cell[1] = t
            
            # for th1 precursor cells, if time value is greater than the
            # time when cell fate was decided (cell[1]), check if fate change occurs
            # by comparing with response time distribution
            if cell[0] == cparams.th1_0_idx and t > cell[1]:
                if n_rnd_1 > cell[2]:
                    cell[2] = cell[2]+(gamma_cdf(t_new-cell[1],alpha_th1,beta_th1)-gamma_cdf(t-cell[1],alpha_th1,beta_th1))
                else:
                    cell[0] = cparams.th1_idx
                    
            # for th2 precursor cells, if time value is greater than the
            # time when cell fate was decided (cell[1]), check if fate change occurs
            # by comparing with response time distribution                    
            if cell[0] == cparams.th2_0_idx and t > cell[1]:
                if n_rnd_1 > cell[2]:
                    cell[2] = cell[2]+(gamma_cdf(t_new-cell[1],alpha_th2,beta_th2)-gamma_cdf(t-cell[1],alpha_th2,beta_th2))
                else:
                    cell[0] = cparams.th2_idx
            
            # update each cell for the next time step
            cells[idx,i+1,:] = cell
            
    return [cells, time]    

def semi_markov_simulation(start, stop, nsteps, ncells, alpha, beta, nstates = 2):
    """
    run a simulation of the th1 th2 models   
    this function uses the a cumulative gamma fct integral to calculate transition probability
    """
    cells = np.zeros((ncells,nsteps,nstates))
    time = np.linspace(start,stop,nsteps)
    numbers_rnd = [np.random.rand() for i in range(ncells)]
    
    for i in range(len(time)-1):
        t = time[i]
        #print t
        cell_j = cells[:,i,:]
        #print cell_j.shape
        for j in range(ncells):
            cell = cell_j[j,:]
            n_rnd = numbers_rnd[j]
            # is there a cell fate switch?
            #print cell.shape
            #print cell                
            if cell[0] == 0:
                p_survive = survival_fct(t,alpha,beta)
                
                if n_rnd > p_survive:
                    cell[0] = 1
                    #print rnd_no, p_survive, t
                    #print t
                    
            cells[j,i+1,:] = cell
            
    return [cells, time]               

def semi_markov_simulation2(start, stop, nsteps, ncells, alpha, beta, nstates = 2):
    """
    run a simulation of the th1 th2 models   
    this function uses the a cumulative gamma fct integral to calculate transition probability
    """
    cells = np.zeros((ncells,nsteps,nstates))    
    time = np.linspace(start,stop,nsteps)
    numbers_rnd = [np.random.rand() for i in range(ncells)]
        
    for i in range(len(time)-1):
        t = time[i]
        t_new = time[i+1]
        #print t
        cell_j = cells[:,i,:]
        #print cell_j.shape
        for j in range(ncells):
            cell = cell_j[j,:]
            n_rnd = numbers_rnd[j]
            # is there a cell fate switch?
            #print cell.shape
            #print cell                
            if cell[0] == 0:
                if n_rnd > cell[1]:
                    cell[1] = cell[1]+(gamma_cdf(t_new, alpha, beta)-gamma_cdf(t, alpha, beta))
                else:
                    cell[0] = 1
                    
            cells[j,i+1,:] = cell
            
    return [cells, time]     

    
def run_simulation2(start, stop, nsteps, ncells, nstates = 3):
    """
    this function uses survival function to calculate transition probability 
    """
    cells = np.zeros((ncells,nsteps,nstates))    
    time = np.linspace(start,stop,nsteps)
        
    for i in range(len(time)-1):
        t = time[i]
        cell_j = cells[:,i,:]
        
        n_th1 = count_cells(cell_j, cell_idx = cparams.th1_idx)
        n_th2 = count_cells(cell_j, cell_idx = cparams.th2_idx)
        
        ifn = cytokine_prod(n_th1, rate = cparams.rate_ifn)
        il4 = cytokine_prod(n_th2, rate = cparams.rate_il4)
        
        p1 = p_th_diff(ifn,il4,1., hill = cparams.hill_1, half_saturation = cparams.half_saturation)
        p2 = p_th_diff(ifn,il4,1., hill = cparams.hill_2, half_saturation = cparams.half_saturation)
        
        probs = p_norm(p1,p2)
        #print cells.shape
        #print cell_j.shape
        for idx, cell in enumerate(cell_j):
            # is there a cell fate switch?
            #print cell.shape
            if cell[0] == cparams.thn_idx and np.random.exponential(cparams.precursor_rate) < t:
                # calculate probabilities which fate to choose
                # assign new cell fate
                cell[0] = draw_fate(probs)
                cell[1] = t
                
            if cell[0] == cparams.th1_0_idx and t > cell[1]:
                if np.random.rand() > survival_fct(t-cell[1],1,1):
                    cell[0] = cparams.th1_idx
                    
            if cell[0] == cparams.th2_0_idx and t > cell[1]:
                if np.random.rand() > survival_fct(t-cell[1],1,1):
                    cell[0] = cparams.th2_idx
            
            cells[idx,i+1,:] = cell
            
    return [cells, time]          
    #plt.plot(time, cell_j[:,0])

