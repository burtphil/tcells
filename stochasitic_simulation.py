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


def cytokine_prod(n_th, rate, base_rate = 1.0):
    """
    calculate cytokine concentration based on number of cells, the production rate and
    a basal cytokine concentration
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

def p_th_diff(conc_ifn,conc_il4,conc_il12, hill, half_saturation, strength = 1.0):
    """
    returns probability of Th1 Th2 differentiation for given cytokine concentrations
    kinetics are hill-like so hill coefficients for cytokines also need to be provided
    """
    assert conc_ifn >= 0, "ifn conc is "+str(conc_ifn)+", concentrations must be non-negative."
    assert conc_il4 >= 0, "il4 conc is "+str(conc_il4)+", concentrations must be non-negative."
    assert conc_il12 >= 0, "il12conc is "+str(conc_il12)+", concentrations must be non-negative."
    
    # watch out, this is tricky because if concentrations are zero and negative feedback (hill coeff < 0)
    # then you divide by zero, best avoid zero concentrations through base cytokine rate
    ifn_prob = (conc_ifn/half_saturation[0])**hill[0]
    il4_prob = (conc_il4/half_saturation[1])**hill[1]
    il12_prob = (conc_il12/half_saturation[2])**hill[2]

    prob_th_diff = strength*ifn_prob*il4_prob*il12_prob
    
    assert prob_th_diff >= 0, "probability is "+str(prob_th_diff)+" must be non-negative."
    
    return prob_th_diff

#==============================================================================
# simulations
#==============================================================================

def run_stochastic_simulation(start, stop, nsteps, ncells, nstates = 3):
    """
    run a simulation of the th1 th2 models   
    this function uses the a cumulative gamma fct integral to calculate transition probability
    """
    cells = np.zeros((ncells,nsteps,nstates))    
    time = np.linspace(start,stop,nsteps)
        
    for i in range(len(time)-1):
        t = time[i]
        cell_j = cells[:,i,:]
        
        n_th1 = count_cells(cell_j, cell_idx = cparams.th1_idx)
        n_th2 = count_cells(cell_j, cell_idx = cparams.th2_idx)
        
        ifn = cytokine_prod(n_th1, rate = cparams.rate_ifn, base_rate = cparams.kd_ifn)
        il4 = cytokine_prod(n_th2, rate = cparams.rate_il4, base_rate = cparams.kd_il4)
        
        p1 = p_th_diff(ifn,il4,cparams.conc_il12, hill = cparams.hill_1, half_saturation = cparams.half_saturation)
        p2 = p_th_diff(ifn,il4,cparams.conc_il12, hill = cparams.hill_2, half_saturation = cparams.half_saturation)
        
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
                if np.random.rand() > cell[2]:
                    cell[2] = cell[2]+gamma_cdf(t-cell[1],cparams.alpha_th1,cparams.beta_th1)
                else:
                    cell[0] = cparams.th1_idx
                    
            if cell[0] == cparams.th2_0_idx and t > cell[1]:
                if np.random.rand() > cell[2]:
                    cell[2] = cell[2]+gamma_cdf(t-cell[1],cparams.alpha_th2,cparams.beta_th2)
                else:
                    cell[0] = cparams.th2_idx
            
            cells[idx,i+1,:] = cell
            
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

