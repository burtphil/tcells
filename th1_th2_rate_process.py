#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Tue Sep 18 12:13:22 2018

@author: burt
"""
import numpy as np
import th1_th2_conceptual_parameters as cparams
from scipy.integrate import odeint
import matplotlib.pyplot as plt

def p_th_diff(conc_ifn,conc_il4,conc_il12, hill, half_saturation):
    """
    returns probability of Th1 Th2 differentiation for given cytokine concentrations
    kinetics are hill-like so hill coefficients for cytokines also need to be provided
    """
    prob_th_diff = (((conc_ifn/half_saturation[0])**hill[0])*
                    ((conc_il4/half_saturation[1])**hill[1])*
                    ((conc_il12/half_saturation[2])**hill[2]))   
    return prob_th_diff

def th_cell_diff(state, t,rate1,rate2,conc_il12, hill_1, hill_2, 
                 rate_ifn, rate_il4, half_saturation):
    
    
    # calculate interferon gamma (ifn) and il4 concentrations based on the number of th1 and th2 cells
    conc_ifn = rate_ifn*state[1]
    conc_il4 = rate_il4*state[2]

    ### calculate initial th1 and th2 populations from naive cells based on branching probabilities
    # naive cells
    th_0 = state[0]
     
    # branching probablities
    prob_th1 = p_th_diff(conc_ifn,conc_il4,conc_il12, hill_1, half_saturation)
    prob_th2 = p_th_diff(conc_ifn,conc_il4,conc_il12, hill_2, half_saturation)
    
    print prob_th1, prob_th2    
    # normalized branching probabilities
    prob_th1_norm = prob_th1 / (prob_th1+prob_th2)
    prob_th2_norm = prob_th2 / (prob_th1+prob_th2)
    th1_0 = prob_th1_norm*th_0
    th2_0 = prob_th2_norm*th_0
       
    # assign th1 states and th2 states from initial vector based on chain length alpha
    th1 = [th1_0, state[1]]
    th2 = [th2_0, state[2]]
      
    dt_th1_0 = -rate1*th1[0]
    dt_th2_0 = -rate2*th2[0]
    dt_th1_1 = rate1*th1[0]
    dt_th2_1 = rate2*th2[0]
    
    dt_th0 = dt_th1_0+dt_th2_0
      
    # return cell states
    dt_state = np.concatenate(([dt_th0],[dt_th1_1],[dt_th2_1]))
        
    return dt_state

initial_cells = cparams.initial_cells
state_0 = [initial_cells,1,1]
t = np.arange(0,6,0.01)

hill_1 = cparams.hill_1
hill_2 = cparams.hill_2
conc_il12 = cparams.conc_il12
rate1 = 1.2
rate2 = 1.
rate_ifn = cparams.rate_ifn
rate_il4 = cparams.rate_il4
half_saturation = cparams.half_saturation

state = odeint(th_cell_diff, state_0, t, args = (rate1, rate2, conc_il12, hill_1, hill_2,
                                                 rate_ifn, rate_il4, half_saturation))
norm = initial_cells/100
    
th0_cells = state[:,0]/norm
th1_cells = state[:,1]/norm
th2_cells = state[:,2]/norm
    
fig, ax = plt.subplots(1,1, figsize = (5,5))
ax.plot(t,th0_cells)
ax.plot(t,th1_cells)
ax.plot(t,th2_cells)
ax.set_yticks([0,50,100])

