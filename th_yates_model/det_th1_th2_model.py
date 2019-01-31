#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Thu Jan 24 11:18:43 2019

@author: burt
equations for det th diff model with prolif adopted from yates et al
"""
import numpy as np
import matplotlib.pyplot as plt

def menten(x, hill, K):
    dummy = (x**hill / (x**hill + K))
    return dummy

def feedback(x, fb, hill, K):
    # positive fb strength should be greater or equal 1
    """ 
    use neg fb values for negative fb
    0 for no feedback (returns 1)
    and positive values for pos fb
    """

    dummy = 1 + fb * menten(x, hill, K)
    assert dummy >= 0
    
    return dummy

def cytokine_prob(cytokines, feedbacks, hill_coeffs, Ks):
        
    prob_cytokines = []
    
    for cyto, fb, hill, K in zip(cytokines, feedbacks, hill_coeffs, Ks):
        prob = feedback(cyto, fb, hill, K)
        prob_cytokines.append(prob)
    
    prob_cytokines = np.prod(np.asarray(prob_cytokines))
    
    return prob_cytokines


def th_cell_diff(state, 
                 t, 
                 alpha_1, 
                 alpha_2, 
                 beta_1, 
                 beta_2, 
                 alpha_prolif, 
                 beta_prolif, 
                 rate_birth, 
                 rate_death,
                 rate_cytos,
                 fb_th1,
                 fb_th2,
                 hill_th1,
                 hill_th2,
                 K_th1,
                 K_th2,
                 conc_il12,
                 ):

    # calculate cytokine concentrations
    th1 = state[:(alpha_1+alpha_prolif)]
    th2 = state[(alpha_1+alpha_prolif):]
    
    th1_cells = np.sum(th1)
    th2_cells = np.sum(th2)

    conc_ifn = rate_cytos[0] * th1_cells
    conc_il4 = rate_cytos[1] * th2_cells    
    cytokines = [conc_ifn, conc_il4, conc_il12]
    
    # branching probablities
    prob_th1 = cytokine_prob(cytokines, fb_th1, hill_th1, K_th1)
    prob_th2 = cytokine_prob(cytokines, fb_th2, hill_th2, K_th1)

    assert prob_th1 > 0, "prob th1="+str(prob_th1)+" prob th2="+str(prob_th2)+" cannot normalize."
    assert prob_th2 > 0
    
    # normalized branching probabilities
    p_1 = prob_th1 / (prob_th1 + prob_th2)
    p_2 = prob_th2 / (prob_th1 + prob_th2)

    cell_flux = [p_1 * rate_birth, p_2 * rate_birth]
         
    th1 = state[:(alpha_1+alpha_prolif)]
    th2 = state[(alpha_1+alpha_prolif):]
      
    dt_th1 = np.zeros_like(th1)
    dt_th2 = np.zeros_like(th2)
        
    th_states = [th1, th2]
    dt_th_states = [dt_th1, dt_th2]
    rate = [beta_1, beta_2]           
    alphas = [alpha_1, alpha_2]

    ### differential equations depending on the number on intermediary states
    # loop over states of th1 and th2 cells
    for i in range(2):
        th_state = th_states[i]
        dt_state = dt_th_states[i]
        r = rate[i]
        alpha = alphas[i]
        flux = cell_flux[i]
        #print flux
            
    # calculate derivatives
        for j in range(len(th_state)):
            if j == 0:
                dt_state[j] = flux - (r + rate_death) * th_state[j]
                
            elif j < alpha:
                dt_state[j] = r * (th_state[j-1] - th_state[j]) - rate_death * th_state[j]
                
            elif j == alpha:
                dt_state[j] = 2 * (r * th_state[j-1] + beta_prolif * th_state[-1]) - (rate_death + beta_prolif) * th_state[j]
            
            else:
                assert j > alpha
                dt_state[j] = beta_prolif * (th_state[j-1] - th_state[j]) - rate_death * th_state[j]

    # assign new number of naive cells based on the number of present naive cells that were designated th1_0 or th2_0
    #dt_th0 = state[0]    
    # return cell states
    dt_state = np.concatenate([dt_th1, dt_th2])
    
    return dt_state