#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Thu Jan 24 11:18:43 2019

@author: burt
equations for det th diff model with prolif adopted from yates et al
"""
import numpy as np
from scipy.integrate import odeint

def feedback2(x, fb, hill, K):
    dummy = (fb * x + K) / (x + K)
    assert dummy >= 0
    return dummy

def cytokine_prob2(cytokines, feedbacks, hill_coeffs, Ks):
        
    prob_cytokines = []
    
    for cyto, fb, hill, K in zip(cytokines, feedbacks, hill_coeffs, Ks):
        prob = feedback2(cyto, fb, hill, K)
        prob_cytokines.append(prob)
    
    prob_cytokines = np.prod(np.asarray(prob_cytokines))
    
    return prob_cytokines

def th_cell_diff2(state, 
                 t, 
                 alpha_1, 
                 alpha_2, 
                 beta_1, 
                 beta_2, 
                 alpha_prolif, 
                 beta_prolif, 
                 rate_birth, 
                 rate_death,
                 rate_death_th0,
                 rate_diff_th0,
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
    th1 = state[1:(alpha_1+alpha_prolif+1)]
    th2 = state[(alpha_1+alpha_prolif+1):]
    
    th1_cells = th1[-alpha_prolif:]
    th1_cells = np.sum(th1_cells)
    
    th2_cells = th2[-alpha_prolif:]
    th2_cells = np.sum(th2_cells)

    conc_ifn = rate_cytos[0] * th1_cells
    conc_il4 = rate_cytos[1] * th2_cells    
    cytokines = [conc_ifn, conc_il4, conc_il12]
    
    # branching probablities
    prob_th1 = cytokine_prob2(cytokines, fb_th1, hill_th1, K_th1)
    prob_th2 = cytokine_prob2(cytokines, fb_th2, hill_th2, K_th1)

    assert prob_th1 > 0, "prob th1="+str(prob_th1)+" prob th2="+str(prob_th2)+" cannot normalize."
    assert prob_th2 > 0
    
    # normalized branching probabilities
    p_1 = prob_th1 / (prob_th1 + prob_th2)
    p_2 = prob_th2 / (prob_th1 + prob_th2)

    probs = [p_1, p_2]
               
    dt_th1 = np.zeros_like(th1)
    dt_th2 = np.zeros_like(th2)
        
    th_states = [th1, th2]
    dt_th_states = [dt_th1, dt_th2]
    rate = [beta_1, beta_2]           
    alphas = [alpha_1, alpha_2]

    ### differential equations depending on the number on intermediary states
    # loop over states of th1 and th2 cells
    dt_th0 = rate_birth - (rate_death_th0 + rate_diff_th0) * state[0]
    
    for i in range(2):
        th_state = th_states[i]
        dt_state = dt_th_states[i]
        r = rate[i]
        alpha = alphas[i]
        p = probs[i]
        #print flux
            
    # calculate derivatives
        for j in range(len(th_state)):
            if j == 0:
                dt_state[j] = p * rate_diff_th0 * state[0] - (r + rate_death) * th_state[j]
                
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
    dt_state = np.concatenate([[dt_th0], dt_th1, dt_th2])
    
    return dt_state

def run_model(
    cells_0,
    time,
    alpha_1, 
    alpha_2, 
    beta_1, 
    beta_2,
    alpha_prolif,
    beta_prolif,
    rate_birth,
    rate_death,
    rate_death_th0,
    rate_diff_th0,
    rate_cytos,
    fb_th1,
    fb_th2,
    hill_th1,
    hill_th2,
    K_th1,
    K_th2,
    conc_il12, 
    ):
    
    no_of_states = int(alpha_1 + alpha_2 + alpha_prolif + alpha_prolif + 1)
    y0 = np.zeros(no_of_states)    
    y0[0] = rate_birth / (rate_death_th0+rate_diff_th0)

    
    state = odeint(th_cell_diff2, 
                   y0, 
                   time, 
                   args = (alpha_1, 
                           alpha_2, 
                           beta_1, 
                           beta_2,
                           alpha_prolif,
                           beta_prolif,
                           rate_birth,
                           rate_death,
                           rate_death_th0,
                           rate_diff_th0,
                           rate_cytos,
                           fb_th1,
                           fb_th2,
                           hill_th1,
                           hill_th2,
                           K_th1,
                           K_th2,
                           conc_il12,
                           )
                   )
    return state

def get_cells(state, alpha_1, alpha_2, alpha_prolif):
    
    th1 = state[:, 1:(alpha_1+alpha_prolif+1)]
    th1_cells = th1[:, -alpha_prolif:]
    
    th2 = state[:, (alpha_1+alpha_prolif+1):]
    th2_cells = th2[:, -alpha_prolif:]
    
    th1_cells = np.sum(th1_cells, axis = 1)
    th2_cells = np.sum(th2_cells, axis = 1)
    
    return th1_cells, th2_cells

def find_nearest(array, value):
    """
    return index of element in array that is closest to value
    """
    array = np.asarray(array)
    idx = (np.abs(array - value)).argmin()
    return idx

def get_tau(time, cells):
    tau = time[find_nearest(cells, cells[-1] / 2.)]
    
    return tau

def vary_param(
    param_name,
    param_arr,
    cells_0,
    time,
    alpha_1, 
    alpha_2, 
    beta_1, 
    beta_2,
    alpha_prolif,
    beta_prolif,
    rate_birth,
    rate_death,
    rate_death_th0,
    rate_diff_th0,
    rate_cytos,
    fb_th1,
    fb_th2,
    hill_th1,
    hill_th2,
    K_th1,
    K_th2,
    conc_il12,               
        ):
    
    th1_endstates = []
    th1_tau = []
    th2_endstates = []
    th2_tau = [] 
    
    
    for i in param_arr:
        if param_name == "rate_death":
            rate_death = i
            
        elif param_name == "chain":
            assert isinstance(i, (int, np.integer))
            assert i > 0
            alpha_1 = i
            beta_1 = float(i)

        elif param_name == "prolif":
            assert isinstance(i, (int, np.integer))
            assert i > 0
            alpha_prolif = i
            beta_prolif = float(i)
            
        elif param_name == "rate_birth":
            rate_birth = i
            
        elif param_name == "fb_th1":
            fb_th1 = list(fb_th1)
            fb_th1[0] = i
        
        elif param_name == "rate_cytos":
            rate_cytos = [i, i]
            
        elif param_name == "rate_diff_th0":
            rate_diff_th0 = i
            
        elif param_name == "rate_death_th0":
            rate_death_th0 = i
        else:
            break
                
        state = run_model(cells_0, 
                          time, 
                          alpha_1, 
                          alpha_2, 
                          beta_1, 
                          beta_2,
                          alpha_prolif,
                          beta_prolif,
                          rate_birth,
                          rate_death,
                          rate_death_th0,
                          rate_diff_th0,
                          rate_cytos,
                          fb_th1,
                          fb_th2,
                          hill_th1,
                          hill_th2,
                          K_th1,
                          K_th2,
                          conc_il12,                                  
                          )
         
        th1_cells = get_cells(state, alpha_1, alpha_2, alpha_prolif)[0]
        th1_endstate = th1_cells[-1]
        th1_endstates.append(th1_endstate)

        tau_th1 = get_tau(time, th1_cells)
        th1_tau.append(tau_th1)
        
        th2_cells = get_cells(state, alpha_1, alpha_2, alpha_prolif)[1]
        th2_endstate = th2_cells[-1]
        th2_endstates.append(th2_endstate)
        
        tau_th2 = get_tau(time, th2_cells)
        th2_tau.append(tau_th2)
        
               
    return [th1_endstates, th2_endstates, th1_tau, th2_tau]