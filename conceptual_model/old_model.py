#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Tue Oct  2 17:51:52 2018

@author: burt
"""

def th_cell_diff(state,t,alpha_1,alpha_2,rate1,rate2,conc_il12, hill_1, hill_2, 
                 rate_ifn, rate_il4, half_saturation, base_production_rate_ifn, base_production_rate_il4, degradation):
        
    # calculate interferon gamma (ifn) and il4 concentrations based on the number of th1 and th2 cells
    assert type(alpha_1) == int, "alpha dtype is "+str(type(alpha_1))+" but alpha must be an integer."
    assert type(alpha_2) == int, "alpha dtype is "+str(type(alpha_2))+" but alpha must be an integer."
        
    base_cytokine_rate = 0.00001
    conc_ifn = rate_ifn*state[alpha_1]+base_cytokine_rate
    conc_il4 = rate_il4*state[-1]+base_cytokine_rate

    ### calculate initial th1 and th2 populations from naive cells based on branching probabilities
    th_0 = state[0]
    assert th_0 > 0, "no initial cells provided"
    
    # branching probablities
    prob_th1 = p_th_diff(conc_ifn,conc_il4,conc_il12, hill_1, half_saturation, base_production_rate_ifn)
    prob_th2 = p_th_diff(conc_ifn,conc_il4,conc_il12, hill_2, half_saturation, base_production_rate_il4)    
    assert prob_th1+prob_th2 > 0, "prob th1="+str(prob_th1)+" prob th2="+str(prob_th2)+" cannot normalize."

    # normalized branching probabilities
    prob_th1_norm = prob_th1/(prob_th1+prob_th2)
    prob_th2_norm = prob_th2/(prob_th1+prob_th2)
    th1_0 = prob_th1_norm*th_0
    th2_0 = prob_th2_norm*th_0


    # assign th1 states and th2 states from initial vector based on chain length alpha
    th1 = np.concatenate(([th1_0],state[1:(alpha_1+1)]))
    th2 = np.concatenate(([th2_0],state[(alpha_1+1):]))      
    dt_th1 = np.ones_like(th1)
    dt_th2 = np.ones_like(th2)        
    th_states = [th1,th2]
    dt_th_states = [dt_th1,dt_th2]
    rate = [rate1,rate2]
        
    ### differential equations depending on the number on intermediary states
    # loop over states of th1 and th2 cells
    for i in range(2):
        th_state = th_states[i]
        dt_state = dt_th_states[i]
        r = rate[i]
            
    # calculate derivatives
        for j in range(len(th_state)):
            if j == 0:
                dt_state[j] = -r*th_state[j]
            elif j != (len(th_state)-1):
                dt_state[j] = r*(th_state[j-1]-th_state[j])
            else:
                dt_state[j] = r*th_state[j-1] - degradation

    # assign new number of naive cells based on the number of present naive cells that were designated th1_0 or th2_0    
    dt_th0 = dt_th1[0]+dt_th2[0]
        
    # return cell states
    dt_state = np.concatenate(([dt_th0],dt_th1[1:],dt_th2[1:]))
    
    assert np.isnan(dt_state).any() == False, "nans detected in dt_state array."
    
    return dt_state  

def run_model(title, parameters, model_name = th_cell_diff, save = True):
    """ 
    run ode model with parameters from params file
    needs shape and rate params as input as well as simulation time
    """
    (alpha_1, alpha_2, rate1, rate2, simulation_time, conc_il12, hill_1, hill_2,
     rate_ifn, rate_il4, half_saturation, base_production_rate_ifn, 
     base_production_rate_il4, initial_cells, degradation) = parameters

    # define initial conditions based on number of intermediary states
    no_of_states = int(alpha_1+alpha_2+1)
    ini_cond = np.zeros(no_of_states)
    ini_cond[0] = initial_cells
    
    state = odeint(model_name, ini_cond, simulation_time, 
                   args =(alpha_1,alpha_2,rate1,rate2,conc_il12, hill_1,hill_2,
                          rate_ifn,rate_il4,half_saturation, base_production_rate_ifn,
                          base_production_rate_il4, degradation))
    
    # save both model simulation and associated parameters
    if save == True:
        parameters = np.asarray(parameters, dtype = object)
        np.savez(title, state = state, parameters = parameters)

    return state