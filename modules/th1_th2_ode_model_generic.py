#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Thu Sep  6 13:11:32 2018

@author: burt
"""
import numpy as np#
from scipy.integrate import odeint

def find_nearest(array, value):
    """
    return index of element in array that is closest to value
    """
    array = np.asarray(array)
    idx = (np.abs(array - value)).argmin()
    return idx

def menten(x, hill, K):
    dummy = (x**hill / (x**hill + K))
    return dummy

def neg_fb(x, fb, hill, K):
    assert hill > 0
    assert 0 <= fb <= 1.
    # negative feedback strength should be between 0 and 1
    return 1 - fb * menten(x, hill, K)

def pos_fb(x, fb, hill, K):
    assert hill > 0
    assert fb >= 0 , "fb strength is "+str(fb)
    # positive fb strength should be greater or equal 1
    return 1 + fb * menten(x, hill, K)

def p_menten(conc_ifn, conc_il4, conc_il12, hill, K, fb_strength):
    """
    returns probability of Th1 Th2 differentiation for given cytokine concentrations
    kinetics are hill-like so hill coefficients for cytokines also need to be provided
    """
    assert conc_ifn >= 0, "ifn conc is "+str(conc_ifn)+", concentrations must be non-negative."
    assert conc_il4 >= 0, "il4 conc is "+str(conc_il4)+", concentrations must be non-negative."
    assert conc_il12 >= 0, "il12conc is "+str(conc_il12)+", concentrations must be non-negative."

    # watch out, this is tricky because if concentrations are zero and negative feedback (hill coeff < 0)
    # then you divide by zero, best avoid zero concentrations through base cytokine rate
    
    cytokine = np.array([conc_ifn, conc_il4, conc_il12])
    prob_cytokine = [1., 1., 1.]
    for i in range(3):
        if fb_strength[i] < 0:            
            prob_cytokine[i] = neg_fb(cytokine[i], hill = hill[i], fb = -fb_strength[i], K = K[i])
            #print prob_cytokine[i]
            
        if fb_strength[i] > 0:    
            prob_cytokine[i] = pos_fb(cytokine[i], hill = hill[i], fb = fb_strength[i], K = K[i])
            
            #print prob_cytokine[i]
            
    #print prob_cytokine   
    prob_th_diff = np.prod(prob_cytokine)

    assert prob_th_diff >= 0, "probability is "+str(prob_th_diff)+"but must be non-negative."
    
    return prob_th_diff

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
   
    
def get_cyto_conc(state, simulation_time, alpha_1, rate_ifn, rate_il4, initial_cells):
    """
    get cytokine concentrations for a performed model simulation
    use to derive probabilities by giving this array to fun get_prob
    """
    th1_cells = state[:,alpha_1+1]
    th2_cells = state[:,-1]
    conc_ifn = rate_ifn*th1_cells
    conc_il4 = rate_il4*th2_cells
    return conc_ifn, conc_il4

def get_prob(conc_ifn_list, conc_il4_list, conc_il12, hill_1, hill_2, K, fb_strength, fun_prob = p_gamma):
    """
    calculate branching probs for a list of cytokine concentrations
    return probabilities
    """
    prob_th1 = [fun_prob(conc_ifn, conc_il4, conc_il12, hill_1, K, fb_strength[0]) for conc_ifn, conc_il4 in zip(conc_ifn_list, conc_il4_list)]
    prob_th2 = [fun_prob(conc_ifn, conc_il4, conc_il12, hill_2, K, fb_strength[1]) for conc_ifn, conc_il4 in zip(conc_ifn_list, conc_il4_list)]
    prob_th1 = np.array(prob_th1)
    prob_th2 = np.array(prob_th2)
    p_1 = prob_th1 / (prob_th1 + prob_th2)
    p_2 = prob_th2 / (prob_th1 + prob_th2)
    return p_1, p_2
    
def th_cell_diff(state,
                 t,
                 alpha_1,
                 alpha_2,
                 beta_1,
                 beta_2,
                 conc_il12,
                 hill_1,
                 hill_2, 
                 fb_strength,
                 rate_ifn,
                 rate_il4,
                 K,
                 degradation,
                 fb_start,
                 fb_end,
                 const_thn = False,
                 fun_probability = p_gamma,
                 ):
        
    # calculate interferon gamma (ifn) and il4 concentrations based on the number of th1 and th2 cells
    assert isinstance(alpha_1, (int, np.integer)), "alpha dtype is "+str(type(alpha_1))+" but alpha must be an integer."
    assert isinstance(alpha_2, (int, np.integer)), "alpha dtype is "+str(type(alpha_2))+" but alpha must be an integer."
        
    conc_ifn = rate_ifn*state[alpha_1+1]
    conc_il4 = rate_il4*state[-1]
        #print conc_il4, conc_ifn

    ### calculate initial th1 and th2 populations from naive cells based on branching probabilities    
    th_0 = state[0]    
    if t < 1:   
        assert th_0 > 0, "no initial cells provided"
    
    # branching probablities
    if fb_start <= t <= fb_end:
        prob_th1 = fun_probability(conc_ifn, conc_il4, conc_il12, hill_1, K, fb_strength[0])
        prob_th2 = fun_probability(conc_ifn, conc_il4, conc_il12, hill_2, K, fb_strength[1])    
        assert prob_th1+prob_th2 > 0, "prob th1="+str(prob_th1)+" prob th2="+str(prob_th2)+" cannot normalize."
    
        # normalized branching probabilities
        p_1 = prob_th1 / (prob_th1 + prob_th2)
        p_2 = prob_th2 / (prob_th1 + prob_th2)
    
    else:
        p_1 = 0.5
        p_2 = 0.5
        
    # assign th1 states and th2 states from initial vector based on chain length alpha
    th1 = state[1:(alpha_1+2)]
    th2 = state[(alpha_1+2):]
      
    dt_th1 = np.ones_like(th1)
    dt_th2 = np.ones_like(th2)
        
    th_states = [th1,th2]
    dt_th_states = [dt_th1,dt_th2]
    rate = [beta_1, beta_2]
    prob = [p_1, p_2]
        
    ### differential equations depending on the number on intermediary states
    # loop over states of th1 and th2 cells
    for i in range(2):
        th_state = th_states[i]
        dt_state = dt_th_states[i]
        r = rate[i]
        p = prob[i]
            
    # calculate derivatives
        for j in range(len(th_state)):
            if j == 0:
                dt_state[j] = p * th_0 - r * th_state[j]
                
            elif j != (len(th_state)-1):
                dt_state[j] = r * (th_state[j-1] - th_state[j])
                
            else:
                dt_state[j] = r * th_state[j-1] - degradation * th_state[j]

    # assign new number of naive cells based on the number of present naive cells that were designated th1_0 or th2_0
    # if a constant Th naive cell pool is assumed (default parameter const_thn = True) then change should be 0
    # because pool should not change   
    th0_influx = 1
    
    if const_thn == False:
        dt_th0 = -th_0
    else:
        dt_th0 = th0_influx - th_0
   # return cell states
    dt_state = np.concatenate(([dt_th0], dt_th1, dt_th2))
    
    assert np.isnan(dt_state).any() == False, "nans detected in dt_state array."
    
    return dt_state  

def run_model(alpha_1,
              alpha_2,
              beta_1,
              beta_2,
              simulation_time,
              conc_il12,
              hill_1,
              hill_2,
              fb_strength,
              rate_ifn,
              rate_il4,
              K,
              initial_cells,
              degradation,
              fb_start,
              fb_end
              ):
    """ 
    run ode model with parameters from params file
    needs shape and rate params as input as well as simulation time
    """
    # define initial conditions based on number of intermediary states
    no_of_states = int(alpha_1 + alpha_2 + 3)
    # note that due to the precursor cell states there are even 2 states for each cell when alpha is 1
    ini_cond = np.zeros(no_of_states)
    ini_cond[0] = initial_cells
    
    #hill_coeff = feedback_dict[feedback_type]
    #hill_1 = hill_coeff[0]
    #hill_2 = hill_coeff[1]  
    state = odeint(th_cell_diff, ini_cond, simulation_time, 
                   args = (alpha_1, alpha_2, beta_1, beta_2, conc_il12, hill_1, hill_2,
                           fb_strength, rate_ifn, rate_il4, K, degradation, fb_start, fb_end))

    return state

def get_readouts(state, alpha, initial_cells, simulation_time):
    """
    from single model simulation, extract final cell states and half time to final state
    """
    stepsize = simulation_time[-1] / (len(simulation_time) - 1)       
    th1_endstate = state[-1, alpha + 1]       
    th2_endstate = state[-1, -1]
    
    th1_halfmax = th1_endstate / 2
    th2_halfmax = th2_endstate / 2        
    th1_tau = find_nearest(state[:, alpha + 1], th1_halfmax) * stepsize
    th2_tau = find_nearest(state[:, -1], th2_halfmax) * stepsize
    
    norm = initial_cells / 100      
    th1_endstate = np.asarray(th1_endstate) / norm
    th2_endstate = np.asarray(th2_endstate) / norm

    return [th1_endstate, th2_endstate, th1_tau, th2_tau]
   
def variable_effect(
         variable_arr,
         variable_name,
         alpha_1,
         alpha_2,
         beta_1,
         beta_2,
         simulation_time,
         conc_il12,
         hill_1,
         hill_2,
         fb_strength,
         rate_ifn,
         rate_il4,
         K,
         initial_cells,
         degradation,
         fb_start,
         fb_end,
         ):
    """
    vary model variable (variable_arr) by providing "variable_name"
    and simulate model for each element in variable_arr
    gives cell final states and half time to final state as output
    """
    readouts = []
    
    for i in variable_arr:

        if variable_name == "IL12":
            conc_il12 = i
        
        if variable_name == "chain":
            assert variable_arr[0] >= 1, "chain needs at least one step"
            alpha_1 = int(i)
            alpha_2 = int(i)
            beta_1 = float(i)
            beta_2 = float(i)
                    
        if variable_name == "feedback_strength_pos_Th1":
            fb_strength[0][0] = i
        
        if variable_name == "feedback_duration":
            fb_end = i
            assert fb_end < simulation_time[-1]

        if variable_name == "feedback_timing":
            fb_start = i
            fb_end = fb_start + 1.
            assert fb_end < simulation_time[-1]
        
        if variable_name == "cytokine_rates":
            rate_ifn = i
            rate_il4 = i
                
        state = run_model(alpha_1,
                  alpha_2,
                  beta_1,
                  beta_2,
                  simulation_time,
                  conc_il12,
                  hill_1,
                  hill_2,
                  fb_strength,
                  rate_ifn,
                  rate_il4,
                  K,
                  initial_cells,
                  degradation,
                  fb_start,
                  fb_end,
                  )
        
        readouts.append(get_readouts(state, alpha_1, initial_cells, simulation_time))
        
    return readouts

def get_cell_arr(readouts, readout_type = "cells"):
    """
    take array from fun "variable effect" and return cell final state or half time to final state as list"
    """
    if readout_type == "cells":
        th1_cells = [item[0] for item in readouts]
        th2_cells = [item[1] for item in readouts]
        return [th1_cells, th2_cells]
    
    if readout_type == "tau":
        tau_th1 = [item[2] for item in readouts]
        tau_th2 = [item[3] for item in readouts]
        return [tau_th1, tau_th2]
    

def assign_fb(fb_type, fb_dict, model_type, params):
    """
    take a set of parameters and change the model type
    possible types: single step model, rtm model, asymmetric model
    takes also feedback dictionary and a string "fb_type" to choose which feedback to simulate
    """
    params = list(params)
    params[6] = fb_dict[fb_type][0]
    params[7] = fb_dict[fb_type][1]
    
    assert model_type == "rate" or "rtm" or "rate_th1_rtm_th2" or "rtm_th1_rate_th2"

    if model_type == "rate":
        alpha_1, alpha_2, beta_1, beta_2 = 1, 1, 1, 1
        
    if model_type == "rtm":
        alpha_1, alpha_2, beta_1, beta_2 = 20, 20, 20, 20
    
    if model_type == "rate_th1_rtm_th2":
        alpha_1, alpha_2, beta_1, beta_2 = 1, 20, 1, 20
    
    if model_type == "rtm_th1_rate_th2":
        alpha_1, alpha_2, beta_1, beta_2 = 20, 1, 20, 1
        
    params[0] = alpha_1
    params[1] = alpha_2
    params[2] = beta_1
    params[3] = beta_2
    
    return params    