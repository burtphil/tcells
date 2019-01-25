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
            hill_1[0] = i
        
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