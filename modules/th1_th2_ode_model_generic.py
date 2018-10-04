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

def p_th_diff(conc_ifn,conc_il4,conc_il12, hill, half_saturation, base_production_rate, strength = 1.):
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

    prob_th_diff = strength*ifn_prob*il4_prob*il12_prob+base_production_rate
    
    assert prob_th_diff >= 0, "probability is "+str(prob_th_diff)+" must be non-negative."
    
    return prob_th_diff

def p_menten(conc_ifn,conc_il4,conc_il12, hill, half_saturation, base_production_rate, strength = 1.):
    """
    returns probability of Th1 Th2 differentiation for given cytokine concentrations
    kinetics are hill-like so hill coefficients for cytokines also need to be provided
    """
    assert conc_ifn >= 0, "ifn conc is "+str(conc_ifn)+", concentrations must be non-negative."
    assert conc_il4 >= 0, "il4 conc is "+str(conc_il4)+", concentrations must be non-negative."
    assert conc_il12 >= 0, "il12conc is "+str(conc_il12)+", concentrations must be non-negative."
    
    # watch out, this is tricky because if concentrations are zero and negative feedback (hill coeff < 0)
    # then you divide by zero, best avoid zero concentrations through base cytokine rate
    cytokine = np.array([conc_ifn,conc_il4, conc_il12])
    hill = np.array(hill)
    half_saturation = np.array(half_saturation)
    menten = cytokine**hill/(half_saturation+cytokine**hill)
    prob_th_diff = strength*np.prod(menten)+base_production_rate
    
    assert prob_th_diff >= 0, "probability is "+str(prob_th_diff)+" must be non-negative."
    
    return prob_th_diff



def th_cell_diff(state,t,alpha_1,alpha_2,rate1,rate2,conc_il12, hill_1, hill_2, 
                 rate_ifn, rate_il4, half_saturation, base_production_rate_ifn, base_production_rate_il4, degradation):
        
    # calculate interferon gamma (ifn) and il4 concentrations based on the number of th1 and th2 cells
    assert type(alpha_1) == int, "alpha dtype is "+str(type(alpha_1))+" but alpha must be an integer."
    assert type(alpha_2) == int, "alpha dtype is "+str(type(alpha_2))+" but alpha must be an integer."
        
    base_cytokine_rate = 0.00001
    conc_ifn = rate_ifn*state[alpha_1+1]+base_cytokine_rate
    conc_il4 = rate_il4*state[-1]+base_cytokine_rate

    ### calculate initial th1 and th2 populations from naive cells based on branching probabilities
    
    th_0 = state[0]
    
    if t<1:   
        assert th_0 > 0, "no initial cells provided or cells"
    
    # branching probablities
    prob_th1 = p_th_diff(conc_ifn,conc_il4,conc_il12, hill_1, half_saturation, base_production_rate_ifn)
    prob_th2 = p_th_diff(conc_ifn,conc_il4,conc_il12, hill_2, half_saturation, base_production_rate_il4)    
    assert prob_th1+prob_th2 > 0, "prob th1="+str(prob_th1)+" prob th2="+str(prob_th2)+" cannot normalize."

    # normalized branching probabilities
    p_1 = prob_th1/(prob_th1+prob_th2)
    p_2 = prob_th2/(prob_th1+prob_th2)
    #print p_1,p_2
    # assign th1 states and th2 states from initial vector based on chain length alpha
    th1 = state[1:(alpha_1+2)]
    th2 = state[(alpha_1+2):]
      
    dt_th1 = np.ones_like(th1)
    dt_th2 = np.ones_like(th2)
        
    th_states = [th1,th2]
    dt_th_states = [dt_th1,dt_th2]
    rate = [rate1,rate2]
    prob = [p_1,p_2]
        
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
                dt_state[j] = p*th_0-r*th_state[j]
            elif j != (len(th_state)-1):
                dt_state[j] = r*(th_state[j-1]-th_state[j])
            else:
                dt_state[j] = r*th_state[j-1]-degradation

    # assign new number of naive cells based on the number of present naive cells that were designated th1_0 or th2_0
    #dt_th0 = state[0]
    dt_th0 = -(p_1+p_2)*th_0
    
   # return cell states
    dt_state = np.concatenate(([dt_th0],dt_th1,dt_th2))
    
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
    no_of_states = int(alpha_1+alpha_2+3)
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

def chain(chain_length, parameters, stepsize = 0.01):
    """
    plot steady state and tau 1/2 dependency on chain length
    watch out for the step size
    """
    chain = np.arange(1,chain_length,1)
    
    mean_th1 = parameters[0]/parameters[2]
    mean_th2 = parameters[1]/parameters[3]
    
    th1_conc = []
    th2_conc = []
    
    th1_tau = []
    th2_tau = []
    
    for i in chain:
        
        i = int(i)
        parameters[0] = i
        parameters[1] = i
        parameters[2] = i/mean_th1
        parameters[3] = i/mean_th2
        
        state = run_model("", parameters, save = False)
        #plot_time_course(state, alpha_1, alpha_2, simulation_time)
    
        # get end states
        th1_endstate = state[-1,int(i)+1]
        th1_conc.append(th1_endstate)
        
        th2_endstate = state[-1,-1]
        th2_conc.append(th2_endstate)
        
        th1_halfmax = th1_endstate/2
        th2_halfmax = th2_endstate/2
        
        th1_tau_idx = find_nearest(state[:,int(i)], th1_halfmax)*stepsize
        th2_tau_idx = find_nearest(state[:,-1], th2_halfmax)*stepsize
        th1_tau.append(th1_tau_idx)
        th2_tau.append(th2_tau_idx)    
        # normalize to initial cell pop
    
    norm = parameters[-2]/100
    th1_conc = np.array(th1_conc)/norm
    th2_conc = np.array(th2_conc)/norm
    
    return [chain, th1_conc, th2_conc, th1_tau, th2_tau, chain_length]

def chain_one(chain_length, parameters, alpha_idx, stepsize = 0.01):
    """
    plot steady state and tau 1/2 dependency on chain length
    watch out for the step size
    alpha_idx takes a tuple, the first index specifies which alpha is varied
    set to 0 for th1 alpha variation and to 1 for th2 alpha variation
    the second index sets the respective other alpha to the index value
    this needs to be an integer
    """
    parameters = list(parameters)
    chain = np.arange(1,chain_length,1)
    
    mean_th1 = parameters[0]/parameters[2]
    mean_th2 = parameters[1]/parameters[3]
    
    th1_conc = []
    th2_conc = []
    
    th1_tau = []
    th2_tau = []
    
    for i in chain:
        
        i = int(i)
        if alpha_idx[0] == 0:            
            parameters[0] = i
            parameters[2] = i/mean_th1
            parameters[1] = alpha_idx[1]
            parameters[3] = alpha_idx[1]/mean_th2

        elif alpha_idx[0] == 1:
            parameters[0] = alpha_idx[1]
            parameters[2] = alpha_idx[1]/mean_th1                        
            parameters[1] = i
            parameters[3] = i/mean_th2
        
        state = run_model("", parameters, save = False)
        #plot_time_course(state, alpha_1, alpha_2, simulation_time)
    
        # get end states
        th1_endstate = state[-1,int(i)+1]
        th1_conc.append(th1_endstate)
        
        th2_endstate = state[-1,-1]
        th2_conc.append(th2_endstate)
        
        th1_halfmax = th1_endstate/2
        th2_halfmax = th2_endstate/2
        
        th1_tau_idx = find_nearest(state[:,int(i)], th1_halfmax)*stepsize
        th2_tau_idx = find_nearest(state[:,-1], th2_halfmax)*stepsize
        th1_tau.append(th1_tau_idx)
        th2_tau.append(th2_tau_idx)    
        # normalize to initial cell pop
    
    norm = parameters[-2]/100
    th1_conc = np.array(th1_conc)/norm
    th2_conc = np.array(th2_conc)/norm
    
    return [chain, th1_conc, th2_conc, th1_tau, th2_tau, chain_length]

def chain_th1(chain_length, parameters, stepsize = 0.01):
    """
    plot steady state and tau 1/2 dependency on chain length
    watch out for the step size
    alpha_idx takes a tuple, the first index specifies which alpha is varied
    set to 0 for th1 alpha variation and to 1 for th2 alpha variation
    the second index sets the respective other alpha to the index value
    this needs to be an integer
    """
    parameters = list(parameters)
    chain = np.arange(1,chain_length,1)
    
    mean_th1 = parameters[0]/parameters[2]
    
    th1_conc = []
    th2_conc = []
    
    th1_tau = []
    th2_tau = []
    
    for i in chain:
        
        i = int(i)
        
        parameters[0] = i
        parameters[2] = i/mean_th1
        #print parameters[1],parameters[3]
        state = run_model("", parameters, save = False)
        #fig, ax = plt.subplots(1,1, figsize = (5,5))
        #ax.plot(parameters[4],state)
        #ax.set_title(str(parameters[0])+", "+str(parameters[1]))
        #plot_time_course(state, alpha_1, alpha_2, simulation_time)
        
        # get end states
        th1_endstate = state[-1,int(i)+1]
        th1_conc.append(th1_endstate)

        th2_endstate = state[-1,-1]
        th2_conc.append(th2_endstate)
        
        #print th1_endstate, th2_endstate
            
    norm = parameters[-2]/100
    th1_conc = np.array(th1_conc)/norm
    th2_conc = np.array(th2_conc)/norm
    
    return [chain, th1_conc, th2_conc, th1_tau, th2_tau, chain_length]

def chain_th2(chain_length, parameters, stepsize = 0.01):
    """
    plot steady state and tau 1/2 dependency on chain length
    watch out for the step size
    alpha_idx takes a tuple, the first index specifies which alpha is varied
    set to 0 for th1 alpha variation and to 1 for th2 alpha variation
    the second index sets the respective other alpha to the index value
    this needs to be an integer
    """
    parameters = list(parameters)
    chain = np.arange(1,chain_length,1)
    
    mean_th2 = parameters[1]/parameters[3]
    
    th1_conc = []
    th2_conc = []
    
    th1_tau = []
    th2_tau = []
    
    for i in chain:
        
        i = int(i)
        parameters[1] = i
        parameters[3] = i/mean_th2
        
        state = run_model("", parameters, save = False)
        #plot_time_course(state, alpha_1, alpha_2, simulation_time)
    
        # get end states
        th1_endstate = state[-1,parameters[0]+1]
        th1_conc.append(th1_endstate)
        
        th2_endstate = state[-1,-1]
        th2_conc.append(th2_endstate)
        
        th1_halfmax = th1_endstate/2
        th2_halfmax = th2_endstate/2
        
        th1_tau_idx = find_nearest(state[:,parameters[0]], th1_halfmax)*stepsize
        th2_tau_idx = find_nearest(state[:,-1], th2_halfmax)*stepsize
        th1_tau.append(th1_tau_idx)
        th2_tau.append(th2_tau_idx)    
        # normalize to initial cell pop
    
    norm = parameters[-2]/100
    th1_conc = np.array(th1_conc)/norm
    th2_conc = np.array(th2_conc)/norm
    
    return [chain, th1_conc, th2_conc, th1_tau, th2_tau, chain_length]

def il12(il12_conc, parameters):
    """
    plot steady state dependency of il12
    """
    th1_conc = []
    th2_conc = []
    th0_conc = []
    
    norm = parameters[-2]/100
    for i in il12_conc:
        
        parameters[5] = i
        
        state = run_model("", parameters)
        
        th0_endstate = state[-1,0]
        th0_conc.append(th0_endstate)
            
        th1_endstate = state[-1,int(parameters[0])+1]
        th1_conc.append(th1_endstate)
        
        th2_endstate = state[-1,-1]
        th2_conc.append(th2_endstate)
        
    th0_conc = np.asarray(th0_conc)/norm
    th1_conc = np.asarray(th1_conc)/norm
    th2_conc = np.asarray(th2_conc)/norm
    
    return [il12_conc, th1_conc, th2_conc]
  