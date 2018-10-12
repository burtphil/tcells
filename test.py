# -*- coding: utf-8 -*-
"""
Created on Fri Oct 12 15:31:31 2018

@author: Philipp
"""
import numpy as np
import matplotlib.pyplot as plt
from scipy.integrate import odeint
import modules.th1_th2_conceptual_parameters as params
from modules.th1_th2_plotting import plot_time_course
cparams = params.parameters

def p_menten(conc_ifn,conc_il4,conc_il12, hill, half_saturation,  strength = 1.):
   
    # watch out, this is tricky because if concentrations are zero and negative feedback (hill coeff < 0)
    # then you divide by zero, best avoid zero concentrations through base cytokine rate
    cytokine = np.array([conc_ifn,conc_il4, conc_il12])
    hill = np.array(hill)
    half_saturation = np.array(half_saturation)
    menten = cytokine**hill/(half_saturation+cytokine**hill)
    prob_th_diff = strength*np.prod(menten)
    
    assert prob_th_diff >= 0, "probability is "+str(prob_th_diff)+" must be non-negative."
    
    return prob_th_diff



def th_cell_diff(state,t,alpha_1,alpha_2,rate1,rate2,conc_il12, hill_1, hill_2, 
                 rate_ifn, rate_il4, half_saturation, degradation, fun_probability = p_menten):
        
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
    prob_th1 = fun_probability(conc_ifn,conc_il4,conc_il12, hill_1, half_saturation)
    prob_th2 = fun_probability(conc_ifn,conc_il4,conc_il12, hill_2, half_saturation)    
    assert prob_th1+prob_th2 > 0, "prob th1="+str(prob_th1)+" prob th2="+str(prob_th2)+" cannot normalize."

    # normalized branching probabilities
    p_1 = prob_th1/(prob_th1+prob_th2)
    p_2 = prob_th2/(prob_th1+prob_th2)
    #print p_1
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
                print th_state[j-1]
                dt_state[j] = r*th_state[j-1]-2000
    
    # assign new number of naive cells based on the number of present naive cells that were designated th1_0 or th2_0
    #dt_th0 = state[0]
    dt_th0 = 0
    
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
     rate_ifn, rate_il4, half_saturation, initial_cells, degradation) = parameters

    # define initial conditions based on number of intermediary states
    no_of_states = int(alpha_1+alpha_2+3)
    ini_cond = np.zeros(no_of_states)
    ini_cond[0] = initial_cells
    
    state = odeint(model_name, ini_cond, simulation_time, 
                   args =(alpha_1,alpha_2,rate1,rate2,conc_il12, hill_1,hill_2,
                          rate_ifn,rate_il4,half_saturation, degradation))
    
    # save both model simulation and associated parameters
    if save == True:
        parameters = np.asarray(parameters, dtype = object)
        np.savez(title, state = state, parameters = parameters)

    return state

state = run_model("test", cparams)

#plot_time_course(state, cparams)

th1_cells = state[:,3]
th2_cells = state[:,-2]

fig, ax = plt.subplots(1,1, figsize = (5,5))
ax.plot(params.simulation_time,th1_cells, "tab:blue",params.simulation_time,th2_cells, "tab:red")
#ax.set_yticks([0,50,100])
#ax.set_ylim([0,100])
ax.set_xlim([params.simulation_time[0],params.simulation_time[-1]])
ax.set_xlabel("Time")
ax.set_ylabel("% cells")
