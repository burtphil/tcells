#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Wed Dec  5 09:39:51 2018

@author: burt
"""
#==============================================================================
# modules
#==============================================================================
import numpy as np
import seaborn as sns
import matplotlib.pyplot as plt
import os
os.chdir("/home/burt/Documents/code/th_cell_differentiation")
import modules.th1_th2_conceptual_parameters as cparams
from scipy.integrate import odeint

#==============================================================================
# choose plotting style and saving path
#==============================================================================
sns.set(context = "talk", style = "ticks")
save_path = "/home/burt/Documents/tcell_project/figures/"

#==============================================================================
# choose type of feedback
#==============================================================================
feedback_dict = cparams.feedback_gamma
feedback_type = "no_fb"

#==============================================================================
# choose feedback duration
#==============================================================================
fb_start = 0
fb_end = cparams.stop
assert fb_end <= cparams.stop

#==============================================================================
# import parameters
#==============================================================================
alpha_1 = 1
alpha_2 = 10
beta_1 = float(alpha_1)
beta_2 = float(alpha_2)
K = cparams.K
degradation = 1.0
rate_ifn = cparams.rate_ifn
rate_il4 = cparams.rate_il4
conc_il12 = cparams.conc_il12
simulation_time = cparams.simulation_time
initial_cells = cparams.initial_cells
fb_strength = feedback_dict[feedback_type]
hill_1 = [1,1,1]
hill_2 = hill_1

rate_params = [alpha_1, alpha_2, beta_1, beta_2, simulation_time, conc_il12,
              hill_1, hill_2, fb_strength, rate_ifn, rate_il4, K,
              initial_cells, degradation, fb_start, fb_end]

rtm_params = [10, 10, 10, 10, simulation_time, conc_il12,
              hill_1, hill_2, fb_strength, rate_ifn, rate_il4, K,
              initial_cells, degradation, fb_start, fb_end]

#==============================================================================
# functions
#==============================================================================
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

def ax_time_course(state, ax, simulation_time, initial_cells, alpha_1, linestyle = "-"):
    
    norm = initial_cells / 100
    th1_cells = state[:,alpha_1] / norm
    th2_cells = state[:,-1] / norm
    
    #print th1_cells[-1],th2_cells[-1]
    #ax.plot(simulation_time, th0_cells, color = "k", linestyle = linestyle)
    ax.plot(simulation_time, th1_cells, color = "tab:blue", linestyle = "--")
    ax.plot(simulation_time, th2_cells, color = "tab:red", linestyle = linestyle)
    #ax.set_yticks([0, 50, 100])
    #ax.set_ylim([0, 100])
    ax.set_xlim([simulation_time[0], simulation_time[-1]])
    ax.set_xlabel("time")
    ax.set_ylabel("% Th cells")

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
                
    conc_ifn = rate_ifn*state[alpha_1]
    conc_il4 = rate_il4*state[-1]
     
    # branching probablities
    prob_th1 = fun_probability(conc_ifn, conc_il4, conc_il12, hill_1, K, fb_strength[0])
    prob_th2 = fun_probability(conc_ifn, conc_il4, conc_il12, hill_2, K, fb_strength[1])    
    assert prob_th1 + prob_th2 > 0, "prob th1="+str(prob_th1)+" prob th2="+str(prob_th2)+" cannot normalize."
    
    # normalized branching probabilities
    p_1 = prob_th1 / (prob_th1 + prob_th2)
    p_2 = prob_th2 / (prob_th1 + prob_th2)

    # assign th1 states and th2 states from initial vector based on chain length alpha
    th_0 = state[0] 
    th1_0 = p_1 * th_0     
    th2_0 = p_2 * th_0
    
    th1 = state[1:(alpha_1+1)]
    th2 = state[(alpha_1+1):]
      
    dt_th1 = np.ones_like(th1)
    dt_th2 = np.ones_like(th2)
    
    th_prec = [th1_0, th2_0]
    th_states = [th1,th2]
    dt_th_states = [dt_th1,dt_th2]
    rate = [beta_1, beta_2]
    #print beta_1
    #print rate
    ### differential equations depending on the number on intermediary states
    # loop over states of th1 and th2 cells
    for i in range(2):
        th_state = th_states[i]
        dt_state = dt_th_states[i]
        r = rate[i]
        th_p = th_prec[i]
            
    # calculate derivatives
        assert len(th_state) >= 1, "len th state is " + str(len(th_state))
            
        if len(th_state) == 1:
            dt_state[0] = r * th_p - degradation * th_state[0]
                
        elif len(th_state) > 1:
            for j in range(len(th_state)):
                if j == 0:
                    # this step is unclear, need to understand if it should
                    # be th_p - r * th_state or r * (th_p - th_state)
                    dt_state[j] = th_p - r * th_state[j]
                    
                elif j != (len(th_state)-1):
                    dt_state[j] = r * (th_state[j-1] - th_state[j])
                    
                else:
                    dt_state[j] = r * th_state[j-1] - degradation * th_state[j]

    # assign new number of naive cells based on the number of present naive cells that were designated th1_0 or th2_0
    # if a constant Th naive cell pool is assumed (default parameter const_thn = True) then change should be 0
    # because pool should not change   
    dt_th0 = 1 - state[0]
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
    no_of_states = int(alpha_1 + alpha_2 + 1)
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

#==============================================================================
# time course simulation
#==============================================================================
state = run_model(*rate_params)
#state_rtm = run_model(*rtm_params)

fig, ax = plt.subplots(figsize = (5,4))

norm = 100
th1_cells = state[:, alpha_1] * norm
th2_cells = state[:, -1] * norm


ax.plot(simulation_time, th1_cells, "tab:blue", linestyle = "-")
ax.plot(simulation_time, th2_cells, "tab:red", linestyle = "-")

ax.legend([r"$\alpha=1$", r"$\alpha=1$"])
ax.set_xlabel("time")
ax.set_ylabel("% Th cells")
plt.tight_layout()