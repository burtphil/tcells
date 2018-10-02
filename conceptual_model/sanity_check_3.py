#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Thu Sep 27 17:47:00 2018

@author: burt
"""
import numpy as np
from scipy.integrate import odeint
import matplotlib.pyplot as plt

def th_cell_diff(state,t,alpha_1,alpha_2,beta_1,beta_2, th0_influx = 0, degradation = 0):
        
    th_0 = state[0]
    th1_0 = 0.5*th_0
    th2_0 = 0.5*th_0

    #print th1_0,th2_0
    # assign th1 states and th2 states from initial vector based on chain length alpha
    th1 = np.concatenate(([th1_0],state[1:(alpha_1+1)]))
    th2 = np.concatenate(([th2_0],state[(alpha_1+1):]))
      
    dt_th1 = np.ones_like(th1)
    dt_th2 = np.ones_like(th2)
        
    th_states = [th1,th2]
    dt_th_states = [dt_th1,dt_th2]
    rate = [beta_1,beta_2]
        
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
                dt_state[j] = r*th_state[j-1]-degradation

    # assign new number of naive cells based on the number of present naive cells that were designated th1_0 or th2_0
    #dt_th0 = state[0]
    dt_th0 = dt_th1[0]+dt_th2[0]+th0_influx
    
    # return cell states
    dt_state = np.concatenate(([dt_th0],dt_th1[1:],dt_th2[1:]))
    
    return dt_state  

def run_model(title, parameters, model_name = th_cell_diff, save = True):
    """ 
    run ode model with parameters from params file
    needs shape and rate params as input as well as simulation time
    """
    alpha_1, alpha_2, beta_1, beta_2, simulation_time = parameters
    # define initial conditions based on number of intermediary states
    no_of_states = int(alpha_1+alpha_2+1)
    ini_cond = np.zeros(no_of_states)
    ini_cond[0] = 10000
    
    state = odeint(model_name, ini_cond, simulation_time, 
                   args =(alpha_1,alpha_2,beta_1,beta_2,))
    
    # save both model simulation and associated parameters
    if save == True:
        parameters = np.asarray(parameters, dtype = object)
        np.savez(title, state = state, parameters = parameters)
    return state



#==============================================================================
# import parameters
#==============================================================================
simulation_time = np.arange(0,10,0.01)
conceptual_params = [1,1.,1,1., simulation_time]
initial_cells= 10000

#==============================================================================
# run time course simulation
#==============================================================================
simulation_name = "sim"

th1_th2_model = run_model(simulation_name, parameters = conceptual_params)

test_simulation = np.load(simulation_name+".npz")
state = test_simulation["state"]
parameters = test_simulation["parameters"]

# plot time course
#alpha_params = list(conceptual_params)
alpha_params = list(conceptual_params)
norm = initial_cells/100

fig, ax = plt.subplots(1,1, figsize = (5,4))

for i in range(1,3):
    alpha_params = list(conceptual_params)
    alpha_params[0] = int(i)
    alpha_params[2] = int(i)
    state = run_model(simulation_name, parameters = alpha_params)
    th0_cells = state[:,0]/norm
    th2_cells = state[:,-1]/norm
    alpha_1 = alpha_params[0]
    th1_cells = state[:,int(alpha_1)]/norm
    #ax.plot(simulation_time, th0_cells, label = r"Th0, $\alpha=$"+str(i))
    ax.plot(simulation_time, th1_cells, label = r"Th1, $\alpha=$"+str(i))
    ax.plot(simulation_time, th2_cells, label = r"Th2, $\alpha_{Th1} =$"+str(i))
    
    ax.set_yticks([0,50,100])
    #ax.set_ylim([0,100])
    ax.set_xlim([simulation_time[0],simulation_time[-1]])
    ax.set_xlabel("Time [h]")
    ax.set_ylabel("% Th cells")
    ax.set_title(r"$\alpha_{Th1} = 1$")
    ax.legend()

plt.tight_layout()


"""
fig, ax = plt.subplots(1,1, figsize = (5,4))

for i in range(1,3):
    alpha_params = list(conceptual_params)
    alpha_params[0] = int(i)
    alpha_params[2] = int(i)
    state = run_model(simulation_name, parameters = alpha_params)
    th0_cells = state[:,0]
    th2_cells = state[:,-1]
    th1_cells = state[:,int(i)]
    
    #th1_cells = th1_cells/(th1_cells+th2_cells+th0_cells)
    #th2_cells = th2_cells/(th1_cells+th2_cells+th0_cells)
    #ax.plot(simulation_time, th0_cells)
    ax.plot(simulation_time, th1_cells, label = "Th1, alpha="+str(i))
    ax.plot(simulation_time, th2_cells, label = "Th2, alpha Th1 ="+str(i))
    
    #ax.set_yticks([0,50,100])
    #ax.set_ylim([0,100])
    ax.set_xlim([simulation_time[0],simulation_time[-1]])
    ax.set_xlabel("Time [h]")
    ax.set_ylabel("% Th cells")
    ax.legend()

plt.tight_layout()
"""