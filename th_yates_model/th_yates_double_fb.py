#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Thu Jan 24 11:20:42 2019

@author: burt
analysis of the th1 th2 diff model (deterministic) with prolif from yates et al
corresponding parameters stored in "det_model_params.py"
and model stored in "det_th1_th2_model.py"
here, I analyze a th1 th2 double feedback system and vary only one feedback
make sure to check that params are adjusted accordingly
the question here is if the feedback system is more robust
if we look at rtm process compared with rate process
also, finding is that with same means, rate and rtm process differ
"""
#==============================================================================
# import modules
#==============================================================================
from scipy.integrate import odeint
import numpy as np
import matplotlib.pyplot as plt
# own modules
from det_th1_th2_model import th_cell_diff
import th_diff_model_parameters as params
import seaborn as sns

#==============================================================================
# choose plotting style and saving path
#==============================================================================
sns.set(context = "talk", style = "ticks")
save_path = "/home/burt/Documents/tcell_project/figures/"

#==============================================================================
# functions
#==============================================================================
def remove_spines(ax):
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)
    ax.get_xaxis().tick_bottom()
    ax.get_yaxis().tick_left()

def adjust_lim(ax):
    ax.set_xlim(xmin = 0)
    ax.set_ylim(ymin = 0)

def get_cells(state, alpha_1, alpha_2, alpha_prolif):
    
    th1_cells = state[:, alpha_1:(alpha_1+alpha_prolif)]
    th2_cells = state[:, (alpha_1+alpha_prolif+alpha_2):]
    
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
    rate_cytos,
    fb_th1,
    fb_th2,
    hill_th1,
    hill_th2,
    K_th1,
    K_th2,
    conc_il12, 
    ):
    
    no_of_states = int(alpha_1 + alpha_2 + alpha_prolif + alpha_prolif)
    y0 = np.zeros(no_of_states)    
    y0[0] = cells_0[0]
    y0[alpha_1+alpha_prolif] = cells_0[1]
    
    state = odeint(th_cell_diff, 
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
        
        elif param_name == "rate_cyto":
            rate_cytos = [i, i]
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
#==============================================================================
# set up parameters
#==============================================================================
time = params.time
cells_0 = params.cells_0

alpha_1 = params.alpha_1
alpha_2 = params.alpha_2

beta_1 = params.beta_1
beta_2 = params.beta_2

alpha_prolif = params.alpha_prolif
beta_prolif = params.beta_prolif

rate_birth = params.rate_birth
rate_death = params.rate_death

rate_cytos = params.rate_cytos
fb_th1 = params.fb_th1
fb_th2 = params.fb_th2
hill_th1 = params.hill_th1
hill_th2 = params.hill_th2
K_th1 = params.K_th1
K_th2 = params.K_th2
conc_il12 = params.conc_il12

#assert fb_th1[0] != 0
#assert fb_th1[1] != 0
#assert fb_th2[0] != 0
#assert fb_th2[1] != 0

rtm_params = [cells_0,
              time,
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
              conc_il12,]

#==============================================================================
# run model w/o any modifications time course
#==============================================================================

# simple time course
state = run_model(*rtm_params)

th1_cells, th2_cells = get_cells(state, alpha_1, alpha_2, alpha_prolif)

fig, ax = plt.subplots(figsize = (5,4))
ax.plot(time, th1_cells)
ax.plot(time, th2_cells, "tab:red", linestyle = "--")
ax.set_xlabel("time")
ax.set_ylabel("n cells")
ax.set_ylim(ymin = 0)
ax.set_xlim(0, time[-1])
plt.tight_layout()

#==============================================================================
# vary alpha diff
#==============================================================================
param_name = "chain"
param_arr = np.arange(1,11,1)

rtm_params_chain = list(rtm_params)

rtm_arr = vary_param(param_name, param_arr, *rtm_params_chain)

fig, ax = plt.subplots(1,2, figsize = (9,4))

ax[0].scatter(param_arr, rtm_arr[0], c = "tab:blue")
ax[0].scatter(param_arr, rtm_arr[1], c = "tab:red")
ax[1].scatter(param_arr, rtm_arr[2], c = "tab:blue")
ax[1].scatter(param_arr, rtm_arr[3], c = "tab:red")

ax[0].set_xlabel(r"$\alpha$ differentiation")
ax[0].set_ylabel(r"% Th1 cells after "+str(round(time[-1],1))+ " hrs")
ax[1].set_xlabel(r"$\alpha$ differentiation")
ax[1].set_ylabel(r"$\tau_{1/2}$")

plt.tight_layout()

#==============================================================================
# vary alpha prolif
#==============================================================================
param_name = "prolif"
param_arr = np.arange(1,11,1)

rtm_params_prolif = list(rtm_params)
rtm_arr = vary_param(param_name, param_arr, *rtm_params_prolif)

fig, ax = plt.subplots(1,2, figsize = (9,4))

ax[0].scatter(param_arr, rtm_arr[0], c = "tab:blue")
ax[0].scatter(param_arr, rtm_arr[1], c = "tab:red")
ax[1].scatter(param_arr, rtm_arr[2], c = "tab:blue")
ax[1].scatter(param_arr, rtm_arr[3], c = "tab:red")

ax[0].set_xlabel(r"$\alpha$ proliferation")
ax[0].set_ylabel(r"% Th1 cells after "+str(round(time[-1],1))+ " hrs")
ax[1].set_xlabel(r"$\alpha$ proliferation")
ax[1].set_ylabel(r"$\tau_{1/2}$")

plt.tight_layout()



#==============================================================================
# death rate
#==============================================================================
param_name = "rate_death"
param_arr = np.arange(0.5,1,0.02)

rtm_params_death = list(rtm_params)
rtm_arr = vary_param(param_name, param_arr, *rtm_params_death)

fig, ax = plt.subplots(1,2, figsize = (9,4))

ax[0].plot(param_arr, rtm_arr[0], c = "tab:blue")
ax[0].plot(param_arr, rtm_arr[1], c = "tab:red", linestyle = "--")
ax[1].plot(param_arr, rtm_arr[2], c = "tab:blue")
ax[1].plot(param_arr, rtm_arr[3], c = "tab:red", linestyle = "--")


ax[0].set_xlabel("death rate d")
ax[0].set_ylabel(r"% Th1 cells after "+str(round(time[-1],1))+ " hrs")
ax[1].set_xlabel("death rate d")
ax[1].set_ylabel(r"$\tau_{1/2}$")

plt.tight_layout()

#==============================================================================
# birth rate
#==============================================================================
param_name = "rate_birth"
param_arr = np.arange(0.1, 1, 0.1)

rtm_params_birth = list(rtm_params)
rtm_arr = vary_param(param_name, param_arr, *rtm_params_birth)

fig, ax = plt.subplots(1,2, figsize = (9,4))

ax[0].plot(param_arr, rtm_arr[0], c = "tab:blue")
ax[0].plot(param_arr, rtm_arr[1], c = "tab:red", linestyle = "--")
ax[1].plot(param_arr, rtm_arr[2], c = "tab:blue")
ax[1].plot(param_arr, rtm_arr[3], c = "tab:red", linestyle = "--")


ax[0].set_xlabel("birth rate a")
ax[0].set_ylabel(r"% Th1 cells after "+str(round(time[-1],1))+ " hrs")
ax[1].set_xlabel("birth rate a")
ax[1].set_ylabel(r"$\tau_{1/2}$")

plt.tight_layout()

#==============================================================================
# vary the feedback strength
#==============================================================================
param_name = "fb_th1"
param_arr = np.arange(-1,20,0.1)

rtm_params_fb = list(rtm_params)
rtm_arr = vary_param(param_name, param_arr, *rtm_params_fb)

fig, ax = plt.subplots(1,2, figsize = (9,4))

ax[0].plot(param_arr, rtm_arr[0], c = "tab:blue")
ax[0].plot(param_arr, rtm_arr[1], c = "tab:red", linestyle = "--")
ax[1].plot(param_arr, rtm_arr[2], c = "tab:blue")
ax[1].plot(param_arr, rtm_arr[3], c = "tab:red", linestyle = "--")

#ax[0].set_xticks[-1,0,1]
ax[0].set_xlabel("fb strength")
ax[0].set_ylabel(r"% Th1 cells after "+str(round(time[-1],1))+ " hrs")
ax[1].set_xlabel("fb strength")
ax[1].set_ylabel(r"$\tau_{1/2}$")
#ax[1].set_xticks[-1,0,1]
plt.tight_layout()

#==============================================================================
# cytokine rate
#==============================================================================
param_name = "rate_cyto"
param_arr = np.arange(1, 3, 0.1)

rtm_params_cyto = list(rtm_params)
rtm_arr = vary_param(param_name, param_arr, *rtm_params_cyto)

fig, ax = plt.subplots(1,2, figsize = (9,4))

ax[0].plot(param_arr, rtm_arr[0], c = "tab:blue")
ax[0].plot(param_arr, rtm_arr[1], c = "tab:red", linestyle = "--")
ax[1].plot(param_arr, rtm_arr[2], c = "tab:blue")
ax[1].plot(param_arr, rtm_arr[3], c = "tab:red", linestyle = "--")


ax[0].set_xlabel("cytokine secretion rate r")
ax[0].set_ylabel(r"% Th1 cells after "+str(round(time[-1],1))+ " hrs")
ax[1].set_xlabel("cytokine secretion rate r")
ax[1].set_ylabel(r"$\tau_{1/2}$")

plt.tight_layout()
