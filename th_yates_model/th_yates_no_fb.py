#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Thu Jan 24 11:20:42 2019

@author: burt
analysis of the th1 th2 diff model (deterministic) with prolif from yates et al
corresponding parameters stored in "det_model_params.py"
and model stored in "det_th1_th2_model.py"
note that this model has no feedback
note that this only works if th2 is seen as th1 with alpha = 1 (rate process)
this is possible because circuits are disconnected
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
save_path = "/home/burt/Documents/phd_summary/20190128/"

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

params = [cells_0,
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

rate_params = [cells_0,
              time,
              1,
              alpha_2,
              1.,
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
state = run_model(*params)

th1_cells, th2_cells = get_cells(state, alpha_1, alpha_2, alpha_prolif)

#fig, ax = plt.subplots(figsize = (5, 4))
#ax.plot(time, th1_cells)
#ax.plot(time, th2_cells, "tab:blue", linestyle = "--")
#ax.set_xlabel("time")
#ax.set_ylabel("n cells")
#ax.set_ylim(ymin = 0)
#ax.set_xlim(0, time[-1])
#plt.tight_layout()

#fig.savefig("one_celltype_time.pdf", bbox_inches = "tight")

#==============================================================================
# plot with tau and max readouts and effect of alpha diff
#==============================================================================
tau = get_tau(time, th1_cells)
th1_max = th1_cells[-1]
th1_halfmax = th1_max / 2.

alpha_arr = np.arange(1,11,1)

params_chain = list(params)
param_name = "chain"

alpha_sim = vary_param(param_name, alpha_arr, *params_chain)

fig, axes = plt.subplots(1, 3, figsize = (14,4))

ax = axes[0]
ax1 = axes[1]
ax2 = axes[2]

ax.plot(time, th1_cells)
ax.plot(time, th2_cells, "tab:blue", linestyle = "--")
ax.hlines(th1_max, 0, time[-1], linestyles = ":", colors = "tab:blue")
ax.hlines(th1_halfmax, 0, tau, linestyles = ":", colors = "tab:blue")
ax.vlines(tau, 0, th1_halfmax, linestyles = ":", colors = "tab:blue")

ax.set_xlabel("time")
ax.set_ylabel(r"$n_\mathrm{cells}$")

ax.set_xlim(0, time[-1])
ax.set_ylim(ymin = 0)

ax.set_yticks([0, th1_halfmax, th1_max])
s = r"$\frac{1}{2}\mathrm{Th1}_\mathrm{max}$"
s2 = r"$\mathrm{Th1}_\mathrm{max}$"
ax.set_yticklabels(["0", s, s2])
ax.set_xticks([0, tau, time[-1]])
ax.set_xticklabels(["0", r"$\tau_{1/2}$", str(round(time[-1],0))])

ax1.scatter(alpha_arr, alpha_sim[0])
ax1.set_xlabel(r"$\alpha$ differentiation")
ax1.set_ylabel(r"$\mathrm{Th1}_\mathrm{max}$")
ax1.set_xticks([0,5,10])

ax2.scatter(alpha_arr, alpha_sim[1])
ax2.set_xlabel(r"$\alpha$ differentiation")
ax2.set_ylabel(r"$\tau_{1/2}$")
ax2.set_xticks([0,5,10])
#for ax in axes: adjust_lim(ax)

ax.legend(["rtm", "rate process"], loc = 1)

plt.tight_layout()

fig.savefig(save_path+"one_celltype_diff.pdf", bbox_inches = "tight")
#==============================================================================
# vary alpha prolif
#==============================================================================
param_name = "prolif"
param_arr = np.arange(1,11,1)

params_prolif = list(params)
rtm_arr = vary_param(param_name, param_arr, *params_prolif)

fig, ax = plt.subplots(1,2, figsize = (9,4))

ax[0].scatter(param_arr, rtm_arr[0], c = "tab:blue", label = "rtm")
ax[1].scatter(param_arr, rtm_arr[2], c = "tab:blue")
ax[0].scatter(param_arr, rtm_arr[1], c = "tab:blue", marker = "x", label = "rate process")
ax[1].scatter(param_arr, rtm_arr[3], c = "tab:blue", marker = "x")


ax[0].set_xlabel(r"$\alpha$ proliferation")
ax[0].set_ylabel(r"% Th1 cells after "+str(round(time[-1],1))+ " hrs")
ax[1].set_xlabel(r"$\alpha$ proliferation")
ax[1].set_ylabel(r"$\tau_{1/2}$")
ax[0].legend()
plt.tight_layout()

fig.savefig(save_path+"one_celltype_prolif.pdf", bbox_inches = "tight")
#==============================================================================
# death rate
#==============================================================================
param_name = "rate_death"
param_arr = np.arange(0.5,1,0.02)

rtm_params_death = list(params)
rtm_arr = vary_param(param_name, param_arr, *rtm_params_death)

fig, ax = plt.subplots(1,2, figsize = (9,4))

ax[0].plot(param_arr, rtm_arr[0], c = "tab:blue", label = "rtm")
ax[0].plot(param_arr, rtm_arr[1], c = "tab:blue", linestyle = "--", label = "rate process")
ax[1].plot(param_arr, rtm_arr[2], c = "tab:blue")
ax[1].plot(param_arr, rtm_arr[3], c = "tab:blue", linestyle = "--")


ax[0].set_xlabel("death rate d")
ax[0].set_ylabel(r"% Th1 cells after "+str(round(time[-1],1))+ " hrs")
ax[1].set_xlabel("death rate d")
ax[1].set_ylabel(r"$\tau_{1/2}$")

plt.tight_layout()

fig.savefig(save_path+"one_celltype_death.pdf", bbox_inches = "tight")
#==============================================================================
# birth rate
#==============================================================================
param_name = "rate_birth"
param_arr = np.arange(0.1, 3, 0.1)

rtm_params_birth = list(params)
rtm_arr = vary_param(param_name, param_arr, *rtm_params_birth)

fig, ax = plt.subplots(1,2, figsize = (9,4))

ax[0].plot(param_arr, rtm_arr[0], c = "tab:blue")
ax[0].plot(param_arr, rtm_arr[1], c = "tab:blue", linestyle = "--")
ax[1].plot(param_arr, rtm_arr[2], c = "tab:blue")
ax[1].plot(param_arr, rtm_arr[3], c = "tab:blue", linestyle = "--")

ax[0].set_xlabel("birth rate a")
ax[0].set_ylabel(r"% Th1 cells after "+str(round(time[-1],1))+ " hrs")
ax[1].set_xlabel("birth rate a")
ax[1].set_ylabel(r"$\tau_{1/2}$")

plt.tight_layout()

fig.savefig(save_path+"one_celltype_birth.pdf", bbox_inches = "tight")
#==============================================================================
# now feedback comes into play
#==============================================================================
param_name = "fb_th1"
param_arr = np.arange(-1,10,0.1)

rtm_params_fb = list(params)
rate_params_fb = list(rate_params)

rtm_arr = vary_param(param_name, param_arr, *rtm_params_fb)
rate_arr = vary_param(param_name, param_arr, *rate_params_fb)

fig, ax = plt.subplots(1,2, figsize = (9,4))

ax[0].plot(param_arr, rtm_arr[0], c = "tab:blue")
ax[0].plot(param_arr, rate_arr[0], c = "tab:blue", linestyle = "--")
ax[1].plot(param_arr, rtm_arr[2], c = "tab:blue")
ax[1].plot(param_arr, rate_arr[2], c = "tab:blue", linestyle = "--")

ax[0].set_xlabel("feedback strength")
ax[0].set_ylabel(r"% Th1 cells after "+str(round(time[-1],1))+ " hrs")
ax[1].set_xlabel("feedback strength")
ax[1].set_ylabel(r"$\tau_{1/2}$")

plt.tight_layout()

fig.savefig(save_path+"one_celltype_feedback.pdf", bbox_inches = "tight")