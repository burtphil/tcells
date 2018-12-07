#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Fri Nov 30 11:26:19 2018

@author: burt
phd summary report from november 30.th 2018
should be able to reproduce all plots shown in the report together with the
scripts th1_th2_ode_model_generic.py
th1_th2_plotting.py
gamma_dist.py
only exception: stochastic simulations, for that, look into 
"""

#==============================================================================
# import modules
#==============================================================================
import numpy as np
import seaborn as sns
import matplotlib.pyplot as plt
import os
os.chdir("/home/burt/Documents/code/th_cell_differentiation")
import modules.th1_th2_conceptual_parameters as cparams
from modules.th1_th2_ode_model_generic import run_model
from modules.th1_th2_ode_model_generic import variable_effect
from modules.th1_th2_ode_model_generic import get_cell_arr
from modules.th1_th2_ode_model_generic import get_cyto_conc
from modules.th1_th2_ode_model_generic import get_prob
from modules.th1_th2_plotting import ax_time_course
from modules.th1_th2_plotting import ax_var_effect
from modules.th1_th2_plotting import ax_chain

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

rtm_params = [20, 20, 20, 20, simulation_time, conc_il12,
              hill_1, hill_2, fb_strength, rate_ifn, rate_il4, K,
              initial_cells, degradation, fb_start, fb_end]

#==============================================================================
# time course simulation
#==============================================================================
state = run_model(*rate_params)
state_rtm = run_model(*rtm_params)

norm = 100
th1_cells = state[:, alpha_1+1] * norm
th2_cells = state[:, -1] * norm

simulation_time2 = simulation_time - (1-np.exp(-simulation_time))

fig, ax = plt.subplots(1, 1, figsize = (5, 4))
ax.plot(simulation_time, th1_cells, "tab:blue", linestyle = "--", label = r"$\alpha=1$")
ax.plot(simulation_time, th2_cells, "tab:red", linestyle = "-", label = r"$\alpha=10$")
ax.set_xlabel("regular simulation time")
plt.tight_layout()


fig, ax = plt.subplots(1, 1, figsize = (5, 4))
ax.plot(simulation_time2, th1_cells, "tab:blue", linestyle = "--", label = r"$\alpha=1$")
ax.plot(simulation_time2, th2_cells, "tab:red", linestyle = "-", label = r"$\alpha=10$")
ax.set_xlabel("transformed simulation time")
plt.tight_layout()


fig, ax = plt.subplots(1, 2, figsize = (9, 4))
ax[0].plot(simulation_time, th1_cells, "tab:blue", linestyle = "--", label = r"$\alpha=1$")
ax[0].plot(simulation_time, th2_cells, "tab:red", linestyle = "-", label = r"$\alpha=10$")
ax[0].set_xlabel("regular simulation time")

ax[1].plot(simulation_time2, th1_cells, "tab:blue", linestyle = "--", label = r"$\alpha=1$")
ax[1].plot(simulation_time2, th2_cells, "tab:red", linestyle = "-", label = r"$\alpha=10$")
ax[1].set_xlabel("transformed simulation time")

for a in ax:
    a.legend()
    a.set_ylabel("Th cells")

plt.tight_layout()

"""
fig, ax = plt.subplots()
ax_time_course(state, 
               ax,
               simulation_time = simulation_time, 
               alpha_1 = alpha_1, 
               initial_cells = initial_cells, 
               linestyle = "--")

ax_time_course(state_rtm, 
               ax,
               simulation_time = simulation_time, 
               alpha_1 = 20, 
               initial_cells = initial_cells, 
               linestyle = "-")
ax.set_title("time course")

plt.tight_layout()

#==============================================================================
# IL12
#==============================================================================
il12_params = list(rate_params)
il12_arr = np.linspace(0,2,50)
il12_readouts = variable_effect(il12_arr, "IL12", *il12_params)
il12_cells = get_cell_arr(il12_readouts)

il12_params_rtm = list(rtm_params)
il12_arr_rtm = np.linspace(0,2,50)
il12_readouts_rtm = variable_effect(il12_arr_rtm, "IL12", *il12_params_rtm)
il12_cells_rtm = get_cell_arr(il12_readouts_rtm)

fig, ax = plt.subplots()
ax_var_effect(il12_arr, il12_cells, ax, linestyle = "--", plot_both = False)
ax_var_effect(il12_arr, il12_cells_rtm, ax, plot_both = False)
ax.set_xlabel("IL12")
ax.set_ylabel("%Th cells")
plt.tight_layout()

#==============================================================================
# chain
#==============================================================================
chain_params = list(rate_params)
chain_arr = np.arange(1,20,1)
chain_readouts = variable_effect(chain_arr, "chain", *chain_params)
chain_cells = get_cell_arr(chain_readouts)

fig, ax = plt.subplots(1,1, figsize = (5,4))
ax_chain(chain_arr, chain_cells, ax, plot_both = False, xsteps = 5)
ax.set_ylabel("% Th1 cells at final state")
ax.set_ylabel("%Th cells")
plt.tight_layout()
#fig.savefig(save_path+"chain_pos_fb.svg", bbox_inches = "tight")

#==============================================================================
# fb strength
#==============================================================================
fb_params = list(rate_params)
fb_arr = np.linspace(0,10,100)
fb_readouts = variable_effect(fb_arr, "feedback_strength_pos_Th1", *fb_params)
fb_cells = get_cell_arr(fb_readouts)

fb_params_rtm = list(rtm_params)
fb_arr_rtm = np.linspace(0,10,100)
fb_readouts_rtm = variable_effect(fb_arr_rtm, "feedback_strength_pos_Th1", *fb_params_rtm)
fb_cells_rtm = get_cell_arr(fb_readouts_rtm)


fig, ax = plt.subplots()
ax_var_effect(fb_arr, fb_cells, ax, linestyle = "--")
ax_var_effect(fb_arr_rtm, fb_cells_rtm, ax)
ax.set_xscale('log')
ax.set_xticks([1, 10])
ax.set_xlim(0,10)
#ax.set_ylim(0,100)

ax.set_xlabel("feedback strength")
ax.set_ylabel("%Th cells")
plt.tight_layout()
#==============================================================================
# cytokine rates
#==============================================================================
cyto_params = list(rate_params)
cyto_arr = np.linspace(0,10,100)
cyto_readouts = variable_effect(cyto_arr, "cytokine_rates", *cyto_params)
cyto_cells = get_cell_arr(cyto_readouts)

cyto_params_rtm = list(rtm_params)
cyto_arr_rtm = np.linspace(0,10,100)
cyto_readouts_rtm = variable_effect(cyto_arr_rtm, "cytokine_rates", *cyto_params_rtm)
cyto_cells_rtm = get_cell_arr(cyto_readouts_rtm)


fig, ax = plt.subplots()
ax_var_effect(cyto_arr, cyto_cells, ax, linestyle = "--")
ax_var_effect(cyto_arr_rtm, cyto_cells_rtm, ax)
ax.set_xlabel("rate ifng / il4")
ax.set_ylabel("%Th cells")
plt.tight_layout()
#==============================================================================
# probabilities
#==============================================================================
rate_cytokines = get_cyto_conc(state, simulation_time, alpha_1, rate_ifn, rate_il4, initial_cells)
probabilities = get_prob(rate_cytokines[0], rate_cytokines[1], conc_il12, hill_1, hill_2, K, fb_strength)

fig,ax = plt.subplots()
ax.plot(simulation_time, probabilities[0], "tab:blue")
ax.plot(simulation_time, probabilities[1], "tab:red")
ax.set_xlabel("time")
ax.set_ylabel("probability")
plt.tight_layout()
"""