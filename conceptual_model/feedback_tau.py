#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Tue Oct  2 17:32:22 2018

@author: burt
"""

import numpy as np
import os

os.chdir("/home/burt/Documents/code/th_cell_differentiation")

### run_model takes default parameters stored in th1_th2_parameters.py
from modules.th1_th2_ode_model_generic import run_model, il12, feedback_strength
from modules.th1_th2_plotting import plot_time_course, plot_il12, ax_time_course, ax_il12, ax_feedback_tau
import modules.th1_th2_conceptual_parameters
import matplotlib.pyplot as plt
import seaborn as sns
save_path = "/home/burt/Documents/tcell_project/figures/"
#==============================================================================
# import parameters
#==============================================================================
rtm_params = modules.th1_th2_conceptual_parameters.rtm_params
rate_params = modules.th1_th2_conceptual_parameters.rate_params
#==============================================================================
# run time for rtm and for rate process
#==============================================================================
simulation_name = ""
rate_sim = run_model(simulation_name, parameters = rate_params)
rtm_sim = run_model(simulation_name, parameters = rtm_params)

fig, ax = plt.subplots(1,1, figsize = (5,4))

ax_time_course(rate_sim, ax, rate_params)
ax_time_course(rtm_sim, ax, rtm_params, linestyle = "--")
#ax.legend([r"Th0","Th1","Th2"])
ax.set_xlabel("time")
plt.tight_layout()
sns.set_style("ticks")

#fig.savefig("th1_th2.svg", bbox_inches = "tight")
# plot time course
#plot_time_course(th1_th2_model, parameters)
#==============================================================================
# run time course for asymmetric case (change params in conc param file for that)
#==============================================================================
new_params1 = list(rtm_params)
new_params2 = list(rtm_params)
new_params1[0] = 1
new_params1[2] = 1.
new_params2[1] = 1
new_params2[3] = 1.

run1 = run_model(simulation_name, parameters = new_params1)
run2 = run_model(simulation_name, parameters = new_params2)

fig, ax = plt.subplots(1,3, figsize = (12,4))
ax = ax.flatten()
ax_time_course(rate_sim, ax[0], rate_params)
ax_time_course(rtm_sim, ax[0], rtm_params, linestyle = "--")
#ax.legend([r"Th0","Th1","Th2"])
ax_time_course(run1, ax[1], new_params1)
ax_time_course(run2, ax[2], new_params2)
#ax.legend([r"Th0","Th1","Th2"])
plt.tight_layout()
sns.set_style("ticks")
#fig.savefig(save_path+"th1_th2_mixed_model.svg", bbox_inches = "tight")

il12_conc = np.arange(0,2,0.01)

il12_rate_arr = il12(il12_conc, rate_params)
il12_rtm_arr = il12(il12_conc, rtm_params)

il12_conc_tau = np.arange(0,0.2,0.001)

fb_rate_arr = feedback_strength(il12_conc_tau, rate_params, cytokine_type = "IL12")
fb_rtm_arr = feedback_strength(il12_conc_tau, rtm_params, cytokine_type = "IL12")



fig, ax = plt.subplots(1,1, figsize = (5,4))
ax_il12(il12_rate_arr, ax)
ax_il12(il12_rtm_arr, ax, linestyle = "--")
ax.set_ylabel("%Th cells in steady state")
plt.tight_layout()

fig, ax = plt.subplots(1,1, figsize = (5,4))
ax_feedback_tau(fb_rate_arr, ax)
ax_feedback_tau(fb_rtm_arr, ax, linestyle = "--")
#ax.set_ylabel("%Th cells in steady state")
plt.tight_layout()