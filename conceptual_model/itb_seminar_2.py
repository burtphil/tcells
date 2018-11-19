#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Thu Nov 15 11:42:41 2018

@author: burt

feedback strength vs tau and steady state
"""
import numpy as np
import os

os.chdir("/home/burt/Documents/code/th_cell_differentiation")

### run_model takes default parameters stored in th1_th2_parameters.py
from modules.th1_th2_ode_model_generic import il12, feedback_strength
from modules.th1_th2_plotting import ax_il12, ax_feedback_tau
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
# perturbations
#==============================================================================
fig_name = "pos_fb"

il12_conc = np.arange(0,2,0.01)
il12_conc_tau = np.arange(0,0.01,0.001)

il12_rate_arr = feedback_strength(il12_conc_tau, rate_params, cytokine_type = "IFNG")
il12_rtm_arr = feedback_strength(il12_conc_tau, rtm_params, cytokine_type = "IFNG")


fig, ax = plt.subplots(1,1, figsize = (5,4))
ax_il12(il12_rate_arr, ax)
ax_il12(il12_rtm_arr, ax, linestyle = "--")
ax.set_ylabel("%Th cells in steady state")
ax.set_xlabel("IFNg [a.u.]")
plt.tight_layout()
sns.set_style("ticks")
#fig.savefig(save_path+fig_name+".svg", bbox_inches = "tight")

fb_rate_arr = feedback_strength(il12_conc_tau, rate_params, cytokine_type = "IFNG")
fb_rtm_arr = feedback_strength(il12_conc_tau, rtm_params, cytokine_type = "IFNG")


fig, ax = plt.subplots(1,1, figsize = (5,4))
ax_feedback_tau(fb_rate_arr, ax)
ax_feedback_tau(fb_rtm_arr, ax, linestyle = "--")
#ax.set_ylabel("%Th cells in steady state")
plt.tight_layout()
sns.set_style("ticks")
ax.legend(["Th1","Th2"])
#fig.savefig(save_path+fig_name+"_tau.svg", bbox_inches = "tight")