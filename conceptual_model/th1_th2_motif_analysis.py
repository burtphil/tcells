#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Wed Sep 19 14:42:32 2018

@author: burt
very basic model T can convert to T1 or T2 and I change the feedback
"""
import os

windows_path= "C:/Users/Philipp/Documents/tcells/"
linux_path= "/home/burt/Documents/code/th_cell_differentiation"
os.chdir(linux_path)

### run_model takes default parameters stored in th1_th2_parameters.py
from modules.th1_th2_ode_model_generic import chain, il12, run_model, chain_one
from modules.th1_th2_plotting import plot_chain, ax_chain, subplot_tau, plot_il12,plot_time_course, ax_time_course, ax_il12
import matplotlib.pyplot as plt
import numpy as np
import modules.th1_th2_conceptual_parameters as cparams
import matplotlib.gridspec as gridspec

save_path = "/home/burt/Documents/figures/"

#==============================================================================
# model parameters
#==============================================================================

conceptual_params = cparams.parameters

#==============================================================================
# case: pos feedback on Th1
#==============================================================================
params = list(conceptual_params)

# define hill parameters for pos feedback on Th1

hill_1 = [1,0,1]
hill_2 = [0,0,0]

params[6] = hill_1
params[7] = hill_2

chain_length = 20
#==============================================================================
# run time course and IL12 effect simulation
#==============================================================================
alpha_1_params = list(params)
alpha_2_params = list(params)
il12_alpha_1_params = list(params)
il12_alpha_2_params = list(params)
chain_both_params = list(params)
chain_th1_params = list(params)
chain_th2_params = list(params)

alpha_1_params[0],il12_alpha_1_params[0] = 1,1
alpha_1_params[1],il12_alpha_1_params[1] = 1,1

alpha_2_params[0], il12_alpha_2_params[1] = 2,2
alpha_2_params[1], il12_alpha_2_params[1] = 2,2

simulation_name = "sim"
alpha_1 = run_model(simulation_name, parameters = alpha_1_params)

simulation_name = "sim"
alpha_2 = run_model(simulation_name, parameters = alpha_2_params)

il12_conc = np.linspace(0,2,100)

#==============================================================================
# plotting
#==============================================================================
"""
fig = plt.figure(figsize = (9,7))
gs = gridspec.GridSpec(ncols=2, nrows=2)
ax1 = fig.add_subplot(gs[0, 0])
ax2 = fig.add_subplot(gs[0, 1])
ax3 = fig.add_subplot(gs[1, 0])
ax4 = fig.add_subplot(gs[1, 1])

ax_time_course(state = alpha_1, ax = ax1, parameters = alpha_1_params)
ax_time_course(state = alpha_2, ax = ax2, parameters = alpha_2_params)
ax_il12(il12(il12_conc, parameters = il12_alpha_1_params), ax = ax3)
ax_il12(il12(il12_conc, parameters = il12_alpha_2_params), ax = ax4)
ax1.set_title(r"$\alpha=1$")
ax3.set_title(r"$\alpha=1$")
ax2.set_title(r"$\alpha=2$")
ax4.set_title(r"$\alpha=2$")
ax2.legend(["Th0","Th1","Th2"])
plt.tight_layout()
"""
fig, ax = plt.subplots(1,3, figsize = (12,4))
ax_chain(chain_array = chain(chain_length = chain_length, parameters = chain_both_params), ax = ax[0])
ax[0].set_title(r"$\alpha_{Th1}=\alpha_{Th2}$")
ax[0].legend(["Th1","Th2"])
alpha_idx_th1 = [0,1]
ax_chain(chain_array = chain_one(chain_length = chain_length, alpha_idx = alpha_idx_th1, parameters = chain_th1_params), ax = ax[1])
ax[1].set_title(r"$\alpha_2=$"+str(alpha_idx_th1[1]))

alpha_idx_th2 = [1,1]
ax_chain(chain_array = chain_one(chain_length = chain_length, alpha_idx = alpha_idx_th2, parameters = chain_th2_params), ax = ax[2])
ax[2].set_title(r"$\alpha_1=$"+str(alpha_idx_th1[1]))

plt.tight_layout()

fig, ax = plt.subplots(1,3, figsize = (12,4))
ax_chain(chain_array = chain(chain_length = chain_length, parameters = chain_both_params), ax = ax[0])
ax[0].set_title(r"$\alpha_{Th1}=\alpha_{Th2}$")
ax[0].legend(["Th1","Th2"])
alpha_idx_th1 = [0,2]
ax_chain(chain_array = chain_one(chain_length = chain_length, alpha_idx = alpha_idx_th1, parameters = chain_th1_params), ax = ax[1])
ax[1].set_title(r"$\alpha_{Th2}=$"+str(alpha_idx_th1[1]))

alpha_idx_th2 = [1,2]
ax_chain(chain_array = chain_one(chain_length = chain_length, alpha_idx = alpha_idx_th2, parameters = chain_th2_params), ax = ax[2])
ax[2].set_title(r"$\alpha_{Th1}=$"+str(alpha_idx_th1[1]))

plt.tight_layout()
#plot_chain(chain(chain_length, chain_params))

#plot_il12(il12(il_12, il_12_params), xlabel ="IL12")

#==============================================================================
# chain effect constant IL12 different means
#==============================================================================
"""
chain_params = list(params)
#plot_chain(chain(25.,chain_params))

mean_diff = [0.5,0.75,1.0,1.25,1.5]
mean_diff_params = list(params)

fig, axes = plt.subplots(1,len(mean_diff), figsize = (15,4))

for i, ax_chain in zip(mean_diff, axes):
    mean_diff_params[3] = mean_diff_params[1]/i
    chain_array = chain(chain_length, mean_diff_params, stepsize = cparams.stepsize)
    subplot_tau(chain_array, ax_chain)
    ax_chain.set_title("mean diff ="+str(i))
axes[0].set_ylabel(r"$\tau_{1/2}$")
axes[2].set_xlabel(r"chain length $\alpha$")
plt.tight_layout()
#fig.savefig(save_path+"mean_diff_tau_il12_1.svg", bbox_inches = "tight", dpi = 1200)

#==============================================================================
#  chain effect constant mean different IL12
#==============================================================================
"""