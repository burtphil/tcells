#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Fri Sep 14 11:27:20 2018

@author: burt
plt chain effects
"""
import numpy as np
import os

os.chdir("/home/burt/Documents/scripts/th_cell_differentiation")

### run_model takes default parameters stored in th1_th2_parameters.py
from th1_th2_ode_model_generic import chain, il12
from th1_th2_plotting import plot_chain, subplot_chain, subplot_tau, plot_il12
import th1_th2_conceptual_parameters as cparams
import matplotlib.pyplot as plt
save_path = "/home/burt/Documents/figures/"
#==============================================================================
# model parameters
#==============================================================================

conceptual_params = cparams.parameters
stepsize = cparams.stepsize
chain_length = 25

#==============================================================================
# IL12 effect
#==============================================================================
mean_diff = [0.8,0.9,1.0,1.1,1.2]
mean_diff_params = list(conceptual_params)

fig, axes = plt.subplots(1,len(mean_diff), figsize = (15,4))

for i, ax_chain in zip(mean_diff, axes):
    mean_diff_params[3] = mean_diff_params[1]/i
    chain_array = chain(chain_length, mean_diff_params, stepsize = stepsize)
    subplot_chain(chain_array, ax_chain)
    ax_chain.set_title("mean diff ="+str(i))
axes[0].set_ylabel(r"$\tau_{1/2}$")
axes[2].set_xlabel(r"chain length $\alpha$")
plt.tight_layout()
#fig.savefig(save_path+"mean_diff_tau_il12_1.svg", bbox_inches = "tight", dpi = 1200)


"""
il12_concentrations = np.linspace(0,2,100)
il12_params = list(conceptual_params)
plot_il12(il12(il12_concentrations, il12_params), factor_x_axis = 1, xlabel = "IL12 [a.u.]", save = "il12_conc_model_mean_diff_1")

#==============================================================================
# chain effects
#==============================================================================
# plot chain effect
chain_params = list(conceptual_params)
plot_chain(chain(chain_length, parameters = chain_params, stepsize = stepsize), save = "chain_conc_model_il12_1_mean_diff_1")


# plot effect of IL12 for fixed mean differene
# set mean difference factor to 1.5!


il_12 = [0.5,0.75,1.0,1.25,1.5]

il_12_params = list(conceptual_params)

fig, axes = plt.subplots(1,len(il_12), figsize = (15,4))

for i, ax_chain in zip(il_12, axes):
    il_12_params[5] = i
    chain_array = chain(chain_length, il_12_params, stepsize = stepsize)
    subplot_tau(chain_array, ax_chain)
    ax_chain.set_title("IL12 ="+str(i))
axes[0].set_ylabel(r"$\tau_{1/2}$")
axes[2].set_xlabel(r"chain length $\alpha$")
plt.tight_layout()
#fig.savefig(save_path+"il12_tau_mean_diff_1.svg", bbox_inches = "tight", dpi = 1200)


mean_diff = [0.5,0.75,1.0,1.25,1.5]
mean_diff_params = list(conceptual_params)

fig, axes = plt.subplots(1,len(mean_diff), figsize = (15,4))

for i, ax_chain in zip(mean_diff, axes):
    mean_diff_params[3] = mean_diff_params[1]/i
    chain_array = chain(chain_length, mean_diff_params, stepsize = stepsize)
    subplot_tau(chain_array, ax_chain)
    ax_chain.set_title("mean diff ="+str(i))
axes[0].set_ylabel(r"$\tau_{1/2}$")
axes[2].set_xlabel(r"chain length $\alpha$")
plt.tight_layout()
#fig.savefig(save_path+"mean_diff_tau_il12_1.svg", bbox_inches = "tight", dpi = 1200)

mean_diff = [0.5,0.75,1.0,1.25,1.5]
mean_diff_params = list(conceptual_params)


fig, axes = plt.subplots(1,len(mean_diff), figsize = (15,4))

for i, ax_chain in zip(mean_diff, axes):
    mean_diff_params[3] = mean_diff_params[1]/i
    chain_array = chain(chain_length, mean_diff_params, stepsize = stepsize)
    subplot_chain(chain_array, ax_chain)
    ax_chain.set_title("mean diff ="+str(i))
axes[0].set_ylabel("% Th cells after 100 hrs")
axes[2].set_xlabel(r"chain length $\alpha$")
plt.tight_layout()
#fig.savefig(save_path+"mean_diff_chain_length_il12_1.svg", bbox_inches = "tight", dpi = 1200)

il_12 = [0.5,0.75,1.0,1.25,1.5]
il_12_params = list(conceptual_params)

fig, axes = plt.subplots(1,len(il_12), figsize = (15,4))

for i, ax_chain in zip(il_12, axes):
    il_12_params[5] = i
    chain_array = chain(chain_length, il_12_params, stepsize = stepsize)
    subplot_chain(chain_array, ax_chain)
    ax_chain.set_title("IL12 ="+str(i))
axes[0].set_ylabel("% Th cells after 100 hrs")
axes[2].set_xlabel(r"chain length $\alpha$")
plt.tight_layout()
#fig.savefig(save_path+"il12_chain_length_mean_diff_1.svg", bbox_inches = "tight", dpi = 1200)


il_12 = [0.2,0.3,0.4,0.5,0.7]
il_12_params = list(conceptual_params)
il_12_params[3] = il_12_params[1]/1.5

fig, axes = plt.subplots(1,len(il_12), figsize = (15,4))

for i, ax_chain in zip(il_12, axes):
    il_12_params[5] = i
    chain_array = chain(chain_length, il_12_params, stepsize = stepsize)
    subplot_chain(chain_array, ax_chain)
    ax_chain.set_title("IL12 ="+str(i))
axes[0].set_ylabel("% Th cells after 100 hrs")
axes[2].set_xlabel(r"chain length $\alpha$")
plt.tight_layout()
#fig.savefig(save_path+"il12_chain_length_mean_diff_1.5.pdf", bbox_inches = "tight", dpi = 1200)


mean_diff = [0.5,0.75,1,1.25,1.5]
il_12 = mean_diff
params = list(conceptual_params)

fig, axes = plt.subplots(len(il_12), len(mean_diff), figsize = (20,20))
for i,x in enumerate(mean_diff):
    ax_mean_diff = axes[i]
    params[3] = params[1]/x
    for j, y in enumerate(il_12):
        params[5] = y
        chain_array = chain(chain_length, params, stepsize = stepsize)
        ax_il12 = ax_mean_diff[j]
        subplot_chain(chain_array, ax_il12, labels = "off")
        #ax_il12.set_title("IL12="+str(y)+", mean diff="+str(x), loc = "right")
        if i == 4:
            ax_il12.set_xticks([0, 12.5,25])
            if j == 2:
                ax_il12.set_xlabel(r"chain length $\alpha$")
            
        if j == 0:
            ax_il12.set_yticks([0, 50, 100])
            if i == 2:
                ax_il12.set_ylabel("% Th cells after 100 hrs")
            
plt.tight_layout()
#fig. savefig(save_path+"il12_mean_diff_chain_effect.svg", bbox_inches = "tight", dpi=1200)


mean_diff = [0.5,0.75,1,1.25,1.5]
il_12 = mean_diff
params = list(conceptual_params)

fig, axes = plt.subplots(len(il_12), len(mean_diff), figsize = (20,20))
for i,x in enumerate(mean_diff):
    ax_mean_diff = axes[i]
    params[3] = params[1]/x
    for j, y in enumerate(il_12):
        params[5] = y
        chain_array = chain(chain_length, params, stepsize = stepsize)
        ax_il12 = ax_mean_diff[j]
        subplot_tau(chain_array, ax_il12)
        #ax_il12.set_title("IL12="+str(y)+", mean diff="+str(x), loc = "right")
        if i == 4:
            #ax_il12.set_xticks([0, 12.5,25])
            if j == 2:
                ax_il12.set_xlabel(r"chain length $\alpha$")
            
        if j == 0:
            #ax_il12.set_yticks([0, 20, 40])
            if i == 2:
                ax_il12.set_ylabel(r"$\tau_{1/2}$")
            
plt.tight_layout()
fig. savefig(save_path+"il12_mean_diff_tau_effect.pdf", bbox_inches = "tight", dpi=1200)
"""