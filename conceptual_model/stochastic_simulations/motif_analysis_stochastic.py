#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Tue Nov 13 10:09:58 2018

@author: burt
"""
import os
linux_path= "/home/burt/Documents/code/th_cell_differentiation"
os.chdir(linux_path)

import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
import matplotlib.image as mpimg
import modules.th1_th2_conceptual_parameters as cparams
from modules.th1_th2_ode_model_generic import chain, chain_th1, chain_th2
from modules.th1_th2_plotting import ax_chain

def plt_chain_stochastic(chain_arr, ax, alpha = 1, label = False):
    th1_endstates, th2_endstates = chain_arr
    chain_arr = range(1,len(th1_endstates)+1)
    
    if label == False:
        ax.scatter(chain_arr, th1_endstates, c = "tab:blue", alpha = alpha)
        ax.scatter(chain_arr, th2_endstates, c = "tab:red", alpha = alpha)
    else:
        ax.scatter(chain_arr, th1_endstates, c = "tab:blue", alpha = alpha, label = "Th1")
        ax.scatter(chain_arr, th2_endstates, c = "tab:red", alpha = alpha, label = "Th2")
        
    ax.set_ylim([0,1])

def update_params(parameters, hill_1, hill_2):
    parameters = list(parameters)
    parameters[6] = hill_1
    parameters[7] = hill_2
    
    return parameters

stored_simulations_path = "/home/burt/Documents/tcell_project/saved_simulations/"
image_path = "/home/burt/Documents/code/th_cell_differentiation/images/"
save_image_path = "/home/burt/Documents/tcell_project/figures/"

chain_length = 10

hill_coeff_1 = [
    [1,0,1],
    [0,0,1],
    [1,-1,1],
    [0,-1,1],
    [1,0,1],
    [1,0,1],
    [0,0,1],
    ]

hill_coeff_2 = [
    [0,0,0],
    [-1,0,0],
    [-1,1,0],
    [-1,1,0],
    [-1,1,0],
    [-1,0,0],
    [-1,1,0],
    ]

fig_names = [
    "pos_feedback",
    "neg_feedback",
    "pos_neg_feedback",
    "neg_th1_neg_pos_th2",
    "pos_th1_neg_pos_th2",
    "pos_th1_neg_th2",
    "pos_th2_neg_th2",
    ]

chain_array = [
    "chain_arr_both",
    "chain_arr_th2",
    "chain_arr_th1",
    ]

title_array = [
    r"$\alpha_{Th1}=\alpha_{Th2}$",
    r"$\alpha_{Th1}=1$",
    r"$\alpha_{Th2}=1$",
    ]

# load array with stochastic simulations
stochastic_chain_simulations = [np.load(stored_simulations_path+fig_name+"_stoc.npz") for fig_name in fig_names]

# generate deterministic simulations
params = cparams.parameters

params_list = [update_params(params, hill_1, hill_2) for hill_1, hill_2 in zip(hill_coeff_1, hill_coeff_2)]
chain_both_params = list(params_list)
chain_th1_params = list(params_list)
chain_th2_params = list(params_list)

#chain_arr_both = [chain(chain_length, param) for param in chain_both_params]
#chain_arr_th1 = [chain_th1(chain_length, param) for param in chain_th1_params]
#chain_arr_th2 = [chain_th2(chain_length, param) for param in chain_th2_params]

for fig_name, stoc_sim in zip(fig_names, stochastic_chain_simulations):
    fig, axes = plt.subplots(1,4, figsize = (16,4))
    img=mpimg.imread(image_path+fig_name+".png")
    axes[0].imshow(img)
    axes[0].axis('off')
    
    for idx in range(1,4):
        arr = stoc_sim[chain_array[idx-1]]
        ax = axes[idx]
        ax.set_xlabel(r"chain length $\alpha$")
        ax.set_ylabel("% Th cells")
        ax.set_title(title_array[idx-1])
        
        for idy,sim in enumerate(arr):
            if idy == 0:
                plt_chain_stochastic(sim, ax, alpha = 0.5, label = True)
            else:
                plt_chain_stochastic(sim, ax, alpha = 0.5)
    
    plt.legend()
    plt.tight_layout()
    
    #fig.savefig(save_image_path+fig_name+"_stoc.pdf", bbox_inches = "tight", dpi = 1200)   

"""
# use when you want to plot all in one figure
for axes, fig_name in zip(all_axes, fig_names):
    img=mpimg.imread(image_path+fig_name+".png")
    axes[0].imshow(img)
    axes[0].axis('off')
    axes[0].set_title(fig_name)
    
    for stoc_sim in stochastic_chain_simulations:    
        for idx in range(1,4):
            arr = stoc_sim[chain_array[idx-1]]
            ax = axes[idx]
            ax.set_xlabel(r"chain length $\alpha$")
            ax.set_ylabel("% Th cells")
            ax.set_title(chain_array[idx-1])
            
            for idy,sim in enumerate(arr):
                if idy == 0:
                    plt_chain_stochastic(sim, ax, alpha = 0.1, label = True)
                else:
                    plt_chain_stochastic(sim, ax, alpha = 0.1)
        
plt.legend()
plt.tight_layout()   
"""