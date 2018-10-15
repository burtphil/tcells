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
import matplotlib.image as mpimg
### run_model takes default parameters stored in th1_th2_parameters.py
from modules.th1_th2_ode_model_generic import chain, il12, run_model, chain_one, chain_th1, chain_th2
from modules.th1_th2_plotting import plot_chain, ax_chain, ax_tau, plot_il12,plot_time_course, ax_time_course, ax_il12
import matplotlib.pyplot as plt
import numpy as np
import modules.th1_th2_conceptual_parameters as cparams

def p_th_diff(conc_ifn,conc_il4,conc_il12, hill, half_saturation, strength = 1.0):
    """
    returns probability of Th1 Th2 differentiation for given cytokine concentrations
    kinetics are hill-like so hill coefficients for cytokines also need to be provided
    """
    
    # watch out, this is tricky because if concentrations are zero and negative feedback (hill coeff < 0)
    # then you divide by zero, best avoid zero concentrations through base cytokine rate
    ifn_prob = (conc_ifn/half_saturation[0])**hill[0]
    il4_prob = (conc_il4/half_saturation[1])**hill[1]
    il12_prob = (conc_il12/half_saturation[2])**hill[2]

    prob_th_diff = strength*ifn_prob*il4_prob*il12_prob
   
    return prob_th_diff

#==============================================================================
# model parameters
#==============================================================================
chain_length = 10

# hill coefficients for th1 cells

fig_names = ["pos_feedback",
            "neg_feedback",
            "pos_neg_feedback",
            "neg_th1_neg_pos_th2",
            "pos_th1_neg_pos_th2",
            "pos_th1_neg_th2",
            "pos_th2_neg_th2"]

hill_1 = [[1,0,1],
          [0,0,1],
          [1,-1,1],
          [0,-1,1],
          [1,0,1],
          [1,0,1],
          [0,0,1]]

hill_2 = [[0,0,0],
          [-1,0,0],
          [-1,1,0],
          [-1,1,0],
          [-1,1,0],
          [-1,0,0],
          [-1,1,0]]

for coeff_1, coeff_2, fig_name in zip(hill_1,hill_2,fig_names):
    #==============================================================================
    # parameter settings
    #==============================================================================
    params = list(cparams.parameters)
    params[6]=coeff_1
    # hill coefficients for th2 cells
    params[7]=coeff_2
    
    chain_both_params = list(params)
    chain_th1_params = list(params)
    chain_th2_params = list(params)
    alpha_1_params = list(params)
    
    tau_both_params = list(params)
    tau_th1_params = list(params)
    tau_th2_params = list(params)
    
    alpha1 = 9
    alpha2 = 1
    alpha_1_params[0],alpha_1_params[2]=alpha1,alpha1
    alpha_1_params[1],alpha_1_params[3]=alpha2,alpha2
    
    simulation_name = "sim"
    alpha_1 = run_model(simulation_name, parameters = alpha_1_params)
    
    #==============================================================================
    # plotting
    #==============================================================================
    fig, ax = plt.subplots(3,3, figsize = (12,10))
    
    ax = ax.flatten()
    img=mpimg.imread('images/'+fig_name+".png")
    ax[0].imshow(img)
    ax[0].axis('off')
    
    ax_time_course(state = alpha_1, ax = ax[1], parameters = alpha_1_params)
    ax[1].set_title(r"$\alpha_{Th1}=$"+str(alpha_1_params[0])+r", $\alpha_{Th2}=$"+str(alpha_1_params[1]))
    ax[1].legend(["Thn","Th1","Th2"])
    ax[1].set_xticks([0,cparams.stop/2.,cparams.stop])
    ### calculate probabilities based on state array (output from odeint)
    
    state = alpha_1
    half_saturation = [1.,1.,1.]
    conc_il12 = 1.
    conc_ifn = cparams.rate_ifn*state[:,alpha1+1]+half_saturation[0]
    conc_il4 = cparams.rate_il4*state[:,-1]+half_saturation[1]
    ### calculate initial th1 and th2 populations from naive cells based on branching probabilities
    th_0 = state[0]
    # branching probablities
    
    prob_th1 = p_th_diff(conc_ifn,conc_il4,conc_il12, params[6], half_saturation)
    prob_th2 = p_th_diff(conc_ifn,conc_il4,conc_il12, params[7], half_saturation)    
    # normalized branching probabilities
    p_1 = prob_th1/(prob_th1+prob_th2)
    p_2 = prob_th2/(prob_th1+prob_th2)
    ax[2].plot(cparams.simulation_time, p_1)
    ax[2].plot(cparams.simulation_time, p_2, "r")
    ax[2].legend(["$p_{Th1}$","$p_{Th2}$"])
    ax[2].set_xlabel("time")
    ax[2].set_ylabel("probability")
    ax[2].set_xticks([0,cparams.stop/2.,cparams.stop])
    ax[2].set_yticks([0,0.5,1])
    ax[2].set_ylim([0,1])
    ax_chain(chain_array = chain(chain_length = chain_length, parameters = chain_both_params), ax = ax[3])
    ax[3].set_title(r"$\alpha_{Th1}=\alpha_{Th2}$")
    ax[3].legend(["Th1","Th2"])
    # vary alpha1 (Th1), set alpha Th2 to 1
    alpha_idx_th2 = [1,1]
    ax_chain(chain_array = chain_th2(chain_length = chain_length, parameters = chain_th2_params), ax = ax[4])
    ax[4].set_title(r"$\alpha_{Th1}=$"+str(alpha_idx_th2[1]))
    
    alpha_idx_th1 = [0,1]
    ax_chain(chain_array = chain_th1(chain_length = chain_length, parameters = chain_th1_params), ax = ax[5])
    ax[5].set_title(r"$\alpha_{Th2}=$"+str(alpha_idx_th1[1]))
    
    for i in range(3,6):
        ax[i].set_ylabel("% Th cells after "+str(cparams.stop)+" hrs")
    
    ax_tau(chain_array = chain(chain_length = chain_length, parameters = tau_both_params), ax = ax[6])
    ax_tau(chain_array = chain_th2(chain_length = chain_length, parameters = tau_th2_params), ax = ax[7])
    ax_tau(chain_array = chain_th1(chain_length = chain_length, parameters = tau_th1_params), ax = ax[8])
    
    for i in range(6,9):
        ax[i].set_ylim([0,3])
        ax[i].set_yticks([0,1.5,3])
    
    plt.tight_layout()
    
    path = "/home/burt/Documents/tcell_project/figures/model_simulations/conceptual_model/motif_analysis/new_motifs/"
    fig.savefig(path+fig_name+"_motif.pdf", bbox_inches = "tight", dpi = 1200)
    plt.close(fig)
"""
# vary alpha1 (Th2), set alpha Th1 to 1

# chain effect
fig, ax = plt.subplots(1,3, figsize = (12,4))
ax_time_course(state = alpha_1, ax = ax[0], parameters = alpha_1_params)
ax[0].set_title("alpha1=2, alpha2=1")
ax[0].legend(["Th0","Th1","Th2"])
ax_time_course(state = alpha_2, ax = ax[1], parameters = alpha_2_params)
ax[1].set_title("alpha1=1, alpha2=2")
alpha_idx_th2 = [1,1]
ax_chain(chain_array = chain_th2(chain_length = chain_length, parameters = chain_th2_params), ax = ax[2])
ax[2].set_title(r"$\alpha_{Th1}=$"+str(alpha_idx_th2[1]))

plt.tight_layout()

#==============================================================================
# chain effect constant IL12 different means
#==============================================================================

chain_params = list(params)
#plot_chain(chain(25.,chain_params))

mean_diff = [0.5,0.75,1.0,1.25,1.5]
mean_diff_params = list(params)

fig, axes = plt.subplots(1,len(mean_diff), figsize = (15,4))

for i, ax in zip(mean_diff, axes):
    mean_diff_params[3] = mean_diff_params[1]/i
    chain_array = chain(chain_length, mean_diff_params, stepsize = cparams.stepsize)
    subplot_tau(chain_array, ax)
    ax.set_title("mean diff ="+str(i))
axes[0].set_ylabel(r"$\tau_{1/2}$")
axes[2].set_xlabel(r"chain length $\alpha$")
plt.tight_layout()
#fig.savefig(save_path+"mean_diff_tau_il12_1.svg", bbox_inches = "tight", dpi = 1200)
"""
#==============================================================================
#  chain effect constant mean different IL12
#==============================================================================
#==============================================================================
# save output
#==============================================================================

#save_path = "/home/burt/Documents/tcell_project/figures/model_simulations/conceptual_model/motif_analysis/"
#pdf = matplotlib.backends.backend_pdf.PdfPages(save_path+motif_name+".pdf")
#for fig in xrange(1, plt.gcf().number + 1):
#    pdf.savefig(fig)
    #plt.close(fig)
#pdf.close()
