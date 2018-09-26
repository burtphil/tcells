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
from modules.th1_th2_ode_model_generic import chain, il12, run_model, chain_one, chain_th1, chain_th2
from modules.th1_th2_plotting import plot_chain, ax_chain, subplot_tau, plot_il12,plot_time_course, ax_time_course, ax_il12
import matplotlib.pyplot as plt
import numpy as np
import modules.th1_th2_conceptual_parameters as cparams
import matplotlib.gridspec as gridspec
import matplotlib.backends.backend_pdf

#==============================================================================
# model parameters
#==============================================================================
#motif_name = "neg_pos_feedback_th1"
conceptual_params = cparams.parameters

chain_length = 10

hill1 = [[0,0,1],[1,0,1],[0,-1,1],[1,-1,1]]
hill2 = [[0,0,0],[0,1,0],[-1,0,0],[-1,1,0]]

hill_coefficients = [[x,y] for x in hill1 for y in hill2]

fb_th1 = ["no_th1","pos_th1","neg_th1","pos_neg_th1"]
fb_th2 = ["no_th2","pos_th2","neg_th2","pos_neg_th2"]
motif_names = [x+"_"+y for x in fb_th1 for y in fb_th2]



#==============================================================================
# control run
#==============================================================================
#simulation_name = "sim"
#state = run_model(simulation_name, parameters = params)
#plot_time_course(state,params)
    
#==============================================================================
# case: pos feedback on Th1
#==============================================================================
for hill_coeff, motif_name in zip(hill_coefficients, motif_names):
    params = list(conceptual_params)
    params[6]=hill_coeff[0]
    params[7]=hill_coeff[1]
    #==============================================================================
    # parameter settings
    #==============================================================================
    
    alpha_1_params = list(params)
    alpha_2_params = list(params)
    il12_alpha_1_params = list(params)
    il12_alpha_2_params = list(params)
    chain_both_params = list(params)
    chain_th1_params = list(params)
    chain_th2_params = list(params)
    
    # settings with alpha 1 = 2 and alpha 2 = 2
    
    # replace original parameters alpha1,alpha2,beta1,beta2 with above defined settings
    alpha_1_params[0], alpha_1_params[1],alpha_1_params[2],alpha_1_params[3] = 3,1,3,1
    alpha_2_params[0], alpha_2_params[1],alpha_2_params[2],alpha_2_params[3] = 1,3,1,3
    
    il12_alpha_1_params = list(alpha_1_params)
    il12_alpha_2_params = list(alpha_2_params)
    
    # il12 settings
    il12_conc = np.linspace(0,2,50)
    
    #==============================================================================
    # run alpha 1 and alpha 2 simulations
    #==============================================================================
    simulation_name = "sim"
    alpha_1 = run_model(simulation_name, parameters = alpha_1_params)
    
    simulation_name = "sim"
    alpha_2 = run_model(simulation_name, parameters = alpha_2_params)
    
    
    
    #==============================================================================
    # plotting
    #==============================================================================
    ### time course and IL12 effect
    
    fig4 = plt.figure(figsize = (9,7))
    gs = gridspec.GridSpec(ncols=2, nrows=2)
    ax1 = fig4.add_subplot(gs[0, 0])
    ax2 = fig4.add_subplot(gs[0, 1])
    ax3 = fig4.add_subplot(gs[1, 0])
    ax4 = fig4.add_subplot(gs[1, 1])
    
    ax_time_course(state = alpha_1, ax = ax1, parameters = alpha_1_params)
    ax_time_course(state = alpha_2, ax = ax2, parameters = alpha_2_params)
    ax_il12(il12(il12_conc, parameters = il12_alpha_1_params), ax = ax3)
    ax_il12(il12(il12_conc, parameters = il12_alpha_2_params), ax = ax4)
    ax2.set_title(r"$\alpha_1=1, \alpha_2=3$")
    ax4.set_title(r"$\alpha_1=1, \alpha_2=3$")
    ax1.set_title(r"$\alpha_1=3, \alpha_2=1$")
    ax3.set_title(r"$\alpha_1=3, \alpha_2=1$")
    ax2.legend(["Th0","Th1","Th2"])
    plt.tight_layout()
    
    # chain effect
    
    fig, ax = plt.subplots(1,3, figsize = (12,4))
    ax_chain(chain_array = chain(chain_length = chain_length, parameters = chain_both_params), ax = ax[0])
    ax[0].set_title(r"$\alpha_{Th1}=\alpha_{Th2}$")
    ax[0].legend(["Th1","Th2"])
    # vary alpha1 (Th1), set alpha Th2 to 1
    
    alpha_idx_th2 = [1,1]
    ax_chain(chain_array = chain_th2(chain_length = chain_length, parameters = chain_th2_params), ax = ax[1])
    ax[1].set_title(r"$\alpha_{Th1}=$"+str(alpha_idx_th2[1]))
    
    alpha_idx_th1 = [0,1]
    ax_chain(chain_array = chain_th1(chain_length = chain_length, parameters = chain_th1_params), ax = ax[2])
    ax[2].set_title(r"$\alpha_{Th2}=$"+str(alpha_idx_th1[1]))
    plt.tight_layout()
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
    """
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
    
    #==============================================================================
    #  chain effect constant mean different IL12
    #==============================================================================
    #==============================================================================
    # save output
    #==============================================================================
    
    save_path = "/home/burt/Documents/tcell_project/figures/model_simulations/conceptual_model/motif_analysis/automated/"
    pdf = matplotlib.backends.backend_pdf.PdfPages(save_path+motif_name+".pdf")
    for fig in xrange(1, plt.gcf().number + 1):
        pdf.savefig(fig)
        plt.close(fig)
    pdf.close()
