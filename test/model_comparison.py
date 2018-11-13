#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Fri Oct 19 14:45:18 2018

@author: burt
compare stochastic and ode model simulation
"""
import numpy as np
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
import os
os.chdir("/home/burt/Documents/code/th_cell_differentiation")
import modules.th1_th2_conceptual_parameters as cparams
from modules.th1_th2_ode_model_generic import run_model
from stochasitic_simulation import run_stochastic_simulation
from modules.th1_th2_plotting import ax_time_course
#==============================================================================
# import parameters
#==============================================================================
model_params = cparams.parameters

#==============================================================================
# run time course simulation
#==============================================================================
simulation_name = "sim"

th1_th2_model = run_model(simulation_name, parameters = model_params)

test_simulation = np.load(simulation_name+".npz")
state = test_simulation["state"]
parameters = test_simulation["parameters"]

   

# plot time course
#plot_time_course(th1_th2_model, parameters)

#==============================================================================
# perform and plot multiple simulations
#==============================================================================
def get_cells(cells, time):
    all_cells = cells[:,:,0]
        
    # make some dummy lists
    naive_cells = []
    th1_cells = []
    th2_cells = []
    
    # for each time step, check how many of the cells belong to each subtype
    for t in range(len(time)):
        x = all_cells[:,t]
        naive_cells.append(len(x[x==cparams.thn_idx]))
        th1_cells.append(len(x[x==cparams.th1_idx]))
        th2_cells.append(len(x[x==cparams.th2_idx]))
    
    return [naive_cells,th1_cells,th2_cells]

def plt_stochastic_model():
    """
    deprecated
    """
    fig, ax = plt.subplots(1,2,figsize = (12,5))
    
    for i in range(1):
        cells, time = run_stochastic_simulation(start = cparams.start, stop = cparams.stop, nsteps = cparams.nsteps, ncells = cparams.ncells)
        
        # get all cells for each time step but only the cell state (not probability or tau)
        
        naive_cells,th1_cells,th2_cells = get_cells(cells,time)
        
        if i == 0:
            ax[0].plot(time,naive_cells, "k", label = "Thn")
            ax[0].plot(time,th1_cells, "tab:blue", label = "Th1")
            ax[0].plot(time,th2_cells, "r", label = "Th2")
        else:
            ax[0].plot(time,naive_cells, "k")
            ax[0].plot(time,th1_cells, "tab:blue")
            ax[0].plot(time,th2_cells, "r")
    
    ax[0].set_xlabel("time")
    ax[0].set_ylabel("cells")
    ax[0].legend()
    ax[0].set_xlim([cparams.start,cparams.stop])
    ax[0].set_ylim([0,cparams.initial_cells])
    ax[0].set_yticks([0,int(cparams.initial_cells/2),cparams.initial_cells])
    #ax[0].set_title("stochastic simulation")
    
    
    ax_time_course(state = state, ax = ax[1], parameters = model_params)
    ax[1].set_xlim([cparams.start,cparams.stop])
    ax[1].set_xlabel("time")
    
    ax[1].set_yticks([0,50,100])
    #ax[1].set_title("ODE simulation")
    fig.suptitle(r"$\alpha_1=10,\alpha_2=1$, 10000 cells, stochastic sim (left) vs ODE sim (right)", fontsize = 12)
    plt.tight_layout()


def plt_stochastic_model2(cells, time, ax):            
    naive_cells,th1_cells,th2_cells = get_cells(cells, time)     
    ax.plot(time, naive_cells, "k")
    ax.plot(time, th1_cells, "tab:blue")
    ax.plot(time, th2_cells, "r")
    
#fig, ax = plt.subplots(1,1,figsize =(5,5))
#for __ in range(1):
#    cells, time = run_stochastic_simulation(start = cparams.start, stop = cparams.stop, nsteps = cparams.nsteps, ncells = cparams.ncells)
#    plt_stochastic_model2(cells, time, ax)


def chain_stochastic(chain_length,
                     hill_1 = cparams.hill_1,
                     hill_2 = cparams.hill_2,
                     alpha_var = "both"):
    """
    varies chain length and saves th1 and th2 concentration at end point
    alpha_var denotes which chains should be varied, use both for both th1 and th2
    chains, "th1" for "th1" chain variation and "th2" for th2 chain variation
    """
    th1_endstates = []
    th2_endstates = []
    
    alpha_th1 = cparams.alpha_th1
    alpha_th2 = cparams.alpha_th2
    beta_th1 = cparams.beta_th1
    beta_th2 = cparams.beta_th2
    
    for i in range(1, chain_length):
        
        if alpha_var == "both":
            alpha_th1 = i
            alpha_th2 = i
            beta_th1 = i/1.
            beta_th2 = i/1.            
        elif alpha_var == "th1":
            alpha_th1 = i
            beta_th1 = i/1.            
        else:
            alpha_th2 = i
            beta_th2 = i/1.
            
        cells, time = run_stochastic_simulation(
            start = cparams.start,
            stop = cparams.stop,
            nsteps = cparams.nsteps,
            ncells = cparams.ncells,
            alpha_th1 = alpha_th1,
            alpha_th2 = alpha_th2,
            beta_th1 = beta_th1,
            beta_th2 = beta_th2,
            hill_1 = hill_1,
            hill_2 = hill_2,
            )
    
        naive_cells,th1_cells,th2_cells = get_cells(cells, time)
        
        th1_endstate = th1_cells[-1]/cparams.initial_cells
        th2_endstate = th2_cells[-1]/cparams.initial_cells        
        th1_endstates.append(th1_endstate)
        th2_endstates.append(th2_endstate)
    
    return [th1_endstates, th2_endstates]

def plt_chain_stochastic(chain_arr, ax, alpha = 1, label = False):
    th1_endstates, th2_endstates = chain_arr
    chain = range(1,len(th1_endstates)+1)
    
    if label == False:
        ax.scatter(chain, th1_endstates, c = "tab:blue", alpha = alpha)
        ax.scatter(chain, th2_endstates, c = "tab:red", alpha = alpha)
    else:
        ax.scatter(chain, th1_endstates, c = "tab:blue", alpha = alpha, label = "Th1")
        ax.scatter(chain, th2_endstates, c = "tab:red", alpha = alpha, label = "Th2")
        
    ax.set_ylim([0,1])


fig_names = [
    "pos_feedback",
    "neg_feedback",
    "pos_neg_feedback",
    "neg_th1_neg_pos_th2",
    "pos_th1_neg_pos_th2",
    "pos_th1_neg_th2",
    "pos_th2_neg_th2",
    ]


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

CHAIN_LENGTH = cparams.chain_length_stoc
SAVE_PATH = "/home/burt/Documents/tcell_project/saved_simulations/"

for hill_1, hill_2, fig_name in zip(hill_coeff_1, hill_coeff_2, fig_names):
    
    chain_arr_both = [chain_stochastic(CHAIN_LENGTH, hill_1, hill_2) for i in range(cparams.nsim)]
    chain_arr_th1 = [chain_stochastic(CHAIN_LENGTH, hill_1, hill_2, alpha_var = "th1") for i in range(cparams.nsim)]
    chain_arr_th2 = [chain_stochastic(CHAIN_LENGTH, hill_1, hill_2, alpha_var = "th2") for i in range(cparams.nsim)]
    np.savez(SAVE_PATH+fig_name+"_stoc",
             chain_arr_both = chain_arr_both,
             chain_arr_th1 = chain_arr_th1,
             chain_arr_th2 = chain_arr_th2,
             )
    
#nsim = cparams.nsim
#fig, ax = plt.subplots(1,1,figsize =(7,5))

#for i in range(nsim):
#    chain_length = 10
#    chain_arr = chain_stochastic(chain_length)
#    
#    if i == 0:
#        plt_chain_stochastic(chain_arr, ax, alpha = 1./nsim, label = True)
#    else:
#        plt_chain_stochastic(chain_arr, ax, alpha = 1./nsim)
        
#ax.set_xlabel(r"chain length $\alpha$")
#ax.set_ylabel("% Th cells after "+str(cparams.stop)+" h")

#plt.legend()
#plt.tight_layout()

#save_path = "/home/burt/Documents/tcell_project/figures/"
#fig.savefig(save_path+"chain_stochastic_pos_fb_th1.pdf", bbox_inches = "tight", dpi = 1200)