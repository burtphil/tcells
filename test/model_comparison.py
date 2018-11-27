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
feedback_dict = cparams.feedback_new
feedback_type = "no_fb"

alpha_1 = cparams.alpha_th1
alpha_2 = cparams.alpha_th2
beta_1 = cparams.beta_th2
beta_2 = cparams.beta_th2
half_saturation = cparams.half_saturation
degradation = cparams.degradation
rate_ifn = cparams.rate_ifn
rate_il4 = cparams.rate_il4
conc_il12 = cparams.conc_il12
simulation_time = cparams.simulation_time
initial_cells = cparams.initial_cells
fb_ifn = cparams.fb_ifn
fb_il4 = cparams.fb_il4
fb_il12 = cparams.fb_il12
hill_1 = feedback_dict[feedback_type][0]
hill_2 = feedback_dict[feedback_type][1]
fb_start = 0
fb_end = 6
assert fb_end <= simulation_time[-1]+0.1

rate_params = [alpha_1, alpha_2, beta_1, beta_2, simulation_time, conc_il12,
              hill_1, hill_2, rate_ifn, rate_il4, half_saturation,
              initial_cells, degradation, fb_ifn, fb_il4, fb_il12, fb_start, fb_end]

rtm_params = [20, 20, 20, 20, simulation_time, conc_il12,
              hill_1, hill_2, rate_ifn, rate_il4, half_saturation,
              initial_cells, degradation, fb_ifn, fb_il4, fb_il12, fb_start, fb_end]

#==============================================================================
# run time course simulation
#==============================================================================
th1_th2_model = run_model(*rate_params)
rtm_simu = run_model(*rtm_params)
# plot time course
#plot_time_course(th1_th2_model, parameters)
#fig, ax = plt.subplots()
#ax_time_course(th1_th2_model, ax = ax, simulation_time = simulation_time, initial_cells = initial_cells, alpha_1 = alpha_1)
#ax_time_course(rtm_simu, ax = ax, simulation_time = simulation_time, initial_cells = initial_cells, alpha_1 = 20)
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
        cells, time = run_stochastic_simulation(start = cparams.start,
                                                stop = cparams.stop,
                                                nsteps = cparams.nsteps,
                                                ncells = cparams.ncells)
        
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
    #ax.plot(time, naive_cells, "k")
    ax.plot(time, th1_cells, "tab:blue")
    ax.plot(time, th2_cells, "r")

stoc_simus = [run_stochastic_simulation(start = cparams.start,stop = cparams.stop,nsteps = cparams.nsteps,ncells = cparams.ncells) for _ in range(10)]
th1_th2_cells = [get_cells(simu[0], simu[1])[1:] for simu in stoc_simus]

label = ["Th1","Th2"]
labels = [label for _ in range(10)]
flat_labels = [item for sublist in labels for item in sublist]
flat_cells = [item for sublist in th1_th2_cells for item in sublist]
df = pd.DataFrame.from_records(flat_cells).transpose()
df.columns = flat_labels

df["time"] = simulation_time
df_long = df.melt(id_vars = ["time"])

stoc_simus_rtm = [run_stochastic_simulation(start = cparams.start,stop = cparams.stop,nsteps = cparams.nsteps,ncells = cparams.ncells, alpha_th1 = 20, alpha_th2 = 20, beta_th1 = 20, beta_th2 = 20) for _ in range(10)]
th1_th2_cells_rtm = [get_cells(simu[0], simu[1])[1:] for simu in stoc_simus_rtm]
flat_cells_rtm = [item for sublist in th1_th2_cells_rtm for item in sublist]
df_rtm = pd.DataFrame.from_records(flat_cells_rtm).transpose()
df_rtm.columns = flat_labels

df_rtm["time"] = simulation_time
df_long_rtm = df_rtm.melt(id_vars = ["time"])

fig, ax = plt.subplots(1, 2, figsize =(10,4))
stoc_plot = sns.lineplot(x = "time", y = "value",  data = df_long, hue = "variable", ci = "sd", ax = ax[0], palette = ["tab:blue", "tab:red"])
### for this I set the probabilities constant to pTh1 =0.7 and pTh2 = 0.3
for i in range(2):
    stoc_plot.lines[i].set_linestyle("--")

stoc_plot_rtm = sns.lineplot(x = "time", y = "value", data = df_long_rtm, hue = "variable", ci = "sd", ax = ax[0], palette = ["tab:blue", "tab:red"])

ax[0].set_ylim(0,1000)
ax[0].set_yticks([0,500,1000])
ax[0].set_ylabel(r"$n_{Tcells}$")
ax[0].set_xlim(0,6)

ax_time_course(th1_th2_model, ax = ax[1], simulation_time = simulation_time, initial_cells = initial_cells, alpha_1 = alpha_1, linestyle = "--")
ax_time_course(rtm_simu, ax = ax[1], simulation_time = simulation_time, initial_cells = initial_cells, alpha_1 = 20)
ax[1].set_xticks([0,2,4,6])
plt.tight_layout()

save_path = "/home/burt/Documents/tcell_project/figures/"
#fig.savefig(save_path+"model_comparison.svg", bbox_inches = "tight")

def chain_stochastic(chain_length,
                     hill_1 = cparams.stoc_hill_1,
                     hill_2 = cparams.stoc_hill_2,
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

"""
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
"""