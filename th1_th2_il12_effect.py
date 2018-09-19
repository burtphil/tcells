#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Fri Sep  7 15:26:58 2018

@author: burt
"""

import numpy as np
import os
import matplotlib.pyplot as plt
from th1_th2_plotting import plot_time_course

os.chdir("/home/burt/Documents/scripts/th_cell_differentiation")

### import parameters from parameter
#import th1_th2_parameters as params
from th1_th2_ode_model_generic import run_model

######## plotting params
SMALL_SIZE = 15
MEDIUM_SIZE = 20

plt.rc('font', size=MEDIUM_SIZE)          # controls default text sizes
plt.rc('axes', titlesize=SMALL_SIZE)     # fontsize of the axes title
plt.rc('axes', labelsize=MEDIUM_SIZE)    # fontsize of the x and y labels
plt.rc('xtick', labelsize=SMALL_SIZE)    # fontsize of the tick labels
plt.rc('ytick', labelsize=SMALL_SIZE)    # fontsize of the tick labels
plt.rc('legend', fontsize=15)    # legend fontsize
plt.rc('figure', titlesize=MEDIUM_SIZE)

### import fit parameters
alpha_1, alpha_2, rate1, rate2 = np.load('th1_th2_gamma_fit_params.npy')

### simulation time
starttime = 0
endtime = 100
stepsize = 0.01
simulation_time = np.arange(starttime,endtime,stepsize)


def find_nearest(array, value):
    """
    return index of element in array that is closest to value
    """
    array = np.asarray(array)
    idx = (np.abs(array - value)).argmin()
    return idx
##############################################################################
########################### effect of IL12 ###################################
##############################################################################
"""
il12_concentrations = np.linspace(0,6*10**(-12),100)

th1_conc = []
th2_conc = []
th0_conc = []

for i in il12_concentrations:
    
    state = run_model(alpha_1, alpha_2, rate1, rate2, simulation_time, conc_il12=i)
    
    th0_endstate = state[-1,0]
    th0_conc.append(th0_endstate)
        
    th1_endstate = state[-1,int(alpha_1)]
    th1_conc.append(th1_endstate)
    
    th2_endstate = state[-1,-1]
    th2_conc.append(th2_endstate)
    
th0_conc = np.asarray(th0_conc)
th1_conc = np.asarray(th1_conc)
th2_conc = np.asarray(th2_conc)

all_tcells = th0_conc+th1_conc+th2_conc

th0_conc = (th0_conc / all_tcells)*100
th1_conc = (th1_conc / all_tcells)*100
th2_conc = (th2_conc / all_tcells)*100

il12_conc_pm = il12_concentrations*(10**(12))


fig, ax = plt.subplots(1,1,figsize=(5,5))

#ax.plot(il12_conc_pm, th0_conc)
ax.plot(il12_conc_pm, th1_conc, "tab:blue")
ax.plot(il12_conc_pm, th2_conc, "tab:red")
ax.set_xlabel("IL12 [pM]")
ax.set_yticks([0,50,100])
#ax.set_title("Effect of IL-12 conc. on Th cell balance after 3 hours")
ax.set_ylim([0,100])
ax.set_xlim([il12_conc_pm[0],il12_conc_pm[-1]])

ax.set_ylabel("% Th cells")
plt.legend(["Th1","Th2"])
plt.tight_layout()
"""
##############################################################################
########################### effect of chain length ###########################
##############################################################################
chain_length = 25
chain = np.arange(1,chain_length,1)

mean_th1 = alpha_1 / rate1
mean_th2 = alpha_2 /rate2

th1_conc = []
th2_conc = []

th1_tau = []
th2_tau = []

for i in chain:
    
    alpha_1 = i
    alpha_2 = i
    rate1 = alpha_1/mean_th1
    rate2 = alpha_2/mean_th2
    state = run_model(alpha_1, alpha_2, rate1, rate2, simulation_time)
    #plot_time_course(state, alpha_1, alpha_2, simulation_time)

    # get end states
    th1_endstate = state[-1,int(alpha_1)]
    th1_conc.append(th1_endstate)
    
    th2_endstate = state[-1,-1]
    th2_conc.append(th2_endstate)
    
    th1_halfmax = th1_endstate/2
    th2_halfmax = th2_endstate/2
    
    th1_tau_idx = find_nearest(state[:,int(alpha_1)], th1_halfmax)*stepsize
    th2_tau_idx = find_nearest(state[:,-1], th2_halfmax)*stepsize
    th1_tau.append(th1_tau_idx)
    th2_tau.append(th2_tau_idx)    
    # normalize to initial cell pop
    
th1_conc = np.array(th1_conc)/100
th2_conc = np.array(th2_conc)/100

fig, ax = plt.subplots(1,1, figsize = (5,5))
ax.scatter(chain, th1_conc, c= "tab:blue")
ax.scatter(chain, th2_conc, c= "tab:red")
ax.set_xlabel(r"chain length $\alpha$")
ax.set_ylabel("% Th cells after 100 hrs")
ax.set_ylim([0,100])
ax.set_xlim([0, chain_length])
ax.set_yticks([0,50,100])

plt.tight_layout()
fig.savefig("th1_th2_steady_state_vs_alpha.svg", bbox_inches = "tight", dpi = 1200)

fig, ax = plt.subplots(1,1, figsize = (5,5))
ax.scatter(chain, th1_tau, c= "tab:blue")
ax.scatter(chain, th2_tau, c= "tab:red")
ax.set_xlabel(r"chain length $\alpha$")
ax.set_ylabel(r"$\tau_{1/2}$")
ax.set_ylim([0,50])
ax.set_xlim([0, chain_length])
ax.set_yticks([0,25,50])
plt.tight_layout()

fig.savefig("th1_th2_tau_vs_alpha.svg", bbox_inches = "tight", dpi = 1200)
