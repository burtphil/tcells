#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Thu Jan 24 11:20:42 2019

@author: burt
analysis of the th1 th2 diff model (deterministic) with prolif from yates et al
corresponding parameters stored in "det_model_params.py"
and model stored in "det_th1_th2_model.py"
here, I analyze the feedback of the system
"""
#==============================================================================
# import modules
#==============================================================================
from scipy.integrate import odeint
import numpy as np
import matplotlib.pyplot as plt
# own modules
from det_th1_th2_model import th_cell_diff
import th_diff_model_parameters as params
import seaborn as sns

#==============================================================================
# choose plotting style and saving path
#==============================================================================
sns.set(context = "talk", style = "ticks")
save_path = "/home/burt/Documents/tcell_project/figures/"

#==============================================================================
# functions
#==============================================================================
def remove_spines(ax):
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)
    ax.get_xaxis().tick_bottom()
    ax.get_yaxis().tick_left()

def adjust_lim(ax):
    ax.set_xlim(xmin = 0)
    ax.set_ylim(ymin = 0)

def get_cells(state, alpha_1, alpha_2, alpha_prolif):
    
    th1_cells = state[:, alpha_1:(alpha_1+alpha_prolif)]
    th2_cells = state[:, (alpha_1+alpha_prolif+alpha_2):]
    
    th1_cells = np.sum(th1_cells, axis = 1)
    th2_cells = np.sum(th2_cells, axis = 1)
    
    return th1_cells, th2_cells

def find_nearest(array, value):
    """
    return index of element in array that is closest to value
    """
    array = np.asarray(array)
    idx = (np.abs(array - value)).argmin()
    return idx

def get_tau(time, cells):
    tau = time[find_nearest(cells, cells[-1] / 2.)]
    
    return tau
#==============================================================================
# set up parameters
#==============================================================================
time = params.time
y0 = params.y0

alpha_1 = params.alpha_1
alpha_2 = params.alpha_2

beta_1 = params.beta_1
beta_2 = params.beta_2

alpha_prolif = params.alpha_prolif
beta_prolif = params.beta_prolif

rate_birth = params.rate_birth
rate_death = params.rate_death

rate_cytos = params.rate_cytos
fb_th1 = [10, 0, 0]

fb_th2 = params.fb_th2
hill_th1 = params.hill_th1
hill_th2 = params.hill_th2
K_th1 = params.K_th1
K_th2 = params.K_th2
conc_il12 = params.conc_il12
#==============================================================================
# run model w/o any modifications
#==============================================================================

# simple time course
state = odeint(th_cell_diff, 
               y0, 
               time, 
               args = (alpha_1, 
                       alpha_2, 
                       beta_1, 
                       beta_2,
                       alpha_prolif,
                       beta_prolif,
                       rate_birth,
                       rate_death,
                       rate_cytos,
                       fb_th1,
                       fb_th2,
                       hill_th1,
                       hill_th2,
                       K_th1,
                       K_th2,
                       conc_il12,
                       )
               )


th1_cells, th2_cells = get_cells(state, alpha_1, alpha_2, alpha_prolif)

fig, ax = plt.subplots()
ax.plot(time, th1_cells)
ax.plot(time, th2_cells)
plt.tight_layout()

#==============================================================================
# get endstates
#==============================================================================
th1_endstates = []
alpha_arr = np.arange(1,11,1)
th1_tau = []

for alpha in alpha_arr:
    
    state = odeint(th_cell_diff, 
                   y0, 
                   time, 
                   args = (alpha, 
                           alpha_2, 
                           float(alpha), 
                           beta_2,
                           alpha_prolif,
                           beta_prolif,
                           rate_birth,
                           rate_death,
                           rate_cytos,
                           fb_th1,
                           fb_th2,
                           hill_th1,
                           hill_th2,
                           K_th1,
                           K_th2,
                           conc_il12,
                           )
                   )
                   
    th1_cells = get_cells(state, alpha, alpha_2, alpha_prolif)[0]
    th1_endstate = th1_cells[-1]
    th1_endstates.append(th1_endstate)
    
    tau = get_tau(time, th1_cells)
    th1_tau.append(tau)

#==============================================================================
# plot with tau and max readouts
#==============================================================================
tau = get_tau(time, th1_cells)
th1_max = th1_cells[-1]
th1_halfmax = th1_max / 2.

fig, axes = plt.subplots(1, 3, figsize = (14,4))

ax = axes[0]
ax1 = axes[1]
ax2 = axes[2]

ax.plot(time, th1_cells)
ax.hlines(th1_max, 0, time[-1], linestyles = ":", colors = "tab:blue")
ax.hlines(th1_halfmax, 0, tau, linestyles = ":", colors = "tab:blue")
ax.vlines(tau, 0, th1_halfmax, linestyles = ":", colors = "tab:blue")
ax.text(tau+0.05,th1_halfmax / 2., r"$\tau_{1/2}$")

s = r"$1/2\ \mathrm{Th1}_\mathrm{max}$"
s2 = r"$\mathrm{Th1}_\mathrm{max}$"
ax.text(tau / 8,th1_halfmax + 0.001, s)
ax.text(tau / 8,th1_max + 0.001, s2)

ax.set_xlabel("time")
ax.set_ylabel("n_cells")

remove_spines(ax)
adjust_lim(ax)

ax1.scatter(alpha_arr, th1_endstates)
ax1.set_xlabel(r"chain length $\alpha$")
ax1.set_ylabel(r"$\mathrm{Th1}_\mathrm{max}$")
ax1.set_xticks([0,5,10])

ax2.scatter(alpha_arr, th1_tau)
ax2.set_xlabel(r"chain length $\alpha$")
ax2.set_ylabel(r"$\tau_{1/2}$")
ax2.set_xticks([0,5,10])
#for ax in axes: adjust_lim(ax)

plt.tight_layout()

#==============================================================================
# vary the feedback strength
#==============================================================================

fb_th1_c = np.asarray(fb_th1).copy()
fb_arr = np.arange(0, 20, 0.1)

th1_endstates = []
th1_tau = []

for fb in fb_arr:
    fb_th1_c[0] = fb
    
    state = odeint(th_cell_diff, 
                   y0, 
                   time, 
                   args = (alpha_1, 
                           alpha_2, 
                           beta_1, 
                           beta_2,
                           alpha_prolif,
                           beta_prolif,
                           rate_birth,
                           rate_death,
                           rate_cytos,
                           fb_th1_c,
                           fb_th2,
                           hill_th1,
                           hill_th2,
                           K_th1,
                           K_th2,
                           conc_il12,
                           )
                   )
                   
    th1_cells = get_cells(state, alpha_1, alpha_2, alpha_prolif)[0]
    th1_endstate = th1_cells[-1]
    th1_endstates.append(th1_endstate)
    
    tau = get_tau(time, th1_cells)
    th1_tau.append(tau)

# do the same but make it a rate process
th1_endstates_rate = []
th1_tau_rate = []

for fb in fb_arr:
    fb_th1_c[0] = fb
    
    state = odeint(th_cell_diff, 
                   y0, 
                   time, 
                   args = (1, 
                           alpha_2, 
                           1., 
                           beta_2,
                           alpha_prolif,
                           beta_prolif,
                           rate_birth,
                           rate_death,
                           rate_cytos,
                           fb_th1_c,
                           fb_th2,
                           hill_th1,
                           hill_th2,
                           K_th1,
                           K_th2,
                           conc_il12,
                           )
                   )
                   
    th1_cells = get_cells(state, 1, alpha_2, alpha_prolif)[0]
    th1_endstate = th1_cells[-1]
    th1_endstates_rate.append(th1_endstate)
    
    tau = get_tau(time, th1_cells)
    th1_tau_rate.append(tau)
    
fig, ax = plt.subplots(1,2, figsize = (9,4))

ax[0].plot(fb_arr, th1_endstates,)
ax[0].plot(fb_arr, th1_endstates_rate, linestyle = "--")
ax[0].set_xlabel("fb strength")
ax[0].set_ylabel(r"% Th1 cells after "+str(round(time[-1],1))+ " hrs")

ax[1].plot(fb_arr, th1_tau)
ax[1].plot(fb_arr, th1_tau_rate, linestyle = "--")
ax[1].set_xlabel("fb strength")
ax[1].set_ylabel(r"$\tau_{1/2}$")

plt.tight_layout()