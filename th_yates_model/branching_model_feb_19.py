#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Fri Feb  1 11:28:36 2019

@author: burt
"""
import numpy as np
import matplotlib.pyplot as plt
# own modules
from det_th1_th2_model import run_model
from det_th1_th2_model import get_tau
from det_th1_th2_model import get_cells
from det_th1_th2_model import vary_param
import th_diff_model_parameters as params
import seaborn as sns

#==============================================================================
# choose plotting style and saving path
#==============================================================================
sns.set(context = "talk", style = "ticks")
save_path = "/home/burt/Documents/phd_summary/20190204/mutual_inhibition/"
#==============================================================================
# functions
#==============================================================================
def adjust_lim(ax, arr):
    ax.set_xlim(0, arr[-1])
    
def run_analysis(arr_dict, params):
    for key, val in arr_dict.iteritems():
        new_params = list(params)
        arr = vary_param(key, val[0], *new_params)
              
        fig, ax = plt.subplots(1,2, figsize = (9,4))
        
        ylabels = ["Th max", r"$\tau_{1/2}$"]
    
        for i in range(len(ylabels)):
            ax[i].set_ylabel(ylabels[i])
            ax[i].set_xlabel(val[1])
            
        if key == "chain":
            
            ax[0].scatter(val[0], arr[0], c = "tab:blue")
            ax[0].scatter(val[0], arr[1], c = "tab:red")
            ax[1].scatter(val[0], arr[2], c = "tab:blue")
            ax[1].scatter(val[0], arr[3], c = "tab:red")

        elif key == "prolif":
            ax[0].scatter(val[0], arr[0], c = "tab:blue")
            ax[0].scatter(val[0], arr[1], c = "tab:red")
            ax[1].scatter(val[0], arr[2], c = "tab:blue")
            ax[1].scatter(val[0], arr[3], c = "tab:red")
        
        elif key == "fb_th1":
            ax[0].plot(val[0], arr[0], c = "tab:blue")
            ax[0].plot(val[0], arr[1], c = "tab:red", linestyle = "--")
            ax[1].plot(val[0], arr[2], c = "tab:blue")
            ax[1].plot(val[0], arr[3], c = "tab:red", linestyle = "--")            
            ax[0].set_xscale('log')
            ax[1].set_xscale('log')
            ax[0].set_xticks([1, 10])
            ax[0].set_xlim(0,10)        
            ax[1].set_xticks([1, 10])
            ax[1].set_xlim(0,10)
            
        else:    
            ax[0].plot(val[0], arr[0], c = "tab:blue")
            ax[0].plot(val[0], arr[1], c = "tab:red", linestyle = "--")
            ax[1].plot(val[0], arr[2], c = "tab:blue")
            ax[1].plot(val[0], arr[3], c = "tab:red", linestyle = "--")
            
            #adjust_lim(ax[0], val)
            #adjust_lim(ax[1], val)
            
        plt.tight_layout()
        fig.savefig(save_path+val[1]+".pdf", bbox_inches = "tight")
#==============================================================================
# parameters
#==============================================================================
time = params.time
cells_0 = params.cells_0

alpha_1 = params.alpha_1
alpha_2 = params.alpha_2

beta_1 = params.beta_1
beta_2 = params.beta_2

alpha_prolif = params.alpha_prolif
beta_prolif = params.beta_prolif

rate_birth = params.rate_birth
rate_death = params.rate_death
rate_diff_th0 = params.rate_diff_th0
rate_death_th0 = params.rate_death_th0

rate_cytos = params.rate_cytos
fb_th1 = params.fb_th1
fb_th2 = params.fb_th2
hill_th1 = params.hill_th1
hill_th2 = params.hill_th2
K_th1 = params.K_th1
K_th2 = params.K_th2
conc_il12 = params.conc_il12

params = [cells_0,
          time,
          alpha_1,
          alpha_2,
          beta_1,
          beta_2,
          alpha_prolif,
          beta_prolif,
          rate_birth,
          rate_death,
          rate_death_th0,
          rate_diff_th0,
          rate_cytos,
          fb_th1,
          fb_th2,
          hill_th1,
          hill_th2,
          K_th1,
          K_th2,
          conc_il12,]

#==============================================================================
# output
#==============================================================================
# simple time course
state = run_model(*params)

th1_cells, th2_cells = get_cells(state, alpha_1, alpha_2, alpha_prolif)

fig, ax = plt.subplots(figsize = (5,4))
ax.plot(time, th1_cells, label = "rtm")
ax.plot(time, th2_cells, "tab:red", linestyle = "--", label = "rate process")
ax.legend(loc = 5)
ax.set_xlabel("time")
ax.set_ylabel("Th cells")
adjust_lim(ax, time)
ax.set_ylim(ymin = 0)
plt.tight_layout()
#fig.savefig(save_path+"semi_quant.pdf", bbox_inches = "tight")
#fig, ax = plt.subplots()
#ax.plot(time, state)
#ax.legend(["Th0", "Th1p", "Th1_diff", "Th2p", "Th2diff"])

# make a dictionary with all the parameters and corresponding arrays to set up
arr_dict = {
        "chain" : [np.arange(1,11,1), r"$\alpha$ diff"],
        "prolif" : [np.arange(1,11,1), r"$\alpha$ prolif"],
        "rate_birth" : [np.arange(0.1,5, 0.1), "rate influx"],       
        "rate_death" : [np.arange(3.0,5,0.1), "rate death th diff"],
        "rate_death_th0" : [np.arange(0.1,5,0.1), "rate death th0"],
        "rate_diff_th0" : [np.arange(0.1,5,0.1), "rate diff th0"],
        "fb_th1" : [np.arange(0, 10.1,0.1), "feedback fold-change"],
        "rate_cytos" : [np.arange(0,100,0.1), "rate cytokine secretion"],
        }

dict2 = {
        "fb_th1" : [np.arange(0, 10.1,0.1), "feedback fold-change"]
        }
      
run_analysis(arr_dict, params)
#run_analysis(dict2, params)

#==============================================================================
# feedback sensibility
#==============================================================================
param_name = "fb_th1"
param_arr = np.arange(0,10.1,0.1)

params_fb_rate = [cells_0,
                  time,
                  1,
                  1,
                  1,
                  1,
                  alpha_prolif,
                  beta_prolif,
                  rate_birth,
                  rate_death,
                  rate_death_th0,
                  rate_diff_th0,
                  rate_cytos,
                  fb_th1,
                  fb_th2,
                  hill_th1,
                  hill_th2,
                  K_th1,
                  K_th2,
                  conc_il12,]

params_fb_rtm = [cells_0,
  time,
  10,
  10,
  10,
  10,
  alpha_prolif,
  beta_prolif,
  rate_birth,
  rate_death,
  rate_death_th0,
  rate_diff_th0,
  rate_cytos,
  fb_th1,
  fb_th2,
  hill_th1,
  hill_th2,
  K_th1,
  K_th2,
  conc_il12,]


rate_arr = vary_param(param_name, param_arr, *params_fb_rate)
rtm_arr = vary_param(param_name, param_arr, *params_fb_rtm)

"""
fig, axes = plt.subplots(1,2, figsize = (9,4))

ax = axes[0]
ax1 = axes[1]

ax.plot(param_arr, rtm_arr[0], c = "tab:blue", label = "rtm")
ax.plot(param_arr, rate_arr[0], c = "tab:blue", linestyle = "--", label = "rate process")
ax.set_xscale('log')
ax.set_xticks([1, 10])
ax.set_xlim(0,10)        
ax.set_ylabel("Th max")
ax.set_xlabel("feedback fold-change")
ax.legend(loc = 1)
th1_rate = rate_arr[0]
th1_rate_0 = th1_rate[10]

th1_rtm = rtm_arr[0]
th1_rtm_0 = th1_rtm[10]

th1_rate_rel = (th1_rate - th1_rate_0) / th1_rate_0
th1_rtm_rel = (th1_rtm - th1_rtm_0) / th1_rtm_0

ax1.plot(param_arr, th1_rtm_rel)
ax1.plot(param_arr, th1_rate_rel, "tab:blue", linestyle = "--") 
ax1.set_xscale('log')
ax1.set_xticks([1, 10])
ax1.set_xlim(0,10)
ax1.set_xlabel("feedback fold-change")
ax1.set_ylabel("relative feedback sensitivity")
plt.tight_layout()  
#fig.savefig(save_path+"fb_robustness.pdf", bbox_inches = "tight")
"""