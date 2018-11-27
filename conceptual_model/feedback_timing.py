#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Wed Nov 21 17:12:05 2018

@author: burt
"""

#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Wed Nov 14 17:40:39 2018

@author: burt

time course simulations
"""
import numpy as np
import os
import matplotlib.image as mpimg
os.chdir("/home/burt/Documents/code/th_cell_differentiation")

### run_model takes default parameters stored in th1_th2_parameters.py
from modules.th1_th2_ode_model_generic import run_model, feedback_timing, feedback_duration
from modules.th1_th2_plotting import ax_il12, ax_chain
import modules.th1_th2_conceptual_parameters as params
import matplotlib.pyplot as plt
import seaborn as sns

#==============================================================================
# define params to store figure
#==============================================================================
sns.set(context = "talk", style = "ticks")
save_path = "/home/burt/Documents/tcell_project/figures/"
image_path = "/home/burt/Documents/code/th_cell_differentiation/images/"

#==============================================================================
# define plotting params
#==============================================================================

#==============================================================================
# load parameters
#==============================================================================
feedback_dict = params.feedback_new
feedback_type = "pos_th1"

alpha_1 = params.alpha_th1
alpha_2 = params.alpha_th2
beta_1 = params.beta_th2
beta_2 = params.beta_th2
half_saturation = params.half_saturation
degradation = params.degradation
rate_ifn = params.rate_ifn
rate_il4 = params.rate_il4
conc_il12 = params.conc_il12
simulation_time = params.simulation_time
initial_cells = params.initial_cells
fb_ifn = params.fb_ifn
fb_il4 = params.fb_il4
fb_il12 = params.fb_il12
hill_1 = feedback_dict[feedback_type][0]
hill_2 = feedback_dict[feedback_type][1]
fb_start = 0
window = 0.1
fb_end = 1
assert fb_end <= simulation_time[-1]

rate_params = [alpha_1, alpha_2, beta_1, beta_2, simulation_time, conc_il12,
              hill_1, hill_2, rate_ifn, rate_il4, half_saturation,
              initial_cells, degradation, fb_ifn, fb_il4, fb_il12, fb_start, fb_end]


#==============================================================================
# functions
#==============================================================================
def ax_time_course(state, ax, simulation_time, initial_cells, alpha_1, linestyle = "-"):
    
    norm = initial_cells / 100
    th1_cells = state[:,alpha_1+1] / norm
    th2_cells = state[:,-1] / norm
    
    #print th1_cells[-1],th2_cells[-1]
    #ax.plot(simulation_time, th0_cells, color = "k", linestyle = linestyle)
    ax.plot(simulation_time, th1_cells, color = "tab:blue", linestyle = linestyle)
    #ax.plot(simulation_time, th2_cells, color = "tab:red", linestyle = linestyle)
    ax.set_yticks([0, 50, 100])
    ax.set_ylim([0, 100])
    ax.set_xlim([simulation_time[0], simulation_time[-1]])
    ax.set_xlabel("time")
    ax.set_ylabel("% Th cells")

def get_hill(feedback_type, fb_dict = feedback_dict):
    return feedback_dict[feedback_type]

def assign_fb(feedback, model_type, params = rate_params):
    params = list(params)
    params[6] = get_hill(feedback)[0]
    params[7] = get_hill(feedback)[1]
    
    assert model_type == "rate" or "rtm" or "rate_th1_rtm_th2" or "rtm_th1_rate_th2"

    if model_type == "rate":
        alpha_1, alpha_2, beta_1, beta_2 = 1, 1, 1, 1
        
    if model_type == "rtm":
        alpha_1, alpha_2, beta_1, beta_2 = 20, 20, 20, 20
    
    if model_type == "rate_th1_rtm_th2":
        alpha_1, alpha_2, beta_1, beta_2 = 1, 10, 1, 10
    
    if model_type == "rtm_th1_rate_th2":
        alpha_1, alpha_2, beta_1, beta_2 = 10, 1, 10, 1
        
    params[0] = alpha_1
    params[1] = alpha_2
    params[2] = beta_1
    params[3] = beta_2
    
    return params

def assign_window(params, window):
    params = list(params)
    params[-1] = window
    return params
#==============================================================================
# run simulations
#==============================================================================
rtm_params = assign_fb(feedback = "pos_th1", model_type = "rtm", params = rate_params)

rate_simu = run_model(*rate_params)
rtm_simu = run_model(*rtm_params)

rate_params2 = list(rate_params)
rate_params2[-1] = 2
rate_params2[-2] = 1

rtm_params2 = list(rtm_params)
rtm_params2[-1] = 2
rtm_params2[-2] = 1

rate_simu2 = run_model(*rate_params2)
rtm_simu2 = run_model(*rtm_params2)


fig, ax = plt.subplots(1,2, figsize = (10,4))
ax_time_course(rate_simu, 
               ax[0], 
               simulation_time = simulation_time, 
               alpha_1 = 1, 
               initial_cells = initial_cells, 
               linestyle = "--")

ax_time_course(rtm_simu, 
               ax[0], 
               simulation_time = simulation_time, 
               alpha_1 = 20, 
               initial_cells = initial_cells)


ax_time_course(rate_simu2, 
               ax[1], 
               simulation_time = simulation_time, 
               alpha_1 = 1, 
               initial_cells = initial_cells, 
               linestyle = "--")

ax_time_course(rtm_simu2, 
               ax[1], 
               simulation_time = simulation_time, 
               alpha_1 = 20, 
               initial_cells = initial_cells)

ax[0].set_xticks([0,2,4,6])
ax[0].axvspan(0, 1, alpha=0.5, color='gray')
ax[1].axvspan(1, 2, alpha=0.5, color='gray')
plt.tight_layout()

# change ax_time_course function before to include th2 cellss
#fig.savefig(save_path+"fb_onset_time_course.svg", bbox_inches = "tight")


fb_start_arr = np.linspace(0,2,100)
windows = [0.25,0.5,1.0]
fb_durations = np.linspace(0,3,100)



timing_params_rtm = [assign_window(rtm_params, w) for w in windows]
timing_params_rate = [assign_window(rate_params, w) for w in windows]

rtm_timing_simu = [feedback_timing(fb_start_arr, *p) for p in timing_params_rtm]
rate_timing_simu = [feedback_timing(fb_start_arr, *p) for p in timing_params_rate]

diff = []

for i in range(len(rtm_timing_simu)):
    rate_simu = rate_timing_simu[i][1]
    rtm_simu = rtm_timing_simu[i][1]
    th1_diff = rate_simu - rtm_simu
    diff.append(th1_diff)
    
fig, ax = plt.subplots()
ax.plot(fb_start_arr, diff[0], label = "fb duration = 0.25")
ax.plot(fb_start_arr, diff[1], label = "fb duration = 0.5")
ax.plot(fb_start_arr, diff[2], label = "fb duration = 1.0")
ax.set_xlabel("feedback onset time")
ax.set_ylabel(r"Th1 difference")   
ax.set_ylim(-2,15)
ax.legend()
ax.set_yticks([0, 7.5, 15])
plt.tight_layout()
#fig.savefig(save_path+"feedback_onset_diff.svg", bbox_inches = "tight")

fig, axes = plt.subplots(1, 3, figsize = (12,4))
axes = axes.flatten()

for ax, rtm_simu, rate_simu, w in zip(axes, rtm_timing_simu, rate_timing_simu, windows):
    
    ax_il12(rate_simu, ax, linestyle = "--")
    ax_il12(rtm_simu, ax)
    ax.set_title(str(w))

plt.tight_layout()


fig, ax = plt.subplots(1, 1, figsize = (5,4))

for rtm_simu, rate_simu, w in zip(rtm_timing_simu, rate_timing_simu, windows):
    y1 = rate_simu[1]
    y2 = rtm_simu[1]
    ax.plot(fb_start_arr, y1, linestyle = "--")
    ax.plot(fb_start_arr, y2)    
    
plt.tight_layout()

fb_dur_rate = list(rate_params)
fb_dur_rtm = list(rtm_params)
fb_duration_arr = np.linspace(0,2,100)

fb_dur_simu_rate = feedback_duration(fb_duration_arr, *fb_dur_rate)
fb_dur_simu_rtm = feedback_duration(fb_duration_arr, *fb_dur_rtm)


fig, ax = plt.subplots(1,2, figsize = (10,4))

ax_il12(fb_dur_simu_rate, ax[0], plot_both = False, linestyle = "--")
ax_il12(fb_dur_simu_rtm, ax[0], plot_both = False)
ax[0].set_xlabel("feedback duration")
ax[0].set_ylim(49,80)
ax[0].set_yticks([50, 60, 70, 80])

y1 = rate_timing_simu[1][1]
y2 = rtm_timing_simu[1][1]
ax[1].plot(fb_start_arr, y1, "tab:blue",  linestyle = "--")
ax[1].plot(fb_start_arr, y2, "tab:blue")
ax[1].set_xlabel("feedback onset time")
ax[1].set_ylabel("% Th1 cells")  
ax[1].set_ylim(49,80)
ax[1].set_xlim(0,2)
#ax.set_ylim(25,75)
plt.tight_layout()
fig.savefig(save_path+"feedback_timings.svg", bbox_inches = "tight")
