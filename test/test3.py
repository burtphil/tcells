#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Thu Nov 22 11:52:55 2018

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
from modules.th1_th2_ode_model_generic import run_model, il12, get_readouts, hilli_sym, hilli_pos_th1
from modules.th1_th2_plotting import ax_time_course, ax_il12
import modules.th1_th2_conceptual_parameters as params
import matplotlib.pyplot as plt
import seaborn as sns

#==============================================================================
# define params to store figure
#==============================================================================
sns.set(context = "talk", style = "ticks")
save_path = "/home/burt/Documents/tcell_project/figures/"
image_path = "/home/burt/Documents/code/th_cell_differentiation/images/"

fig_names = ["pos_feedback",
            "neg_feedback",
            "pos_neg_feedback",
            "neg_th1_neg_pos_th2",
            "pos_th1_neg_pos_th2",
            "pos_th1_neg_th2",
            "pos_th2_neg_th2"]
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
fb_start = params.fb_start
fb_end = params.fb_end

rate_params = [alpha_1, alpha_2, beta_1, beta_2, simulation_time, conc_il12,
              hill_1, hill_2, rate_ifn, rate_il4, half_saturation,
              initial_cells, degradation, fb_ifn, fb_il4, fb_il12, fb_start, fb_end]

#==============================================================================
# functions
#==============================================================================
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
        alpha_1, alpha_2, beta_1, beta_2 = 1, 20, 1, 20
    
    if model_type == "rtm_th1_rate_th2":
        alpha_1, alpha_2, beta_1, beta_2 = 20, 1, 20, 1
        
    params[0] = alpha_1
    params[1] = alpha_2
    params[2] = beta_1
    params[3] = beta_2
    
    return params
#==============================================================================
# initialize models
#==============================================================================
model_types = ["rate", "rtm", "rtm_th1_rate_th2", "rate_th1_rtm_th2"]

#==============================================================================
# run simulations and get readouts
#==============================================================================
rtm_params = assign_fb(feedback = "pos_th1", model_type = "rtm", params = rate_params)

rate_simu = run_model(*rate_params)
rtm_simu = run_model(*rtm_params)
#==============================================================================
# load images
#==============================================================================
no_fb_img = mpimg.imread(image_path + "no_fb.png")
pos_fb_img = mpimg.imread(image_path + "pos_feedback.png")

#==============================================================================
# compare pos fb and no fb for rtm and rate model
#==============================================================================

fig, ax = plt.subplots(1,1, figsize = (8,7))


ax_time_course(rate_simu, 
               ax, 
               simulation_time = simulation_time, 
               alpha_1 = 1, 
               initial_cells = initial_cells, 
               linestyle = "--")

ax_time_course(rtm_simu, 
               ax, 
               simulation_time = simulation_time, 
               alpha_1 = 20, 
               initial_cells = initial_cells)


#ax.legend([r"Th0","Th1","Th2"])
ax.set_xlabel("time")
plt.tight_layout()

fig, ax = plt.subplots(1,1, figsize = (5,4))
xlabel = "fb strength pos fb on Th1"
rate_params = assign_fb(feedback = "pos_th1", model_type = "rate", params = rate_params)
rtm_params = assign_fb(feedback = "pos_th1", model_type = "rtm", params = rate_params)

il12_arr = np.linspace(0, 10, 100)

rate_il12 = hilli_pos_th1(il12_arr, *rate_params)
rtm_il12 = hilli_pos_th1(il12_arr, *rtm_params)
ax_il12(rate_il12, ax, linestyle = "--")
ax_il12(rtm_il12, ax)
ax.set_title("pos fb")
ax.set_xlabel(xlabel)
plt.tight_layout()