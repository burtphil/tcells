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
from modules.th1_th2_ode_model_generic import run_model
from modules.th1_th2_plotting import ax_time_course, ax_il12, ax_chain
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
fb_start = 1
fb_end = 2

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
        alpha_1, alpha_2, beta_1, beta_2 = 1, 10, 1, 10
    
    if model_type == "rtm_th1_rate_th2":
        alpha_1, alpha_2, beta_1, beta_2 = 10, 1, 10, 1
        
    params[0] = alpha_1
    params[1] = alpha_2
    params[2] = beta_1
    params[3] = beta_2
    
    return params

#==============================================================================
# run simulations
#==============================================================================
rtm_params = assign_fb(feedback = "pos_th1", model_type = "rtm", params = rate_params)

rate_simu = run_model(*rate_params)
rtm_simu = run_model(*rtm_params)

fig, ax = plt.subplots(1,1, figsize = (5,4))
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
plt.tight_layout()