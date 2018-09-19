#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Wed Sep 19 14:42:32 2018

@author: burt
very basic model T can convert to T1 or T2 and I change the feedback
"""
import os

os.chdir("/home/burt/Documents/code/th_cell_differentiation")

### run_model takes default parameters stored in th1_th2_parameters.py
from th1_th2_ode_model_generic import chain, il12, run_model
from th1_th2_plotting import plot_chain, subplot_chain, subplot_tau, plot_il12,plot_time_course
import matplotlib.pyplot as plt
import numpy as np
import th1_th2_conceptual_parameters as cparams
save_path = "/home/burt/Documents/figures/"

#==============================================================================
# model parameters
#==============================================================================

conceptual_params = cparams.parameters

#==============================================================================
# case: pos feedback on Th1
#==============================================================================
params = list(conceptual_params)

# define hill parameters for pos feedback on Th1

hill_1 = [1,0,1]
hill_2 = [0,1,0]

params[6] = hill_1
params[7] = hill_2

#==============================================================================
# run time course simulation
#==============================================================================
simulation_name = "sim"

th1_th2_model = run_model(simulation_name, parameters = params)
plot_time_course(th1_th2_model, params)

#plt.plot(cparams.simulation_time,th1_th2_model)

il_12 = np.linspace(0,2,100)

#plot_il12(il12(il_12, params), xlabel ="IL12")

