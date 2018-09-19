#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Fri Sep 14 11:27:12 2018

@author: burt
chain_effects
"""

import numpy as np
import os

os.chdir("/home/burt/Documents/scripts/th_cell_differentiation")

### run_model takes default parameters stored in th1_th2_parameters.py
from th1_th2_ode_model_generic import il12
from th1_th2_plotting import  plot_il12
import th1_th2_parameters as params
import th1_th2_conceptual_parameters as cparams

#==============================================================================
# import parameters
#==============================================================================
model_params = params.parameters
conceptual_params = cparams.parameters

mean_th1 = params.alpha_th1/params.beta_th1
mean_th2 = params.alpha_th2/params.beta_th2

#==============================================================================
# IL12 effect
#==============================================================================

# plot IL12 effect
il12_concentrations = np.linspace(0,2,100)
plot_il12(il12(il12_concentrations, conceptual_params), factor_x_axis = 1, xlabel = "IL12 [a.u.]", save = "il12_conc_model_mean_diff_1")
