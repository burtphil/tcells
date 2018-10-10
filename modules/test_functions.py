#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Wed Oct 10 09:41:13 2018

@author: burt

test script
"""

#==============================================================================
# test conceptual model
# first perform simulation
#==============================================================================
from th1_th2_ode_model_generic import chain, il12, run_model, chain_one, chain_th1, chain_th2
from th1_th2_plotting import plot_chain, ax_chain, subplot_tau, plot_il12,plot_time_course, ax_time_course, ax_il12
import matplotlib.pyplot as plt
import numpy as np
import th1_th2_conceptual_parameters as cparams
import th1_th2_parameters as params

#==============================================================================
# time courses 
#==============================================================================
data_params = params.parameters
conc_params = cparams.parameters

simulation_name = "sim"
th1_th2_model = run_model(simulation_name, parameters = data_params)

test_simulation = np.load(simulation_name+".npz")
state = test_simulation["state"]
parameters = test_simulation["parameters"]

# plot time course
plot_time_course(th1_th2_model, parameters)

simulation_name = "sim"
th1_th2_model = run_model(simulation_name, parameters = conc_params)

test_simulation = np.load(simulation_name+".npz")
state = test_simulation["state"]
parameters = test_simulation["parameters"]


# plot time course
plot_time_course(th1_th2_model, parameters)

#==============================================================================
# IL12 concentration 
#==============================================================================
data_il12 = np.linspace(0,5*10**(-12),100)
conc_il12 = np.linspace(0,2,100)

plot_il12(il12(data_il12, data_params))
plot_il12(il12(conc_il12, conc_params))

#==============================================================================
# Chain effects 
#==============================================================================
chain_length = 10

data_chain = chain(chain_length, data_params)
conc_chain = chain(chain_length, conc_params)
plot_chain(data_chain)
plot_chain(conc_chain)

