#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Thu Sep 13 11:48:04 2018

@author: burt
parameters for conceptual model simulation with equal means and all
"""
import numpy as np
#hill coefficients for probability of th1 differentiation
k_th1_ifn = 1
k_th1_il4 = 0
k_th1_il12 = 1

hill_1 = [k_th1_ifn, k_th1_il4, k_th1_il12]

#hill coefficients for probability of th2 differentiation
k_th2_ifn = 0
k_th2_il4 = 1
k_th2_il12 = 0

hill_2 = [k_th2_ifn, k_th2_il4, k_th2_il12]

# extracellular il12 concentration
conc_il12 = 1.0

#production rates cytokines
rate_ifn = 1.
rate_il4 = 1.

# half saturation constants
kd_ifn = 1.
kd_il4 = 1.
kd_il12 = 1.

half_saturation = [kd_ifn, kd_il4, kd_il12]

base_production_rate_il4 = 0.
base_production_rate_ifn = 0.1
### at some point I need to change this to cell densities
initial_cells = 10.

#
mean_th1 = 1.
mean_th2 = 1.
alpha_th1 = 1
alpha_th2 = 1
beta_th1 = alpha_th1/mean_th1
beta_th2 = alpha_th2/mean_th2

# simulation time

start = 0
stop = 6
stepsize = .01
simulation_time = np.arange(start, stop, stepsize)

parameters = [alpha_th1, alpha_th2, beta_th1, beta_th2, simulation_time,
              conc_il12, hill_1, hill_2, rate_ifn, rate_il4, half_saturation,
              base_production_rate_ifn, base_production_rate_il4, initial_cells]