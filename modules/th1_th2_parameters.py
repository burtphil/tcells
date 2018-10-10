#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Thu Sep  6 12:56:08 2018

@author: burt
"""
import numpy as np
#hill coefficients for probability of th1 differentiation
k_th1_ifn = 1
k_th1_il4 = -1
k_th1_il12 = 1

hill_1 = [k_th1_ifn, k_th1_il4, k_th1_il12]

#hill coefficients for probability of th2 differentiation
k_th2_ifn = -1
k_th2_il4 = 1
k_th2_il12 = 0

hill_2 = [k_th2_ifn, k_th2_il4, k_th2_il12]

# extracellular il12 concentration
conc_il12 = 5*10**(-12)

#production rates cytokines
rate_ifn = 10**(-15)
rate_il4 = 10**(-15)

# half saturation constants
kd_ifn = 10**(-12)
kd_il4 = 10**(-12)
kd_il12 = 10**(-12)

half_saturation = [kd_ifn, kd_il4, kd_il12]

### at some point I need to change this to cell densities
initial_cells = 10000.

### import fit parameters
alpha_th1, alpha_th2, beta_th1, beta_th2 = np.load('/home/burt/Documents/code/th_cell_differentiation/th1_th2_gamma_fit_params.npy')

alpha_th1 = int(alpha_th1)
alpha_th2 = int(alpha_th2)
# simulation time
start = 0
stop = 100
stepsize = .1
simulation_time = np.arange(start, stop, stepsize)

degradation = 0

parameters = [alpha_th1, alpha_th2, beta_th1, beta_th2, simulation_time,
              conc_il12, hill_1, hill_2, rate_ifn, rate_il4, half_saturation,
              initial_cells, degradation]