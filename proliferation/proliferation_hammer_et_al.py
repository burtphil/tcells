#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Wed Sep 12 10:06:23 2018

@author: burt

implementation of Hammer et al proliferation model
with direct analytical solution for division states > 1
"""
import numpy as np
from scipy.special import gamma, gammainc
from scipy.integrate import odeint
import matplotlib.pyplot as plt
import os
windows_path= "C:/Users/Philipp/Documents/tcells/modules"
linux_path = "/home/burt/Documents/code/th_cell_differentiation/modules"
os.chdir(linux_path)
import prolif_params as param
from scipy.interpolate import interp1d

def gamma_dist(t, alpha, beta):
    return np.exp(-beta*t)*(beta**alpha)*(t**(alpha-1))/gamma(alpha)

def gamma_cdf(t, alpha, beta):
    dummy = t*beta
    return gammainc(alpha,dummy)

def N_0(t, non_dividing_cells, d, alpha, beta):
    """
    analytical function N_0(t) for cell numbers that have undergone 0 divisions
    """
    cell0 = (1-gamma_cdf(t, alpha, beta))+non_dividing_cells*np.exp(-d*t)
    return cell0

def N_i(t,i, d, div_t, n_1):
    """
    analytical function N_i(t) for cell numbers that have undergone i divisions
    depends on number of cells that have undergone 1 division (n_1)
    """
    return (2*np.exp(-d*div_t))**(i-1) * n_1(t-(i-1)*div_t)
 
def dN_1_dt(cell_numbers, t, alpha, beta, div_t, d):
    """
    ODE for number of cells that made 1 division
    note that the 2* is not the same in de Boer et al (2006)
    """
    dn1_dt = (2*gamma_dist(t, alpha, beta)-
              gamma_dist(t-div_t, alpha, beta)*
               np.exp(-d*div_t)-
               d*cell_numbers)
    return dn1_dt


simulation_time_n1 = np.arange(0,250,0.01)
state = odeint(dN_1_dt, param.initial_cells, simulation_time_n1,
               args = (param.alpha,param.beta, param.div_t, param.d))

# n_1 is a function object!
n_1 = interp1d(simulation_time_n1, state, axis = 0, kind='cubic', fill_value = 0, bounds_error = False)

#plt.plot(simulation_time_n1, n_1(simulation_time_n1))
### calculate cells in N1 generation through integration

###


fig, ax = plt.subplots(1,1, figsize = (5,3.5))

for i in range(7):
    
    if i == 0:
        ax.plot(param.simulation_time, n_1(param.simulation_time), label = "1. Division")
    else:
        n_i = N_i(t=param.simulation_time, i = i, d= param.d, div_t = param.div_t,
              n_1 = n_1)
        ax.plot(param.simulation_time, n_i, label = str(i+1)+ ". Division")
        ax.legend()
    ax.set_xlabel("time (h)")
    ax.set_ylabel("cells")

plt.tight_layout()

