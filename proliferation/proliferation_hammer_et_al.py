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
linux_path = "home/burt/documents/"
os.chdir(windows_path)
import prolif_params as param

def gamma_dist(t, alpha, beta):
    return np.exp(-beta*t)*(beta**alpha)*(t**(alpha-1))/gamma(alpha)

def gamma_cdf(t, alpha, beta):
    dummy = t*beta
    return gammainc(alpha,dummy)

def N_0(t, non_dividing_cells,death_rate, alpha, beta):
    """
    analytical function N_0(t) for cell numbers that have undergone 0 divisions
    """
    cell0 = (1-gamma_cdf(t, alpha, beta))+non_dividing_cells*np.exp(-death_rate*t)
    return cell0

def N_i(t,state_idx, alpha, beta, death_rate, division_time, n_1):
    """
    analytical function N_i(t) for cell numbers that have undergone i divisions
    depends on number of cells that have undergone 1 division (n_1)
    """
    n1_indices =(t-(state_idx-1)*division_time)*100
    n1_indices = n1_indices.astype(int)
    
    return ((2*np.exp(-death_rate*division_time))**(state_idx-1)*n_1[n1_indices])
 
def dN_1_dt(cell_numbers, t, alpha, beta, division_time, death_rate):
    """
    ODE for number of cells that made 1 division
    note that the 2* is not the same in de Boer et al (2006)
    """
    dn1_dt = (2*gamma_dist(t, alpha, beta)-
              gamma_dist(t-division_time, alpha, beta)*
               np.exp(-death_rate*division_time)-
               death_rate*cell_numbers)
    return dn1_dt


### calculate cells in N1 generation through integration
simulation_time_n1 = np.arange(0,600,0.01)
state = odeint(dN_1_dt, param.initial_cells, simulation_time_n1, args = (param.alpha,
                                                      param.beta,
                                                      param.div_t, param.d))

###
fig, ax = plt.subplots(3,3, figsize = (12,12))
ax = ax.flatten()
for i in range(9):
    
    n_i = N_i(t=param.simulation_time, state_idx = i, alpha = param.alpha,
          beta = param.beta, death_rate= param.d, division_time = param.div_t,
          n_1 = state)
    ax[i].plot(param.simulation_time, n_i, "k", label = str(i)+ ". Division")
    ax[i].legend()
plt.tight_layout()
