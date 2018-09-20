#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Wed Sep 12 11:12:39 2018

@author: burt
implementation of Hammer et al prolif model with derivative of
analytical solution embedded in ODE system
"""
import numpy as np
from scipy.special import gamma, gammainc
from scipy.integrate import odeint
import matplotlib.pyplot as plt
import prolif_params as param

def gamma_dist(t, alpha, beta):
    return np.exp(-beta*t)*(beta**alpha)*(t**(alpha-1))/gamma(alpha)

def gamma_cdf(t, alpha, beta):
    dummy = t*beta
    return gammainc(alpha,dummy)

def dN_1_dt(cell_numbers, t, alpha, beta, division_time, death_rate):
    """
    derivative of cells that have undergonge 1 division
    note that here hammer et al differ from de boer with the coefficient 2 before first gamma dist
    but it does not effect outcome
    """
    dn1_dt = (gamma_dist(t, alpha, beta)-
              gamma_dist(t-division_time, alpha, beta)*
               np.exp(-death_rate*division_time)-
               death_rate*cell_numbers[1])
    return dn1_dt

def cells_per_state(state_idx, cell_numbers, t, alpha, beta, death_rate, division_time):
    """
    ODE for dN_i/dt as derivation from analytical solution from Hammer et al
    """
    t = t-(state_idx-1)*division_time    
    dummy = (2*np.exp(-death_rate*division_time)**(state_idx-1)*
             dN_1_dt(cell_numbers, t, alpha, beta, division_time, death_rate))
    
    return dummy

def proliferate2(cell_numbers, t, alpha, beta, death_rate, division_time, non_dividing_cells):
    """
    proliferation model
    """        
    dn0_dt = (-gamma_dist(t, alpha, beta)-
              death_rate*non_dividing_cells*np.exp(-death_rate*t))
    
    dn1_dt = dN_1_dt(cell_numbers, t, alpha, beta, division_time, death_rate)
    cells = [cells_per_state(state_idx, cell_numbers, t, alpha, beta, 
                             death_rate, division_time) for state_idx in range(2,10)]
    
    all_cells = np.concatenate(([dn0_dt],[dn1_dt],cells))
    return all_cells


plt.plot(param.simulation_time, gamma_dist(param.simulation_time, param.alpha_first_division, param.beta_first_division))

cell_states = np.zeros(10)
cell_states[0] = 1

state = odeint(proliferate2, cell_states, param.simulation_time,
               args =(param.alpha_first_division, param.beta_first_division, param.death_rate, 
                      param.division_time, param.non_dividing_cells,))


fig, ax = plt.subplots(3,4, figsize = (12,9))
ax = ax.flatten()
for i in range(len(cell_states)):
    ax[i].plot(param.simulation_time, state[:,i], label = str(i)+ ". Division")
    ax[i].legend()
plt.tight_layout()
