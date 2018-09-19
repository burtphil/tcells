#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Tue Sep 11 15:06:24 2018

@author: burt
proliferation model implemented according to de Boer et al (2006) eq. (13)
"""
import numpy as np
from scipy.special import gamma, gammainc
from scipy.integrate import odeint
import matplotlib.pyplot as plt

def gamma_dist(t, alpha, beta):
    return np.exp(-beta*t)*(beta**alpha)*(t**(alpha-1))/gamma(alpha)

def gamma_cdf(t, alpha, beta):
    dummy = t*beta
    return gammainc(alpha,dummy)

#plt.plot(simulation_time,state[:,1])

def N_0(t, non_dividing_cells,death_rate, alpha, beta):
    """
    analytical function N_0(t) for cell numbers that have undergone 0 divisions
    """
    cell0 = (1-gamma_cdf(t, alpha, beta))+non_dividing_cells*np.exp(-death_rate*t)
    return cell0

def dN_1_dt(cell_numbers, t, alpha, beta, division_time, death_rate):
    dn1_dt = (gamma_dist(t, alpha, beta)-
              gamma_dist(t-division_time, alpha, beta)*
               np.exp(-death_rate*division_time)-
               death_rate*cell_numbers)
    return dn1_dt
   
def cells_per_state(state_idx, cell_numbers, t, alpha, beta, death_rate, division_time):
    """
    ODE for dN_i/dt from de Boer et al 2006 eq. (13)
    """    
    dummy = (2**(state_idx-1)*gamma_dist(t-(state_idx-1)*division_time, alpha, beta)*
             np.exp(-death_rate*division_time)-
             2**(state_idx-1)*gamma_dist(t-state_idx*division_time, alpha, beta)*
             np.exp(-death_rate*division_time)-
             death_rate*cell_numbers[state_idx])
    
    return dummy

def proliferate(cell_numbers, t, alpha, beta, death_rate, division_time, non_dividing_cells):
    """
    proliferation model
    """        
    dn0_dt = (-gamma_dist(t, alpha, beta)-
              death_rate*non_dividing_cells*np.exp(-death_rate*t))
    
    dn1_dt = (gamma_dist(t, alpha, beta)-
              gamma_dist(t-division_time, alpha, beta)*
               np.exp(-death_rate*division_time)-
               death_rate*cell_numbers[1])
    
    cells = [cells_per_state(state_idx, cell_numbers, t, alpha, beta, 
                             death_rate, division_time) for state_idx in range(2,10)]
    
    all_cells = np.concatenate(([dn0_dt],[dn1_dt],cells))
    return all_cells

### set up params
initial_cells = 1.
death_rate = 1/24
division_time = 10.
non_dividing_cells = 0.05/np.exp(-death_rate*(6*24))

start = 0
stop = 100
stepsize = 0.1

simulation_time = np.arange(start, stop, stepsize)
cell_states = [initial_cells, 0]

alpha_first_division = 10.
mean_first_division = 50.
beta_first_division = alpha_first_division/ mean_first_division

#plt.plot(simulation_time, gamma_dist(simulation_time, alpha_first_division, beta_first_division))

cell_states = np.zeros(10)
cell_states[0] = 1

state = odeint(proliferate, cell_states, simulation_time,
               args =(alpha_first_division, beta_first_division, death_rate, 
                      division_time, non_dividing_cells,))


fig, ax = plt.subplots(3,4, figsize = (12,9))
ax = ax.flatten()
for i in range(len(cell_states)):
    ax[i].plot(simulation_time, state[:,i], label = str(i)+ ". Division")
    ax[i].legend()
plt.tight_layout()
