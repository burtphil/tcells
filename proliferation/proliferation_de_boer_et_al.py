#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Tue Sep 11 15:06:24 2018

@author: burt
proliferation model implemented according to de Boer et al (2006) eq. (13)
"""
import os
import numpy as np
from scipy.special import gamma, gammainc
from scipy.integrate import odeint
import matplotlib.pyplot as plt

windows_path= "C:/Users/Philipp/Documents/tcells/modules"
linux_path = "home/burt/documents/"
os.chdir(windows_path)
import prolif_params as param

def gamma_dist(t, alpha, beta, factor = 1):
    dummy= factor*np.exp(-beta*t)*(beta**alpha)*(t**(alpha-1))/gamma(alpha)
    
    return dummy

def gamma_cdf(t, alpha, beta):
    dummy = t*beta
    return gammainc(alpha,dummy)

#plt.plot(simulation_time,state[:,1])
   
def cells_per_state(state_idx, cell_numbers, t, alpha, beta, death_rate, division_time):
    """
    ODE for dN_i/dt from de Boer et al 2006 eq. (13)
    """    
    dummy = (2**(state_idx-1)*gamma_dist(t-(state_idx-1)*division_time, alpha, beta)*
             np.exp(-(state_idx-1)*death_rate*division_time)-
             2**(state_idx-1)*gamma_dist(t-state_idx*division_time, alpha, beta)*
             np.exp(-state_idx*death_rate*division_time)-
             death_rate*cell_numbers[state_idx])
    
    return dummy

def proliferate(cell_numbers, t, alpha, beta, death_rate, division_time, non_dividing_cells, n):
    """
    proliferation model
    """        
    dn0_dt = (-gamma_dist(t, alpha, beta)-
              death_rate*non_dividing_cells*cell_numbers[0])
    
    dn1_dt = (gamma_dist(t, alpha, beta)-
              gamma_dist(t-division_time, alpha, beta)*
               np.exp(-death_rate*division_time)-
               death_rate*cell_numbers[1])
    
    cells = [cells_per_state(state_idx, cell_numbers, t, alpha, beta, 
                             death_rate, division_time) for state_idx in range(2,n)]
    
    all_cells = np.concatenate(([dn0_dt],[dn1_dt],cells))
    return all_cells

def proliferate2(cell_numbers, t, alpha, beta, death_rate, division_time):
    """
    proliferation model
    """        
    dn1_dt = (gamma_dist(t, alpha, beta)-
              gamma_dist(t-division_time, alpha, beta)*
               np.exp(-death_rate*division_time)-
               death_rate*cell_numbers[0])
    
    dn2_dt = ((2*gamma_dist(t-division_time, alpha, beta)*np.exp(-death_rate*division_time))-
              (2*gamma_dist(t-(2*division_time), alpha, beta)* np.exp(-2*death_rate*division_time))-
              (death_rate*cell_numbers[1]))
    
    return [dn1_dt,dn2_dt]


state = odeint(proliferate2, [0.,0.], param.simulation_time,
               args =(param.alpha, param.beta, param.d, 
                      param.div_t,))

fig, axes = plt.subplots(1,2, figsize = (10,5))
axes[0].plot(param.simulation_time, state[:,0])
axes[1].plot(param.simulation_time, state[:,1])

def first_gen(x, t, alpha, beta, death_rate, division_time, non_dividing_cells):
    dn1_dt = (gamma_dist(t, alpha, beta)-
              gamma_dist(t-division_time, alpha, beta)*
               np.exp(-death_rate*division_time)-
               death_rate*x)
    
    return dn1_dt



x0 = 0
t = np.linspace(0,100,10000)

#plt.plot(t,gamma_dist(t, param.alpha_first_division, param.beta_first_division, factor = 10000))
#state = odeint(first_gen, x0, t, args = (param.alpha_first_division, param.beta_first_division, param.death_rate, 
#                      param.division_time, param.non_dividing_cells,))

#plt.plot(t,state)
"""
n = 5
cell_states = np.zeros(n)

state = odeint(proliferate, cell_states, param.simulation_time,
               args =(param.alpha_first_division, param.beta_first_division, param.death_rate, 
                      param.division_time, param.non_dividing_cells,n,))


fig, ax = plt.subplots(1,len(cell_states), figsize = (15,5))
ax = ax.flatten()
for i in range(len(cell_states)):
    ax[i].plot(param.simulation_time, state[:,i], label = str(i)+ ". Division")
    ax[i].legend()
plt.tight_layout()
"""