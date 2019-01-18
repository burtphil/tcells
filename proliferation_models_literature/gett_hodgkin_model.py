#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Thu Sep 27 17:47:00 2018

@author: burt
hodkin model reproduction of de Boer et al 2006
"""
import numpy as np
from scipy.integrate import odeint
import matplotlib.pyplot as plt
from scipy.special import gamma
from scipy.interpolate import interp1d
import seaborn as sns


#==============================================================================
# choose plotting style and saving path
#==============================================================================
sns.set(context = "talk", style = "ticks")
save_path = "/home/burt/Documents/tcell_project/figures/"
#==============================================================================
# import parameters
#==============================================================================
mu = 60.32
sigma = 0.27
C = 10408
delta = -1.1
t_div = 10.39
rate_death = 0.008

mu = 30.9
sigma = 0.68
C = 7920
delta = 35.8
t_div = 9.49
rate_death = 0.09

time = np.linspace(0,110,1000)
#==============================================================================
# functions
#==============================================================================
def H(t):
    dum = 1
    if t < 0 : 
        dum = 0
    return dum

def prec_dist(t, C, sigma, mu, delta):
    
    if t > delta:
        state = (C / (np.sqrt(2*np.pi)*sigma*(t-delta)))*np.exp(-((np.log(t-delta)-np.log(mu))**2) / (2*sigma*sigma))
    else:
        state = 0
        
    return state

def N_i(time, i, d, div_t, n_1):
    """
    analytical function N_i(t) for cell numbers that have undergone i divisions
    depends on number of cells that have undergone 1 division (n_1)
    """
    scale = ((2*np.exp(-d*div_t))**(i-1))
    state = [scale * n_1(t-((i-1)*div_t))*H(t-((i-1)*div_t)) for t in time]
    return state

cells = [prec_dist(t, C, sigma, mu, delta) for t in time]

fig, ax = plt.subplots()
ax.plot(time, cells, label = "R(t)")
ax.legend()
#ax.set_ylim(0,270)
#ax.set_yticks(np.arange(0,300,50))
#ax.set_xticks(np.arange(0,120,20))
ax.set_xlim(0,time[-1])
plt.tight_layout()
  
def th_cell_diff(state, t, C, sigma, mu, delta, t_div, rate_death):
    
    dt_state = (prec_dist(t, C, sigma, mu, delta)
    - (prec_dist(t-t_div, C, sigma, mu, delta) * np.exp(-rate_death * t_div)*H(t-t_div))
    - rate_death * state)
                
    dt_state = dt_state
    
    return dt_state  

y0 = 0
state = odeint(th_cell_diff, y0, time, args = (C, sigma, mu, delta, t_div, rate_death))
#state = state / 1000

#fig, ax = plt.subplots()
#ax.plot(time, state)
#plt.ylim(0,8)
#plt.tight_layout()


# n_1 is a function object!
n_1 = interp1d(time, state, axis = 0, kind='cubic', fill_value = 0, bounds_error = False)

#plt.plot(simulation_time_n1, n_1(simulation_time_n1))
### calculate cells in N1 generation through integration

###
colors = ["k","tab:red", "tab:green", "tab:blue"]

fig, ax = plt.subplots(1,1, figsize = (5,3.5))

for i in range(1,4):
    
    if i == 1:
        ax.plot(time, n_1(time) / 1000, label = "n = 1", c = "k")
    else:
        n_i = N_i(time=time, i = i, d= rate_death, div_t = t_div,
              n_1 = n_1)
        ax.plot(time, np.asarray(n_i) / 1000., label = "n = "+str(i), c = colors[i-1])
        ax.legend()
    ax.set_xlabel("time (h)")
    ax.set_ylabel("cells")

#ax.set_yticks(np.arange(0,10,1))
#ax.set_xticks([0,20,40,60,80, 100])
#ax.set_ylim(0,10)
plt.tight_layout()
