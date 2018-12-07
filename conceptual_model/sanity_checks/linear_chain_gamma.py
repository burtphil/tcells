#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Thu Sep 27 17:47:00 2018

@author: burt
compare ode linear chain and gamma distribution
"""
import numpy as np
from scipy.special import gamma, gammainc
from scipy.integrate import odeint
import matplotlib.pyplot as plt


def gamma_cdf(t, alpha, beta):
    dummy = t*beta
    return gammainc(alpha,dummy)


def th_cell_diff(state, t, alpha, beta):
    
    rate = beta
    dt_state = np.zeros_like(state)
    
    for j in range(len(state)):
        if j == 0:
            dt_state[j] = -rate*state[j]
        elif j != (len(state)-1):
            dt_state[j] = rate*(state[j-1]-state[j])
        else:
            dt_state[j] = rate*state[j-1]

    return dt_state

start=0
stop=5
stepsize=0.001

initial_cells=1.
t=np.arange(start,stop,stepsize)


fig, axes = plt.subplots(1,2, figsize = (8,4))

ax = axes[0]
ax1 = axes[1]

for i in range(1,5):
    
    #vary alpha and set beta = alpha to keep mean constant
    alpha = i
    beta = i
    
    #initialize states depending on the chain length (1 extra state for naive cells)    
    no_of_states = int(alpha+1)
    y0 = np.zeros(no_of_states)
    y0[0] = initial_cells

    state = odeint(th_cell_diff, y0, t, args=(alpha,beta))
    
    #index th0 and th1 cells (omit intermediary states)
    th0_cells = state[:,0]
    th1_cells = state[:,-1]
    #ax[0].plot(t, th0_cells)
    ax.plot(t, th1_cells, label = r"$\alpha=$"+str(i))
    ax.set_yticks([0,0.5,1.])
    ax1.set_yticks([0,0.5,1.])
    ax1.set_xticks([0,1,2,3,4,5])
    ax.set_xticks([0,1,2,3,4,5])
    ax.set_ylim([0,1])
    ax.set_xlim([t[0],t[-1]])
    ax1.set_ylim([0,1])
    ax1.set_xlim([t[0],t[-1]])
    ax.set_xlabel("t (h)")
    ax.set_ylabel("% Th cells")
    ax.legend()
    ax.set_title(r"linear chain ODE, $\lambda_i=\beta$")
    
    ax1.plot(t,gamma_cdf(t, alpha, beta))
    ax1.set_xlabel("t (h)")
    ax1.set_ylabel(r"$\gamma(t,\alpha,\beta)$")
    ax1.set_title(r"$\alpha=\beta$")

plt.tight_layout()