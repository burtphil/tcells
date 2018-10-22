#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Mon Oct 22 11:58:27 2018

@author: burt
"""

import numpy as np
import seaborn as sns
import matplotlib.pyplot as plt
import os
os.chdir("/home/burt/Documents/code/th_cell_differentiation")
import modules.th1_th2_conceptual_parameters as cparams
from scipy.special import gammainc

from stochasitic_simulation import semi_markov_simulation, semi_markov_simulation2
def gamma_cdf(t, alpha, beta):
    dummy = t*beta
    return gammainc(alpha,dummy)   

def survival_fct(t, alpha, beta):
    return 1-gamma_cdf(t,alpha,beta)
#==============================================================================
# import parameters
#==============================================================================
model_params = cparams.parameters

#==============================================================================
# perform and plot multiple simulations
#==============================================================================
fig, axes = plt.subplots(1,2,figsize = (5,5))

axes = axes.flatten()
ax = axes[0]
ax1 = axes[1]

for i in range(1,4):
    cells, time = semi_markov_simulation2(start = cparams.start, 
                                         stop = cparams.stop, 
                                         nsteps = cparams.nsteps, 
                                         ncells = cparams.ncells,
                                         alpha = i,
                                         beta = i)
    all_cells = cells[:,:,0]

    #naive_cells = []
    th1_cells = []

    for t in range(len(time)):
        x = all_cells[:,t]
        #naive_cells.append(len(x[x==0]))
        th1_cells.append(len(x[x==1]))
    

    ax.plot(time,th1_cells, label = r"$\alpha=$"+str(i))
    ax1.plot(time, gamma_cdf(time, i, i), label = r"$\alpha=$"+str(i))

ax.set_xlabel("time")
ax.set_ylabel("cells")
ax.set_xlim([cparams.start,cparams.stop])
ax.legend()
plt.tight_layout()

"""
def semi_markov_simulation3(start=0, stop=2, stepsize=0.0001, alpha = 1, beta = 1):

    cell=0
    a = np.random.rand()
    while start < stop:
        t = start
        if cell==0:           
            if a > survival_fct(t, alpha, beta):                
                cell = 1
                print t

        start = start+stepsize    
        
semi_markov_simulation3()
"""