#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Wed Oct 17 17:10:28 2018

@author: burt
"""

import numpy as np
import matplotlib.pyplot as plt
from scipy.special import gammainc

def gamma_cdf(t, alpha, beta):
    dummy = t*beta
    return gammainc(alpha,dummy)      

def survival_fct(t, alpha, beta):
    return 1-gamma_cdf(t,alpha,beta)

def draw_fate(probs):     
    if np.random.rand() < probs[0]:
        return 1
    else:
        return 2

def prob_fate():
    dummy = [0.2,0.8]
    return dummy

def draw_time():
    time = np.random.exponential(1.)
    return time

def cell_fate():
    probs = prob_fate()
    fate = draw_fate(probs)
    time = draw_time()
    
    return [fate,time]

# set timestep
#time = 0
#timestep = .01
#simulation_time = 2

# cell vector, first entry is state, second is time in state, p is prob to change state


### assign states 
# 0 represents Thn
# 1 represents Th1_0
# 2 represents Th2_0
# 3 represents Th1
# 4 represents Th2

# simulation for one cell
def run_simulation(start = 0, stop = 0.5, nsteps = 100, nstates = 3,ncells = 100):
    
    cells = np.zeros((ncells,nsteps,nstates))
    
    
    time = np.linspace(start,stop,nsteps)
    
    for cell_j in cells:    
        for i in range(len(time)-1):
            t = time[i]
            cell = cell_j[i]
            # is there a cell fate switch?
            if cell[0] == 0 and np.random.exponential(1.) < t:
                # calculate probabilities which fate to choose
                probs = prob_fate()
                # assign new cell fate
                cell[0] = draw_fate(probs)
                cell[1] = t
                
            if cell[0] == 1 and t > cell[1]:
                if np.random.rand() > cell[2]:
                    cell[2] = cell[2] + gamma_cdf(t-cell[1],1,1)
                else:
                    cell[0] = 3
                    
            if cell[0] == 2 and t > cell[1]:
                if np.random.rand() > cell[2]:
                    cell[2] = cell[2] + gamma_cdf(t-cell[1],1,1)
                else:
                    cell[0] = 4
            
            cell_j[i+1] = cell
            
    return [cells, time]   

def run_simulation2(start = 0, stop = 0.5, nsteps = 100, nstates = 3,ncells = 100):
    
    cells = np.zeros((ncells,nsteps,nstates))    
    time = np.linspace(start,stop,nsteps)
    
    
    for i in range(len(time)-1):
        t = time[i]
        cell_j=cells[i]
        for cell in cell_j:
            # is there a cell fate switch?
            if cell[0] == 0 and np.random.exponential(1.) < t:
                # calculate probabilities which fate to choose
                probs = prob_fate()
                # assign new cell fate
                cell[0] = draw_fate(probs)
                cell[1] = t
                
            if cell[0] == 1 and t > cell[1]:
                if np.random.rand() > survival_fct(t-cell[1],1,1):
                    cell[0] = 3
                    
            if cell[0] == 2 and t > cell[1]:
                if np.random.rand() > survival_fct(t-cell[1],1,1):
                    cell[0] = 4
            
            cell_j[i+1] = cell
            
    return [cells, time]          
    #plt.plot(time, cell_j[:,0])
#
fig, ax = plt.subplots(1,1,figsize = (8,5))

for i in range(100):
    cells, time = run_simulation()
    all_cells = cells[:,:,0]

    naive_cells = []
    th1_cells = []
    th2_cells = []

    for t in range(len(time)):
        x = all_cells[:,t]
        naive_cells.append(len(x[x==0]))
        th1_cells.append(len(x[x==3]))
        th2_cells.append(len(x[x==4]))

    
    if i == 0:
        ax.plot(time,naive_cells, "k", label = "Thn")
        ax.plot(time,th1_cells, "tab:blue", label = "Th1")
        ax.plot(time,th2_cells, "r", label = "Th2")
    else:
        ax.plot(time,naive_cells, "k")
        ax.plot(time,th1_cells, "tab:blue")
        ax.plot(time,th2_cells, "r")

ax.set_xlabel("time")
ax.set_ylabel("cells")
ax.legend()
"""
cells = np.zeros((3,10))

for cell in cells:
    
    while time < simulation_time:

        # is there a cell fate switch?
        if cell[0] == 0 and np.random.exponential(1.) < time:
            # calculate probabilities which fate to choose
            probs = prob_fate()
            # assign new cell fate
            cell[0] = draw_fate(probs)
            cell[1] = time
            
        if cell[0] == 1 and time > cell[1]:
            if np.random.rand() > cell[2]:
                cell[2] = cell[2] + gamma_cdf(time-cell[1],1,1)
            else:
                cell[0] = 3
                
        if cell[0] == 2 and time > cell[1]:
            if np.random.rand() > cell[2]:
                cell[2] = cell[2] + gamma_cdf(time-cell[1],1,1)
            else:
                cell[0] = 4
        

             
        time = time + timestep

naive_cells = cells[cells[:,0] == 0]
th1_prec = cells[cells[:,0] == 1]
th2_prec = cells[cells[:,0] == 2]
th1_cells = cells[cells[:,0] == 3]
th2_cells = cells[cells[:,0] == 4]


while time < simulation_time:
    
    naive_cells = cells[cells[:,0] == 0]
    th1_prec = cells[cells[:,0] == 1]
    th2_prec = cells[cells[:,0] == 2]
    th1_cells = cells[cells[:,0] == 3]
    th2_cells = cells[cells[:,0] == 4]

    
    # is there a cell fate switch?
    if cells[:,0] == 0 and np.random.exponential(1.) < time:
        # calculate probabilities which fate to choose
        probs = prob_fate()
        # assign new cell fate
        cell[0] = draw_fate(probs)
        cell[1] = time
        
    if cell[0] == 1 and time > cell[1]:
        if np.random.rand() > cell[2]:
            cell[2] = cell[2] + gamma_cdf(time-cell[1],1,1)
        else:
            cell[0] = 3
            
    if cell[0] == 2 and time > cell[1]:
        if np.random.rand() > cell[2]:
            cell[2] = cell[2] + gamma_cdf(time-cell[1],1,1)
        else:
            cell[0] = 4
    
          
    time = time + timestep
"""