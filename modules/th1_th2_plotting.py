#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Mon Sep 10 09:26:10 2018

@author: burt
"""
import numpy as np

def ax_time_course(state, ax, simulation_time, initial_cells, alpha_1, linestyle = "-"):
    
    norm = initial_cells / 100
    th1_cells = state[:,alpha_1+1] / norm
    th2_cells = state[:,-1] / norm
    
    #print th1_cells[-1],th2_cells[-1]
    #ax.plot(simulation_time, th0_cells, color = "k", linestyle = linestyle)
    ax.plot(simulation_time, th1_cells, color = "tab:blue", linestyle = linestyle)
    ax.plot(simulation_time, th2_cells, color = "tab:red", linestyle = linestyle)
    ax.set_yticks([0, 50, 100])
    ax.set_ylim([0, 100])
    ax.set_xlim([simulation_time[0], simulation_time[-1]])
    ax.set_xlabel("time")
    ax.set_ylabel("% Th cells")
       
def ax_var_effect(x_var, y_val, ax, linestyle = "-", plot_both = True):
    ax.plot(x_var, y_val[0], "tab:blue", linestyle = linestyle)    
    if plot_both == True:
        ax.plot(x_var, y_val[1], "tab:red", linestyle = linestyle)
    
    ax.set_xlim(x_var[0], x_var[-1])
    
def ax_chain(x_var, y_val, ax, xsteps = 2, linestyle = "-", plot_both = True):
    ax.scatter(x_var, y_val[0], c = "tab:blue")    
    if plot_both == True:
        ax.scatter(x_var, y_val[1], c = "tab:red")

    ax.set_xlabel(r"chain length $\alpha$")
    ax.set_xlim([0, x_var[-1]+1])
    ax.set_xticks(np.arange(0, x_var[-1] + xsteps, xsteps))  
    