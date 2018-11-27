#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Mon Sep 10 09:26:10 2018

@author: burt
"""
import matplotlib.pyplot as plt
import numpy as np
######## plotting params
SMALL_SIZE = 15
MEDIUM_SIZE = 20

save_path = "/home/burt/Documents/figures/"
### plot time course th1 th2
        
def ax_time_course(state, ax, simulation_time, initial_cells, alpha_1, linestyle = "-"):
    
    norm = initial_cells / 100
    th0_cells = state[:,0] / norm
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
       
def ax_chain(chain_array, ax, output = "steady_state", steps_x_ticks = 2):
    """
    takes an array produced by function chain
    can produce a plot of chain length vs cell steady state (choose option output = steady_state)
    alternatively, use output = tau to plot chain length vs tau1/2
    returns axes object
    """
    
    th1_conc = chain_array[0]
    th2_conc = chain_array[1]
    th1_tau = chain_array[2]
    th2_tau = chain_array[3]
    chain_length = chain_array[4]
    chain = [i for i in range(1, chain_length)]
    
    if output == "steady_state":
        scatter_output = [th1_conc, th2_conc]
        ylabel = "% Th cells"
    else:
        scatter_output = [th1_tau, th2_tau]
        ylabel = r"$\tau_{1/2}$"
        
    ax.scatter(chain, scatter_output[0], c = "tab:blue")
    ax.scatter(chain, scatter_output[1], c = "tab:red")
    ax.set_xlabel(r"chain length $\alpha$")
    ax.set_ylabel(ylabel)
    ax.set_xlim([0, chain_length])
    ax.set_xticks(np.arange(0, chain_length + steps_x_ticks, steps_x_ticks))

def ax_il12(il12_array, ax, plot_both = True, output = "steady_state", factor_x_axis = 1.0, linestyle = "-"):
    """
    return axes object for il12 dependency
    """
    il12_conc = il12_array[0]
    il12_conc_pm = il12_conc*factor_x_axis
    th1_conc = il12_array[1]
    th2_conc = il12_array[2]
    th1_tau = il12_array[3]
    th2_tau = il12_array[4]
    
    if output == "steady_state":
        yval = [th1_conc, th2_conc]
        ylabel = "% Th cells"
        ax.set_ylim([0,100])
        ax.set_yticks([0,50,100])
    else:
        yval = [th1_tau, th2_tau]
        ylabel = r"$\tau_{1/2}$"
        
    ax.plot(il12_conc_pm, yval[0], "tab:blue", linestyle = linestyle)
    if plot_both == True:
        ax.plot(il12_conc_pm, yval[1], "tab:red", linestyle = linestyle)
    ax.set_ylabel(ylabel)
    
    #ax.set_title("Effect of IL-12 conc. on Th cell balance after 3 hours")
    
    ax.set_xlim([il12_conc_pm[0],il12_conc_pm[-1]])
    
    #ax.set_ylabel("% Th cells after")
def ax_feedback_tau(cytokine_conc_arr, ax, xlabel = "IFNg [a.u.]", ylabel = r"$\tau_{1/2}$", linestyle = "-"):
    x_values = cytokine_conc_arr[0]
    #th1_conc = chain_array[1]
    #th2_conc = chain_array[2]
    th1_tau = cytokine_conc_arr[3]
    th2_tau = cytokine_conc_arr[4]   
    ax.plot(x_values, th1_tau, "tab:blue", linestyle = linestyle)
    ax.plot(x_values, th2_tau, "tab:red", linestyle = linestyle)
    ax.set_xlabel(xlabel)
    ax.set_ylabel(ylabel)
    #ax[1].set_ylim([0,50])
    ax.set_xlim([x_values[0], x_values[-1]])
