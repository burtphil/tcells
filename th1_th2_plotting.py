#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Mon Sep 10 09:26:10 2018

@author: burt
"""
import matplotlib.pyplot as plt

######## plotting params
SMALL_SIZE = 15
MEDIUM_SIZE = 20

plt.rc('font', size=MEDIUM_SIZE)          # controls default text sizes
plt.rc('axes', titlesize=SMALL_SIZE)     # fontsize of the axes title
plt.rc('axes', labelsize=MEDIUM_SIZE)    # fontsize of the x and y labels
plt.rc('xtick', labelsize=SMALL_SIZE)    # fontsize of the tick labels
plt.rc('ytick', labelsize=SMALL_SIZE)    # fontsize of the tick labels
plt.rc('legend', fontsize=15)    # legend fontsize
plt.rc('figure', titlesize=MEDIUM_SIZE)

save_path = "/home/burt/Documents/figures/"
### plot time course th1 th2

def plot_time_course(th1_th2_model, parameters, save = "no", title = True, title_name = ""):
    
    alpha_1, alpha_2, rate1, rate2, simulation_time, conc_il12, hill_1, hill_2, rate_ifn, rate_il4, half_saturation, initial_cells = parameters
    norm = initial_cells/100
    th0_cells = th1_th2_model[:,0]/norm
    th1_cells = th1_th2_model[:,int(alpha_1)]/norm
    th2_cells = th1_th2_model[:,-1]/norm
    
    fig, ax = plt.subplots(1,1, figsize = (5,5))
    ax.plot(simulation_time,th0_cells, "k",
        simulation_time,th1_cells, "tab:blue",
        simulation_time,th2_cells, "tab:red")
    ax.set_yticks([0,50,100])
    ax.set_ylim([0,100])
    ax.set_xlim([simulation_time[0],simulation_time[-1]])
    ax.set_xlabel("Time [h]")
    ax.set_ylabel("% Th cells")
    
    if title == True:
        ax.set_title(r"$\alpha_1$ = "+str(alpha_1)+r" $\alpha_2$ = "+str(alpha_2)+" "+title_name)
    
    plt.legend(["Th0","Th1","Th2"], loc = "right")
    plt.tight_layout()
    
    if save != "no":
        fig.savefig(save_path+save+".svg",bbox_inches="tight", dpi=1200)
        

def plot_chain(chain_array, save = "no"):
    chain = chain_array[0]
    th1_conc = chain_array[1]
    th2_conc = chain_array[2]
    th1_tau = chain_array[3]
    th2_tau = chain_array[4]
    chain_length = chain_array[-1]
    
    fig, ax = plt.subplots(1,2, figsize = (10,4))
    ax[0].scatter(chain, th1_conc, c= "tab:blue")
    ax[0].scatter(chain, th2_conc, c= "tab:red")
    ax[0].set_xlabel(r"chain length $\alpha$")
    ax[0].set_ylabel("% Th cells after 100 hrs")
    ax[0].set_ylim([0,100])
    ax[0].set_xlim([0, chain_length])
    ax[0].set_yticks([0,50,100])
    
    ax[1].scatter(chain, th1_tau, c= "tab:blue")
    ax[1].scatter(chain, th2_tau, c= "tab:red")
    ax[1].set_xlabel(r"chain length $\alpha$")
    ax[1].set_ylabel(r"$\tau_{1/2}$")
    #ax[1].set_ylim([0,50])
    ax[1].set_xlim([0, chain_length])
    #ax[1].set_yticks([0,25,50])
    plt.tight_layout()
    
    if save != "no":
        fig.savefig(save_path+save+".svg", bbox_inches = "tight", dpi = 1200)

def subplot_chain(chain_array, ax, labels = "on"):
    chain = chain_array[0]
    th1_conc = chain_array[1]
    th2_conc = chain_array[2]
    #th1_tau = chain_array[3]
    #th2_tau = chain_array[4]
    chain_length = chain_array[-1]
    
    ax.scatter(chain, th1_conc, c= "tab:blue")
    ax.scatter(chain, th2_conc, c= "tab:red")
    #ax.set_xlabel(r"chain length $\alpha$")
    #ax.set_ylabel("% Th cells after 100 hrs")
    ax.set_ylim([0,100])
    ax.set_xlim([0, chain_length])
    if labels == "on":
        ax.set_yticks([0,50,100])
    else:
        ax.set_yticks([])
        ax.set_xticks([])
        
def subplot_tau(chain_array, ax):
    chain = chain_array[0]
    #th1_conc = chain_array[1]
    #th2_conc = chain_array[2]
    th1_tau = chain_array[3]
    th2_tau = chain_array[4]
    chain_length = chain_array[-1]
    
    ax.scatter(chain, th1_tau, c= "tab:blue")
    ax.scatter(chain, th2_tau, c= "tab:red")
    #ax.set_xlabel(r"chain length $\alpha$")
    #ax.set_ylabel(r"$\tau_{1/2}$")
    #ax[1].set_ylim([0,50])
    ax.set_xlim([0, chain_length])

def plot_il12(il12_array, xlabel = "IL-12 [pm]", factor_x_axis = 10**12, save = "no"):
    
    il12_conc = il12_array[0]
    il12_conc_pm = il12_conc*factor_x_axis
    th1_conc = il12_array[1]
    th2_conc = il12_array[2]
    
    fig, ax = plt.subplots(1,1,figsize=(5,5))
    ax.plot(il12_conc_pm, th1_conc, "tab:blue")
    ax.plot(il12_conc_pm, th2_conc, "tab:red")
    ax.set_xlabel(xlabel)
    ax.set_yticks([0,50,100])
    #ax.set_title("Effect of IL-12 conc. on Th cell balance after 3 hours")
    ax.set_ylim([0,100])
    ax.set_xlim([il12_conc_pm[0],il12_conc_pm[-1]])
    
    ax.set_ylabel("% Th cells")
    plt.legend(["Th1","Th2"])
    plt.tight_layout()  
    
    if save != "no":
        fig.savefig(save_path+save+".svg", bbox_inches = "tight", dpi = 1200)