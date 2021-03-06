#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Thu Sep 13 11:48:04 2018

@author: burt
parameters for conceptual model simulation with equal means and all
"""
import numpy as np
#==============================================================================
# model params
#==============================================================================
pos = 2.
neg = -0.5
neut = 0


feedback_menten = {
        "no_fb" : [[neut,neut,pos], [neut,neut,neut]],
        "pos_neg_fb" : [[pos,neg,pos], [neg,pos,neut]],
        "pos_th1" : [[pos,neut,pos], [neut,neut,neut]],
        "pos_th2" : [[neut,neut,pos], [neut,pos,neut]], 
        "neg_th1" : [[neut,neg,pos],[neut,neut,neut]],
        "neg_th2" : [[neut,neut,pos], [neg,neut,neut]],
        "pos_th2_neg_th2" : [[neut,neut,pos],[neg,pos,neut]]
        }

pos = 10.
neg = 0.1
neut = 1.


feedback_gamma = {
        "no_fb" : [[neut,neut,pos], [neut,neut,neut]],
        "pos_neg_fb" : [[pos,neg,pos], [neg,pos,neut]],
        "pos_th1" : [[pos,neut,pos], [neut,neut,neut]],
        "pos_th2" : [[neut,neut,pos], [neut,pos,neut]], 
        "neg_th1" : [[neut,neg,pos],[neut,neut,neut]],
        "neg_th2" : [[neut,neut,pos], [neg,neut,neut]],
        "pos_th2_neg_th2" : [[neut,neut,pos],[neg,pos,neut]]
        }
# extracellular il12 concentration
conc_il12 = 0

# half saturation constants
kd_ifn = 1.0
kd_il4 = 1.0
kd_il12 = 1.0

K = [kd_ifn, kd_il4, kd_il12]

#production rates cytokines
# these rates have little effect over wide range of parameters but sometimes,
# curves look more smooth if they are set higher, then probabilities go up faster
rate_ifn = 3*kd_ifn
rate_il4 = 3*kd_il4

### at some point I need to change this to cell densities
initial_cells = 1.0

#
mean_th1 = 1.
mean_th2 = 1.
alpha_th1 = 1
alpha_th2 = 1
beta_th1 = alpha_th1/mean_th1
beta_th2 = alpha_th2/mean_th2

alpha_th1_rtm = 20
alpha_th2_rtm = 20
beta_th1_rtm = 20.
beta_th2_rtm = 20.

degradation = 0
#==============================================================================
# simulation time
#==============================================================================
start = 0
stop = 10
# watch out, method chain also takes stepsize as independent argument
stepsize = .01
simulation_time = np.arange(start, stop, stepsize)

fb_start = 0
fb_end = stop
#==============================================================================
# stochastic parameters 
#==============================================================================
chain_length_stoc = 10
thn_idx = 0
th1_0_idx = 1
th2_0_idx = 2
th1_idx = 3
th2_idx = 4

precursor_rate = 1.

nsteps = stop*100
ncells = int(1000)
nsim = 10

stoc_hill_1 = [1,1,1]
stoc_hill_2 = [1,1,1]