#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Fri Oct 19 14:45:18 2018

@author: burt
compare stochastic and ode model simulation
"""
import numpy as np
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
import os
os.chdir("/home/burt/Documents/code/th_cell_differentiation")
import modules.th1_th2_conceptual_parameters as cparams
from modules.th1_th2_ode_model_generic import run_model
from stochasitic_simulation import run_stochastic_simulation
from modules.th1_th2_plotting import ax_time_course
#==============================================================================
# import parameters
#==============================================================================
model_params = cparams.parameters

#==============================================================================
# run time course simulation
#==============================================================================
simulation_name = "sim"

th1_th2_model = run_model(simulation_name, parameters = model_params)

test_simulation = np.load(simulation_name+".npz")
state = test_simulation["state"]
parameters = test_simulation["parameters"]


# plot time course
#plot_time_course(th1_th2_model, parameters)

#==============================================================================
# perform and plot multiple simulations
#==============================================================================
fig, ax = plt.subplots(1,2,figsize = (12,5))

n_simulations = cparams.nsim
for i in range(n_simulations):
    cells, time = run_stochastic_simulation(start = cparams.start, stop = cparams.stop, nsteps = cparams.nsteps, ncells = cparams.ncells)
    all_cells = cells[:,:,0]

    naive_cells = []
    th1_cells = []
    th2_cells = []

    for t in range(len(time)):
        x = all_cells[:,t]
        naive_cells.append(len(x[x==cparams.thn_idx]))
        th1_cells.append(len(x[x==cparams.th1_idx]))
        th2_cells.append(len(x[x==cparams.th2_idx]))
    
    if i == 0:
        ax[0].plot(time,naive_cells, "k", label = "Thn")
        ax[0].plot(time,th1_cells, "tab:blue", label = "Th1")
        ax[0].plot(time,th2_cells, "r", label = "Th2")
    else:
        ax[0].plot(time,naive_cells, "k")
        ax[0].plot(time,th1_cells, "tab:blue")
        ax[0].plot(time,th2_cells, "r")

ax[0].set_xlabel("time")
ax[0].set_ylabel("cells")
ax[0].legend()
ax[0].set_xlim([cparams.start,cparams.stop])


ax_time_course(state = state, ax = ax[1], parameters = model_params)
ax[1].set_xlim([cparams.start,cparams.stop])
ax[0].set_yticks([0,25,50,75,100])
ax[1].set_yticks([0,25,50,75,100])
fig.suptitle(r"$\alpha_1=1,\alpha_2=1$, all parameters set equal for ODE and stochastic simulation", fontsize = 12)
plt.tight_layout()

cells, time = run_stochastic_simulation(start = cparams.start, stop = cparams.stop, nsteps = cparams.nsteps, ncells = cparams.ncells)
all_cells = cells[:,:,0]

df = pd.DataFrame(data = all_cells)
df = df.T
df["time"] = time

df = pd.melt(df, id_vars = ["time"], var_name = "cell_no", value_name = "cell_type")

df2 = df[["time","cell_type"]]
df2 =df2.groupby(['time', 'cell_type']).size().reset_index(name='counts')


naive_cells = df2[df2["cell_type"]==0]
th1_cells = df2[df2["cell_type"]==3]
th2_cells = df2[df2["cell_type"]==4]

new_df = naive_cells
new_df = new_df.append(th1_cells)
new_df = new_df.append(th2_cells)

sns.relplot(x="time", y="counts", hue = "cell_type", kind="line", data=new_df)