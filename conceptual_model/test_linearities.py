#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Mon Dec  3 16:35:53 2018

@author: burt
"""

#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Fri Nov 30 11:26:19 2018

@author: burt
phd summary report from november 30.th 2018
should be able to reproduce all plots shown in the report together with the
scripts th1_th2_ode_model_generic.py
th1_th2_plotting.py
gamma_dist.py
only exception: stochastic simulations, for that, look into 
"""

#==============================================================================
# import modules
#==============================================================================
import numpy as np
import seaborn as sns
import matplotlib.pyplot as plt
import os
os.chdir("/home/burt/Documents/code/th_cell_differentiation")
import modules.th1_th2_conceptual_parameters as cparams
from modules.th1_th2_ode_model_generic import run_model
from modules.th1_th2_ode_model_generic import variable_effect
from modules.th1_th2_ode_model_generic import get_cell_arr
from modules.th1_th2_ode_model_generic import get_cyto_conc
from modules.th1_th2_ode_model_generic import get_prob
from modules.th1_th2_plotting import ax_time_course
from modules.th1_th2_plotting import ax_var_effect
from modules.th1_th2_plotting import ax_chain

#==============================================================================
# choose plotting style and saving path
#==============================================================================
sns.set(context = "talk", style = "ticks")
save_path = "/home/burt/Documents/tcell_project/figures/"

#==============================================================================
# choose type of feedback
#==============================================================================
feedback_dict = cparams.feedback_menten
feedback_type = "pos_th1"

#==============================================================================
# choose feedback duration
#==============================================================================
fb_start = 0
fb_end = cparams.stop
assert fb_end <= cparams.stop

#==============================================================================
# import parameters
#==============================================================================
alpha_1 = cparams.alpha_th1
alpha_2 = cparams.alpha_th2
beta_1 = cparams.beta_th2
beta_2 = cparams.beta_th2
K = cparams.K
degradation = cparams.degradation
rate_ifn = cparams.rate_ifn
rate_il4 = cparams.rate_il4
conc_il12 = cparams.conc_il12
simulation_time = cparams.simulation_time
initial_cells = cparams.initial_cells
fb_strength = feedback_dict[feedback_type]
hill_1 = [1,1,1]
hill_2 = hill_1

rate_params = [alpha_1, alpha_2, beta_1, beta_2, simulation_time, conc_il12,
              hill_1, hill_2, fb_strength, rate_ifn, rate_il4, K,
              initial_cells, degradation, fb_start, fb_end]

rtm_params = [20, 20, 20, 20, simulation_time, conc_il12,
              hill_1, hill_2, fb_strength, rate_ifn, rate_il4, K,
              initial_cells, degradation, fb_start, fb_end]

hill_1 = [5,5,5]
hill_2 = hill_1

hill5_params = [alpha_1, alpha_2, beta_1, beta_2, simulation_time, conc_il12,
              hill_1, hill_2, fb_strength, rate_ifn, rate_il4, K,
              initial_cells, degradation, fb_start, fb_end]

hill5_rtm_params = [20, 20, 20, 20, simulation_time, conc_il12,
              hill_1, hill_2, fb_strength, rate_ifn, rate_il4, K,
              initial_cells, degradation, fb_start, fb_end]

hill_1 = [10,10,10]
hill_2 = hill_1

hill10_params = [alpha_1, alpha_2, beta_1, beta_2, simulation_time, conc_il12,
              hill_1, hill_2, fb_strength, rate_ifn, rate_il4, K,
              initial_cells, degradation, fb_start, fb_end]

hill10_rtm_params = [20, 20, 20, 20, simulation_time, conc_il12,
              hill_1, hill_2, fb_strength, rate_ifn, rate_il4, K,
              initial_cells, degradation, fb_start, fb_end]

#==============================================================================
# chain
#==============================================================================
chain_params = list(rate_params)
chain5_params = list(hill5_params)
chain10_params = list(hill10_params)

fig, ax = plt.subplots(1,1, figsize = (5,4))
chain_arr = np.arange(1,20,1)

chain_readouts = variable_effect(chain_arr, "chain", *chain_params)
chain_cells = get_cell_arr(chain_readouts)
ax_chain(chain_arr, chain_cells, ax, plot_both = False, xsteps = 5)

chain_readouts = variable_effect(chain_arr, "chain", *chain5_params)
chain_cells = get_cell_arr(chain_readouts)
ax_chain(chain_arr, chain_cells, ax, plot_both = False, xsteps = 5)

chain_readouts = variable_effect(chain_arr, "chain", *chain10_params)
chain_cells = get_cell_arr(chain_readouts)
ax_chain(chain_arr, chain_cells, ax, plot_both = False, xsteps = 5)


ax.set_ylabel("% Th1 cells at final state")
ax.set_ylabel("%Th cells")
plt.tight_layout()
#fig.savefig(save_path+"chain_pos_fb.svg", bbox_inches = "tight")

#==============================================================================
# fb strength
#==============================================================================
fb_params = list(rate_params)
fb_arr = np.linspace(0,10,100)
fb_readouts = variable_effect(fb_arr, "feedback_strength_pos_Th1", *fb_params)
fb_cells = get_cell_arr(fb_readouts)

fb_params_rtm = list(rtm_params)
fb_arr_rtm = np.linspace(0,10,100)
fb_readouts_rtm = variable_effect(fb_arr_rtm, "feedback_strength_pos_Th1", *fb_params_rtm)
fb_cells_rtm = get_cell_arr(fb_readouts_rtm)


fig, ax = plt.subplots()
ax_var_effect(fb_arr, fb_cells, ax, linestyle = "--")
ax_var_effect(fb_arr_rtm, fb_cells_rtm, ax)
ax.set_xscale('log')
ax.set_xticks([1, 10])
ax.set_xlim(0,10)
#ax.set_ylim(0,100)

ax.set_xlabel("feedback strength")
ax.set_ylabel("%Th cells")
plt.tight_layout()