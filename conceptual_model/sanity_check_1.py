import numpy as np
import os

os.chdir("/home/burt/Documents/code/th_cell_differentiation")

### run_model takes default parameters stored in th1_th2_parameters.py
from modules.th1_th2_ode_model_generic import run_model, il12
from modules.th1_th2_plotting import plot_time_course, plot_il12, ax_time_course, ax_il12
import modules.th1_th2_parameters as params
import modules.th1_th2_conceptual_parameters as cparams
import matplotlib.pyplot as plt
#==============================================================================
# import parameters
#==============================================================================
model_params = params.parameters
conceptual_params = cparams.parameters

#==============================================================================
# run time course simulation
#==============================================================================
simulation_name = "sim"

th1_th2_model = run_model(simulation_name, parameters = conceptual_params)

test_simulation = np.load(simulation_name+".npz")
state = test_simulation["state"]
parameters = test_simulation["parameters"]

# plot time course
#alpha_params = list(conceptual_params)
alpha_params = list(conceptual_params)
norm = cparams.initial_cells/100

fig, ax = plt.subplots(1,1, figsize = (5,4))

for i in range(1,5):
    alpha_params = list(conceptual_params)
    alpha_params[0] = int(i)
    alpha_params[2] = int(i)
    state = run_model(simulation_name, parameters = alpha_params)
    th0_cells = state[:,0]/norm
    alpha_1 = alpha_params[0]
    th1_cells = state[:,int(alpha_1)]/norm
    ax.plot(cparams.simulation_time, th0_cells)
    ax.plot(cparams.simulation_time, th1_cells, label = "alpha="+str(i))
    ax.set_yticks([0,50,100])
    ax.set_ylim([0,100])
    ax.set_xlim([cparams.simulation_time[0],cparams.simulation_time[-1]])
    ax.set_xlabel("Time [h]")
    ax.set_ylabel("% Th cells")
    ax.legend()

plt.tight_layout()
fig.savefig("conceptual_model_sanity.pdf", dpi=1200, bbox_inches = "tight")  

