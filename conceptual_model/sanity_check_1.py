import numpy as np
import os

os.chdir("/home/burt/Documents/code/th_cell_differentiation")

### run_model takes default parameters stored in th1_th2_parameters.py
from modules.th1_th2_ode_model_generic import run_model2, il12
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

th1_th2_model = run_model2(simulation_name, parameters = conceptual_params)

test_simulation = np.load(simulation_name+".npz")
state = test_simulation["state"]
parameters = test_simulation["parameters"]

#plt.plot(cparams.simulation_time,state[:,-1])

# plot time course
#alpha_params = list(conceptual_params)
alpha_params = list(conceptual_params)
norm = cparams.initial_cells/100

fig, ax = plt.subplots(1,1, figsize = (5,4))

for i in range(1,10):
    alpha_params = list(conceptual_params)
    alpha_params[0] = int(i)
    alpha_params[1] = int(i)
    alpha_params[2] = int(i)
    alpha_params[3] = int(i)
    state = run_model2(simulation_name, parameters = alpha_params)
    #th0_cells = state[:,0]/norm
    alpha_1 = alpha_params[0]
    th1_cells = state[:,alpha_1+1]/norm
    th2_cells = state[:,-1]/norm
    #ax.plot(cparams.simulation_time, th2_cells, label = r"th2, $\alpha= 1$")
    ax.plot(cparams.simulation_time, th1_cells, label = r"$\alpha=$"+str(i))
    ax.set_yticks([0,50,100])
    ax.set_ylim([0,100])
    ax.set_xlim([cparams.simulation_time[0],cparams.simulation_time[-1]])
    ax.set_xlabel("Time [h]")
    ax.set_ylabel("% Th cells")
    ax.legend()
    ax.set_title(r"$p_{Th1}=p_{Th2}, \alpha_{Th1}=\alpha_{Th2}$")

plt.tight_layout()
#fig.savefig("conceptual_model_sanity.pdf", dpi=1200, bbox_inches = "tight")  


fig, ax = plt.subplots(1,1, figsize = (5,4))

for i in range(1,5):
    alpha_params = list(conceptual_params)
    alpha_params[0] = int(i)
    alpha_params[1] = int(1)
    alpha_params[2] = i
    alpha_params[3] = 1.0
    state = run_model2(simulation_name, parameters = alpha_params)
    #th0_cells = state[:,0]/norm
    alpha_1 = alpha_params[0]
    th1_cells = state[:,alpha_1+1]/norm
    th2_cells = state[:,-1]/norm
    ax.plot(cparams.simulation_time, th2_cells, label = r"th2, $\alpha= 1$")
    ax.plot(cparams.simulation_time, th1_cells, label = r"$\alpha_1=$"+str(i))
    ax.set_yticks([0,50,100])
    ax.set_ylim([0,100])
    ax.set_xlim([cparams.simulation_time[0],cparams.simulation_time[-1]])
    ax.set_xlabel("Time [h]")
    ax.set_ylabel("% Th cells")
    ax.legend()
    ax.set_title(r"$p_{Th1}=p_{Th2}, \alpha_2=1$")

plt.tight_layout()
