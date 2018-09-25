import numpy as np
import os

os.chdir("/home/burt/Documents/code/th_cell_differentiation")

### run_model takes default parameters stored in th1_th2_parameters.py
from modules.th1_th2_ode_model_generic import run_model, il12
from modules.th1_th2_plotting import plot_time_course, plot_il12
import modules.th1_th2_parameters as params
import modules.th1_th2_conceptual_parameters as cparams

#==============================================================================
# import parameters
#==============================================================================
model_params = params.parameters
conceptual_params = cparams.parameters

#==============================================================================
# run time course simulation
#==============================================================================
simulation_name = "sim"

th1_th2_model = run_model(simulation_name, parameters = model_params)

test_simulation = np.load(simulation_name+".npz")
state = test_simulation["state"]
parameters = test_simulation["parameters"]

# plot time course
plot_time_course(th1_th2_model, parameters)

il_12_conc = np.linspace(0,10**(-11),100)

plot_il12(il12(il_12_conc, parameters = parameters))