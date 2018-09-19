import numpy as np
import os

os.chdir("/home/burt/Documents/scripts/th_cell_differentiation")

### run_model takes default parameters stored in th1_th2_parameters.py
from th1_th2_ode_model_generic import run_model, chain
from th1_th2_plotting import plot_time_course,plot_chain
import th1_th2_parameters as params
import th1_th2_conceptual_parameters as cparams

#==============================================================================
# import parameters
#==============================================================================
model_params = params.parameters
conceptual_params = cparams.parameters

mean_th1 = params.alpha_th1/params.beta_th1
mean_th2 = params.alpha_th2/params.beta_th2

#==============================================================================
# run time course simulation
#==============================================================================
simulation_name = "sim"

th1_th2_model = run_model(simulation_name, parameters = conceptual_params)

test_simulation = np.load(simulation_name+".npz")
state = test_simulation["state"]
parameters = test_simulation["parameters"]

# plot time course
plot_time_course(th1_th2_model, parameters, save = "time_course_conc_model_il12_1")

chain_length = 25.

plot_chain(chain(chain_length, conceptual_params, stepsize = cparams.stepsize))
