import numpy as np
import os

os.chdir("/home/burt/Documents/code/th_cell_differentiation")

### run_model takes default parameters stored in th1_th2_parameters.py
from modules.th1_th2_ode_model_generic import run_model, il12, th_cell_diff
from modules.th1_th2_plotting import plot_time_course, plot_il12, ax_time_course, ax_il12
import modules.th1_th2_parameters as params
import matplotlib.pyplot as plt
from scipy.integrate import odeint
#==============================================================================
# import parameters
#==============================================================================
model_params = params.parameters

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

il12_conc = np.linspace(0,5*10**(-12),100)

#plot_il12(il12(il_12_conc, parameters = parameters))

fig, ax = plt.subplots(1,2, figsize = (10,5))

ax_time_course(state, ax[0], parameters)
ax_il12(il12(il12_conc, parameters), ax[1], xlabel = "IL-12 [pm]", factor_x_axis = 10**12)
ax[1].set_ylabel("% Th cells after 100 h")
ax[0].legend([r"Th0","Th1","Th2"])
plt.tight_layout()

fig.savefig("/home/burt/Documents/tcell_project/figures/model_simulations/data_model.svg", bbox_inches = "tight", dpi = 1200)
