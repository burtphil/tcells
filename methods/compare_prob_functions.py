#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Thu Nov 22 11:52:55 2018

@author: burt
"""

#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
@author: burt
compare probability functions
"""

import numpy as np
import os
import matplotlib.image as mpimg
os.chdir("/home/burt/Documents/code/th_cell_differentiation")

### run_model takes default parameters stored in th1_th2_parameters.py
from modules.th1_th2_ode_model_generic import run_model, il12, get_readouts, hilli_sym, hilli_pos_th1
from modules.th1_th2_plotting import ax_time_course, ax_il12
import modules.th1_th2_conceptual_parameters as params
import matplotlib.pyplot as plt
import seaborn as sns

#==============================================================================
# define params to store figure
#==============================================================================
sns.set(context = "talk", style = "ticks")
save_path = "/home/burt/Documents/tcell_project/figures/"
image_path = "/home/burt/Documents/code/th_cell_differentiation/images/"

fig_names = ["pos_feedback",
            "neg_feedback",
            "pos_neg_feedback",
            "neg_th1_neg_pos_th2",
            "pos_th1_neg_pos_th2",
            "pos_th1_neg_th2",
            "pos_th2_neg_th2"]
#==============================================================================
# define plotting params
#==============================================================================
def menten(x, p = 1, K = 1):
    dummy = (x**p / (x**p + K))
    return dummy

def neg_fb(x, fb = 0.5, p = 1, K = 1):
    # negative feedback strength should be between 0 and 1
    return 1 - fb * menten(x, p, K)

def pos_fb(x, fb = 2, p = 1, K = 1):
    # positive fb strength should be greater or equal 1
    return 1 + fb * menten(x, p, K)

def normalize(a, b):
    return a / (a + b)

x = np.arange(0,10,0.01)

p1 = pos_fb(x)
p2 = 1

p1_norm = normalize(p1, p2)
p2_norm = normalize(p2, p1)

plt.plot(x, p1_norm)
plt.plot(x, p2_norm)