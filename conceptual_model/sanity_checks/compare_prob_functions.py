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
import matplotlib.pyplot as plt
import seaborn as sns

#==============================================================================
# define params to store figure
#==============================================================================
sns.set(context = "talk", style = "ticks")
save_path = "/home/burt/Documents/tcell_project/figures/"

def menten(x, p, K):
    dummy = (x**p / (x**p + K))
    return dummy

def neg_fb(x, fb = 1.0, p = 0, K = 1):
    # negative feedback strength should be between 0 and 1
    return 1 - fb * menten(x, p, K)

def pos_fb(x, fb = 0.5, p = 1, K = 1):
    # positive fb strength should be greater or equal 1
    return 1 + fb * menten(x, p, K)

def normalize(a, b):
    return a / (a + b)

x = np.arange(0,10,0.01)

p1 = pos_fb(x)
p2 = 1

p1_norm = normalize(p1, p2)
p2_norm = normalize(p2, p1)

plt.plot(x, p1_norm, "tab:blue")
plt.plot(x, p2_norm, "tab:red")