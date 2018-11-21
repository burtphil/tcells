#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Thu Nov 15 12:21:53 2018

@author: burt
"""

import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
import pandas as pd

save_path = "/home/burt/Documents/tcell_project/figures/"

x = [
     "both_rate",
     "both_rtm",
     "fb_rate",
     "fb_rtm",]

Th1 = [
       90,
       75,
       90,
       75,
       ]

Th2 = [
       10,
       25,
       10,
       25,
       ]

data = {"x" : x, "Th1" : Th1, "Th2" : Th2}

df = pd.DataFrame(data)
df = pd.melt(df, id_vars=['x'], value_vars=['Th1', 'Th2'])

barplot = sns.barplot(data = df, x = "x", y = "value", hue = "variable", palette = ["tab:blue", "tab:red"])
sns.set(style="white", context="talk")
barplot.set_xlabel("% Th cells in steady state")
plt.tight_layout()

fig = barplot.get_figure()
#fig.savefig(save_path+"fb_effect_barplot.svg", bbox_inches = "tight")
"""

#==============================================================================
# plot the "faster feedback wins" situation 
#==============================================================================

y = [
     "equal_dynamics",
     "unequal_dynamics",
     ]

Th1 = [
       50,
       85,
       ]

Th2 = [
       50,
       15,
       ]

data = {"y" : y, "Th1" : Th1, "Th2" : Th2}

df = pd.DataFrame(data)
df = pd.melt(df, id_vars=['y'], value_vars=['Th1', 'Th2'])

barplot = sns.barplot(data = df, y = "y", x = "value", hue = "variable", palette = ["tab:blue", "tab:red"])
sns.set(style="white", context="talk")
barplot.set_xlabel("% Th cells in steady state")
plt.tight_layout()


fig = barplot.get_figure()
#fig.savefig(save_path+"symmetry_effect_barplot.svg", bbox_inches = "tight")
"""