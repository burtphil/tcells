#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Fri Dec 14 15:38:12 2018

@author: burt
"""

from scipy.integrate import odeint
import numpy as np
import matplotlib.pyplot as plt
a = 1
b = 2

y0 = 1

def ode(state, t):
    dt_state = a- b * state
    return dt_state

time = np.linspace(0,5,100)
state = odeint(ode, y0, time)

plt.plot(time, state)

#def exp_cdf(t, mean):
#    return 1 - np.exp(-t / mean)

#plt.plot(time, exp_cdf(time, 1/100.))
#plt.plot(time, exp_cdf(time, 1))