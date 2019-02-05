#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Fri Jan 25 15:39:54 2019

@author: burt
script to visualize probabilities
"""

import numpy as np
import matplotlib.pyplot as plt

def menten(x, hill, K):
    dummy = (x**hill / (x**hill + K))
    return dummy

def feedback(x, fb, hill, K):
    # positive fb strength should be greater or equal 1
    """ 
    use neg fb values for negative fb
    0 for no feedback (returns 1)
    and positive values for pos fb
    """

    dummy = 1 + fb * menten(x, hill, K)
    #assert dummy >= 0
    
    return dummy

def gamma(x, fb, hill, K):
    dummy = (x * fb + K) / (x + K)
    return dummy

cells = np.arange(0,5, 0.01)

prob = feedback(cells, fb = 10, hill = 2., K = 1.)
prob = prob / (1 + prob)

#fig, ax = plt.subplots()
#ax.plot(cells, prob)

prob = gamma(cells, fb = 0.1, hill = 2., K = 1.)
prob = prob / (1 + prob)

fig, ax = plt.subplots()
ax.plot(cells, prob)