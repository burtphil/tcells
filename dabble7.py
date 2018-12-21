#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Fri Dec 21 10:30:19 2018

@author: burt
"""

import numpy as np

a = np.array([1,2,3])

def fun(a):
    a = a.view()
    return a*2

b = fun(a)

a = np.zeros((3,2))

a[1,0] = 1
a[2,0] = 1


b = a[:,0] == 1
c = a[b,:]

d = a[:,0]
indices = np.nonzero(b)[0]

# where are the th0 cells in cell array
#th0_bool = cells[:,0] == 0
#th0_cells = cells[th0_bool, :]
#th0_ids = np.nonzero(cells[:,0] == 0)[0]

#bool_2 = th0_cells[:,1] < th0_cells[:,2]

#th0_cells[bool_2, 3] = th0_cells[bool_2, 3] + fun(th0_cells[bool_2, 4]) 

