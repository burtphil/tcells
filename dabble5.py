# -*- coding: utf-8 -*-
"""
Created on Tue Dec 18 09:00:49 2018

@author: Philipp
"""
import numpy as np

def make_cells(ncells, nsteps, nstates):
    cells = np.zeros((ncells, nsteps, nstates))
    for i in range(cells.shape[0]):        
        # random number for th0 to th_eff transition
        cells[i, :, 1] = np.random.rand()  
    return cells

def add_cell(cells):
    """
    insert a new naive cell into cell vector if there are dead cells available
    otherwise append a new naive cell to cell vector
    """
    ncells, nsteps, nstates = cells.shape

    # initialize new cell
    new_cell = make_cells(1, nsteps, nstates)    
    dead_cell_idx = 1
    # check if there are any dead cells
    c = cells[:,:,0] == dead_cell_idx   
    if c.any():
        # if yes, insert new cell into dead cell vector
        d = np.where(cells[:,:,0] == dead_cell_idx)
        # get idx of chell
        d_cell = d[0][0]
        # get time step when this cell fate first occured
        d_time = d[1][0]
        #print new_cell
        cells[d_cell,d_time+1:,:] = new_cell        
    # else append new cell
    else:
        cells = np.concatenate((cells, new_cell))
    
    return cells

cells = np.zeros((3, 6, 2))

cells[2,2:,0] = 1

c = cells[:,2,0] == 1

d = np.where(cells[:,:,0] == 1)
d_cell = d[0][0]
d_time = d[1][0]

# ich will nicht wissen seit wann die zelle tot ist
# sondern seit wann

i = 2
# are there dead cells for time step i?
c = cells[:,i,0] == 1

# ich muss nur das aktuelle cell fate vererben, weil sich das dann
# für jeden Zeitschritt reproduziert
#  okay, what if there are no dead cells?
# I can just initialize 50 dummy dead cells