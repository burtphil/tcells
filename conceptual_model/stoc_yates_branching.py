#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Fri Jan 18 11:41:58 2019

@author: burt
prolif model yates as stochastic process
with branching
"""

import numpy as np
import matplotlib.pyplot as plt
from scipy.special import gammainc
import seaborn as sns
from scipy.integrate import odeint
import pandas as pd
sns.set(context = "talk", style = "ticks")

#==============================================================================
# params
#==============================================================================
start = 0
stop = 2
nsteps = 12000
simulation_time = np.linspace(start, stop, nsteps)

# stoc params
nstates = 6
nsim = 100

# diff rate
alpha_1 = 10
alpha_2 = 1
beta_1 = float(alpha_1)
beta_2 = float(alpha_2)

# prolif rate
alpha_prolif = 20
beta_prolif = float(alpha_prolif)

# other rates
rate_birth = 1.0
rate_death = 1.0


# other params
# one th1 precursor, one th2 precursor
ncells = 2

stoc_params = [start, 
               stop, 
               nsteps, 
               ncells, 
               nstates, 
               rate_birth, 
               alpha_1, 
               beta_1, 
               alpha_2,
               beta_2,
               alpha_prolif,
               beta_prolif,
               rate_death,
               ]

prolif_params = [simulation_time, 
                     alpha_1, 
                     alpha_2, 
                     beta_1, 
                     beta_2, 
                     alpha_prolif,
                     beta_prolif,
                     rate_birth, 
                     rate_death,
                ]
#==============================================================================
# functions
#==============================================================================
def gamma_cdf(t, alpha, beta):
    dummy = t*beta
    return gammainc(alpha,dummy)      


def make_cells(ncells, nsteps, nstates):
    """
    returns an array of cells for given number of
    desired cells, steps and cell states
    initializes random numbers probability
    of cell state transition
    """
    cells = np.zeros((ncells, nsteps, nstates))
    for i in range(ncells):        
        # random number for th0 to th_eff transition
        cells[i, :, 2] = np.random.rand()  
        cells[i, :, 4] = np.random.rand()
    return cells

def exp_cdf(t, rate):
    return 1 - np.exp(-t * rate)


def prob_simple(x):
    y = x / (x + 1)
    return y

def stoc_simulation(start, 
                    stop, 
                    nsteps, 
                    ncells, 
                    nstates, 
                    rate_birth, 
                    alpha_th1_diff, 
                    beta_th1_diff, 
                    alpha_th2_diff, 
                    beta_th2_diff,                    
                    alpha_prolif, 
                    beta_prolif, 
                    rate_death):
    """
    run a simulation of the th1 th2 models   
    this function uses the a cumulative gamma fct integral to calculate transition probability
    """
    # define state indices for the nstates array
    cell_state = 0
    prob_state_change = 1
    rnd_change = 2
    t_last_change = 3
    rnd_death = 4
    prob_death = 5

    # define cell_state_idx state values
    th1_prec_cell = 0
    th1_cell = 1
    th2_prec_cell = 2
    th2_cell = 3
    
    dead_cell = 4

    # initialize some dummy dead cells t
    cells = make_cells(ncells+1, nsteps, nstates)
    cells[0,:, cell_state] = dead_cell
    cells[1,:, cell_state] = th1_prec_cell
    cells[2,:, cell_state] = th2_prec_cell
    
    # initialize one th1_cell
    time = np.linspace(start, stop, nsteps)

    p_th1_birth = 0
    t0_th1 = 0
    rnd_th1_birth = np.random.rand()    

    p_th2_birth = 0
    t0_th2 = 0
    rnd_th2_birth = np.random.rand()    

    # loop over each time point
    for i in range(len(time)-1):
        t = time[i]
        t_new = time[i+1]

        # get all cells for current time step
        cell_j = cells[:,i,:]

        # cell number should be number of alive cells not total number              
        cell_number = cell_j.shape[0]
        #print cell_number, t
        counter_th1 = 0
        counter_th2 = 0
        # loop over each cell for current time point
           
        for j in range(cell_number):
            
            cell = cell_j[j,:] 

            #check if cell differentiates
            if cell[cell_state] == th1_prec_cell:
                                             
                if cell[rnd_change] > cell[prob_state_change]:
                    #print cell[prob_state_change]
                    cell[prob_state_change] = (cell[prob_state_change]+
                        (gamma_cdf(t_new-cell[t_last_change], alpha_th1_diff, beta_th1_diff)-
                         gamma_cdf(t-cell[t_last_change], alpha_th1_diff, beta_th1_diff)))
                    #print cell[prob_state_change]
                else:
                    cell[cell_state] = th1_cell
                    cell[t_last_change] = t
                    cell[prob_state_change] = 0
                    cell[rnd_change] = np.random.rand()
                    cell[rnd_death] = np.random.rand()
                    cell[prob_death] = 0
                    
                    counter_th1 += 1

                if cell[rnd_death] > cell[prob_death]:
                    
                    cell[prob_death] = (cell[prob_death]+
                        (exp_cdf(t_new-cell[t_last_change], rate_death)-
                         exp_cdf(t-cell[t_last_change], rate_death)))
                else:
                    cell[cell_state] = dead_cell

            #check if cell differentiates
            if cell[cell_state] == th1_cell:
                                             
                if cell[rnd_change] > cell[prob_state_change]:
                    #print cell[prob_state_change]
                    cell[prob_state_change] = (cell[prob_state_change]+
                        (gamma_cdf(t_new-cell[t_last_change], alpha_prolif, beta_prolif)-
                         gamma_cdf(t-cell[t_last_change], alpha_prolif, beta_prolif)))
                    #print cell[prob_state_change]
                else:
                    cell[cell_state] = th1_cell
                    cell[t_last_change] = t
                    cell[prob_state_change] = 0
                    cell[rnd_change] = np.random.rand()
                    cell[rnd_death] = np.random.rand()
                    cell[prob_death] = 0
                    
                    counter_th1 += 1

                if cell[rnd_death] > cell[prob_death]:
                    
                    cell[prob_death] = (cell[prob_death]+
                        (exp_cdf(t_new-cell[t_last_change], rate_death)-
                         exp_cdf(t-cell[t_last_change], rate_death)))
                else:
                    cell[cell_state] = dead_cell

            if cell[cell_state] == th2_prec_cell:
                                             
                if cell[rnd_change] > cell[prob_state_change]:
                    #print cell[prob_state_change]
                    cell[prob_state_change] = (cell[prob_state_change]+
                        (gamma_cdf(t_new-cell[t_last_change], alpha_th2_diff, beta_th2_diff)-
                         gamma_cdf(t-cell[t_last_change], alpha_th2_diff, beta_th2_diff)))
                    #print cell[prob_state_change]
                else:
                    cell[cell_state] = th2_cell
                    cell[t_last_change] = t
                    cell[prob_state_change] = 0
                    cell[rnd_change] = np.random.rand()
                    cell[rnd_death] = np.random.rand()
                    cell[prob_death] = 0
                    
                    counter_th2 += 1

                if cell[rnd_death] > cell[prob_death]:
                    
                    cell[prob_death] = (cell[prob_death]+
                        (exp_cdf(t_new-cell[t_last_change], rate_death)-
                         exp_cdf(t-cell[t_last_change], rate_death)))
                else:
                    cell[cell_state] = dead_cell
                    
            #check if cell differentiates
            if cell[cell_state] == th2_cell:
                                             
                if cell[rnd_change] > cell[prob_state_change]:
                    #print cell[prob_state_change]
                    cell[prob_state_change] = (cell[prob_state_change]+
                        (gamma_cdf(t_new-cell[t_last_change], alpha_prolif, beta_prolif)-
                         gamma_cdf(t-cell[t_last_change], alpha_prolif, beta_prolif)))
                    #print cell[prob_state_change]
                else:
                    cell[cell_state] = th2_cell
                    cell[t_last_change] = t
                    cell[prob_state_change] = 0
                    cell[rnd_change] = np.random.rand()
                    cell[rnd_death] = np.random.rand()
                    cell[prob_death] = 0
                    
                    counter_th2 += 1
                    
                if cell[rnd_death] > cell[prob_death]:
                    
                    cell[prob_death] = (cell[prob_death]+
                        (exp_cdf(t_new-cell[t_last_change], rate_death)-
                         exp_cdf(t-cell[t_last_change], rate_death)))
                else:
                    cell[cell_state] = dead_cell
                    
            cells[j,i+1,:] = cell

        # is there a birth? if not, cumulatively add p_birth
        if rnd_th1_birth > p_th1_birth:              
            p_th1_birth = p_th1_birth + (exp_cdf(t_new-t0_th1, 0.8 * rate_birth)-
                                 exp_cdf(t-t0_th1, 0.8 * rate_birth))
        
        # if there is a birth, draw new rnd number 
        # and make a new th0 cell out of a dead cell
        else:       

            p_th1_birth = 0
            t0_th1 = t
            rnd_th1_birth = np.random.rand()            
            #search for a dead cell to transform             
            new_th1_prec_cell = make_cells(1, nsteps, nstates)
            new_th1_prec_cell[:,i+1, t_last_change] = t 
            cells = np.concatenate((cells, new_th1_prec_cell)) 

        if rnd_th2_birth > p_th2_birth:              
            p_th2_birth = p_th2_birth + (exp_cdf(t_new-t0_th2, 0.2 * rate_birth)-
                                 exp_cdf(t-t0_th2, 0.2 * rate_birth))
        
        # if there is a birth, draw new rnd number 
        # and make a new th0 cell out of a dead cell
        else:       

            p_th2_birth = 0
            t0_th2 = t
            rnd_th2_birth = np.random.rand()            
            #search for a dead cell to transform             
            new_th2_prec_cell = make_cells(1, nsteps, nstates)
            new_th2_prec_cell[:,i+1, t_last_change] = t
            new_th2_prec_cell[:,i+1, cell_state] = th2_prec_cell
            cells = np.concatenate((cells, new_th2_prec_cell)) 

            
        new_cells = make_cells(counter_th1, nsteps, nstates)
        new_cells[:,i+1, t_last_change] = t
        new_cells[:,i+1, cell_state] = th1_cell
        cells = np.concatenate((cells, new_cells))        

        new_cells = make_cells(counter_th2, nsteps, nstates)
        new_cells[:,i+1, t_last_change] = t
        new_cells[:,i+1, cell_state] = th2_cell
        cells = np.concatenate((cells, new_cells))
        
    return [cells, time]     


def get_cells(cells, time):
    """
    input: cells and time, which is output from the function: run_stochastic_simulation
    use this for the model with precursor state times included
    """
    all_cells = cells[:,:, 0]
          
    th1_prec_cells = np.sum(all_cells == 0, axis = 0)
    th1_cells = np.sum(all_cells == 1, axis = 0)
    th2_prec_cells = np.sum(all_cells == 2, axis = 0)
    th2_cells = np.sum(all_cells == 3, axis = 0)    
    dead_cells = np.sum(all_cells == 4, axis = 0)
   
    return th1_prec_cells, th1_cells, th2_prec_cells, th2_cells, dead_cells

# =============================================================================
# run stochastic simulation
# =============================================================================

simulation = [stoc_simulation(*stoc_params) for i in range(nsim)]
thn_th1_cells = [get_cells(*simu)[:-1] for simu in simulation]
#==============================================================================
# 
#==============================================================================
def th_cell_diff(state, t, alpha_1, alpha_2, beta_1, beta_2, alpha_prolif, beta_prolif, rate_birth, rate_death):
        
    th1 = state[:(alpha_1+alpha_prolif)]
    th2 = state[(alpha_1+alpha_prolif):]
      
    dt_th1 = np.zeros_like(th1)
    dt_th2 = np.zeros_like(th2)
        
    th_states = [th1, th2]
    dt_th_states = [dt_th1, dt_th2]
    rate = [beta_1, beta_2]
    
    p_1 = 0.5
    p_2 = 0.5
    
    cell_flux = [p_1 * rate_birth, p_2 * rate_birth]
    
    alphas = [alpha_1, alpha_2]
    rate = [beta_1, beta_2]
    ### differential equations depending on the number on intermediary states
    # loop over states of th1 and th2 cells
    for i in range(2):
        th_state = th_states[i]
        dt_state = dt_th_states[i]
        r = rate[i]
        alpha = alphas[i]
        flux = cell_flux[i]
        #print flux
            
    # calculate derivatives
        for j in range(len(th_state)):
            if j == 0:
                dt_state[j] = flux - (r + rate_death) * th_state[j]
                
            elif j < alpha:
                dt_state[j] = r * (th_state[j-1] - th_state[j]) - rate_death * th_state[j]
                
            elif j == alpha:
                dt_state[j] = 2 * (r * th_state[j-1] + beta_prolif * th_state[-1]) - (rate_death + beta_prolif) * th_state[j]
            
            else:
                assert j > alpha
                dt_state[j] = beta_prolif * (th_state[j-1] - th_state[j]) - rate_death * th_state[j]

    # assign new number of naive cells based on the number of present naive cells that were designated th1_0 or th2_0
    #dt_th0 = state[0]    
    # return cell states
    dt_state = np.concatenate([dt_th1, dt_th2])
    
    return dt_state  

def run_model(simulation_time, 
              alpha_1, 
              alpha_2, 
              beta_1, 
              beta_2, 
              alpha_prolif,
              beta_prolif,
              rate_birth, 
              rate_death,
              ):
    """ 
    run ode model with parameters from params file
    needs shape and rate params as input as well as simulation time
    """

    # define initial conditions based on number of intermediary states
    no_of_states = int(alpha_1 + alpha_2 + alpha_prolif + alpha_prolif)
    ini_cond = np.zeros(no_of_states)
    
    initial_cells = 1
    ini_cond[0] = initial_cells
    ini_cond[alpha_1+alpha_prolif] = initial_cells

    
    state = odeint(th_cell_diff, 
                   ini_cond, 
                   simulation_time, 
                   args =(alpha_1, 
                          alpha_2, 
                          beta_1, 
                          beta_2, 
                          alpha_prolif,
                          beta_prolif,
                          rate_birth, 
                          rate_death)
                   )
    
    return state

state = run_model(*prolif_params)

# get cells from the first generation
th1_cells = state[:, alpha_1:(alpha_1+alpha_prolif)]
th2_cells = state[:, (alpha_1+alpha_prolif+alpha_2):]


th1_all_cells = np.sum(th1_cells, axis = 1)
th2_all_cells = np.sum(th2_cells, axis = 1)

fig, ax = plt.subplots()
ax.plot(simulation_time, th1_all_cells)
ax.plot(simulation_time, th2_all_cells)
# =============================================================================
# plot stochastic simulation
# =============================================================================

# make dummy labels for pandas df
label = ["Th1 pre","Th1", "Th2 pre","Th2"]
labels = [label for _ in range(nsim)]
flat_labels = [item for sublist in labels for item in sublist]

# make a df for single step stoc simulation (adjusted time)
flat_cells = [item for sublist in thn_th1_cells for item in sublist]
df = pd.DataFrame.from_records(flat_cells).transpose()
df.columns = flat_labels

df["time"] = simulation_time
df_long = df.melt(id_vars = ["time"])

df_diff = df_long.loc[(df_long["variable"] == "Th1") | (df_long["variable"] == "Th2")]

palette = ["tab:blue", "tab:orange"]

fig, ax = plt.subplots(1, 1, figsize =(5,4))
stoc_plot = sns.lineplot(x = "time", 
                         y = "value", 
                         data = df_diff,
                         hue = "variable",
                         ci = "sd",
                         ax = ax,
                         palette = palette,
                         legend = False)
ax.set_ylabel("$n_{cells}$")
ax.set_xlim(0, simulation_time[-1])

ax.plot(simulation_time, th1_all_cells, c = "tab:blue", linestyle = "--")
ax.plot(simulation_time, th2_all_cells, c = "tab:orange", linestyle = "--")
ax.legend(["Th1 stoc", "Th2 stoc", "Th1 det", "Th2 det"])
#ax.legend(label)
#ax.set_title(r"$\rightarrow$ Thn $\rightarrow$ Th diff $\rightarrow$ $\emptyset$")
plt.tight_layout()
fig.savefig("th1_th2_prolif.pdf", bbox_inches = "tight")