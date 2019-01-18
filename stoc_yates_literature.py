#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Fri Jan 18 11:41:58 2019

@author: burt
prolif model yates as stochastic process
"""

#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Fri Dec 14 10:39:00 2018

@author: burt
--> thn --> th1 --> death
also, th1 cells try to divide but this needs to be fixed
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
stop = 3
nsteps = 12000
simulation_time = np.linspace(start, stop, nsteps)

# stoc params
nstates = 6
nsim = 100
n_dummies = 100

# diff rate
alpha_1 = 10
alpha_2 = 1
beta_1 = float(alpha_1)
beta_2 = float(alpha_2)

# prolif rate
alpha_prolif = 10
beta_prolif = float(alpha_prolif)

# other rates
rate_birth = 0.5
rate_death = 2.0


y0 = np.zeros(alpha_1)
y0[0] = 1

# other params
ncells = 1

stoc_params = [start, 
               stop, 
               nsteps, 
               ncells, 
               nstates, 
               rate_birth, 
               alpha_1, 
               beta_1, 
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


def stoc_simulation(start, stop, nsteps, ncells, nstates, rate_birth, alpha_diff, beta_diff, alpha_prolif, beta_prolif, rate_death):
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
    th0_cell = 0
    th1_cell = 1
    dead_cell = 2


    # initialize some dummy dead cells t
    cells = make_cells(ncells, nsteps, nstates)
    # initialize one th1_cell
    time = np.linspace(start, stop, nsteps)

    p_birth = 0
    t0 = 0
    rnd_birth = np.random.rand()    

    # loop over each time point
    for i in range(len(time)-1):
        t = time[i]
        t_new = time[i+1]

        # get all cells for current time step
        cell_j = cells[:,i,:]

        # cell number should be number of alive cells not total number              
        cell_number = cell_j.shape[0]
        #print cell_number, t
        counter = 0      
        # loop over each cell for current time point
   
        
        for j in range(cell_number):
            
            cell = cell_j[j,:] 

            #check if cell differentiates
            if cell[cell_state] == th0_cell:
                                             
                if cell[rnd_change] > cell[prob_state_change]:
                    #print cell[prob_state_change]
                    cell[prob_state_change] = (cell[prob_state_change]+
                        (gamma_cdf(t_new-cell[t_last_change], alpha_diff, beta_diff)-
                         gamma_cdf(t-cell[t_last_change], alpha_diff, beta_diff)))
                    #print cell[prob_state_change]
                else:
                    cell[cell_state] = th1_cell
                    cell[t_last_change] = t
                    cell[prob_state_change] = 0
                    cell[rnd_change] = np.random.rand()
                    cell[rnd_death] = np.random.rand()
                    
                    counter += 1

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
                    
                    counter += 1
                    
                if cell[rnd_death] > cell[prob_death]:
                    
                    cell[prob_death] = (cell[prob_death]+
                        (exp_cdf(t_new-cell[t_last_change], rate_death)-
                         exp_cdf(t-cell[t_last_change], rate_death)))
                else:
                    cell[cell_state] = dead_cell


            cells[j,i+1,:] = cell

        # is there a birth? if not, cumulatively add p_birth
        if rnd_birth > p_birth:              
            p_birth = p_birth + (exp_cdf(t_new-t0, rate_birth)-
                                 exp_cdf(t-t0, rate_birth))
        
        # if there is a birth, draw new rnd number 
        # and make a new th0 cell out of a dead cell
        else:       

            p_birth = 0
            t0 = t
            rnd_birth = np.random.rand()            
            #search for a dead cell to transform             
            new_th0_cell = make_cells(1, nsteps, nstates)
            new_th0_cell[:,i+1, t_last_change] = t 
            cells = np.concatenate((cells, new_th0_cell)) 
            
        new_cells = make_cells(counter, nsteps, nstates)
        new_cells[:,i+1, t_last_change] = t
        new_cells[:,i+1, cell_state] = th1_cell
        cells = np.concatenate((cells, new_cells))        
        
    return [cells, time]     


def get_cells(cells, time):
    """
    input: cells and time, which is output from the function: run_stochastic_simulation
    use this for the model with precursor state times included
    """
    all_cells = cells[:,:, 0]
          
    th0_cells = np.sum(all_cells == 0, axis = 0)
    th1_cells = np.sum(all_cells == 1, axis = 0)
    dead_cells = np.sum(all_cells == 2, axis = 0)
   
    return th0_cells, th1_cells, dead_cells

# =============================================================================
# run stochastic simulation
# =============================================================================

simulation = [stoc_simulation(*stoc_params) for i in range(nsim)]
thn_th1_cells = [get_cells(*simu)[:2] for simu in simulation]
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
    
    p_1 = 1.0
    p_2 = 0
    
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
                dt_state[j] = flux - r * th_state[j]
                
            elif j < alpha:
                dt_state[j] = r * (th_state[j-1] - th_state[j])
                
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
    ini_cond[alpha_1+alpha_prolif] = 0

    
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
label = ["Thn","Th1"]
labels = [label for _ in range(nsim)]
flat_labels = [item for sublist in labels for item in sublist]

# make a df for single step stoc simulation (adjusted time)
flat_cells = [item for sublist in thn_th1_cells for item in sublist]
df = pd.DataFrame.from_records(flat_cells).transpose()
df.columns = flat_labels

df["time"] = simulation_time
df_long = df.melt(id_vars = ["time"])

fig, ax = plt.subplots(1, 1, figsize =(5,4))
stoc_plot = sns.lineplot(x = "time", 
                         y = "value", 
                         data = df_long,
                         hue = "variable",
                         ci = "sd",
                         ax = ax,
                         legend = False)
ax.set_ylabel("$n_{cells}$")
ax.set_xlim(0, simulation_time[-1])

ax.plot(simulation_time, th1_all_cells, c = "tab:orange", linestyle = "--")
#ax.legend(label)
#ax.set_title(r"$\rightarrow$ Thn $\rightarrow$ Th diff $\rightarrow$ $\emptyset$")
plt.tight_layout()