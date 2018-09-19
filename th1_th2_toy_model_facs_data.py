#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Thu Aug  9 14:44:26 2018

@author: burt
"""

from scipy.special import gamma, gammainc
from scipy.optimize import curve_fit
import numpy as np
import matplotlib.pyplot as plt
from scipy.integrate import odeint


SMALL_SIZE = 15
MEDIUM_SIZE = 22

plt.rc('font', size=MEDIUM_SIZE)          # controls default text sizes
plt.rc('axes', titlesize=SMALL_SIZE)     # fontsize of the axes title
plt.rc('axes', labelsize=MEDIUM_SIZE)    # fontsize of the x and y labels
plt.rc('xtick', labelsize=SMALL_SIZE)    # fontsize of the tick labels
plt.rc('ytick', labelsize=SMALL_SIZE)    # fontsize of the tick labels
plt.rc('legend', fontsize=15)    # legend fontsize
plt.rc('figure', titlesize=MEDIUM_SIZE)

def gamma_dist(t, alpha, beta):
    return np.exp(-beta*t)*(beta**alpha)*(t**(alpha-1))/gamma(alpha)

def gamma_cdf(t, alpha, beta):
    dummy = t*beta
    return gammainc(alpha,dummy)



### set up exp data and dummy time array
### deprecated because this is restimulation data
"""
t = np.linspace(0,8,100)

exp_time_points = np.array([0,0.5,1,1.5,2,3,3.5,4,5,6,7], dtype = np.float64)

pos_cells_ifn = np.array([0,0,15,35,45,52,53,51,45,38,27],dtype = np.float64)
pos_cells_il4 = np.array([0,0,0,6,13.5,16,15.5,14.5,15,13.5,13.5],dtype = np.float64)


#only look at first n hours
n = 7
exp_time_points = exp_time_points[:7]
pos_cells_ifn = pos_cells_ifn[:7]
pos_cells_il4 = pos_cells_il4[:7]

### normalize data for pdf (deprecated)
pos_cells_ifn_norm = pos_cells_ifn/np.max(pos_cells_ifn)
pos_cells_il4_norm = pos_cells_il4/np.max(pos_cells_il4)


il4_fit, il4_pcov = curve_fit(gamma_cdf,exp_time_points,pos_cells_il4_norm, bounds=([2.,-np.inf], np.inf))
ifn_fit, ifn_pcov = curve_fit(gamma_cdf,exp_time_points,pos_cells_ifn_norm, bounds=([2.,-np.inf], np.inf))

### fit data to gamma dist - this is deprecated!
#il4_fit, il4_pcov = curve_fit(gamma_dist,exp_time_points,pos_cells_il4_norm, bounds=([2.,-np.inf], np.inf))

#ifn_fit, ifn_pcov = curve_fit(gamma_dist,exp_time_points,pos_cells_ifn_norm, bounds=([2.,-np.inf], np.inf))
"""

def norm(arr):
    arr = arr-np.min(arr)
    arr = arr/np.max(arr)
    return arr


gata3_time = np.array([0,3,6,12,24,35,48,73])
tbet_time = np.array([0,3,6,12,35,48,73,96])


gata_3_exp1 = np.array([1,2,1.7,5,5.2,8.5,10.4,9.9])
tbet_exp1 = np.array([1,1.2,2.5,2.6,8.75,6.8,12.4,9.7])

tbet_norm = norm(tbet_exp1)
gata3_norm = norm(gata_3_exp1)

"""
fig, ax = plt.subplots(1,2, figsize=(5,5))
ax[0].plot(gata3_time,gata_3_exp1)
ax[1].plot(tbet_time,tbet_exp1)
"""
tbet_fit, tbet_pcov = curve_fit(gamma_cdf,tbet_time,tbet_norm, bounds=([2.,-np.inf], np.inf))
gata3_fit, gata3_pcov = curve_fit(gamma_cdf,gata3_time,gata3_norm, bounds=([2.,-np.inf], np.inf))

### fit data to gamma dist - this is deprecated!
#il4_fit, il4_pcov = curve_fit(gamma_dist,exp_time_points,pos_cells_il4_norm, bounds=([2.,-np.inf], np.inf))

#ifn_fit, ifn_pcov = curve_fit(gamma_dist,exp_time_points,pos_cells_ifn_norm, bounds=([2.,-np.inf], np.inf))

ylabel_tbet = "tbet positive cells"
ylabel_gata3 = "gata3 positive cells"
xlabel = "Time [h]"
#xlim = [0,4]
t = np.linspace(0,98,98)

#### plot with correct alpha fits
fig, ax = plt.subplots(1,1, figsize = (5,5))

ax.scatter(tbet_time,tbet_norm, label = "Th1")
ax.scatter(gata3_time,gata3_norm, c = "tab:red", label = "Th2")
ax.plot(t,gamma_cdf(t, tbet_fit[0],tbet_fit[1]), c = "tab:blue", linestyle = "-", label = "gamma cdf fit, alpha = "+str(round(tbet_fit[0],2)))
ax.plot(t,gamma_cdf(t, gata3_fit[0],gata3_fit[1]), c = "tab:red", linestyle = "-", label = "gamma cdf fit, alpha = "+str(round(gata3_fit[0],2)))
ax.legend()
ax.set_yticks([0,0.5,1])
ax.set_yticks([0,0.5,1])
plt.tight_layout()


fig, ax = plt.subplots(1,1, figsize = (5,5))

ax.scatter(tbet_time,tbet_norm, label = "T-bet", s= 70)
ax.scatter(gata3_time,gata3_norm, c = "tab:red", label = "GATA3", s = 70)
ax.plot(t,gamma_cdf(t, 2,tbet_fit[1]), c = "tab:blue", linestyle = "-", label = "Model fit")
ax.plot(t,gamma_cdf(t, 2,gata3_fit[1]), c = "tab:red", linestyle = "-", label = "Model fit")
ax.legend()
ax.set_xlabel(xlabel)
ax.set_ylabel("% of maximum")
ax.set_yticks([0,0.5,1])
ax.set_yticks([0,0.5,1])
ax.set_ylim([0,1.05])
ax.set_xlim([0,100])
plt.tight_layout()
#fig.savefig("gamma_fit_th1_th2_facs_data.svg",dpi=1200, bbox_inches = "tight")

def p_th_diff(conc_ifn,conc_il4,conc_il12, hill, conc_cyto = 10**(-12)):
    return ((conc_ifn/conc_cyto)**hill[0])*((conc_il4/conc_cyto)**hill[1])*((conc_il12/conc_cyto)**hill[2])


rate_1 = tbet_fit[1]
rate_2 = gata3_fit[1]


mean_1= 2 / rate_1
mean_2= 2 / rate_2
#conc_il12 = 5

def th_cell_diff(state, t, conc_il12, hill_1, hill_2, rate_ifn, rate_il4):
        ### purpose:simulate Hong et al 2008 model for neuropora clock


        conc_ifn = rate_ifn*state[2]
        conc_il4 = rate_il4*state[4]
        ### define state vector
        th_0 = state[0]
        
        prob_th1 = p_th_diff(conc_ifn,conc_il4,conc_il12, hill_1)
        prob_th2 = p_th_diff(conc_ifn,conc_il4,conc_il12, hill_2)
        
        prob_th1_norm = prob_th1 / (prob_th1+prob_th2)
        prob_th2_norm = prob_th2 / (prob_th1+prob_th2)
        #print prob_th2_norm
        th1_0 = prob_th1_norm*th_0
        th2_0 = prob_th2_norm*th_0
        
        th1_1 = state[1]
        th1_2 = state[2]
        th2_1 = state[3]
        th2_2 = state[4]


        
        
        ###  ODEs Hong et al 2008
        ### letzter summand unklar bei dtfrqmrna
        
        dt_th1_0 = -rate_1*th1_0
        dt_th1_1 = rate_1*(th1_0-th1_1)       
        dt_th1_2 = rate_1*(th1_1)
        
        dt_th2_0 = -rate_2*th2_0
        dt_th2_1 = rate_2*(th2_0-th2_1)   
        dt_th2_2 = rate_2*(th2_1)
        dt_th0 = dt_th1_0+dt_th2_0
        ### derivatives
        
        
        dt_state = [dt_th0,
                    dt_th1_1,
                    dt_th1_2,
                    dt_th2_1,
                    dt_th2_2]
        
        return dt_state
    
def th_cell_diff2(state, t, conc_il12, hill_1, hill_2, rate_ifn, rate_il4):
        ### purpose:simulate Hong et al 2008 model for neuropora clock


        conc_ifn = rate_ifn*state[1]
        conc_il4 = rate_il4*state[2]
        ### define state vector
        th_0 = state[0]
        
        prob_th1 = p_th_diff(conc_ifn,conc_il4,conc_il12, hill_1)
        prob_th2 = p_th_diff(conc_ifn,conc_il4,conc_il12, hill_2)
        
        prob_th1_norm = prob_th1 / (prob_th1+prob_th2)
        prob_th2_norm = prob_th2 / (prob_th1+prob_th2)
        print prob_th2_norm
        th1_0 = prob_th1_norm*th_0
        th2_0 = prob_th2_norm*th_0
        
        dt_th_0 = (-rate_1*th1_0)+(-rate_2*th2_0)       
        dt_th1 = rate_1*(th1_0)
        dt_th2 = rate_2*(th2_0)
        
        
        dt_state = [dt_th_0,
                    dt_th1,
                    dt_th2]
        
        return dt_state
               
### set initial state and time vector
### set initial conditions for each ODE

states = np.ones(5)
states[0] = 10000

### function to run model with parameters
def run_model(model = th_cell_diff, ini_cond = states, t = np.arange(0,100,0.01), conc_il12=3*10**(-12),
              hill_1 = [1,-1,1], hill_2 = [-1,1,0], rate_ifn=10**(-15), rate_il4=10**(-15)):
    
    state = odeint(model,ini_cond,t, args =(conc_il12, hill_1, hill_2, rate_ifn, rate_il4,))

    return state

state = run_model()


th0_cells = state[:,0]/100
th1_cells = state[:,2]/100
th2_cells = state[:,4]/100

th_cells = [th0_cells,th1_cells, th2_cells]
th_cell_sum = th0_cells+th1_cells+th2_cells


fix, ax = plt.subplots(1,1, figsize = (5,5))
ax.plot(t, th0_cells)
ax.plot(t, th1_cells)
ax.plot(t, th2_cells)


th_cells = (th_cells/th_cell_sum)*100

t = np.arange(0,100,0.01)

fig, ax = plt.subplots(1,1, figsize = (5,5))
ax.plot(t,th_cells[0], "k",
        t,th_cells[1], "tab:blue",
        t,th_cells[2], "tab:red")
ax.set_yticks([0,50,100])
ax.set_ylim([0,100])
ax.set_xlim([0,100])
ax.set_xlabel("Time [h]")
ax.set_ylabel("% Th cells")
#ax.set_title("[IL12] = 1 pM")
plt.legend(["Th0","Th1","Th2"], loc = "right")
plt.tight_layout()

fig, ax = plt.subplots(1,1, figsize = (5,5))
ax.plot(t,th_cells[2], "tab:red")
ax.set_yticks([0,50,100])
ax.set_ylim([0,100])
ax.set_xlim([0,100])
ax.set_xlabel("Time [h]")
ax.set_ylabel("% Th cells")
#ax.set_title("[IL12] = 1 pM")
plt.legend(["Th0","Th1","Th2"], loc = "right")
plt.tight_layout()

#fig.savefig("th_model_time_course.svg",bbox_inches="tight", dpi=1200)
"""
inis = np.ones(3)
inis[0] = 10000
state = run_model(model= th_cell_diff2, ini_cond = inis, hill_1 = [1,0,1], hill_2 = [0,1,0])

plt.plot(t, state)

#### effect of IL12 on th1 end state
il12_concentrations = np.linspace(0,6*10**(-12),100)

th1_conc = []
th2_conc = []
th0_conc = []

for i in il12_concentrations:
    
    state = run_model(t = np.arange(0,100,0.01), conc_il12=i)
    
    th1_endstate = state[-1,2]
    th1_conc.append(th1_endstate)
    
    th2_endstate = state[-1,4]
    th2_conc.append(th2_endstate)
    
    th0_endstate = state[-1,0]
    th0_conc.append(th0_endstate)
    
th0_conc = np.asarray(th0_conc)
th1_conc = np.asarray(th1_conc)
th2_conc = np.asarray(th2_conc)

all_tcells = th0_conc+th1_conc+th2_conc

th0_conc = (th0_conc / all_tcells)*100
th1_conc = (th1_conc / all_tcells)*100
th2_conc = (th2_conc / all_tcells)*100

il12_conc_pm = il12_concentrations*(10**(12))

fig, ax = plt.subplots(1,1,figsize=(5,5))

#ax.plot(il12_conc_pm, th0_conc)
ax.plot(il12_conc_pm, th1_conc, "tab:blue")
ax.plot(il12_conc_pm, th2_conc, "tab:red")
ax.set_xlabel("IL12 [pM]")
ax.set_yticks([0,50,100])
ax.set_xlim([0,6])
#ax.set_title("Effect of IL-12 conc. on Th cell balance after 3 hours")
ax.set_ylim([0,100])
#ax.set_xlim([0,0.1])

ax.set_ylabel("% Th cells")
plt.legend(["Th1","Th2"])
plt.tight_layout()
#fig.savefig("th_model_il12_effect.svg",bbox_inches="tight", dpi=1200)
"""
