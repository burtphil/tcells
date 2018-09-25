###############################################################################
###################### import modules  ########################################
###############################################################################

from scipy.special import gamma, gammainc
from scipy.optimize import curve_fit
import numpy as np
import matplotlib.pyplot as plt
import os

os.chdir("/home/burt/Documents/code/th_cell_differentiation")
###############################################################################
###################### set up plotting params #################################
###############################################################################

SMALL_SIZE = 10
MEDIUM_SIZE = 22

plt.rc('font', size=MEDIUM_SIZE)          # controls default text sizes
plt.rc('axes', titlesize=SMALL_SIZE)     # fontsize of the axes title
plt.rc('axes', labelsize=MEDIUM_SIZE)    # fontsize of the x and y labels
plt.rc('xtick', labelsize=SMALL_SIZE)    # fontsize of the tick labels
plt.rc('ytick', labelsize=SMALL_SIZE)    # fontsize of the tick labels
plt.rc('legend', fontsize=12)    # legend fontsize
plt.rc('figure', titlesize=MEDIUM_SIZE)

###############################################################################
###################### define functions #######################################
###############################################################################

def gamma_dist(t, alpha, beta):
    return np.exp(-beta*t)*(beta**alpha)*(t**(alpha-1))/gamma(alpha)

def gamma_cdf(t, alpha, beta):
    dummy = t*beta
    return gammainc(alpha,dummy)

def norm(arr):
    arr = arr-np.min(arr)
    arr = arr/np.max(arr)
    return arr


###############################################################################
###################### set up data  ###########################################
###############################################################################


gata3_time = np.array([0,3,6,12,24,35,48,73])
tbet_time = np.array([0,3,6,12,35,48,73,96])


gata_3_exp1 = np.array([1,2,1.7,5,5.2,8.5,10.4,9.9])
tbet_exp1 = np.array([1,1.2,2.5,2.6,8.75,6.8,12.4,9.7])

# should there be more data I need to adjust
exp_times = [gata3_time,tbet_time]

###############################################################################
###################### gamma fit  #############################################
###############################################################################

tbet_norm = norm(tbet_exp1)
gata3_norm = norm(gata_3_exp1)

tbet_fit, tbet_pcov = curve_fit(gamma_cdf,tbet_time,tbet_norm, bounds=([2.,-np.inf], np.inf))
gata3_fit, gata3_pcov = curve_fit(gamma_cdf,gata3_time,gata3_norm, bounds=([2.,-np.inf], np.inf))


###############################################################################
########## write fit params to file for model simulation  #####################
###############################################################################

np.save('th1_th2_gamma_fit_params', np.array([tbet_fit[0],gata3_fit[0],tbet_fit[1], gata3_fit[1]]))

###############################################################################
###################### visualization  #########################################
###############################################################################

### I might need to make this more generic once other data comes in
ylabel_tbet = "tbet positive cells"
ylabel_gata3 = "gata3 positive cells"
xlabel = "Time [h]"

# make the simulation time an even dividor of 10 that is bigger than the maximum
# entry of the exp time data
simulation_time = np.max(exp_times)-(np.max(exp_times)%10)+10

xlim = [0,simulation_time]

# I might adjust this to the maximum time of the fit
t = np.arange(0,simulation_time,0.01)

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
ax.set_xlim(xlim)
plt.tight_layout()
#fig.savefig("gamma_fit_th1_th2_facs_data.svg",dpi=1200, bbox_inches = "tight")

