#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Tue Nov 27 10:20:29 2018

@author: burt
"""

from scipy.special import gamma

import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
save_path = "/home/burt/Documents/tcell_project/figures/"
sns.set(context = "talk", style = "ticks")


###############################################################################
###################### define functions #######################################
###############################################################################

def gamma_dist(t, alpha, beta):
    return np.exp(-beta*t)*(beta**alpha)*(t**(alpha-1))/gamma(alpha)

times = np.linspace(0,3,100)

rate = gamma_dist(times,1,1)
alpha_10 = gamma_dist(times, 10, 10)
alpha_30 = gamma_dist(times, 30, 30)

fig, ax = plt.subplots(1,1, figsize = (5,4))

ax.plot(times, rate, label = r"$\alpha = 1$")
ax.plot(times, alpha_10, label = r"$\alpha = 10$")
ax.plot(times, alpha_30, label = r"$\alpha = 30$")
ax.set_xlim(0,3)
ax.set_xlabel("time")
ax.legend()
ax.set_ylabel("density")
ax.set_yticks([0,1,2,3])
ax.set_ylim(0,3)
sns.despine()
plt.tight_layout()
#fig.savefig(save_path+"gamma_dist.svg", bbox_inches = "tight")

alpha_20 = gamma_dist(times, 20, 20)
fig, ax = plt.subplots(1,1, figsize = (5,4))
ax.plot(times, rate, "k--", label = "single-step model")
ax.plot(times, alpha_20, "k", label = "response-time model")
ax.set_xlim(0,2)
ax.set_xlabel("time")
ax.legend()
ax.set_ylabel("density")
ax.set_yticks([0,1,2,3])
ax.set_ylim(0,3)
sns.despine()
plt.tight_layout()
fig.savefig(save_path+"rtm_vs_rate_distribution.svg", bbox_inches = "tight")

fig, axes = plt.subplots(3,1, figsize = (5,12))
axes[0].plot(times, rate, "k--")
axes[0].plot(times, alpha_20, "k")

axes[1].plot(times, rate, "tab:blue", linestyle = "--")
axes[1].plot(times, alpha_20, "tab:red")

axes[2].plot(times, rate, "tab:red", linestyle = "--")
axes[2].plot(times, alpha_20, "tab:blue")
for ax in axes:
    ax.set_xlim(0,2)
    ax.set_xlabel("time")
    ax.set_ylabel("density")
    ax.set_yticks([0,1,2,3])
    ax.set_ylim(0,2)
    sns.despine()
plt.tight_layout()
fig.savefig(save_path+"symmetric_branch_distributions.svg", bbox_inches = "tight")
