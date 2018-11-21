#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Wed Nov 21 11:07:16 2018

@author: burt
"""

import numpy as np

def p_new(conc_ifn, conc_il4, conc_il12, hill, half_saturation, fb_ifn, fb_il4, fb_il12):

    #assert conc_ifn >= 0, "ifn conc is "+str(conc_ifn)+", concentrations must be non-negative."
    #assert conc_il4 >= 0, "il4 conc is "+str(conc_il4)+", concentrations must be non-negative."
    #assert conc_il12 >= 0, "il12conc is "+str(conc_il12)+", concentrations must be non-negative."
    
    # watch out, this is tricky because if concentrations are zero and negative feedback (hill coeff < 0)
    # then you divide by zero, best avoid zero concentrations through base cytokine rate
    ifn_prob = (fb_ifn * conc_ifn**hill[0] + 1) / (conc_ifn**hill[0] + 1)
    #assert ifn_prob.all() >= 0
    il4_prob = (fb_il4 * (conc_il4**hill[1]) + 1) / ((conc_il4**hill[1]) + 1)
    #assert il4_prob.all() >= 0, "probability is "+str(il4_prob)
    il12_prob = 1.
    
    prob_th_diff = ifn_prob*il4_prob*il12_prob
    #assert prob_th_diff.all() >= 0, "probability is "+str(prob_th_diff)+" must be non-negative."
    return prob_th_diff


conc_il4 = np.linspace(0,600, 600)

conc_ifn = 1.
conc_il12 = 1.
hill = [1,0,0]
fb_ifn = 1.
fb_il4 = 1.
fb_il12 =1.


a = p_new(conc_ifn = conc_ifn, conc_il4 = conc_il4, conc_il12 = conc_il12, hill = hill, half_saturation = [0,0,0], fb_ifn = fb_ifn, fb_il4 = fb_il4, fb_il12 = fb_il12)