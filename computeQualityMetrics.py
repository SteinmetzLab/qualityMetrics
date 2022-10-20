# -*- coding: utf-8 -*-
"""
Created on Thu Oct 20 09:37:19 2022

@author: Noam Roth
"""

from computeNoiseCutoff import *
from slidingRefractory import *


def computeQualityMetrics(spiketimes, amps):
    '''
    

    Parameters
    ----------
    spiketimes : spike times for the neuron across the entire recording
    amps: amplitudes (or template amplitudes) for the neuron across the entire recording

    Returns
    -------
    pass : binary variable which states whether your neuron passes (1) or fails (0)


    '''
    
    
    nc_pass, cutoff_value, first_low_quantile = noise_cutoff(amps)
    