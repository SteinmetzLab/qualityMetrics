# -*- coding: utf-8 -*-
"""
Created on Thu Oct 20 09:37:19 2022

@author: Noam Roth
"""

from computeNoiseCutoff import *
from slidingRefractory import *


def computeQualityMetrics(spikeTimes, spikeAmps):
    '''

    This code computes the full set of 3 quality metrics below as used in the IBL and Steinmetz Lab.
    Metrics:
        1) Sliding Refractory Period metric (slidingRefractory)
        2) Noise Cutoff metric (noiseCutoff)
        3) median amplitude threshold

    Parameters
    ----------
    spikeTimes : spike times for the neuron across the entire recording
    spikeAmps: amplitudes (or template amplitudes) for the neuron across the entire recording

    Returns
    -------
    pass : binary variable which states whether your neuron passes (1) or fails (0)

    '''
    

    #1) Compute RP values:
    maxConfidenceAt10Cont, minContWith90Confidence, timeOfLowestCont, nSpikesBelow2, confMatrix, cont, rp, nACG, firingRate, secondsElapsed = slidingRP(spikeTimes, params=None)

    # Also compute whether RP pass: #TODO: add nSpikesBelow2 condition
    if minContWith90Confidence < 10
        rp_pass = 1
    else:
        rp_pass = 0

    #2) Compute NC values and whether neuron passes NC:
    nc_pass, cutoff_value, first_low_quantile = noise_cutoff(spikeAmps)

    #3) Compute median amplitude: #TODO: make this optional in params above
    medAmp = np.nanmedian(spikeAmps)
    if medAmp > medAmpThresh:
        ma_pass = 1
    else:
        ma_pass = 0


    #Now compute whether neuron passes all 3 metrics, return as boolean
    neuron_pass = bool(rp_pass and nc_pass and ma_pass)

    return neuron_pass, rp_pass, nc_pass, ma_pass

