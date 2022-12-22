# -*- coding: utf-8 -*-
"""
Created on Thu Oct 20 09:37:19 2022

@author: Noam Roth
"""

from computeNoiseCutoff import *
from slidingRP.metrics import slidingRP
import mpl_scatter_density

def computeQualityMetrics(spikeTimes, spikeAmps, medAmpThresh = 50, verbose = True):
    '''

    This code computes the full set of 3 quality metrics below as used in the IBL and Steinmetz Lab.
    Metrics:
        1) Sliding Refractory Period metric (slidingRefractory)
        2) Noise Cutoff metric (noiseCutoff)
        3) median amplitude threshold (default threshold is 50 uV)

    Parameters
    ----------
    spikeTimes : spike times for the neuron across the entire recording
    spikeAmps: amplitudes (or template amplitudes) for the neuron across the entire recording

    Returns
    -------
    pass : binary variable which states whether your neuron passes (1) or fails (0)

    '''
    

    if verbose:
        print("Computing metrics...")

    #1) Compute RP values:
    maxConfidenceAt10Cont, minContWith90Confidence, timeOfLowestCont, nSpikesBelow2, confMatrix, cont, rp, nACG, firingRate, secondsElapsed = slidingRP(spikeTimes, params=None)

    # Also compute whether RP pass: #TODO: add nSpikesBelow2 condition
    if minContWith90Confidence < 10:
        rp_pass = True
    else:
        rp_pass = False

    #2) Compute NC values and whether neuron passes NC:
    nc_pass, cutoff_value, first_low_quantile = noise_cutoff(spikeAmps)

    #3) Compute median amplitude: #TODO: make this optional in params above
    medAmp = np.nanmedian(spikeAmps) * 1e6  #median amplitude in microvolts

    if (medAmp > medAmpThresh).any:
        ma_pass = True
    else:
        ma_pass = False


    #Now compute whether neuron passes all 3 metrics, return as boolean
    neuron_pass = bool(rp_pass and nc_pass and ma_pass)

    return neuron_pass, rp_pass, nc_pass, ma_pass


def computeQualityMetrics_All(spikeTimes, spikeAmps, spikeClusters, medAmpThresh=50, verbose = True):
    '''

    This code takes in spike times and amplitudes for a population of recorded neurons and
     computes the full set of 3 quality metrics below as used in the IBL and Steinmetz Lab
    (see 'computeQUalityMetrics')


    Parameters
    ----------
    spikeTimes : spike times for the neurons in the population
    spikeAmps: amplitudes (or template amplitudes) for the neurons in the population
    spikeClusters: np.array of size (spikeTimes) with the id of the cluster corresponding to that spike

    Returns
    -------
    metrics : dict which includes the outcome of quality metrics
        (overall_pass is True if a neuron passes, False if a neuron fails)

    '''

    # initialize metrics as dict
    metrics = {}
    metrics['cidx'] = []
    metrics['overall_pass'] = []
    metrics['rp_pass'] = []
    metrics['nc_pass'] = []
    metrics['medAmp_pass'] = []


    cids = np.unique(spikeClusters)

    if verbose:
        print("Computing metrics for %d clusters \n" % len(cids))

    for cidx in range(len(cids)):
        st = spikeTimes[spikeClusters == cids[cidx]]
        sa = spikeAmps[spikeClusters == cids[cidx]]
        neuron_pass, rp_pass, nc_pass, ma_pass = computeQualityMetrics(st, sa, verbose = False )

        metrics['cidx'].append(cids[cidx])
        metrics['overall_pass'].append(neuron_pass)
        metrics['rp_pass'].append(rp_pass)
        metrics['nc_pass'].append(nc_pass)
        metrics['medAmp_pass'].append(ma_pass)

    return metrics






def plotQualityMetrics(spikeTimes, spikeAmps):
    plotSlidingRP(spikeTimes)
    plotNoiseCutoff(spikeAmps,spikeTimes)



def using_mpl_scatter_density(fig, ax1, ax2, x, y):
    density = ax1.scatter_density(x, y, cmap=white_viridis)
    fig.colorbar(density, ax=ax2, label='Number of points per pixel')


def plotNoiseCutoff(amps, spikeTimes, n_bins=100, percent_threshold=0.10):
    '''
    Parameters
    ----------
    amps : ndarray_like
        The amplitudes (in uV) of the spikes. (spikes.amps)

    spikeTimes : ndarray_like
        The spike times. (spikes.times)

    '''

    # Set up colormap and scatter density for noise_cutoff plots
    white_viridis = LinearSegmentedColormap.from_list('white_viridis', [
        (0, '#ffffff'),
        (1e-20, '#440053'),
        (0.2, '#404388'),
        (0.4, '#2a788e'),
        (0.6, '#21a784'),
        (0.8, '#78d151'),
        (1, '#fde624'),
    ], N=256)

    nc_pass, nc_value, first_bin_height = noise_cutoff(amps)

    fig = plt.figure(figsize=(10, 5))
    ax0 = fig.add_subplot(121, projection='scatter_density')
    ax1 = fig.add_subplot(122)

    # Add scatter plot
    using_mpl_scatter_density(fig, ax0, ax1, spikeTimes, amps)
    ax0.set_ylim([0, max(amps)])
    ax0.set_xlabel('Time (s)')
    ax0.set_ylabel('Template amplitude (KS)')
    ax0.set_yticklabels([])

    # Add histogram plot
    ax1.plot()
    bins_list = np.linspace(0, np.max(amps), n_bins)  # list of bins to compute the amplitude histogram for plotting
    n, bins, patches = ax1.hist(amps, bins_list, color='#440053', orientation='horizontal')
    peak_bin_height = np.max(n)
    percent_label = np.round(first_bin_height / peak_bin_height, 2) * 100
    ax1.axvline(x=percent_threshold * peak_bin_height)
    ax1.axes.yaxis.set_ticklabels([])
    ax1.set_ylim([0, max(amps)])
    ax1.set_xlabel('Count')

    # Add titles with metric value and percent of peak; red = fail and green = pass
    if ~nc_pass:
        ax0.set_title('Cutoff metric value: ' + str(round(nc_value, 2)), color='red')
        ax1.set_title('Low bin: ' + str(round(percent_label, 2)) + '{}% of peak ', color='red')
    else:
        ax0.set_title('Cutoff metric value: ' + str(round(nc_value, 2)), color='green')
        ax1.set_title('Low bin: ' + str(round(percent_label, 2)) + '{}% of peak ', color='green')

    fig.show()


def plotSlidingRP(spikeTimes, params=None, symmetrize=True):
    '''


    Parameters
    ----------
    spikeTimes : numpy.ndarray
        array of spike times (ms)
    params : dict
        params.binSizeCorr : bin size for ACG, usually set to 1/sampleRate (s)
        params.sampleRate : sample rate of the recording (Hz)

    Returns
    -------
    None.

    '''

    if params is None:
        clusterlabel = False
    else:
        clusterlabel = params['clusterLabel']

    [maxConfidenceAt10Cont, minContWith90Confidence, timeOfLowestCont,
     nSpikesBelow2, confMatrix, cont, rp, nACG,
     firingRate, xx] = slidingRP(spikeTimes)

    fig, axs = plt.subplots(1, 2, figsize=(10, 5))

    ax = axs[0]
    acg = nACG[0:len(rp)]  # only use values of the acg for which the metric was computed

    if symmetrize:
        rp_plot = np.concatenate((-1 * np.flip(rp), rp)) * 1000
        acg_plot = np.concatenate((np.flip(acg), acg))
        ax.plot(rp_plot, acg_plot, color='k')
        ax.set_ylim(bottom=0)
        # ax.bar(rp_plot, acg_plot, width = np.diff(rp)[0]*1000, color = 'k', edgecolor = (1, 0, 0, 0)) #TODO width??
    else:
        rp_plot = rp * 1000
        acg_plot = acg
        ax.bar(rp_plot, acg_plot, width=np.diff(rp_plot)[0] * 1000, color='k', edgecolor=(1, 0, 0, 0))  # TODO width??

    ax.set_xlabel('Time from spike (ms)')
    ax.set_ylabel('ACG count (spks)')
    if clusterlabel:
        t1 = ('Cluster #%d: FR=%.2f' % (params['cidx'][0], firingRate))
    else:
        t1 = ('FR=%.2f' % firingRate)

    ax.spines['right'].set_visible(False)
    ax.spines['top'].set_visible(False)

    # Plot confidence matrix for slidingRP
    ax = axs[1]
    c = ax.imshow(confMatrix, extent=[rp[0] * 1000, rp[-1] * 1000, cont[0], cont[-1]], aspect='auto', vmin=0, vmax=100,
                  origin='lower')
    ax.spines['right'].set_visible(False)
    ax.spines['top'].set_visible(False)
    ax.spines['left'].set_visible(False)
    ax.spines['bottom'].set_visible(False)
    cbar = fig.colorbar(c, ax=ax, location='right')
    cbar.set_label('Confidence (%)')
    ax.invert_yaxis()
    ax.plot([rp[0] * 1000, rp[-1] * 1000], [10, 10], 'r', linewidth=1)

    # compute and plot 90%contour line
    if ~np.isnan(timeOfLowestCont):
        ax.plot(timeOfLowestCont * 1000 * np.array([1, 1]), [cont[0], cont[-1]], 'r', linewidth=1)

        # compute the conf=90 contour
        # zeropad confMatrix
        z = np.zeros((np.shape(confMatrix)[0] + 1, np.shape(confMatrix)[1]))
        z[1:, :] = confMatrix

        ii = np.argmax(z > 90, 0).astype(float)
        ii[ii == 0] = np.nan
        contContour = np.empty(np.shape(ii));
        contContour[:] = np.nan
        contContour[~np.isnan(ii)] = cont[(ii[~np.isnan(ii)] - 1).astype(int)]
        ax.plot(rp * 1000, contContour, 'r', linewidth=2)

    ax.set_xlabel('Time from spike (ms)')
    ax.set_xlim([0, 10])
    ax.set_ylabel('Contamination (%)')
    t2 = ('max conf = %.2f%%, min cont = %.1f%%, time = %.2f ms' % (
    maxConfidenceAt10Cont, minContWith90Confidence, timeOfLowestCont * 1000))

    # Add title with appropriate color
    if minContWith90Confidence >= 10:
        axs[0].set_title(t1, color='r')
        axs[1].set_title(t2, color='r')
    elif nSpikesBelow2 == 0:
        axs[0].set_title(t1, color='b')
        axs[1].set_title(t2, color='b')
    elif np.isnan(minContWith90Confidence):
        axs[0].set_title(t1, color='r')
        axs[1].set_title(t2, color='r')
    else:
        axs[0].set_title(t1, color='g')
        axs[1].set_title(t2, color='g')

    fig.show()




