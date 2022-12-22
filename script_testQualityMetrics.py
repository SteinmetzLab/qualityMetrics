
#To ensure qualityMetrics is in your system path:
'''
import sys
sys.path.append('filepath of qualityMetrics')
'''

import sys
sys.path.append(r'C:\Users\noamroth\Documents\GitHub\SteinmetzLab\qualityMetrics')
#%%
#imports for downloading sample IBL data and computing Quality metrics
import numpy as np
from one.api import ONE
one = ONE()
from brainbox.io.one import SpikeSortingLoader
from computeQualityMetrics import *
#%%

#for testing, load one IBL example neuron
one = ONE()
#list of pids to check
pids = ['ce397420-3cd2-4a55-8fd1-5e28321981f4',
        'ce397420-3cd2-4a55-8fd1-5e28321981f4',
        'ce397420-3cd2-4a55-8fd1-5e28321981f4']
cids = [410,
        659,
        811]

ind = 2
pid = pids[ind]
cid = cids[ind]

sl = SpikeSortingLoader(pid=pid, one=one)
spikes, clusters, channels = sl.load_spike_sorting()
# clusters = sl.merge_clusters(spikes, clusters, channels)

clusterInd = np.where(spikes.clusters == cid)[0]
spikeAmps =spikes.amps[clusterInd]
spikeTimes = spikes.times[clusterInd]

#%%
#test the three main functions of computeQualityMetrics


#test computeQualityMetrics by computing QC for one neuron
neuron_pass, rp_pass, nc_pass, ma_pass = computeQualityMetrics(spikeTimes, spikeAmps)

#test computeQualityMetrics_All for a 'population' of one neuron
spikeClusters = np.zeros(np.shape(spikeTimes))
metrics = computeQualityMetrics_All(spikeTimes,spikeAmps,spikeClusters)

#test plotQualityMetrics
plotQualityMetrics(spikeTimes,spikeAmps)






