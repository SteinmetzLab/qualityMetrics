

#%%
from one.api import ONE
one = ONE()
from brainbox.io.one import SpikeSortingLoader

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
amps =spikes.amps[clusterInd]
st = spikes.times[clusterInd]

#%%

neuron_pass, rp_pass, nc_pass, ma_pass = computeQualityMetrics(spikeTimes, spikeAmps)

plotQualityMetrics(spikeTimes,spikeAmps)
