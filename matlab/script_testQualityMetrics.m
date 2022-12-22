%% paths
githubDir = 'C:\Users\Steinmetz Lab User\Documents\GitHub\SteinmetzLab';
addpath(genpath(fullfile(githubDir, 'npy-matlab')))
addpath(genpath(fullfile(githubDir, 'qualityMetrics')))
addpath(genpath(fullfile(githubDir, 'slidingRefractory')))
addpath(genpath(fullfile(githubDir, 'noiseCutoff')))

addpath(genpath(fullfile('C:\Users\Steinmetz Lab User\Documents\MATLAB', 'spikes')))


%% load sample data

allst = readNPY('E:\Hopkins_CortexLab\spike_times.npy'); 
allamps = readNPY('E:\Hopkins_CortexLab\amplitudes.npy'); 
clu = readNPY('E:\Hopkins_CortexLab\spike_clusters.npy'); 

allst = double(allst)/30000;
cids = unique(clu);

%%
cidx = 0;
spikeTimes = allst(clu==cidx);
spikeAmps = allamps(clu==cidx);


% test the three main functions of computeQualityMetrics
 
% test computeQualityMetrics by computing QC for one neuron
[neuron_pass, rp_pass, nc_pass, ma_pass] = computeQualityMetrics(spikeTimes, spikeAmps);

% test plotQualityMetrics for one neuron
plotQualityMetrics(spikeTimes,spikeAmps);

% test computeQualityMetrics_All for the full population of neurons
%(Note: this can take a few minutes to run)
spikeClusters = clu;%np.zros(np.shape(spikeTimes))
% metrics = computeQualityMetrics_All(allst,allamps,spikeClusters);


