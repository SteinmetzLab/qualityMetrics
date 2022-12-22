%% paths
githubDir = 'C:\Users\Steinmetz Lab User\Documents\GitHub\';
addpath(genpath(fullfile(githubDir, 'npy-matlab')))
addpath(genpath(fullfile(githubDir, 'slidingRefractory')))
addpath(genpath(fullfile(githubDir, 'SteinmetzLab\qualityMetrics')))
addpath(genpath(fullfile(githubDir, 'noiseCutoff')))

addpath(genpath(fullfile('C:\Users\Steinmetz Lab User\Documents\MATLAB', 'spikes')))


%% load

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
neuron_pass, rp_pass, nc_pass, ma_pass = computeQualityMetrics(spikeTimes, spikeAmps)

% test computeQualityMetrics_All for a 'population' of one neuron
% spikeClusters = np.zeros(np.shape(spikeTimes)) metrics = computeQualityMetrics_All(spikeTimes,spikeAmps,spikeClusters)
% 
% #test plotQualityMetrics
% plotQualityMetrics(spikeTimes,spikeAmps)