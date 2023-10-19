%% paths
% githubDir = 'C:\Users\Steinmetz Lab User\Documents\GitHub\SteinmetzLab';
addpath(genpath(fullfile(githubDir, 'npy-matlab')))
addpath(genpath(fullfile(githubDir, 'qualityMetrics')))
addpath(genpath(fullfile(githubDir, 'slidingRefractory')))
addpath(genpath(fullfile(githubDir, 'noiseCutoff')))
addpath(genpath(fullfile(githubDir, 'spikes')))

dataRoot = 'E:\Hopkins_CortexLab\';

%% load sample data

spikeTimes_samps = readNPY(fullfile(dataRoot, 'spike_times.npy')); 
clu = readNPY(fullfile(dataRoot, 'spike_clusters.npy')); 
stTemps = readNPY(fullfile(dataRoot, 'spike_templates.npy')); 
temps = readNPY(fullfile(dataRoot, 'templates.npy')); 
tempScAmps = readNPY(fullfile(dataRoot, 'amplitudes.npy')); % NOT in units of uV

allst = double(spikeTimes_samps)/30000;
cids = unique(clu);


