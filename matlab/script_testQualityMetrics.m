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


%% run Sliding RP with spike times


%% run noise cutoff after re-scaling amplitudes


%% load raw mean waveforms to get true amplitude in uV of each cluster


%% combine metrics into table

% table columns: 
% - clusterID
% - rp_pass (whether RP test passes)
% - rp_minContamination (at 90% confidence)
% - rp_maxConfidence (at 10% contamination)
% - rp_rpTime (time of min Contamination)
% - nc_pass (whether noise cutoff passes)
% - nc_value (test statistic)
% - amp_pass (whether amp > 50 uV)
% - amp_value (median peak-peak amplitude of the unit in uV)