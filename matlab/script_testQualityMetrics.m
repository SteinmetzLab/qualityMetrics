%% paths
% githubDir = 'C:\Users\Steinmetz Lab User\Documents\GitHub\SteinmetzLab';
addpath(genpath(fullfile(githubDir, 'npy-matlab')))
addpath(genpath(fullfile(githubDir, 'qualityMetrics')))
addpath(genpath(fullfile(githubDir, 'slidingRefractory')))
addpath(genpath(fullfile(githubDir, 'noiseCutoff')))

addpath(genpath(fullfile('C:\Users\Steinmetz Lab User\Documents\MATLAB', 'spikes')))

dataRoot = 'E:\Hopkins_CortexLab\';

%% load sample data

spikeTimes_samps = readNPY(fullfile(dataRoot, 'spike_times.npy')); 
clu = readNPY(fullfile(dataRoot, 'spike_clusters.npy')); 
tempScAmps = readNPY(fullfile(dataRoot, 'amplitudes.npy')); % NOT in units of uV

allst = double(spikeTimes_samps)/30000;
cids = unique(clu);


%%
cidx = 0;
spikeTimes = allst(clu==cidx);
spikeAmps = tempScAmps(clu==cidx);


% test the three main functions of computeQualityMetrics
 
% test computeQualityMetrics by computing QC for one neuron
[neuron_pass, rp_pass, nc_pass] = computeQualityMetrics(spikeTimes, spikeAmps);

% test plotQualityMetrics for one neuron
plotQualityMetrics(spikeTimes,spikeAmps);

% test computeQualityMetrics_All for the full population of neurons
%(Note: this can take a few minutes to run)
spikeClusters = clu;
metrics = computeQualityMetrics_All(allst,allamps,spikeClusters);

%%

gain = 2.34; % 2.34 uV/bit for NP 1.0, 0.76 uV/bit for NP 2.0

nWF = 50; % number of waveforms to load for each cluster

gwfparams.dataDir = dataRoot;    % KiloSort/Phy output folder
d = dir(fullfile(dataRoot, '*.ap.bin')); 
gwfparams.fileName = d(1).name;         % .dat file containing the raw 
gwfparams.dataType = 'int16';            % Data type of .dat file (this should be BP filtered)
gwfparams.nCh = 385;                      % Number of channels that were streamed to disk in .dat file
gwfparams.wfWin = [-40 41];              % Number of samples before and after spiketime to include in waveform
gwfparams.nWf = nWF;                    % Number of waveforms per unit to pull out
gwfparams.spikeTimes =    spikeTimes_samps; % Vector of cluster spike times (in samples) same length as .spikeClusters
gwfparams.spikeClusters = clu; % Vector of cluster IDs (Phy nomenclature)   same length as .spikeTimes

% cidx = 2; 
% gwfparams.spikeTimes =    spikeTimes_samps(clu==cidx); % Vector of cluster spike times (in samples) same length as .spikeClusters
% gwfparams.spikeClusters = clu(clu==cidx); % Vector of cluster IDs (Phy nomenclature)   same length as .spikeTimes
figure; 
wf = getWaveForms(gwfparams);


wfs = wf.waveForms; 
wfs = wfs-mean(mean(mean(wfs,1),2), 4);
wfs = wfs-median(wfs,3); 
wfsMean = squeeze(mean(wfs,1)); 
wfChanAmps = squeeze(max(wfsMean,[],2))-squeeze(min(wfsMean,[],2));
wfAmps = max(wfChanAmps,[],2);
