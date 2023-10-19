%script with plots used to test waveform extraction for
%script_testQualityMetrics.m



%% paths
githubDir = 'C:\Users\noamroth\Documents\GitHub';

addpath(genpath(fullfile(githubDir, 'npy-matlab')))
addpath(genpath(fullfile(githubDir, 'spikes')))
addpath(genpath(fullfile(githubDir, 'Steinmetzlab\qualityMetrics')))
addpath(genpath(fullfile(githubDir, 'Steinmetzlab\slidingRefractory')))
addpath(genpath(fullfile(githubDir, 'SteinmetzLab\noiseCutoff')))
%%

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

nWF = 5; % number of waveforms to load for each cluster

gwfparams.dataDir = dataRoot;    % KiloSort/Phy output folder
d = dir(fullfile(dataRoot, '*.ap*bin')); 
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

%%

wfs = wf.waveForms; 
%sanity check plots
figure(1)
imagesc(squeeze(wfs(1,1,:,:)))
title('raw wfs, one unit and one waveform (y axis chans, x axis time)')
figure(2)
imagesc(squeeze(wfs(1,:,133,:))) %picked the channel based on what looked like the max waveform in the previous plot
title('raw wfs, one unit and all waveforms, one channel')

meanActivityChannels = nanmean(nanmean(nanmean(wfs,1),2), 4);
figure(3)
imagesc(squeeze(meanActivityChannels)) %picked the channel based on what looked like the max waveform in the previous plot
title('mean activity per channel (across wfs and units and time)')



wfs = wfs-nanmean(nanmean(nanmean(wfs,1),2), 4); %subtract mean across units and waveforms and time points (mean activity per channel)
figure(4)
imagesc(squeeze(wfs(1,1,:,:)))
title('raw wfs, one unit and one waveform, subtracting mean activity per channel')

wfs = wfs-nanmedian(wfs,3); %subtract the activity of the unit, waveform, timepoint 
figure(5)
imagesc(squeeze(wfs(1,1,:,:)))
title('raw wfs, one unit and one waveform, subtracting median value across channels')



wfsMean = squeeze(nanmean(wfs,2)); 
figure(6)
imagesc(squeeze(wfsMean(1,:,:)))
title('wfsmean for one unit averaged across waveforms')

%%
wfChanAmps = squeeze(nanmax(wfsMean,[],2))-squeeze(nanmin(wfsMean,[],2));
wfAmps = nanmax(wfChanAmps,[],2); %median with 5 wfs is 66.4098, with 20 is 53.32