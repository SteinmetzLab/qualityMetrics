%% paths
githubDir = 'C:\Users\noamroth\Documents\GitHub';

addpath(genpath(fullfile(githubDir, 'npy-matlab'))) % https://github.com/kwikteam/npy-matlab
addpath(genpath(fullfile(githubDir, 'qualityMetrics'))) % https://github.com/SteinmetzLab/qualityMetrics
addpath(genpath(fullfile(githubDir, 'spikes'))) % https://github.com/cortex-lab/spikes


dataRoot = 'E:\Hopkins_CortexLab\';

%% load sample data

spikeTimes_samps = readNPY(fullfile(dataRoot, 'spike_times.npy')); 
clu = readNPY(fullfile(dataRoot, 'spike_clusters.npy')); 
stTemps = readNPY(fullfile(dataRoot, 'spike_templates.npy')); 
temps = readNPY(fullfile(dataRoot, 'templates.npy')); 
tempScAmps = readNPY(fullfile(dataRoot, 'amplitudes.npy')); % NOT in units of uV

spikeTimes = double(spikeTimes_samps)/30000;
cids = unique(clu);


%% run Sliding RP metric with spike times

verbose = true;
params = struct();
slidingRPStruct = run_slidingRP_All(spikeTimes, clu,verbose,params);


%% run Noise Cutoff metric after re-scaling amplitudes

% rescale amps to avoid issues with merged clusters; they will still be in arbitrary units (not uV)
% first find peak-to-peak amplitude per channel, then per template
ppAmpPerChan = squeeze(max(temps,[],2)-min(temps, [], 2));
[ppAmpPerTemp, maxChanPerTemp] = max(ppAmpPerChan, [], 2);
% now multiply the tempScAmps by their template's p-p amp
rescAmps = tempScAmps.*ppAmpPerTemp(stTemps+1);

% compute noise cutoff metric on rescaled amps
noiseCutoffStruct = run_noiseCutoff_All(rescAmps,clu,verbose,params);


%% load raw mean waveforms to get true amplitude in uV of each cluster,
% then compute median amplitude across sample waveforms

gain = 2.34; % 2.34 uV/bit for NP 1.0, 0.76 uV/bit for NP 2.0

nWF = 25; % number of waveforms to load for each cluster

gwfparams.dataDir = dataRoot;    % KiloSort/Phy output folder
d = dir(fullfile(dataRoot, '*.ap*bin')); 
gwfparams.fileName = d(1).name;         % .dat file containing the raw 
gwfparams.dataType = 'int16';            % Data type of .dat file (this should be BP filtered)
gwfparams.nCh = 385;                      % Number of channels that were streamed to disk in .dat file
gwfparams.wfWin = [-40 41];              % Number of samples before and after spiketime to include in waveform
gwfparams.nWf = nWF;                    % Number of waveforms per unit to pull out
gwfparams.spikeTimes =    spikeTimes_samps; % Vector of cluster spike times (in samples) same length as .spikeClusters
gwfparams.spikeClusters = clu; % Vector of cluster IDs (Phy nomenclature)   same length as .spikeTimes

wf = getWaveForms(gwfparams);

% compute median amplitude across waveform samples
wfs = wf.waveForms; 
wfs = wfs-nanmean(nanmean(nanmean(wfs,1),2), 4); %subtract mean across units and waveforms and time points (mean activity per channel)
wfs = wfs-nanmedian(wfs,3); %subtract the activity of the unit, waveform, timepoint 
wfsMean = squeeze(nanmean(wfs,2)); % mean across waveforms to compute max channel
wfAmpsPerChan = squeeze(nanmax(wfsMean,[],3))-squeeze(nanmin(wfsMean,[],3)); %peak to peak amplitude per channel
[~,maxChan] = nanmax(wfAmpsPerChan,[],2); %indices of max channel for each unit

%for each unit, take all of its waveforms only on its max channel 
for i = 1:size(wfs,1)
    wfMaxChan(i,:,:) = squeeze(wfs(i,:,maxChan(i),:)); 
end

% compute peak-to-peak amplitude for each sample waveform 
wfAmpMaxChan = squeeze(nanmax(wfMaxChan,[],3))-squeeze(nanmin(wfMaxChan,[],3)); 

%multiply by gain factor: 2.34 for NP 1.0:
wfAmp = wfAmpMaxChan * 2.34; 

% median amplitude across sample waveforms for each unit
amp_value = nanmedian(wfAmp,2);
amp_pass = amp_value > 50; %median amplitude threshold is 50 uV

% append these to a struct (to combine all 3 metrics into a table below)
medianAmpStruct.amp_pass = amp_pass;
medianAmpStruct.amp_value = amp_value;


%% combine metrics into table

% table columns: 
% - clusterID
% - overall_pass (whether all 3 metrics pass)
% - rp_pass (whether RP test passes)
% - rp_minContamination (at 90% confidence)
% - rp_maxConfidence (at 10% contamination)
% - rp_rpTime (time of min Contamination)
% - nc_pass (whether noise cutoff passes)
% - nc_value (test statistic)
% - amp_pass (whether amp > 50 uV)
% - amp_value (median peak-peak amplitude of the unit in uV)

%convert structs to tables and concatenate:
srpTable = struct2table(slidingRPStruct);
ncTable = struct2table(noiseCutoffStruct);
maTable = struct2table(medianAmpStruct);

%remove clusterID column from nctable (redundant with first column of srptable)
ncTable(:,1)=[];

%concatenate tables and add overall_pass column
metricsTable = [srpTable ncTable maTable];
metricsTable.overall_pass = metricsTable{:,2} & metricsTable{:,6} & metricsTable{:,8};
metricsTable = metricsTable(:,[1 width(metricsTable) 2:(width(metricsTable)-1)]);


metricsTable