%% paths
% githubDir = 'C:\Users\Steinmetz Lab User\Documents\GitHub\SteinmetzLab';
addpath(genpath(fullfile(githubDir, 'npy-matlab')))
addpath(genpath(fullfile(githubDir, 'qualityMetrics')))
addpath(genpath(fullfile(githubDir, 'slidingRefractory')))
addpath(genpath(fullfile(githubDir, 'noiseCutoff')))
addpath(genpath(fullfile(githubDir, 'spikes')))
% addpath(genpath(fullfile('C:\Users\Steinmetz Lab User\Documents\MATLAB', 'spikes')))

dataRoot = 'E:\Hopkins_CortexLab\';

%% load sample data

spikeTimes_samps = readNPY(fullfile(dataRoot, 'spike_times.npy')); 
clu = readNPY(fullfile(dataRoot, 'spike_clusters.npy')); 
stTemps = readNPY(fullfile(dataRoot, 'spike_templates.npy')); 
temps = readNPY(fullfile(dataRoot, 'templates.npy')); 
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


%% testing re-scaling of amplitudes.npy with template amps

% calculate the p-p amplitude of each template

ppAmpPerChan = squeeze(max(temps,[],2)-min(temps, [], 2));
[ppAmpPerTemp, maxChanPerTemp] = max(ppAmpPerChan, [], 2);

% multiply the tempScAmps by their template's p-p amp
rescAmps = tempScAmps.*ppAmpPerTemp(stTemps+1);

% load a bunch of raw waveforms
nCh = 385;
d = dir(fullfile(dataRoot, '*.ap.bin')); 
fileName = fullfile(dataRoot,d.name);           
dataTypeNBytes = 2; % determine number of bytes per sample
nSamp = d.bytes/(nCh*dataTypeNBytes);  % Number of samples per channel
wfNSamples = 80;
mmf = memmapfile(fileName, 'Format', {'int16', [nCh nSamp], 'x'});

rndIdx = randi(numel(rescAmps), 1000,1); 

allWF = zeros(numel(rndIdx), 384, 81);
for q = 1:numel(rndIdx)
    thisST = spikeTimes_samps(rndIdx(q)); 
    allWF(q,:,:) = mmf.Data.x(1:nCh-1,thisST-40:thisST+40);
end

%% 
% compute amp of each raw waveform on peak channel
firstSamp = allWF(:,:,1); 
allWF = allWF-median(firstSamp,1); 

fshigh = 300; fs = 30000; 
[b1, a1] = butter(3, fshigh/fs*2, 'high');
for q = 1:numel(rndIdx)
    allWF(q,:,:) = filtfilt(b1, a1, squeeze(allWF(q,:,:))')';
end

allWFampsPerChan = max(allWF, [], 3) - min(allWF, [], 3); 

chanPerSpike = maxChanPerTemp(stTemps(rndIdx)+1); 

ampsPerSpike = allWFampsPerChan(sub2ind([1000,384], 1:1000, chanPerSpike'));

% plot rescaled amps vs true amps
figure; 
subplot(2,2,1); 
cx = corr(ampsPerSpike', rescAmps(rndIdx));
plot(ampsPerSpike, rescAmps(rndIdx), '.')
xlabel('measured amp'); ylabel('resc amp'); 
title(cx)

% plot un-rescaled amps vs true amps
subplot(2,2,2); 
cx = corr(ampsPerSpike', tempScAmps(rndIdx));
plot(ampsPerSpike, tempScAmps(rndIdx), '.'); 
xlabel('measured amp'); ylabel('original temp sc amp'); 
title(cx)

subplot(2,2,3); 
plot(rescAmps(rndIdx), tempScAmps(rndIdx), '.'); 
xlabel('resc amp');ylabel('original temp sc amp');