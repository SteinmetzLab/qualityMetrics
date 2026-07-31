% Master Auto QC script
% 260716 ill eventually need to generalize this but thats ok
% eventually should be you can load a recording, load a time stamp, and get
% the metrics from that slice

%% add repos to path

addpath(genpath('C:\...\github\spikes'))            % https://github.com/cortex-lab/spikes/
addpath(genpath('C:\...\github\npy-matlab'))        % https://github.com/kwikteam/npy-matlab

%% set subjects, probes, recordings

%% load ap and lfp traces

genFilePath = fullfile("C:","spikesort_scratch");
apFilePath  = fullfile(genFilePath,"probe00a_pre","traces_cached_seg0.raw");
lfpFilePath = fullfile(genFilePath,"probe00a_lf_pre","traces_cached_seg0.raw");

apRec   = loadRawRecording(apFilePath);
lfpRec  = loadRawRecording(lfpFilePath);

tspan = [1000 1001];
apTraces = getTraces(apRec, tspan);

lfpTraces = getTraces(lfpRec,tspan);

% recording params

fs      = 30000;                                   % Hz, from params.py
nChan   = 383;                                     % channels in recording.dat
artifact_z = 4850;

%% load kilosort results

ksdir           = fullfile("C:","spikesort_scratch","probe00a_full_ks4/");

% deprecate old loading for loadKSdir 260729
ksr.recordingPath= fullfile(ksdir, 'recording.dat');
ksr.templates    = readNPY(fullfile(ksdir, 'templates.npy'));         % [nClu x nT x nChan]
ksr.chanpos      = readNPY(fullfile(ksdir, 'channel_positions.npy')); % [nChan x 2] = [x z]
ksr.st           = double(readNPY(fullfile(ksdir, 'spike_times.npy')));      % sample indices
ksr.sc           = double(readNPY(fullfile(ksdir, 'spike_clusters.npy')));   % 0-based cluster id/spike
ksr.sp           = double(readNPY(fullfile(ksdir, 'spike_positions.npy')));    % [N x 2] = [x z] um
ksr.amp          = double(readNPY(fullfile(ksdir, 'amplitudes.npy')));         % KS template amplitude   
ksr.metrics      = readNPZ(fullfile(ksdir, 'metrics.npz'));
ksr.winv         = double(readNPY(fullfile(ksdir, 'whitening_mat_inv.npy'))); % 383 x 383

% good and mua unit labels

ksr.lbl = readtable(fullfile(ksdir, 'cluster_group.tsv'), ...
                'FileType','text','Delimiter','\t');

%% qc visualization params

% standardized time window metrics

t0      = 1800;      % start time in seconds
tWindow = 40;   % window in ms
channel = 0;    % where the channel window is at

%% raw AP QC

close all    % start clean so qcSavePDF collects only this section's figures

% raw traces
[apTracePlot, apPlotChan]     = plotRawTraces(apRec, ksr, 'ap', artifact_z, t0,'StartChannel',33,'SaveToggle',false);

% heatmap traces
heatMapTest = wallHeatmap(apRec, 'SaveToggle', false);

% RMS over time 
% turned off for testing 260729
plotBandRMS(apRec,lfpRec,'Band','both', 'SaveToggle', true); 


% power spectra

% collate AP section figures into one PDF
qcSavePDF('ap');

%% raw LFP QC

close all    % start clean so qcSavePDF collects only this section's figures

% raw channel traces
[lfpTracePlot, lfpPlotChan]     = plotRawTraces(lfpRec, ksr, 'lfp', artifact_z, t0, 'manualTimescales', [0.04 0.1 1 10 100], 'StartChannel', 33, 'SaveToggle', false);

% channel vs freq vs power
plotDepthPower(lfpRec, 'SaveToggle', false);

% collate LFP section figures into one PDF
qcSavePDF('lfp');

%% post-KS analysis

close all    % start clean so qcSavePDF collects only this section's figures

% raw traces
plotRawWithSpikes(ksr, fs, nChan, artifact_z, t0, tWindow, 'traceChannel', apPlotChan, ...
    'SaveToggle', true);

% drift map
% hasnt cleared the limiter yet
plotUnitDrift(ksr, fs, artifact_z);

% drift map from spikes toolbox
[spikeTimes, spikeAmps, spikeDepths, spikeSites] = ksDriftmap(ksdir);
figure; 
plotDriftmap(spikeTimes, spikeAmps, spikeDepths);

% raster map
ampDepthScatter(ksr, fs, artifact_z, 'SaveToggle', true);

% pdf/cdf from spikes toolbox

depthBins = 0:40:artifact_z;
ampBins = 0:30:min(max(spikeAmps),800);
recordingDur = ksr.st(end);

[pdfs, cdfs] = computeWFampsOverDepth(spikeAmps, spikeDepths, ampBins, depthBins, recordingDur);
plotWFampCDFscopy(pdfs, cdfs, ampBins, depthBins, 'SaveToggle', true);

% collate post-Kilosort section figures into one PDF
qcSavePDF('postks');

% waveform analysis
%% identify unit types

% columns: cluster_id, KSLabel
isGoodClu   = strcmp(ksr.lbl.KSLabel, 'good');           % logical per cluster row
goodIds     = ksr.lbl.cluster_id(isGoodClu);             % 0-based cluster ids, 260730 need to add +1 when indexing in matlab
muaIds      = ksr.lbl.cluster_id(~isGoodClu);

% 260730 figuring out parsing error
% python is 0-index which is what KS labels from
% matlab 1-index

% so it seems unit 375 which is a good unit is index 376 on KSLabel
% unit 377 is index 378 (good on ksLabel, but plotted wrong (probably used
% index 377 from KSLabel)
% so all units are +1 in KSLabel

% this is bc when selecting units i return the cluster unit and then
% convert back to an index? 

%% pick cluster (ripped from plotTemplateVsRaw), used for testing

ptpAll      = squeeze(max(ksr.templates,[],2) - min(ksr.templates,[],2));     % [nClu x nChan]
goodPTP     = ptpAll(isGoodClu,:);
muaPTP      = ptpAll(~isGoodClu,:);

numUnits = 2;
sampleToggle = 1;

toPlot = [];
toPlot.clusters = zeros(numUnits*2, 1);
toPlot.labels   = strings(4,1);

toPlot.labels(1:numUnits,1)   = "good";
toPlot.labels(numUnits+1:numUnits*2,1)   = "mua";

% error in indexing already at toPlot

if sampleToggle
    % "best" and "worst" good units
    % for testing
    [~, bestGood] = max(max(goodPTP,[],2)); % largest largest PTP
    [~, worstGood] = min(max(goodPTP,[],2)); % smallest largest PTP

    toPlot.clusters(1) = goodIds(bestGood);
    toPlot.clusters(numUnits) = goodIds(worstGood);

    % two muas with largest PTP
    [~, muaExPTP] = maxk(max(muaPTP,[],2),numUnits); % largest largest PTP
    toPlot.clusters(numUnits+1:numUnits*2) = muaIds(muaExPTP); % 0-indexed cluster label 
end

%% 

close all

for i = 1:size(toPlot.clusters,1)
    plotWaveform(ksr, fs, nChan, 'clu',toPlot.clusters(i),'cluLabel',toPlot.labels(i),'SaveToggle',true);
end

% testing
plotWaveform(ksr, fs, nChan, 'clu',toPlot.clusters(1),'cluLabel',toPlot.labels(1), 'showSamples', true);