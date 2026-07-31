function out = plotWaveform(ksResults, fs, nChan, varargin)
% plotTemplateVsRaw  Plot a KS4 template, averaged spike waveform, and unwhitened 
% KS4 template.
%
%   out = plotTemplateVsRaw()               % biggest-amplitude cluster
%   out = plotTemplateVsRaw('clu', 375)     % a specific 0-based cluster id
%   out = plotTemplateVsRaw('clu', 375, 'nSpikes', 1000)
%
% INPUT (name/value)
%   'ksdir'   KS4 output folder      (default C:\spikesort_scratch\probe00a_full_ks4)
%   'clu'     0-based cluster id     (default [] -> largest-amplitude template)
%   'nSpikes' max spikes to average  (default 500; sampled at random if more)
%   'nChan'   # channels in overlay grid, by template amplitude (default 6)
%   'gain_uV' int16 -> microvolts    (default 3.02734375, as in plotWaveforms)
%   'unwhiten' also draw the un-whitened template (default true)
%   'save'    save PNG into autoqc/figures (default false)
%   'seed'    RNG seed for spike sampling (default 0)
%
% OUTPUT (struct out)
%   .clu, .nSpikesUsed, .peakChan(0-based), .meanWf_uV [nT x nChan],
%   .corrPeak (corr of template vs mean on peak channel), .t_ms

    %% input parser

    p = inputParser;
    p.addRequired('ksResults',      @(x) isstruct(x));
    p.addParameter('clu',      [], @(x)isempty(x)||isscalar(x));
    p.addParameter('cluLabel', "good", @isstring);
    p.addParameter('isTest',   0, @islogical)
    p.addParameter('nChan',    6, @isnumeric);
    p.addParameter('gain_uV',  3.02734375, @isnumeric);
    p.addParameter('SaveToggle',     false, @islogical);
    p.addParameter('seed',     0, @isnumeric);
    p.addParameter('showSamples', false, @islogical);   % overlay single-spike traces
    p.addParameter('nSample',  50, @isnumeric);         % # single spikes to overlay
    p.parse(ksResults, varargin{:});
    o = p.Results;

    datPath = ksResults.recordingPath;      % int16, [nChan x nSamp], offset 0
    
    % ---------------------------------------------------------------- KS4 arrays
    [nClu, nT, ~] = size(ksResults.templates);
    z = double(ksResults.chanpos(:,2));
    
    %% select clusters
    
    ptpAll = squeeze(max(ksResults.templates,[],2) - min(ksResults.templates,[],2));     % [nClu x nChan]
    if isempty(o.clu)
        [~, row] = max(max(ptpAll,[],2)); 
        o.clu = row - 1;
    else
        row = o.clu + 1;
        assert(row>=1 && row<=nClu, 'clu %d out of range 0..%d', o.clu, nClu-1);
    end
    Wt = squeeze(ksResults.templates(row,:,:));                  % [nT x nChan] template (a.u.)
    [~, peakChan] = max(ptpAll(row,:));                % 1-based channel
    
    % ---------------------------------------------------------------- spike times for this cluster
    pre  = round((nT-1)/3);                            % samples before trough (matches template centring)
    post = nT-1-pre;                                   % samples after
    times = ksResults.st(ksResults.sc == o.clu);
    tot   = numel(times);
    assert(tot>0, 'cluster %d has no spikes', o.clu);
    % keep snippets fully inside the file
    nSamp = fileSamples(datPath, nChan);
    times = times(times>pre & times<=nSamp-post);
    % random subsample for speed
    rng(o.seed);
    
    %% spike limiter for troubleshooting

    if o.isTest
        nSpikes = 500;
        if numel(times) > nSpikes
            times = times(randperm(numel(times), nSpikes));
        end
    end
    
    %% 

    nUse = numel(times);
    
    % ---------------------------------------------------------------- pull + average raw snippets
    map = memmapfile(datPath, 'Format', {'int16',[nChan nSamp],'x'}, 'Offset', 0);
    acc = zeros(nT, nChan);                         % accumulate [time x chan]
    for k = 1:nUse
        s = times(k);
        snip = double(map.Data.x(:, (s-pre):(s+post))).';   % [nT x nChan]
        snip = snip - mean(snip(1:min(15,nT), :), 1);       % baseline = pre-spike mean
        acc  = acc + snip;
    end
    meanWf = (acc / nUse) * o.gain_uV;                 % [nT x nChan] in microvolts
    clear map
    
    % choose the channels to display: top-ptp channels of the TEMPLATE, depth-sorted
    ptpT = max(Wt,[],1) - min(Wt,[],1);
    [~, topc] = sort(ptpT, 'descend');
    showCh = topc(1:min(o.nChan, numel(topc)));
    [~, zo] = sort(z(showCh), 'descend');              % top of probe first
    showCh = showCh(zo);
    
    t_ms = ((0:nT-1) - pre) / fs * 1000;               % ms, 0 = trough
    
    % unwhitened template
    
    winv = ksResults.winv;
    Wt_uw = Wt * winv;                             % back to data space (shape)
    
    % least-squares scale of template onto the measured mean, per channel shown
    lsScale = @(a,b) (a(:).'*b(:)) / max(a(:).'*a(:), eps);   % scale s: min ||b - s*a||
    
    % =====================================================================  FIGURE
    fig = figure('Color','w','Position',[70 70 1200 760], ...
        'Name', sprintf('template vs raw — clu %d', o.clu));
    tl = tiledlayout(fig, 2, 3, 'TileSpacing','compact', 'Padding','compact');
    title(tl, sprintf(['cluster %d (%s) |  %d spikes averaged  |  peak ch %d @ %.0f um' ...
        '   —   black = measured mean waveform (uV),  red = KS template (scaled)'], ...
        o.clu, o.cluLabel, nUse, peakChan-1, z(peakChan)), 'FontWeight','bold');
    
    % pre-create the tile axes so single-spike traces can be drawn UNDERNEATH the mean
    axList = gobjects(1, numel(showCh));
    for i = 1:numel(showCh)
        axList(i) = nexttile(tl, i); hold(axList(i),'on');
    end
    
    hSample = gobjects(0);   % handle for legend (empty unless samples drawn)
    if o.showSamples
        plotSampleWaveforms(datPath, nChan, times, pre, post, showCh, ...
                            axList, t_ms, o.gain_uV, 'nSample', o.nSample, 'seed', o.seed);
        hSample = findobj(axList(1), 'Type', 'line');   % only the sample line exists yet
    end
    
    for i = 1:numel(showCh)
        ci = showCh(i);
        ax = axList(i); hold(ax,'on');
        raw = meanWf(:,ci);
        hMean = plot(ax, t_ms, raw, 'k', 'LineWidth', 1.8);         % measured, uV
        s1 = lsScale(Wt(:,ci), raw);
        hTpl = plot(ax, t_ms, s1*Wt(:,ci), 'Color',[0.85 0.10 0.10], 'LineWidth',1.4);

        % plot unwhiten template
        s2 = lsScale(Wt_uw(:,ci), raw);
        hTplUw = plot(ax, t_ms, s2*Wt_uw(:,ci), '--', 'Color',[0.10 0.35 0.85], 'LineWidth',1.1);
       
        xline(ax, 0, ':', 'Color',[0.7 0.7 0.7]);
        axis(ax,'tight'); grid(ax,'on'); box(ax,'on');
        ttl = sprintf('ch %d  (%.0f um)', ci-1, z(ci));
        if ci==peakChan, ttl = ['\bf' ttl '  [peak]']; end
        title(ax, ttl);
        if i==1
            legH = [hMean, hTpl, hTplUw];
            legL = {'mean (uV)','template','template unwhitened'};
            if o.showSamples && ~isempty(hSample)
                legH = [hSample(1), legH];
                legL = ['single spikes', legL];
            end
            legend(ax, legH, legL, ...
                    'Location','southeast','FontSize',7,'Box','off');
        end
        if i>numel(showCh)-3, xlabel(ax,'time (ms)'); end
        if mod(i-1,3)==0, ylabel(ax,'uV'); end
    end
    
    % correlation on the peak channel (shape agreement, sign-safe)
    corrPeak = corr(meanWf(:,peakChan), Wt(:,peakChan));
    
    % save
    
    if o.SaveToggle
        figType = fullfile('ks','waveforms');
        figName = sprintf('templateVsRaw_clu%03d.png', o.clu);
        qcFigSave(figType, figName);
    end
    
    
    out = struct('clu',o.clu, 'nSpikesUsed',nUse, 'peakChan',peakChan-1, ...
                 'peakDepth',z(peakChan), 'meanWf_uV',meanWf, 'corrPeak',corrPeak, ...
                 't_ms',t_ms);
    fprintf('clu %d: averaged %d/%d spikes, peak ch %d, template-vs-mean corr = %.3f\n', ...
        o.clu, nUse, tot, peakChan-1, corrPeak);
end

% ---- helper: number of time samples in an int16 [nChan x nSamp] file ------
function n = fileSamples(datPath, nChan)
    d = dir(datPath);
    n = d.bytes / (2*nChan);
    assert(n==floor(n), 'file size not divisible by 2*%d channels', nChan);
end
