function info = plotBandRMS(apRec, lfpRec, varargin)
% PLOTBANDRMS  Time x depth RMS-activity map for the AP and/or LFP bands.
%
% Bins the recording by time and by depth along the probe, computes an RMS
% average of the voltage in each (depth-bin, time-bin) block, and draws it as
% a heatmap: depth on the y-axis, time on the x-axis, RMS voltage (uV) as the
% color (viridis: purple = low activity, yellow = high activity). With the
% default 'both' it stacks the AP band (top) over the LFP band (bottom) on a
% shared time axis so you can compare fast (spiking) and slow (field) activity.
%
% Usage
%   plotBandRMS                                  % probe00a defaults, AP + LFP
%   plotBandRMS('Band','ap')                     % AP band only
%   plotBandRMS('t_start_s',1000,'win_s',60)
%   plotBandRMS('timeBin_s',0.5,'depthBin_um',20)
%
% Name/value options (all optional)
%   'apDir'     AP recording folder (SpikeInterface .raw + binary.json)
%               (default C:\spikesort_scratch\probe00a_pre, 30 kHz, 383 ch)
%   'lfpDir'    LFP recording folder
%               (default C:\spikesort_scratch\probe00a_lf_pre, 2500 Hz, 324 ch)
%   'Band'      'both' (default) | 'ap' | 'lfp'
%   'timeBin_s' time-bin width, seconds (x resolution)        (default 1)
%   'depthBin_um' depth-bin width, um (y resolution). Channels whose sites
%               fall in a bin are pooled before the RMS.       (default 40)
%   'z_max'     mask channels with depth > z_max (top artifact band); Inf off
%               (default 4850, the probe00a convention)
%   'CAR'       subtract across-channel mean per sample first  (default false;
%               probe00a *_pre are already CMR'd)
%   'climPct'   robust color saturation percentiles [lo hi]    (default [2 98])
%   'clim'      explicit color limits [lo hi] uV; overrides climPct. Scalar or
%               [] auto. For 'both', a 2x2 [apLo apHi; lfpLo lfpHi] sets each.
%   'SaveToggle'  save PNG to figures\bandRMS                  (default false)
%
% OUTPUT (struct info) — one field per band drawn ('ap' and/or 'lfp'), each:
%   .depth    [nDepthBin x 1]  depth-bin centers (um, ascending)
%   .t        [1 x nTimeBin]   time-bin centers (s, absolute)
%   .rms      [nDepthBin x nTimeBin] RMS voltage (uV); NaN for empty depth bins
%   .fs, .nKept, .nChan
%
% RMS definition: within each (depth-bin, time-bin) block the per-channel DC
% (mean over that time bin) is removed, then rms = sqrt(mean(x.^2)) is taken
% over ALL samples of ALL channels pooled in the block. Time bins are read and
% reduced one at a time, so long windows stay memory-bounded even at 30 kHz.
%
% No toolboxes required. Uses autoqc infra: loadRawRecording, getTraces,
% recordingGeometry, viridisColors.
%
% See also plotDepthPower, plotLFPWelch, wallHeatmap.

    % ---------------------------------------------------------------- parameters
    p = inputParser;
    % p.addParameter('apDir',  'C:\spikesort_scratch\probe00a_pre',    @(s)ischar(s)||isstring(s));
    % p.addParameter('lfpDir', 'C:\spikesort_scratch\probe00a_lf_pre', @(s)ischar(s)||isstring(s));
    p.addParameter('Band',   'both', @(s) any(strcmpi(s,{'both','ap','lfp'})));
    p.addParameter('timeBin_s',   1,    @isnumeric);
    p.addParameter('depthBin_um', 40,   @isnumeric);
    p.addParameter('z_max',       4850, @isnumeric);
    p.addParameter('CAR',         false, @islogical);
    p.addParameter('climPct',     [2 98], @isnumeric);
    p.addParameter('clim',        [],   @isnumeric);
    p.addParameter('SaveToggle',  false, @islogical);
    p.parse(varargin{:});
    o = p.Results;
    
    bands = lower(o.Band);
    switch bands
        case 'both', bandList = {'ap','lfp'};
        case 'ap',   bandList = {'ap'};
        case 'lfp',  bandList = {'lfp'};
    end
    
    % compute RMS
    
    info = struct();
    for k = 1:numel(bandList)
        b = bandList{k};
        if strcmp(b,'ap')
            info.(b) = computeBandRMS(apRec, o);
        else 
            info.(b) = computeBandRMS(lfpRec, o);
        end
    end
    
    % ---------------------------------------------------------------- color limits
    clim = resolveClim(o.clim, bandList);
    for k = 1:numel(bandList)
        b = bandList{k};
        if isempty(clim.(b))
            cl = prctile_(info.(b).rms(:), o.climPct);
            if ~(numel(cl)==2) || diff(cl) <= 0
                v = info.(b).rms(:); v = v(isfinite(v));
                cl = [min(v) max(v)];
            end
            clim.(b) = cl;
        end
    end
    
    % ---------------------------------------------------------------- plot
    figure('Color','w','Position',[100 100 1000 max(420,360*numel(bandList))]);
    tl = tiledlayout(numel(bandList), 1, 'TileSpacing','compact', 'Padding','compact');
    for k = 1:numel(bandList)
        b = bandList{k};
        S = info.(b);
        ax = nexttile;
        im = imagesc(ax, S.t, S.depth, S.rms, clim.(b));
        set(im, 'AlphaData', isfinite(S.rms));      % empty depth bins -> transparent
        set(ax, 'YDir','normal', 'Color',[0.9 0.9 0.9]);
        colormap(ax, viridisColors(256));           % purple = low, yellow = high
        ylabel(ax, 'depth along probe (\mum)');
        if k == numel(bandList), xlabel(ax, 'time (s)'); end
        cb = colorbar(ax); cb.Label.String = sprintf('%s RMS (\\muV)', upper(b));
        title(ax, sprintf('%s band  |  %d/%d ch, %g \\mum bins, %g s bins  @ %g Hz', ...
            upper(b), S.nKept, S.nChan, o.depthBin_um, o.timeBin_s, S.fs));
        box(ax,'on');
    end
    title(tl, sprintf('probe00a band RMS activity'));

    % save
    
    if o.SaveToggle
        figType = 'ap';
        figName = sprintf('bandRMS_%s.png', bands);
        qcFigSave(figType, figName);
    end
end

%% local functions

function S = computeBandRMS(rec, o)
    % Time x depth RMS map for one band, accumulated one time bin at a time.
    fs  = rec.fs;
    
    % ---- geometry / depth, drop top artifact band ----

    geom = recordingGeometry(rec);
    nChan = rec.nChannels;
    if geom.nChannels < nChan
        nChan = geom.nChannels; 
    end
    
    z    = geom.y(1:nChan);
    keep = find(z <= o.z_max);
    zk   = z(keep);
    
    % ---- depth bins ----
    edges   = min(zk) : o.depthBin_um : max(zk) + o.depthBin_um;
    if numel(edges) < 2, edges = [min(zk) max(zk)+eps]; end
    binOf   = discretize(zk, edges);               % depth-bin index per kept channel
    valid   = ~isnan(binOf);
    keep    = keep(valid); binOf = binOf(valid);
    nDepth  = numel(edges) - 1;
    centers = edges(1:end-1)' + o.depthBin_um/2;
    
    % ---- time bins ----
    t0   = 0;
    % 1 second bins
    % nBin = max(1, floor(o.win_s / o.timeBin_s));
    % nBin = max(1, floor(rec.nSamples / rec.fs));
    nBin = 100;
    timebinSize_s = (rec.nSamples/rec.fs)/nBin;   % timebin size in seconds
    tCenters = t0 + ((1:nBin) - 0.5) * timebinSize_s;

    % ---- accumulate RMS: for each time bin, pool channels within each depth bin ----
    rms = nan(nDepth, nBin);
    for ti = 1:nBin
        ta = t0 + (ti-1)*timebinSize_s;
        tb = ta + timebinSize_s;
        X  = getTraces(rec, [ta tb], keep, 'uV');  % [nKeep x nSamp] uV, double
        if isempty(X) || size(X,2) < 1, continue; end
        if o.CAR, X = X - mean(X, 1); end
        X  = X - mean(X, 2);                        % per-channel DC over the bin
        ss = sum(X.^2, 2);                          % [nKeep x 1] sum of squares
        n  = size(X, 2);                            % samples this bin (same all ch)
        % pool over channels sharing a depth bin: RMS = sqrt(sum ss / (nCh*n))
        ssBin  = accumarray(binOf, ss, [nDepth 1], @sum, NaN);
        cntBin = accumarray(binOf, 1,  [nDepth 1], @sum, 0) * n;
        rms(:, ti) = sqrt(ssBin ./ cntBin);
    end

    S = struct('depth', centers, 't', tCenters, 'rms', rms, ...
           'fs', fs, 'nKept', numel(keep), 'nChan', nChan);
end

% ======================================================================
function clim = resolveClim(c, bandList)
% Normalize the 'clim' option into a per-band struct (empty => auto).
clim = struct();
for k = 1:numel(bandList), clim.(bandList{k}) = []; end
if isempty(c), return; end
if isvector(c) && numel(c) == 2
    for k = 1:numel(bandList), clim.(bandList{k}) = c(:)'; end
elseif isequal(size(c), [2 2]) && numel(bandList) == 2
    clim.(bandList{1}) = c(1,:);
    clim.(bandList{2}) = c(2,:);
else
    warning('plotBandRMS:clim', 'Ignoring clim: expected [lo hi] or 2x2.');
end
end

% ======================================================================
function q = prctile_(x, pcts)
% Percentile without the Statistics Toolbox (linear interpolation).
x = sort(x(~isnan(x(:))));
if isempty(x), q = nan(size(pcts)); return; end
n = numel(x);
q = zeros(size(pcts));
for i = 1:numel(pcts)
    r  = pcts(i)/100 * (n-1) + 1;
    lo = max(floor(r),1); hi = min(ceil(r),n);
    q(i) = x(lo) + (r-lo)*(x(hi)-x(lo));
end
end
