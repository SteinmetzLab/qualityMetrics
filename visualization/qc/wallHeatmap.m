function info = wallHeatmap(rec, varargin)
%WALLHEATMAP Heatmap of Neuropixel channel activity by probe for selected window
%   Color = voltage (µV), rows ordered by depth, diverging colormap centered
%   at 0. Draws one panel per probe per shank. 
% USAGE
%   wallHeatmap(rec)                                     % one probe, all shanks
%   wallHeatmap(rec, 't_start_s',30, 'win_ms',10)
%   wallHeatmap({recA, recB, recC}, 'z_max',4850)        % 3 probes side by side
%   info = wallHeatmap(rec, 'SaveDir', ...
%                      'C:\Users\stein\alex\autoqc\figures\wallHeatmap');
%
% INPUTS
%   rec   a recording struct from loadRawRecording, OR a cell array of them
%         (one per probe). Each must work with getTraces / recordingGeometry.
%
% NAME-VALUE OPTIONS
%   'Geometry'    geom struct from recordingGeometry (or a cell array, one per
%                 rec). Default: computed per rec via recordingGeometry.
%   't_start_s'   window start, seconds                       (default 1000)
%   'win_ms'      window length, milliseconds                 (default 100)
%   'CLimUV'      symmetric color saturation, ±µV. [] => robust auto (percentile
%                 of |voltage|, pooled across panels).         (default 100)
%   'CLimPct'     percentile for the auto color limit          (default 99.5)
%   'z_max'       mask channels with depth y > z_max (µm): they stay on the
%                 depth axis but render as a gray band (out-of-brain / top
%                 artifact band). Scalar (all panels), or a cell array {1 x
%                 nProbe} of scalars/per-shank vectors. Inf = no mask.
%                                                              (default Inf)
%   'ProbeNames'  cellstr of probe names for the headers, one per rec
%                                                    (default: model or 'probeK')
%   'Title'       main figure title.       (default: auto from probe/shank count)
%   'SaveDir'     folder to save a PNG into; '' to skip.       (default '')
%
% OUTPUT (struct info)
%   .vlim         color limit used, µV
%   .t_start_s    window start shown (s)
%   .win_ms       window length (ms)
%   .panels       struct array (metadata only, no traces): probe, probeName,
%                 shank, channels (depth-sorted, 1-based), zRange, nChannels,
%                 nMasked
%
% No toolboxes required (local percentile + colormap fall-backs at bottom).
%
% See also getTraces, recordingGeometry, preprocessTraces, plotRawWithSpikes.

% ---------------------------------------------------------------- parameters
p = inputParser;
p.addParameter('Geometry',   []);
p.addParameter('t_start_s',  1000,  @isnumeric);
p.addParameter('win_ms',     100,   @isnumeric);
p.addParameter('CLimUV',     100,   @(x) isempty(x) || isnumeric(x));
p.addParameter('CLimPct',    99.5,  @isnumeric);
p.addParameter('z_max',      Inf);
p.addParameter('ProbeNames', {},    @iscell);
p.addParameter('SaveToggle', false);
p.parse(varargin{:});
o = p.Results;

% ---- normalize rec + geometry to per-probe cell arrays ----------------
if ~iscell(rec)
    recs = {rec}; 
else 
    recs = rec; 
end
nProbe = numel(recs);

if isempty(o.Geometry)
    geoms = cell(1, nProbe);
    for i = 1:nProbe, geoms{i} = recordingGeometry(recs{i}); end
elseif ~iscell(o.Geometry)
    geoms = {o.Geometry};
else
    geoms = o.Geometry;
end
if numel(geoms) ~= nProbe
    error('wallHeatmap:geom', 'Need one Geometry per probe (%d given, %d probes).', ...
        numel(geoms), nProbe);
end

% probe display names
if isempty(o.ProbeNames)
    names = cell(1, nProbe);
    for i = 1:nProbe
        % 260716 this is weird, renames it oddly
        m = sprintf('probe%02d', i-1);
        names{i} = m;
    end
else
    names = o.ProbeNames;
end

% ---------------------------------------------------------------- gather panels
t0 = o.t_start_s;
t1 = t0 + o.win_ms/1000;

panels = struct('probe',{}, 'probeName',{}, 'shank',{}, 'channels',{}, ...
                'X',{}, 'z',{}, 't_ms',{}, 'nMasked',{});
for i = 1:nProbe
    g  = geoms{i};
    fs = recs{i}.fs;
    shanks = unique(g.shank(:))';                 % one panel per shank
    for si = 1:numel(shanks)
        s  = shanks(si);
        ch = find(g.shank == s);                  % ALL channels on this shank
        if isempty(ch), continue; end
        X = getTraces(recs{i}, [t0 t1], ch, 'uV');      % [nCh x nSamp], µV

        if ~rec.preprocessed
            % preprocess
            [X, ff] = preprocessTraces(X, rec.fs, 'ap', 'doCMR', doCMR);
        else
            ff.passHz       = 300; 
            ff.bandpassType = 'highpass';
            ff.CMR          = true;
        end

        z = g.y(ch);
        [zs, ord] = sort(z);                            % ascending depth from tip
        X  = X(ord, :);
        ch = ch(ord);

        % mask the out-of-brain / top artifact band -> NaN -> gray on the axis
        thr = zmaxFor_(o.z_max, i, si);
        masked = zs > thr;
        X(masked, :) = NaN;

        nS = size(X, 2);
        panels(end+1) = struct( ...  %#ok<AGROW>
            'probe', i, 'probeName', names{i}, 'shank', s, ...
            'channels', ch, 'X', X, 'z', zs, 't_ms', (0:nS-1)/fs*1000, ...
            'nMasked', nnz(masked));
    end
end
if isempty(panels)
    error('wallHeatmap:empty', 'No channels to plot (check geometry).');
end
nPan = numel(panels);

% ---------------------------------------------------------------- color limit
if isempty(o.CLimUV)
    pooled = cell2mat(arrayfun(@(pnl) abs(pnl.X(:)), panels, 'uni', 0)');
    vlim = prctile_(pooled, o.CLimPct);
else
    vlim = o.CLimUV;
end
if ~isfinite(vlim) || vlim <= 0, vlim = 1; end

% shared depth axis so panels are directly comparable
zmin = min(arrayfun(@(pnl) min(pnl.z), panels));
zmax = max(arrayfun(@(pnl) max(pnl.z), panels));

% per-probe colors for the headers / shank labels
probeCol = probeColors_(nProbe);
grayC    = [0.93 0.93 0.93];        % masked (out-of-brain) band color

% ---------------------------------------------------------------- plot
figW = min(2400, 150*nPan + 260);
fig = figure('Color', 'w', 'Position', [60 60 figW 860]);
tl  = tiledlayout(fig, 1, nPan, 'TileSpacing', 'tight', 'Padding', 'compact');
tl.OuterPosition = [0 0 1 0.90];    % leave a top strip for probe headers + title
cmap = blueWhiteRed_(256);          % RdBu_r: blue = negative, white = 0, red = +

axArr = gobjects(1, nPan);
for k = 1:nPan
    ax = nexttile(tl);
    axArr(k) = ax;
    im = imagesc(ax, panels(k).t_ms, panels(k).z, panels(k).X, [-vlim vlim]);
    set(im, 'AlphaData', ~isnan(panels(k).X));   % NaNs show the gray axis color
    set(ax, 'YDir', 'normal', 'Color', grayC);   % depth increases upward
    colormap(ax, cmap);
    ylim(ax, [zmin zmax]);
    xlim(ax, [0 o.win_ms]);
    set(ax, 'TickDir', 'out', 'XTick', []); box(ax, 'off');
    if k == 1
        ylabel(ax, 'depth from probe tip (\mum)   (deep = bottom)');
    else
        set(ax, 'YTickLabel', []);
    end
    pc = probeCol(panels(k).probe, :);
    title(ax, sprintf('sh%d', panels(k).shank), 'Color', pc, 'FontWeight', 'normal');
end

% shared colorbar
cb = colorbar(axArr(end)); cb.Layout.Tile = 'east';
cb.Label.String = 'voltage (\muV)';

% titles / caption
ppStr = ''; 

if exist('ff', 'var')
    ppStr = sprintf('%g Hz high-pass + per-shank CMR', ff.passHz); 
end

mainTitle = sprintf('Wall heatmap — %d probe(s), %d shank(s), %d channels', ...
    nProbe, nPan, sum(arrayfun(@(pnl) numel(pnl.channels), panels)));


capStr = sprintf('time (ms) — %g ms window at t = %g s', o.win_ms, o.t_start_s);
if ~isempty(ppStr), capStr = sprintf('%s,  %s', capStr, ppStr); end
xlabel(tl, capStr);

% per-probe header brackets spanning each probe's shank tiles, then main title
% above them (annotation rather than tl.Title so it clears the headers).
drawnow;   % force tiled layout to compute tile positions
drawProbeHeaders_(fig, axArr, [panels.probe], names, probeCol);
annotation(fig, 'textbox', [0.06 0.955 0.88 0.042], 'String', mainTitle, ...
    'HorizontalAlignment', 'center', 'VerticalAlignment', 'middle', ...
    'EdgeColor', 'none', 'FontWeight', 'bold', 'FontSize', 12, 'FitBoxToText', 'off');

% save
if o.SaveToggle
    figType = 'ap';
    figName = sprintf('wallHeatmap_%gms_t%gs.png', o.win_ms, o.t_start_s);
    qcFigSave(figType, figName);
end

% ---------------------------------------------------------------- report + output
fprintf('Wall heatmap: %d probe(s), %d panel(s), window %g s +%g ms, color +/-%.0f µV\n', ...
    nProbe, nPan, o.t_start_s, o.win_ms, vlim);
for k = 1:nPan
    fprintf('  %-10s sh%d : %3d ch (%d masked), depth %.0f–%.0f µm\n', ...
        panels(k).probeName, panels(k).shank, numel(panels(k).channels), ...
        panels(k).nMasked, min(panels(k).z), max(panels(k).z));
end

info = struct();
info.vlim = vlim;
info.t_start_s = o.t_start_s;
info.win_ms    = o.win_ms;
meta = rmfield(panels, {'X', 't_ms'});
for k = 1:nPan
    meta(k).zRange    = [min(panels(k).z) max(panels(k).z)];
    meta(k).nChannels = numel(panels(k).channels);
end
info.panels = meta;
end

% ======================================================================
function thr = zmaxFor_(zm, iProbe, iShank)
% Resolve the z_max threshold for one (probe, shank). zm may be a scalar, or a
% cell {1 x nProbe} whose entries are scalars or per-shank vectors.
if iscell(zm)
    zi = zm{iProbe};
else
    zi = zm;
end
if isscalar(zi)
    thr = zi;
else
    thr = zi(min(iShank, numel(zi)));
end
end

% ======================================================================
function drawProbeHeaders_(fig, axArr, probeIdx, names, probeCol)
% Draw a colored bracket line + probe label above each probe's group of shank
% tiles. The header band sits in the gap between the shank titles and the main
% title (fixed figure fractions), spanning the x-extent of the probe's tiles.
figPix = getpixelposition(fig);
ny = 0.895;                              % bracket line, figure fraction
probes = unique(probeIdx, 'stable');
for pp = probes
    grp = find(probeIdx == pp);
    L = getpixelposition(axArr(grp(1)),  true);   % [x y w h] in figure pixels
    R = getpixelposition(axArr(grp(end)), true);
    nx0 = L(1) / figPix(3);
    nx1 = (R(1) + R(3)) / figPix(3);
    annotation(fig, 'line', [nx0 nx1], [ny ny], ...
        'Color', probeCol(pp,:), 'LineWidth', 2.5);
    annotation(fig, 'textbox', [nx0, ny+0.004, nx1-nx0, 0.035], ...
        'String', names{pp}, 'Color', probeCol(pp,:), ...
        'HorizontalAlignment', 'center', 'VerticalAlignment', 'bottom', ...
        'EdgeColor', 'none', 'FontWeight', 'bold', 'FitBoxToText', 'off');
end
end

% ======================================================================
function C = probeColors_(n)
% Distinct probe colors (blue, red, green, ... matching the reference figure).
base = [0.122 0.365 0.749;   % blue
        0.839 0.153 0.157;   % red
        0.173 0.549 0.204;   % green
        0.580 0.404 0.741;   % purple
        0.890 0.467 0.114;   % orange
        0.090 0.639 0.678];  % teal
if n <= size(base,1)
    C = base(1:n, :);
else
    C = base(mod(0:n-1, size(base,1)) + 1, :);
end
end

% ======================================================================
function v = prctile_(x, pInPct)
% Linear-interpolated percentile matching MATLAB's prctile convention
% (positions at 100*(i-0.5)/n). NaNs are ignored. No Statistics Toolbox needed.
x = x(:); x = x(~isnan(x));
x = sort(x);
n = numel(x);
if n == 0, v = NaN; return; end
if n == 1, v = x(1); return; end
pos = pInPct/100 * n + 0.5;          % invert (i-0.5)/n
pos = min(max(pos, 1), n);
lo = floor(pos); hi = ceil(pos);
if lo == hi
    v = x(lo);
else
    v = x(lo) + (pos - lo) * (x(hi) - x(lo));
end
end

% ======================================================================
function cmap = blueWhiteRed_(n)
% Diverging blue-white-red colormap (RdBu reversed): blue = negative, white = 0,
% red = positive. Matches the RdBu_r used in the Python original.
if nargin < 1, n = 256; end
blue  = [0.129 0.400 0.674];
white = [1.000 1.000 1.000];
red   = [0.698 0.094 0.168];
m  = floor(n/2);
b2w = [linspace(blue(1),  white(1), m)', ...
       linspace(blue(2),  white(2), m)', ...
       linspace(blue(3),  white(3), m)'];
w2r = [linspace(white(1), red(1), n-m)', ...
       linspace(white(2), red(2), n-m)', ...
       linspace(white(3), red(3), n-m)'];
cmap = [b2w; w2r];
end
