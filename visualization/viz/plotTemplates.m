function out = plotTemplates(varargin)
% plotTemplates  Visualise a Kilosort4 template from templates.npy.
%
%   out = plotTemplates()                 % plots the biggest-amplitude template
%   out = plotTemplates('clu', 47)        % plots cluster/template id 47 (0-based)
%   out = plotTemplates('clu', 47, 'unwhiten', true)
%
% WHAT templates.npy IS
% ---------------------
% Kilosort4 writes one *template* per cluster: the canonical spike waveform it
% fit to every spike in that cluster. The file is a single 3-D array:
%
%       templates.npy  =  [nTemplates x nTimePoints x nChannels]
%                      =  [   525      x     61       x   383   ]   (this sort)
%
%   dim 1 (nTemplates) : one row per cluster. Row (tid+1) in MATLAB is the
%                        template for 0-based cluster id `tid`. For raw KS4
%                        output this index also equals the value stored in
%                        spike_templates.npy for that cluster's spikes.
%   dim 2 (nTimePoints): 61 samples ≈ 2.03 ms at 30 kHz. The waveform is
%                        centred: the spike trough sits ~1/3 of the way in.
%   dim 3 (nChannels)  : ALL 383 AP channels. Templates are stored DENSE here
%                        (not sparse), but a real neuron only shows up on a
%                        handful of nearby channels — every other channel is
%                        ~flat noise. That spatial spread across a few channels
%                        is the "footprint" that lets KS localise the unit.
%
% So a single template W = squeeze(templates(tid+1,:,:)) is a [61 x 383] matrix:
% a little voltage snippet on each channel. Plot those snippets at each
% channel's physical (x, depth) position and you see the footprint travel down
% the probe. That is panel A below.
%
% UNITS / WHITENING
% -----------------
% KS4 templates live in *whitened* space (arbitrary units), which is why the
% y-axis is unitless by default. Pass 'unwhiten', true to right-multiply by
% whitening_mat_inv.npy (channel x channel) to get back to data space — the
% shape is identical, only the per-channel scaling changes. Handy for
% comparing amplitudes to the raw traces; unnecessary just to see the shape.
%
% INPUT (name/value)
%   'ksdir'    KS4 output folder            (default C:\spikesort_scratch\probe00a_full_ks4)
%   'clu'      0-based cluster id to plot   (default [] -> largest-amplitude template)
%   'unwhiten' apply whitening_mat_inv      (default false)
%   'nDraw'    max # channels in footprint  (default 24, the nearest to the peak)
%   'z_max'    hide channels with depth > z_max um (default Inf; 4850 = artifact band)
%   'save'     save PNG into autoqc/figures (default false)
%
% OUTPUT (struct out)
%   .clu, .peakChan, .peakDepth, .ptp (per-channel), .W ([61 x 383]), .t_ms
%
% See also readNPY, viridisColors, plotWaveforms.

p = inputParser;
p.addParameter('ksdir', 'C:\spikesort_scratch\probe00a_full_ks4', @(s)ischar(s)||isstring(s));
p.addParameter('clu',      [], @(x)isempty(x)||isscalar(x));
p.addParameter('unwhiten', false, @islogical);
p.addParameter('nDraw',    24, @isnumeric);
p.addParameter('z_max',    Inf, @isnumeric);
p.addParameter('save',     false, @islogical);
p.parse(varargin{:});
o = p.Results;
o.ksdir = char(o.ksdir);

fs = 30000;   % Hz, from params.py

% ---------------------------------------------------------------- load arrays
templates = readNPY(fullfile(o.ksdir, 'templates.npy'));        % [nClu x nT x nChan]
chanpos   = readNPY(fullfile(o.ksdir, 'channel_positions.npy')); % [nChan x 2] = [x z]
[nClu, nT, nChan] = size(templates);
x = double(chanpos(:,1));            % lateral column position (um), here {0,32}
z = double(chanpos(:,2));            % depth (um)

% optional un-whitening back to data space (per-channel rescale, shape unchanged)
if o.unwhiten
    winv = readNPY(fullfile(o.ksdir, 'whitening_mat_inv.npy'));  % [nChan x nChan]
    % apply to every template: [nT x nChan] * [nChan x nChan]
    for k = 1:nClu
        templates(k,:,:) = squeeze(templates(k,:,:)) * winv;
    end
end

% per-template, per-channel peak-to-peak amplitude -> pick a template to show
ptpAll = squeeze(max(templates,[],2) - min(templates,[],2));    % [nClu x nChan]
if isempty(o.clu)
    [~, row] = max(max(ptpAll,[],2));    % template with the largest footprint
    o.clu = row - 1;                     % report as 0-based id
else
    row = o.clu + 1;                     % 0-based id -> MATLAB row
    assert(row>=1 && row<=nClu, 'clu %d out of range 0..%d', o.clu, nClu-1);
end

W   = squeeze(templates(row,:,:));       % [nT x nChan] : the template
ptp = ptpAll(row,:).';                   % [nChan x 1]  : amplitude per channel
[peakPtp, peakChan] = max(ptp);          % channel carrying the spike
t_ms = ((0:nT-1) - (nT-1)/3) / fs * 1000; % ms, roughly centred on the trough

% channels to draw in the footprint: nearest `nDraw` to the peak channel,
% within the depth mask, that actually carry signal.
alive = ptp > 0.02*peakPtp & z <= o.z_max;
dist  = abs(z - z(peakChan)) + abs(x - x(peakChan));  % city-block, keeps columns
dist(~alive) = Inf;
[~, near] = sort(dist);
drawCh = near(1:min(o.nDraw, sum(isfinite(dist))));

% =====================================================================  FIGURE
fig = figure('Color','w','Position',[80 80 1180 720], ...
    'Name', sprintf('template %d', o.clu));
tl = tiledlayout(fig, 2, 3, 'TileSpacing','compact', 'Padding','compact');
title(tl, sprintf('KS4 template  |  cluster %d  |  peak ch %d @ %.0f um  |  ptp %.2f%s', ...
    o.clu, peakChan-1, z(peakChan), peakPtp, ternary(o.unwhiten,' (unwhitened)','')), ...
    'FontWeight','bold');

% ---- Panel A: spatial footprint on the probe (spans left 2/3, both rows) ----
axA = nexttile(tl, 1, [2 2]); hold(axA,'on');
% scale each snippet so it reads at its channel's (x,z) location
xpitch = min(diff(unique(x(x>=0)))); if isempty(xpitch)||xpitch==0, xpitch = 32; end
zs = sort(unique(z)); dz = diff(zs); zpitch = median(dz(dz>0)); if isempty(zpitch), zpitch = 15; end
xScale = 0.80*xpitch / (t_ms(end)-t_ms(1));     % ms  -> um (fit inside a column)
aScale = 6.0*zpitch / peakPtp;                  % amp -> um (peak spans ~6 rows)
for ci = drawCh(:).'
    isPk = ci==peakChan;
    plot(axA, x(ci) + (t_ms-t_ms(1))*xScale - 0.40*xpitch, ...
              z(ci) + W(:,ci)*aScale, ...
        'Color', ternary(isPk,[0.85 0.10 0.10],[0.30 0.30 0.35]), ...
        'LineWidth', ternary(isPk,1.8,0.9));
    plot(axA, x(ci)-0.40*xpitch, z(ci), '.', 'Color',[0.7 0.7 0.7], 'MarkerSize',4);
end
axis(axA,'tight'); grid(axA,'on'); box(axA,'on');
xlabel(axA,'probe x (um) + time'); ylabel(axA,'depth (um)');
title(axA,'A  spatial footprint (each trace = one channel at its site)');
% scale bars
xr = xlim(axA); yr = ylim(axA);
plot(axA, xr(1)+[2 2], yr(1)+[2 2]+[0 aScale*peakPtp/2], 'k', 'LineWidth',2);
text(axA, xr(1)+4, yr(1)+2+aScale*peakPtp/4, sprintf('%.1f a.u.',peakPtp/2), 'FontSize',8);

% ---- Panel B: depth x time heatmap of the same template --------------------
axB = nexttile(tl, 3);
[~, zord] = sort(z(drawCh));           % channels sorted by depth
Wd = W(:, drawCh(zord)).';             % [nDraw x nT]
imagesc(axB, t_ms, 1:numel(zord), Wd);
set(axB,'YDir','normal'); colormap(axB, divergingMap(256));
cl = max(abs(Wd(:))); if cl==0, cl=1; end; caxis(axB, [-cl cl]);
yticks(axB, 1:numel(zord)); yticklabels(axB, compose('%.0f', z(drawCh(zord))));
xlabel(axB,'time (ms)'); ylabel(axB,'depth (um)');
title(axB,'B  amplitude heatmap'); colorbar(axB);

% ---- Panel C: peak-channel waveform ----------------------------------------
axC = nexttile(tl, 6); hold(axC,'on');
plot(axC, t_ms, W(:,peakChan), 'Color',[0.85 0.10 0.10], 'LineWidth',1.8);
yline(axC, 0, ':', 'Color',[0.6 0.6 0.6]);
[~,it] = min(W(:,peakChan));
plot(axC, t_ms(it), W(it,peakChan), 'kv', 'MarkerFaceColor','k','MarkerSize',5);
axis(axC,'tight'); grid(axC,'on'); box(axC,'on');
xlabel(axC,'time (ms)'); ylabel('amplitude');
title(axC, sprintf('C  peak channel %d', peakChan-1));

% ---------------------------------------------------------------- save / return
if o.save
    figdir = fullfile(fileparts(fileparts(mfilename('fullpath'))), 'figures');
    if ~isfolder(figdir), mkdir(figdir); end
    fn = fullfile(figdir, sprintf('template_clu%03d.png', o.clu));
    exportgraphics(fig, fn, 'Resolution',150);
    fprintf('saved %s\n', fn);
end

out = struct('clu',o.clu, 'peakChan',peakChan-1, 'peakDepth',z(peakChan), ...
             'ptp',ptp, 'W',W, 't_ms',t_ms, 'nTemplates',nClu);
end

% ---- tiny helpers ----------------------------------------------------------
function v = ternary(c,a,b), if c, v=a; else, v=b; end, end

function cmap = divergingMap(n)
% blue -> white -> red, so + and - deflections are visible (sign matters here).
if nargin<1, n=256; end
h = linspace(0,1,n).';
lo = [0.13 0.28 0.75]; mid = [1 1 1]; hi = [0.75 0.12 0.12];
cmap = zeros(n,3);
for k=1:3
    cmap(:,k) = interp1([0 .5 1], [lo(k) mid(k) hi(k)], h, 'linear');
end
end
