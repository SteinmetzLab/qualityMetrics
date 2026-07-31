function output = ampDepthScatter(ksResults, fs, z_max, varargin)
% qcAmpDepthScatter  Post-sort QC scatter: neuron amplitude vs probe depth.
%
% Panel A: colour = mean firing rate (Hz, log colour axis, viridis).
% Panel B: colour = waveform trough-to-peak duration (ms, red->green,
%          red = 0 ms short/narrow-spiking, green = 1 ms long/broad-spiking).
%
% INPUTS
%   ksResults (struct)    = kilosort outputs structure. See autoqc_working
%   t0  (double)    = window start in seconds
% 
% Amplitude / depth / firing rate are imported from the Python side (metrics.npz,
% via readNPZ) -- the validated ptp(templates @ whitening_mat_inv)*gain amplitude,
% depth = channel_positions(peakch,2). Trough-to-peak duration is computed here
% from the un-whitened template on each unit's peak channel (templates.npy).
%
% Note: the z>4850 um band was flagged as artifact from LFP (2026-07-15); this KS4
% AP sort is not yet updated to exclude it, so that band is shaded for reference.

p = inputParser;
p.addParameter('SaveToggle', false);
p.addParameter('dispSumm', false);
p.parse(varargin{:});
o = p.Results;

%% load directories and files

amp    = double(ksResults.metrics.amp_uV(:));
fr     = double(ksResults.metrics.fr(:));
depth  = double(ksResults.metrics.peak_depth(:));
peakch = double(ksResults.metrics.peakch(:));          % 0-based channel index (from Python)
n      = numel(amp);

%% calculate trough-to-peak duration
% --- trough-to-peak duration (ms) from un-whitened template on peak channel ---
dur_ms = zeros(n, 1);
for u = 1:n
    T  = double(squeeze(ksResults.templates(u, :, :)));   % 61 x 383 (ADC units, whitened)
    wf = T * ksResults.winv(:, peakch(u) + 1);            % 61 x 1 un-whitened peak-ch waveform
    [~, ti] = min(wf);                          % trough (global min)
    [~, rel] = max(wf(ti:end));                 % first positive peak after trough
    dur_ms(u) = (rel - 1) / fs * 1000;
end

% columns: cluster_id, KSLabel
isGoodClu   = strcmp(ksResults.lbl.KSLabel, 'good');           % logical per cluster row

%% plot depth/firing rate

% open figure

fig = figure('Color', 'w', 'Position', [80 80 1500 800]);
tl  = tiledlayout(fig, 1, 2, 'TileSpacing', 'compact', 'Padding', 'compact');

xl = [0.9*min(amp) 1.1*max(amp)];
yl = [min(depth)-100 max(depth)+100];

% Panel A -- firing rate (log colour axis, viridis)
axA = nexttile(tl);
drawBase(axA, xl, yl, z_max, 'colored by mean firing rate');

scatter(axA, amp(isGoodClu), depth(isGoodClu), 36, fr(isGoodClu), ...
    'filled', 'MarkerEdgeColor', 'k', 'SizeData', 50, 'LineWidth', 0.3);
scatter(axA, amp(~isGoodClu), depth(~isGoodClu), 36, fr(~isGoodClu), ...
    'filled', 'MarkerEdgeColor', 'k', 'SizeData', 10, 'LineWidth', 0.3);

colormap(axA, viridisMap(256));
set(axA, 'ColorScale', 'log');
clim(axA, [max(min(fr), 1e-3) max(fr)]);
cA = colorbar(axA); cA.Label.String = 'firing rate (Hz)';
text(axA, xl(1)*1.03, z_max+60, 'artifact band (z>4850 \mum)', ...
    'FontSize', 8, 'Color', [.4 .4 .4], 'VerticalAlignment', 'bottom');

%% trough to peak

% Panel B -- trough-to-peak duration (red 0 ms -> green 1 ms)
axB = nexttile(tl);
drawBase(axB, xl, yl, z_max, 'colored by trough-to-peak duration');
scatter(axB, amp(isGoodClu), depth(isGoodClu), 36, dur_ms(isGoodClu), ...
    'filled', 'MarkerEdgeColor', 'k', 'SizeData', 50, 'LineWidth', 0.3);
scatter(axB, amp(~isGoodClu), depth(~isGoodClu), 36, dur_ms(~isGoodClu), ...
    'filled', 'MarkerEdgeColor', 'k', 'SizeData', 10, 'LineWidth', 0.3);
colormap(axB, rdYlGnMap(256));
clim(axB, [0 1]);
cB = colorbar(axB); cB.Label.String = 'trough-to-peak duration (ms)';

title(tl, sprintf('probe00a full KS4 — amplitude vs depth  (%d units)', n), ...
    'FontWeight', 'bold');

%% save and summary

    if o.SaveToggle
        figType = 'ks';
        figName = 'ampDepthScatter';
        qcFigSave(figType, figName);
    end

    % summary statement
    
    if o.dispSumm
        fprintf(['units=%d  amp %.0f-%.0f uV  depth %.0f-%.0f um\n' ...
            'fr %.3f-%.1f Hz  dur %.2f-%.2f ms (median %.2f)  n(z>4850)=%d\n'], ...
            n, min(amp), max(amp), min(depth), max(depth), ...
            min(fr), max(fr), min(dur_ms), max(dur_ms), median(dur_ms), sum(depth > z_max));
    end
end

%% helper functions

function drawBase(ax, xl, yl, artZ, ttl)
    hold(ax, 'on');
    % shaded artifact band (drawn first so points sit on top)
    patch(ax, [xl(1) xl(2) xl(2) xl(1)], [artZ artZ yl(2) yl(2)], [.85 .85 .85], ...
        'EdgeColor', 'none', 'FaceAlpha', 0.6);
    yline(ax, artZ, '--', 'Color', [.5 .5 .5], 'LineWidth', 1);
    % set(ax, 'XScale', 'log', 'XLim', xl, 'YLim', yl);
    set(ax, 'XLim', xl, 'YLim', yl);
    grid(ax, 'on'); set(ax, 'GridAlpha', 0.25, 'Layer', 'top');
    xlabel(ax, 'spike amplitude (\muV)');
    ylabel(ax, 'probe depth (\mum)');
    title(ax, ttl);
end

function cmap = viridisMap(nc)
    % viridis anchor stops (matplotlib), interpolated to nc colours.
    anchors = [ 68   1  84; 72  40 120; 62  74 137; 49 104 142; ...
                38 130 142; 31 158 137; 53 183 121;109 205  89; ...
               180 222  44;253 231  37] / 255;
    cmap = interpMap(anchors, nc);
end

function cmap = rdYlGnMap(nc)
    % RdYlGn (matplotlib) red -> yellow -> green, interpolated to nc colours.
    anchors = [165   0  38;215  48  39;244 109  67;253 174  97; ...
               254 224 139;255 255 191;217 239 139;166 217 106; ...
               102 189  99; 26 152  80;  0 104  55] / 255;
    cmap = interpMap(anchors, nc);
end

function cmap = interpMap(anchors, nc)
    xi = linspace(0, 1, size(anchors, 1));
    xq = linspace(0, 1, nc);
    cmap = [interp1(xi, anchors(:,1), xq)', ...
            interp1(xi, anchors(:,2), xq)', ...
            interp1(xi, anchors(:,3), xq)'];
    cmap = min(max(cmap, 0), 1);
end
