function output = plotUnitDrift(ksResults, fs, z_max, varargin)
% plotUnitDrift  Pre-correction drift map: per-spike depth vs time.
%
% The "smoking gun" drift plot, made from the sorter's own per-spike
% localizations (spike_positions.npy, col 2 = depth in um) BEFORE any DREDGE
% motion correction. Two panels of the same scatter, colored two ways
% (SKILL golden rule #1 -- never certify drift from one coloring alone):
%
%   Panel A: colour = unit (golden-angle hue, ordered by median depth).
%            A continuous depth band that SWITCHES COLOUR over time is one
%            drifting neuron the sorter re-clustered = a drift-split you can see.
%   Panel B: colour = spike amplitude (log, viridis). Residual coherent motion
%            (bands sagging/rising together) shows up most clearly here.
%
% y-axis is physical depth in MICRONS (not channel index or depth-rank) so the
% geometry that drift lives in is undistorted.
%
% INPUTS
%   ksdir  (char/string) = Kilosort4 output dir (spike_times/clusters/positions,
%                          amplitudes, channel_positions .npy). If omitted,
%                          defaults to the probe00a full sort.
%   fs     (double)      = sampling rate (Hz), default 30000.
%   z_max  (double)      = depth (um) above which the band is flagged artifact
%                          (z>4850 um, from LFP 2026-07-15). Shaded for reference;
%                          this AP sort is not yet updated to exclude it.
%
% NAME/VALUE
%   'SaveToggle' (false) = write figures/unitDrift.png at 150 dpi.
%   'dispSumm'   (false) = print a one-line summary.
    
    p = inputParser;
    p.addParameter('downsample', true, @isnumeric);
    p.addParameter('SaveToggle', false);
    p.addParameter('dispSumm',   false);
    p.parse(varargin{:});
    o = p.Results;
        
    %% params
    
    depth = ksResults.sp(:, 2);
    nSpikes = numel(ksResults.st);
    
    %% per-unit colour rank (ordered by median depth, like the gist)
    
    [uc, ~, g] = unique(ksResults.sc);              % g in 1..nU (g is the cluster id for every spike)
    nClusters = numel(uc);
    medDepth = accumarray(g, depth, [], @median);
    [~, ord] = sort(medDepth);
    rankOfUnit = zeros(nClusters, 1);
    rankOfUnit(ord) = 1:nClusters;              % 1 = deepest-... (sorted ascending)
    
    % golden-angle hue so depth-adjacent units separate in colour
    GOLDEN = 0.6180339887;
    hue    = mod((rankOfUnit - 1) * GOLDEN, 1);
    unitRGB = hsv2rgb([hue, 0.80 * ones(nClusters, 1), 0.95 * ones(nClusters, 1)]);   % nU x 3
    
    %% downsample: per-unit random cap (keeps all spikes of sparse units)

    if o.downsample
        maxPoints = 4e5;

        counts   = accumarray(g, 1); % num spikes per unit
        plotFrac = floor((counts(g) ./ nSpikes) * maxPoints); % theoretical number of spikes to plot per cluster according to max points
        keepFrac = plotFrac ./ counts(g); 
        keepA    = rand(nSpikes, 1) < keepFrac;               

        % Panel B: subsample keepA further to a global cap
        idxA = find(keepA);
        if numel(idxA) > maxPoints
            idxB = idxA(randperm(numel(idxA), round(maxPoints)));
        else
            idxB = idxA;
        end
    end
    
    %% figure
    
    fig = figure('Color', 'w', 'Position', [60 60 1700 850]);
    tl  = tiledlayout(fig, 1, 2, 'TileSpacing', 'compact', 'Padding', 'compact');
    
    xl = [0 max(ksResults.st)];
    yl = [min(depth) - 20, max(depth) + 20];
    
    % Panel A -- colour by unit
    axA = nexttile(tl);
    drawBase(axA, xl, yl, z_max, 'colour = unit  (band switching colour = drift-split)');
    cA = unitRGB(g(idxA), :);
    scatter(axA, ksResults.st(idxA), depth(idxA), 2, cA, 'filled', ...
        'MarkerFaceAlpha', 0.55);
    
    % Panel B -- colour by amplitude (log viridis)
    axB = nexttile(tl);
    drawBase(axB, xl, yl, z_max, 'colour = spike amplitude (log)');
    scatter(axB, ksResults.st(idxB), depth(idxB), 2, ksResults.amp(idxB), 'filled', ...
        'MarkerFaceAlpha', 0.55);
    colormap(axB, viridisMap(256));
    set(axB, 'ColorScale', 'log');
    clim(axB, prctile(ksResults.amp(idxB), [1 99]));
    cb = colorbar(axB); cb.Label.String = 'spike amplitude (a.u.)';
    text(axB, xl(2)*0.01, z_max, 'artifact band (z>4850 \mum)', ...
        'FontSize', 8, 'Color', [.4 .4 .4], 'VerticalAlignment', 'bottom');
    
    title(tl, sprintf(['probe00a full KS4 — per-spike depth vs time, PRE-DREDGE  ' ...
        '(%d units, %s spikes)'], nClusters, addComma(nSpikes)), 'FontWeight', 'bold');
    
    %% save & summary
    
    if o.SaveToggle
        figType = 'ks';
        figName = 'unitDrift';
        qcFigSave(figType, figName);
    end
    
    if o.dispSumm
        fprintf(['units=%d  spikes=%d  dur=%.1f s  depth %.0f-%.0f um  ' ...
            'plotted A=%d B=%d  n(z>%g)=%d\n'], ...
            nClusters, nSpikes, max(ksResults.st), min(depth), max(depth), ...
            numel(idxA), numel(idxB), z_max, sum(depth > z_max));
    end
    
    output.nUnits   = nClusters;
    output.nSpikes  = nSpikes;
    output.depthLim = [min(depth) max(depth)];
    output.durSec   = max(ksResults.st);
end

%% helpers

function drawBase(ax, xl, yl, artZ, ttl)
    hold(ax, 'on');
    patch(ax, [xl(1) xl(2) xl(2) xl(1)], [artZ artZ yl(2) yl(2)], [.85 .85 .85], ...
        'EdgeColor', 'none', 'FaceAlpha', 0.6);
    yline(ax, artZ, '--', 'Color', [.5 .5 .5], 'LineWidth', 1);
    set(ax, 'XLim', xl, 'YLim', yl);
    grid(ax, 'on'); set(ax, 'GridAlpha', 0.25, 'Layer', 'top');
    xlabel(ax, 'time (s)');
    ylabel(ax, 'spike depth along shank (\mum from tip)');
    title(ax, ttl);
end

function cmap = viridisMap(nc)
    anchors = [ 68   1  84; 72  40 120; 62  74 137; 49 104 142; ...
                38 130 142; 31 158 137; 53 183 121;109 205  89; ...
               180 222  44;253 231  37] / 255;
    xi = linspace(0, 1, size(anchors, 1));
    xq = linspace(0, 1, nc);
    cmap = [interp1(xi, anchors(:,1), xq)', ...
            interp1(xi, anchors(:,2), xq)', ...
            interp1(xi, anchors(:,3), xq)'];
    cmap = min(max(cmap, 0), 1);
end

function s = addComma(v)
    s = regexprep(sprintf('%d', v), '\d(?=(\d{3})+$)', '$0,');
end