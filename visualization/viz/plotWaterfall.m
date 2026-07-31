function f = plotWaterfall(td, Y, t0, dur, channels, ff, varargin)
%PLOTWATERFALL  Stacked multi-channel raw-voltage "waterfall" for one window.
%
%   f = plotWaterfall(td, Y, t0, dur, channels, ff, Name, Value, ...)
%
% Draws the classic Neuropixels raw view: each channel as a thin voltage
% trace, vertically offset and coloured by depth (viridis), so spikes appear
% as sharp deflections across neighbouring channels.
%
% INPUTS
%   td          (double)    = sample times within the window in seconds
%   Y           (double)    = [nChannels x nSamples] voltage traces in uV
%   t0          (double)    = window start in seconds
%   dur         (double)    = window duration in seconds, e.g. 0.01, 0.1, 1, 10
%   channels    (double)    = 1-based channel indices (contiguous span)
%   ff          (struct)    = filter flags from preprocessTraces (title text)
%
% NAME-VALUE OPTIONS
%   'Figure'    (graphics)  = existing figure/axes handle to draw into (default new)
%
% OUTPUT
%   f   (graphics) = figure handle

    p = inputParser;
    p.addParameter('Figure', []);
    p.parse(varargin{:});
    o = p.Results;
    
    %% params

    nch = numel(channels);

    % plot params

    spacingUV = [];
    scaleBarUV = 200;
    titleTag = [];
    maxPointsPerTrace = 50000;
    lineWidth = 0.4;
    
    %%
    
    % ---- channel spacing & colours -----------------------------------------
    rstd = robustStd(Y);                          % [nch x 1]
    if isempty(spacingUV)
        spacing = 6 * median(rstd(rstd > 0));     % auto: 6x typical channel noise
        if ~isfinite(spacing) || spacing <= 0, spacing = 4*scaleBarUV; end
    else
        spacing = spacingUV;
    end
    offsets = (0:nch-1)' * spacing;               % channel 1 at bottom
    cmap    = viridisColors(nch);                 % row 1 purple (bottom) -> yellow (top)
    
    % time axes
    if dur < 1
        tplot = td * 1000; tunit = 'Time (ms)'; xmax = dur*1000;
    else
        tplot = td;        tunit = 'Time (s)';  xmax = dur;
    end
    
    % figure
    if isempty(o.Figure)
        f = figure('Color', 'w', 'Position', [100 100 900 1100]);
        ax = axes(f);
    elseif isgraphics(o.Figure, 'axes')
        ax = o.Figure; f = ancestor(ax, 'figure');
    else
        f = o.Figure; ax = axes(f);
    end
    hold(ax, 'on');
    
    for i = 1:nch
        [ti, xi] = minmaxEnvelope(tplot, Y(i,:) + offsets(i), maxPointsPerTrace);
        plot(ax, ti, xi, 'Color', cmap(i,:), 'LineWidth', lineWidth);
    end
    
    % scale bar
    xr = xmax * 0.012;                            % small inset from left edge
    yb = -spacing;                                % just below the bottom channel
    plot(ax, [xr xr], [yb yb+scaleBarUV], 'k', 'LineWidth', 2);
    text(ax, xr, yb + scaleBarUV/2, sprintf(' %g µV', scaleBarUV), ...
        'Rotation', 90, 'HorizontalAlignment', 'center', ...
        'VerticalAlignment', 'bottom', 'FontSize', 9);
    
    % cosmetics
    xlim(ax, [0 xmax]);
    ylim(ax, [yb - spacing, offsets(end) + spacing*1.5]);
    set(ax, 'YTick', []);
    xlabel(ax, tunit);
    ylabel(ax, sprintf('Channels by depth (%d contiguous)', nch));
    box(ax, 'off');
    set(ax, 'TickDir', 'out', 'FontSize', 10);
    
    title(ax, buildTitle(titleTag, channels, t0, dur, ff), ...
        'Interpreter', 'none');
    hold(ax, 'off');
end

% ======================================================================
function s = robustStd(Y)
    % Per-row robust std via MAD.
    med = median(Y, 2);
    s = 1.4826 * median(abs(Y - med), 2);
end

% ----------------------------------------------------------------------
function [ti, xi] = minmaxEnvelope(t, x, maxPoints)
    % Downsample a single trace to <= maxPoints while preserving spike extremes:
    % split into bins and keep the min and max sample of each (at their true times).
    n = numel(x);
    if isempty(maxPoints) || maxPoints <= 0 || n <= maxPoints
        ti = t; xi = x; return;
    end
    nbins = max(1, floor(maxPoints/2));
    edges = round(linspace(1, n+1, nbins+1));
    ti = zeros(1, nbins*2); xi = ti; k = 0;
    for b = 1:nbins
        i0 = edges(b); i1 = edges(b+1) - 1;
        if i1 < i0, continue; end
        seg = x(i0:i1);
        [~, imn] = min(seg); [~, imx] = max(seg);
        ii = sort([i0+imn-1, i0+imx-1]);
        ti(k+1:k+2) = t(ii); xi(k+1:k+2) = x(ii); k = k + 2;
    end
    ti = ti(1:k); xi = xi(1:k);
end

% ----------------------------------------------------------------------
function str = buildTitle(tag, channels, t0, dur, ff)
    if dur < 1
        durStr = sprintf('%g ms', dur*1000); 
    else
        durStr = sprintf('%g s', dur); 
    end

    startStr = sprintf('t0 = %g s', t0);

    parts = {};

    if ~isempty(ff.passHz)
        if size(ff.passHz, 2) == 1
            parts{end+1} = sprintf('%g Hz %s', ff.passHz, ff.bandpassType);
        elseif size(ff.passHz, 2) == 2
            parts{end+1} = sprintf('%g\x2013%d Hz %s', ff.passHz(1), ff.passHz(2), ff.bandpassType);
        end
    end
    
    if ~isempty(ff.CMR)
        parts{end+1} = 'CMR'; 
    end
    
    if isempty(parts)
        procStr = 'raw'; 
    else
        procStr = strjoin(parts, ' + '); 
    end
    
    chStr = sprintf('ch %d\x2013%d', channels(1), channels(end));
    head = tag; 
    
    if isempty(head)
        head = chStr; 
    else
        head = sprintf('%s  (%s)', tag, chStr); 
    end

    str = sprintf('%s \x2014 %s, %s, %s', head, startStr, durStr, procStr);
end
