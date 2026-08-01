function info = plotDepthPower(rec, band, varargin)
% PLOTDEPTHPOWER  Depth-resolved LFP power map: channels (depth) on the
% y-axis, frequency on the x-axis, time-averaged power as the color map.
%
% Usage
%   plotLFPSpectrogram                                  % probe00a defaults
%   plotLFPSpectrogram('t_start_s',1000,'win_s',10)
%   plotLFPSpectrogram('Normalize','freq')             % divide out 1/f per band
%
% Name/value options (all optional)
%   'nFreq'     number of log-spaced frequencies              (default 60)
%   'nCyc'      Morlet cycles: scalar, or [lo hi] log-ramped low->high f
%               (more cycles = sharper in frequency, blurrier in time)
%               (default [3 10])
%   'z_max'     mask channels with depth > z_max (top artifact band); Inf off
%               (default 4850, the probe00a convention)
%   'CAR'       subtract across-channel mean before analysis   (default false;
%               probe00a lf_pre is already CMR'd)
%   'Normalize' 'none' (absolute dB) | 'freq' (dB re. median-over-depth at
%               each frequency; cancels the 1/f slope to reveal depth structure)
%               (default 'none')
%   'climPct'   robust color saturation percentiles [lo hi]   (default [2 98])
%   'clim'      explicit color limits [lo hi] in dB; overrides climPct (default [])
%   'SaveToggle'  save PNG to figures\lfpSpectrogram          (default false)
%
% OUTPUT (struct info)
%   .f          [1 x nFreq]  frequencies (Hz)
%   .depth      [nKept x 1]  channel depths shown (µm, ascending)
%   .powerDB    [nKept x nFreq] plotted values (dB)
%   .t_start_s, .win_s, .fs, .nCyc

    %% input parser
    
    p = inputParser;
    p.addRequired('rec',        @(x) isstruct(x));
    p.addRequired('band',       @ischar)
    p.addParameter('nFreq',     60,   @isnumeric);
    p.addParameter('nCyc',      [3 10], @isnumeric);
    p.addParameter('z_max',     4850, @isnumeric);
    p.addParameter('Normalize', 'none', @(s) any(strcmpi(s,{'none','freq'})));
    p.addParameter('climPct',   [2 98], @isnumeric);
    p.addParameter('clim',      [],   @isnumeric);
    p.addParameter('SaveToggle',false, @islogical);
    p.parse(rec, band, varargin{:});
    o = p.Results;
    
    %% parameters
    
    numWindows = 10;

    % create bins for averaging
    % determine bin starts > + window to get bin ends
    % windowStarts = :

    switch band
        case 'ap'
            fMin    = 300;
            fMax    = 5000;
            win_s   = 2;             % window size, seconds

            % graph parameters

            depthPlotXTicks = [300 500 1000 2000 5000];
        case 'lfp'
            fMin    = 1;               % lowest frequency, Hz
            fMax    = 300;             % highest frequency, Hz
            win_s   = 10;             % window size, seconds

            % graph parameters

            depthPlotXTicks = [1 2 4 8 16 32 64 128 256];
    end

    windowEdges = linspace(300, rec.nSamples / rec.fs, 11); % 10 window start times
    
    %% geometry / depth
    
    geom = recordingGeometry(fileparts(rec.filePath));       % .y = depth (µm), aligned to channels

    z = geom.y(1:rec.nChannels);
    keep = z <= o.z_max;                            % drop top artifact band
    [zsorted, ord] = sort(z(keep));                 % ascending depth
    idxKeep = find(keep);
    order   = idxKeep(ord);                         % display-row -> channel index
    nKept = numel(order);

    %% load recording
    
    fs  = rec.fs;                                   % 2500 Hz for probe00a lf_pre
    if fMax > fs/2
        warning('fMax (%g Hz) exceeds Nyquist (%g Hz); clamping.', fMax, fs/2);
        fMax = floor(fs/2) - 1;
    end
    
    %% initialize data structure

    powerDBHold = zeros(nKept, o.nFreq, numWindows);

    %% convolution

    fprintf(['------ \n' ...
        'Beginning Morlet convolution. Progress: '])
    for i = 1:numWindows
        fprintf('%d%% ', i*10)
        t0 = windowEdges(i);
        t1 = t0 + win_s;
        X  = getTraces(rec, [t0 t1], [], 'uV');             % [nChan x nSamp] µV, double

        [nChan, nSamp] = size(X);
        if nSamp < 2
            error('plotLFPSpectrogram:emptyWindow', 'Requested window is empty.');
        end

        if ~rec.preprocessed
            X = X - mean(X, 1);
        end
    
        X = X - mean(X, 2);                                 % remove per-channel DC
        
        % more geometry
        if geom.nChannels ~= nChan
            warning('Geometry has %d channels but data has %d; using min.', geom.nChannels, nChan);
        end

        Xk = X(order, :);

        % frequency bank

        f = logspace(log10(fMin), log10(fMax), o.nFreq);   % [1 x nFreq]
        if isscalar(o.nCyc)
            cyc = repmat(o.nCyc, 1, o.nFreq);
        else
            % log-ramp the cycle count from o.nCyc(1) at fMin to o.nCyc(2) at fMax
            cyc = logspace(log10(o.nCyc(1)), log10(o.nCyc(2)), o.nFreq);
        end
        
        % Morlet convolution
        % For each frequency: build a complex Morlet wavelet, FFT-convolve it against
        % every channel at once ('same' output), take |.|^2, average over time (with
        % the edge-effect region trimmed), giving one power value per (channel, freq).
        powerDB = zeros(nKept, o.nFreq);
        for fi = 1:o.nFreq
            s   = cyc(fi) / (2*pi*f(fi));                % Gaussian SD (s), in seconds
            tw  = -4*s : 1/fs : 4*s;                     % wavelet support: +/- 4 SD
            if numel(tw) < 3, tw = (-1:1)/fs; end
            w   = exp(2i*pi*f(fi)*tw) .* exp(-tw.^2 / (2*s^2));
            w   = w / sqrt(sum(abs(w).^2));              % unit-energy normalization
            nK  = numel(w);
        
            nConv = nSamp + nK - 1;
            nFFT  = 2^nextpow2(nConv);
            Xf    = fft(Xk, nFFT, 2);                    % [nKept x nFFT]
            Wf    = fft(w,  nFFT);                        % [1 x nFFT] (broadcasts)
            C     = ifft(Xf .* Wf, [], 2);               % full convolution
            half  = floor(nK/2);
            C     = C(:, half + (1:nSamp));              % trim to 'same' length
        
            pw = abs(C).^2;                              % instantaneous power
            edge = min(nK, floor(nSamp/2)-1);            % drop wavelet-tainted edges
            edge = max(edge, 0);
            if 2*edge < nSamp
                pw = pw(:, edge+1 : nSamp-edge);
            end
            powerDB(:, fi) = 10*log10(mean(pw, 2) + eps); % time-averaged power -> dB
        end
        
        % normalization
    
        switch lower(o.Normalize)
            case 'freq'
                powerDB = powerDB - median(powerDB, 1);  % dB re. across-depth median per freq
                cbLabel = append(band, ' power (dB re. median depth)');
            otherwise
                cbLabel = append(band, ' power (dB)');
        end
        
        powerDBHold(:, :, i) = powerDB;
    end

    % average

    avPower = mean(powerDBHold, 3);
    
    % color limits
    
    if ~isempty(o.clim)
        cl = o.clim(:)';
    else
        cl = prctile_(avPower(:), o.climPct);
        if diff(cl) <= 0, cl = [min(avPower(:)) max(avPower(:))]; end
    end
    
    % ---------------------------------------------------------------- plot
    figure('Color','w','Position',[100 100 900 760]);
    imagesc(f, zsorted, avPower, cl);
    set(gca, 'YDir', 'normal');                      % depth increases upward
    set(gca, 'XScale', 'log');                       % log frequency axis for LFP
    colormap(viridisColors(256));
    xlim([fMin fMax]);
    xticks(depthPlotXTicks);
    xlabel('frequency (Hz)');
    ylabel('depth along probe (\mum)');
    title(sprintf('probe00a %s power by depth (Morlet)  | %g s window (%d windows averaged) |  %s', ...
          band, win_s, numWindows, o.Normalize));
    cb = colorbar; cb.Label.String = cbLabel;
    box on;
    
    fprintf(['%s depth-spectrum: %d/%d channels (z<=%g), %d freqs %.3g-%.3g Hz, ' ...
             'cyc %.3g-%.3g, %g s window @ %g Hz\n'], ...
        band, nKept, nChan, o.z_max, o.nFreq, fMin, fMax, cyc(1), cyc(end), win_s, fs);
    
    %% output / save

    info = struct('f', f, 'depth', zsorted, 'powerDB', avPower, ...
                  't_start_s', t0, 'win_s', win_s, 'fs', fs, 'nCyc', cyc);

    if o.SaveToggle
        figType = band;
        figName = sprintf('%sSpectrogram_t%g_w%g.png', band, t0, win_s);
        qcFigSave(figType, figName);
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
        r = pcts(i)/100 * (n-1) + 1;      % 1-based fractional rank
        lo = floor(r); hi = ceil(r);
        if lo < 1, lo = 1; end
        if hi > n, hi = n; end
        q(i) = x(lo) + (r-lo)*(x(hi)-x(lo));
    end
end
