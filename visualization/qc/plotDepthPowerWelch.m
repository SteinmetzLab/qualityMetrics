function info = plotDepthPowerWelch(varargin)
% PLOTLFPWELCH  Welch power-spectral-density of the LFP band, as a depth map or
% a classic per-channel PSD plot. Companion to plotLFPSpectrogram (Morlet).
%
% Welch's method: split each channel into overlapping Hann-windowed segments,
% periodogram each, average -> a low-variance PSD estimate (µV^2/Hz). Cheaper
% and smoother than the Morlet convolution; the two should broadly agree, so
% comparing them is itself a QC check.
%
% Two views (option 'Style'):
%   'depth' (default)  channels (depth) on y, frequency on x, PSD (dB) as color
%                      -- same layout as plotLFPSpectrogram for direct A/B.
%   'lines'            classic PSD: frequency (log x) vs PSD (dB/Hz, y), one
%                      thin trace per channel colored by depth, plus the mean.
%
% Usage
%   plotLFPWelch                                   % probe00a defaults, depth map
%   plotLFPWelch('Style','lines')                 % classic spaghetti PSD
%   plotLFPWelch('SegLen_s',4,'win_s',40)         % finer low-freq resolution
%
% Name/value options (all optional)
%   'lfpDir'    LFP recording folder (default C:\spikesort_scratch\probe00a_lf_pre)
%   't_start_s' window start, seconds                         (default 1000)
%   'win_s'     total analysis window, seconds                (default 20)
%   'SegLen_s'  Welch segment length, seconds (sets freq resolution ~1/SegLen)
%                                                             (default 2 -> 0.5 Hz)
%   'Overlap'   fractional segment overlap [0,1)              (default 0.5)
%   'fMin'      lowest  frequency shown, Hz                   (default 1)
%   'fMax'      highest frequency shown, Hz (<= LFP LP ~300)  (default 300)
%   'nFreq'     log-spaced freq bins for the depth map        (default 60)
%   'Style'     'depth' | 'lines'                             (default 'depth')
%   'z_max'     mask channels with depth > z_max; Inf off     (default 4850)
%   'CAR'       subtract across-channel mean first            (default false)
%   'Normalize' 'none' | 'freq' (dB re. median-over-depth per freq; depth map only)
%                                                             (default 'none')
%   'climPct'   robust color saturation percentiles [lo hi]   (default [2 98])
%   'clim'      explicit color limits [lo hi] dB (depth map)  (default [])
%   'SaveToggle'  save PNG to figures\lfpWelch                (default false)
%
% OUTPUT (struct info)
%   .f          frequencies (Hz): log grid ('depth') or native Welch bins ('lines')
%   .depth      [nKept x 1] channel depths shown (µm, ascending)
%   .psdDB      [nKept x nFreq] PSD in dB (dB/Hz)
%   .fNative, .psdNativeDB   native (linear-bin) Welch PSD before log-regrid
%   .t_start_s, .win_s, .fs, .segLen_s
%
% No toolboxes required (own Hann + Welch; dependency-free helpers at bottom).
% Uses autoqc infra: loadRawRecording, getTraces, recordingGeometry, viridisColors.
%
% See also plotLFPSpectrogram, wallHeatmap, preprocessLFP.

% ---------------------------------------------------------------- parameters
p = inputParser;
p.addParameter('lfpDir', 'C:\spikesort_scratch\probe00a_lf_pre', @(s)ischar(s)||isstring(s));
p.addParameter('t_start_s', 1000, @isnumeric);
p.addParameter('win_s',     20,   @isnumeric);
p.addParameter('SegLen_s',  2,    @isnumeric);
p.addParameter('Overlap',   0.5,  @isnumeric);
p.addParameter('fMin',      1,    @isnumeric);
p.addParameter('fMax',      300,  @isnumeric);
p.addParameter('nFreq',     60,   @isnumeric);
p.addParameter('Style',     'depth', @(s) any(strcmpi(s,{'depth','lines'})));
p.addParameter('z_max',     4850, @isnumeric);
p.addParameter('CAR',       false, @islogical);
p.addParameter('Normalize', 'none', @(s) any(strcmpi(s,{'none','freq'})));
p.addParameter('climPct',   [2 98], @isnumeric);
p.addParameter('clim',      [],   @isnumeric);
p.addParameter('SaveToggle',false, @islogical);
p.parse(varargin{:});
o = p.Results;

% ---------------------------------------------------------------- load LFP window
rec = loadRawRecording(char(o.lfpDir));
fs  = rec.fs;
if o.fMax > fs/2
    warning('fMax (%g Hz) exceeds Nyquist (%g Hz); clamping.', o.fMax, fs/2);
    o.fMax = floor(fs/2) - 1;
end

X = getTraces(rec, [o.t_start_s, o.t_start_s + o.win_s], [], 'uV');   % [nChan x nSamp] µV
[nChan, nSamp] = size(X);
if nSamp < 4, error('plotLFPWelch:emptyWindow', 'Requested window is too short.'); end
if o.CAR, X = X - mean(X, 1); end
X = X - mean(X, 2);                              % per-channel DC removal

% ---------------------------------------------------------------- geometry / depth
geom = recordingGeometry(char(o.lfpDir));
if geom.nChannels ~= nChan
    warning('Geometry has %d channels but data has %d; using min.', geom.nChannels, nChan);
end
z = geom.y(1:nChan);
keep = z <= o.z_max;
[zsorted, ord] = sort(z(keep));
idxKeep = find(keep);
order   = idxKeep(ord);
Xk = X(order, :);
nKept = numel(order);

% ---------------------------------------------------------------- Welch PSD
[psd, fN] = welchPSD_(Xk, fs, o.SegLen_s, o.Overlap);   % psd [nKept x nF] µV^2/Hz
% restrict to display band (drop DC; log axis can't show f=0)
sel = fN >= max(o.fMin, fN(2)) & fN <= o.fMax;
fN  = fN(sel);
psd = psd(:, sel);
psdNativeDB = 10*log10(psd + eps);

% ---------------------------------------------------------------- assemble figure
figure('Color','w','Position',[100 100 900 760]);

switch lower(o.Style)
    % ------------------------------------------------------------ depth map
    case 'depth'
        % regrid onto a log-spaced frequency axis (comparable to the Morlet fig)
        fLog  = logspace(log10(o.fMin), log10(o.fMax), o.nFreq);
        fLog  = fLog(fLog >= fN(1) & fLog <= fN(end));
        psdDB = zeros(nKept, numel(fLog));
        for c = 1:nKept
            psdDB(c,:) = interp1(fN, psdNativeDB(c,:), fLog, 'linear');
        end

        switch lower(o.Normalize)
            case 'freq'
                psdDB = psdDB - median(psdDB, 1);
                cbLabel = 'LFP PSD (dB re. median depth)';
            otherwise
                cbLabel = 'LFP PSD (dB/Hz)';
        end

        if ~isempty(o.clim)
            cl = o.clim(:)';
        else
            cl = prctile_(psdDB(:), o.climPct);
            if diff(cl) <= 0, cl = [min(psdDB(:)) max(psdDB(:))]; end
        end

        imagesc(fLog, zsorted, psdDB, cl);
        set(gca, 'YDir', 'normal', 'XScale', 'log');
        colormap(viridisColors(256));
        xlim([o.fMin o.fMax]); xticks([1 2 4 8 16 32 64 128 256]);
        xlabel('frequency (Hz)'); ylabel('depth along probe (\mum)');
        title(sprintf('probe00a LFP PSD by depth (Welch)  |  T=%g s, %g s / %g s seg  |  %s', ...
              o.t_start_s, o.win_s, o.SegLen_s, o.Normalize));
        cb = colorbar; cb.Label.String = cbLabel; box on;
        info.f = fLog; info.psdDB = psdDB;

    % ------------------------------------------------------------ classic PSD lines
    case 'lines'
        cmap = viridisColors(nKept);            % row 1 = shallowest (min depth)
        hold on;
        for c = 1:nKept
            h = plot(fN, psdNativeDB(c,:), '-', 'LineWidth', 0.5);
            h.Color = [cmap(c,:) 0.30];         % depth-colored, semi-transparent
        end
        plot(fN, mean(psdNativeDB,1), 'k-', 'LineWidth', 1.8);   % mean across depth
        set(gca, 'XScale', 'log');
        xlim([o.fMin o.fMax]); xticks([1 2 4 8 16 32 64 128 256]);
        grid on;
        xlabel('frequency (Hz)'); ylabel('LFP PSD (dB/Hz)');
        title(sprintf('probe00a LFP power spectrum (Welch)  |  T=%g s, %g s / %g s seg  |  black = mean', ...
              o.t_start_s, o.win_s, o.SegLen_s));
        colormap(viridisColors(256));
        clim([min(zsorted) max(zsorted)]);
        cb = colorbar; cb.Label.String = 'depth along probe (\mum)';
        info.f = fN; info.psdDB = psdNativeDB;
end

fprintf(['LFP Welch PSD (%s): %d/%d channels (z<=%g), %g s window, %g s segments ' ...
         '(~%.2f Hz res), %.3g-%.3g Hz\n'], ...
    o.Style, nKept, nChan, o.z_max, o.win_s, o.SegLen_s, 1/o.SegLen_s, fN(1), fN(end));

% ---------------------------------------------------------------- output / save
info.depth = zsorted; info.fNative = fN; info.psdNativeDB = psdNativeDB;
info.t_start_s = o.t_start_s; info.win_s = o.win_s; info.fs = fs; info.segLen_s = o.SegLen_s;

    if o.SaveToggle
        figType = 'lfp';
        figName = 'lfpPowerDepthWelch';
        qcFigSave(figType, figName);
    end
end

% ======================================================================
function [pxx, f] = welchPSD_(X, fs, segLen_s, overlap)
% One-sided Welch PSD. X [nChan x nSamp] -> pxx [nChan x nFreq] in units^2/Hz.
% Hann window, `overlap` fractional overlap, averaged periodograms.
nSamp = size(X, 2);
L = max(4, round(segLen_s * fs));
L = min(L, nSamp);
step = max(1, round(L * (1 - overlap)));
starts = 1 : step : (nSamp - L + 1);
if isempty(starts), starts = 1; L = nSamp; end

w = 0.5 * (1 - cos(2*pi*(0:L-1)/(L-1)));        % symmetric Hann
S = sum(w.^2);                                   % window power (for PSD norm)
nOneSided = floor(L/2) + 1;
f = (0:nOneSided-1) * (fs / L);

pxx = zeros(size(X,1), nOneSided);
for s0 = starts
    seg = X(:, s0:s0+L-1) .* w;                  % broadcast window over channels
    Y   = fft(seg, L, 2);
    P   = (abs(Y(:, 1:nOneSided)).^2) / (fs * S);
    P(:, 2:end-1) = 2 * P(:, 2:end-1);           % one-sided (not DC/Nyquist)
    if mod(L,2) == 1, P(:, end) = 2 * P(:, end); end  % odd L has no true Nyquist bin
    pxx = pxx + P;
end
pxx = pxx / numel(starts);
end

% ----------------------------------------------------------------------
function q = prctile_(x, pcts)
% Percentile without the Statistics Toolbox (linear interpolation).
x = sort(x(~isnan(x(:))));
if isempty(x), q = nan(size(pcts)); return; end
n = numel(x); q = zeros(size(pcts));
for i = 1:numel(pcts)
    r = pcts(i)/100 * (n-1) + 1;
    lo = max(1, floor(r)); hi = min(n, ceil(r));
    q(i) = x(lo) + (r-lo)*(x(hi)-x(lo));
end
end
