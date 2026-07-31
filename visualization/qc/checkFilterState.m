function fs_out = checkFilterState(rec, varargin)
% checkFilterState  Estimate the power spectrum to reveal prior filtering.
%
%   info = checkFilterState(rec, Name, Value, ...)
%
% Computes an average power spectral density over a set of channels and a short
% window, then reports the ratio of low-band (< cutoff) to high-band power. A
% recording that has already been high-pass filtered at ~300 Hz shows a steep
% roll-off below 300 Hz (little LFP power) — which means it is ready for AP
% analysis but CANNOT yield LFP by low-pass filtering.
%
% NAME-VALUE OPTIONS
%   'Channels'  channels to average (default 32 spread across the probe)
%   't0'        window start (s)                (default midpoint)
%   'WinSec'    window length (s)               (default 2)
%   'CutoffHz'  low/high split for the ratio    (default 300)
%   'Plot'      draw the PSD                     (default true)
%
% OUTPUT (struct info)
%   .f, .psd          frequency (Hz) and mean PSD (µV^2/Hz)
%   .lowPower/.highPower  band-integrated power below/above cutoff
%   .lowFrac          lowPower / total  (small => already high-passed)
%   .likelyHighpassed logical

p = inputParser;
p.addParameter('Channels', []);
p.addParameter('t0', []);
p.addParameter('WinSec', 2);
p.addParameter('CutoffHz', 300);
p.addParameter('Plot', true);
p.parse(varargin{:});
o = p.Results;

if isempty(o.Channels)
    o.Channels = round(linspace(1, rec.nChannels, min(32, rec.nChannels)));
end
if isempty(o.t0), o.t0 = rec.duration/2; end

X = getTraces(rec, [o.t0, o.t0+o.WinSec], o.Channels, 'uV');   % [nch x nsamp]
X = X - mean(X, 2);                                            % remove DC per channel

nfft = 2^nextpow2(min(size(X,2), round(rec.fs)));   % ~1 s segments
[psd, f] = welchPSD(X, rec.fs, nfft);

df = f(2) - f(1);
low  = f < o.CutoffHz;
high = f >= o.CutoffHz;
lowP  = sum(psd(low))  * df;
highP = sum(psd(high)) * df;
lowFrac = lowP / (lowP + highP);

fs_out = struct('f', f, 'psd', psd, 'lowPower', lowP, 'highPower', highP, ...
                'lowFrac', lowFrac, 'likelyHighpassed', lowFrac < 0.2);

fprintf('Filter-state check (< %g Hz vs above):\n', o.CutoffHz);
fprintf('  low-band power fraction = %.1f%%  =>  %s\n', 100*lowFrac, ...
    ternary(fs_out.likelyHighpassed, ...
        'ALREADY HIGH-PASSED (LFP removed; ready for AP, no LFP here)', ...
        'wideband (LFP present; low-pass would yield real LFP)'));

if o.Plot
    figure('Color','w','Position',[100 100 700 500]);
    loglog(f, psd, 'LineWidth', 1.2); hold on;
    xline(o.CutoffHz, '--', sprintf('%g Hz', o.CutoffHz), 'Color',[0.85 0.1 0.1]);
    grid on; xlabel('frequency (Hz)'); ylabel('PSD (µV^2/Hz)');
    title('Mean power spectrum (reveals prior high-pass)');
    xlim([1 rec.fs/2]);
end
end

% ----------------------------------------------------------------------
function [psd, f] = welchPSD(X, fs, nfft)
% Welch PSD averaged over channels and 50%-overlap segments. Uses pwelch if the
% Signal Processing Toolbox is present, else a base-MATLAB Hann-window Welch.
if exist('pwelch', 'file') == 2
    [psd, f] = pwelch(X', hann(nfft), floor(nfft/2), nfft, fs);
    psd = mean(psd, 2);
    return;
end
w = 0.5 - 0.5*cos(2*pi*(0:nfft-1)'/(nfft-1));   % Hann
U = sum(w.^2);
step = floor(nfft/2);
acc = zeros(nfft/2+1, 1); nseg = 0;
for c = 1:size(X,1)
    x = X(c,:)';
    starts = 1:step:(numel(x)-nfft+1);
    for s = starts
        seg = x(s:s+nfft-1) .* w;
        P = abs(fft(seg)).^2 / (fs*U);
        P = P(1:nfft/2+1); P(2:end-1) = 2*P(2:end-1);
        acc = acc + P; nseg = nseg + 1;
    end
end
psd = acc / max(nseg,1);
f = (0:nfft/2)' * (fs/nfft);
end

% ----------------------------------------------------------------------
function out = ternary(c,a,b), if c, out=a; else, out=b; end
end
