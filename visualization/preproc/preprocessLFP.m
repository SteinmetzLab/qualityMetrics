function [Y, fsLFP] = preprocessLFP(X, fs, varargin)
% preprocessLFP  LFP-band preprocessing: low-pass filter + downsample.
%
%   [Y, fsLFP] = preprocessLFP(X, fs)
%   [Y, fsLFP] = preprocessLFP(X, fs, 'LowpassHz', 300, 'TargetFs', 2500)
%
% Derives the LFP band from a WIDEBAND recording: anti-alias low-pass then
% decimate to a lower sampling rate. On Neuropixels 2.0 there is no separate LF
% stream (lf_sample_frequency_hz = 0), so LFP MUST be obtained this way from the
% wideband/AP data.
%
% INPUTS
%   X   [nChannels x nSamples] WIDEBAND traces in µV.
%   fs  input sampling frequency (Hz).
% OPTIONS
%   'LowpassHz'  low-pass cutoff (default 300)
%   'TargetFs'   output sampling rate after decimation (default 2500)
%   'CAR'        subtract across-channel mean (common average ref) (default false)
%
% OUTPUTS
%   Y      [nChannels x nSamplesDown] LFP-band traces (µV)
%   fsLFP  output sampling rate (Hz)
%
% ⚠ IMPORTANT for probe00a: the cached `probe00a_pre` file is ALREADY 300 Hz
% high-passed — its LFP band has been removed. Running this on it yields
% near-zero garbage. Use checkFilterState() to confirm a recording is wideband
% BEFORE calling this. The wideband source for probe00a is the original .cbin
% (Y:\temp\QBtempdata\probe00a\...ap.cbin), not the pre-processed cache.
%
% See also preprocessAP, checkFilterState.

p = inputParser;
p.addParameter('LowpassHz', 300);
p.addParameter('TargetFs', 2500);
p.addParameter('CAR', false);
p.parse(varargin{:});
o = p.Results;

X = double(X);
if o.CAR
    X = X - mean(X, 1);
end

% integer decimation factor
q = max(1, round(fs / o.TargetFs));
fsLFP = fs / q;

Y = lowpassFilter(X, fs, o.LowpassHz);
if q > 1
    Y = Y(:, 1:q:end);       % decimate (already anti-alias low-passed above)
end
end

% ----------------------------------------------------------------------
function Y = lowpassFilter(X, fs, fc)
% Zero-phase low-pass. Butterworth (Signal Processing Toolbox) if available,
% else a base-MATLAB moving-average (boxcar) low-pass.
try
    [b, a] = butter(4, fc/(fs/2), 'low');
    Y = filtfilt(b, a, X')';
catch
    win = max(3, round(fs/fc));
    Y = movmean(X, win, 2);
end
end
