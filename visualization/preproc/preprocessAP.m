function Y = preprocessAP(X, fs, varargin)
% preprocessAP  AP-band preprocessing: high-pass + common median reference.
%
%   Y = preprocessAP(X, fs)
%   Y = preprocessAP(X, fs, 'HighpassHz', 300, 'CMR', true)
%
% The standard pre-spike-sorting / raw-AP-analysis chain (SKILL Step 2):
% high-pass ~300 Hz then subtract the across-channel median (global CMR).
%
% INPUTS
%   X   [nChannels x nSamples] traces in µV (e.g. from getTraces(...,'uV')).
%   fs  sampling frequency (Hz).
% OPTIONS
%   'HighpassHz'  cutoff (default 300; 0 disables)
%   'CMR'         apply common median reference (default true)
%
% NOTE for probe00a: the cached `probe00a_pre` file is ALREADY high-pass+CMR'd,
% so running this on it is (near) idempotent. This function is what you apply
% when starting from the WIDEBAND raw (the .cbin), not the pre-processed cache.
%
% See also preprocessLFP, preprocessTraces, checkFilterState.

p = inputParser;
p.addParameter('HighpassHz', 300);
p.addParameter('CMR', true);
p.parse(varargin{:});
o = p.Results;

% preprocessTraces already implements HP (Butterworth/movmean) + CMR.
Y = preprocessTraces(X, fs, o.HighpassHz, o.CMR);
end
