function traces = getTraces(rec, tspan, channels, unit)
% GETTRACES Pull a time/channel window out of a memory-mapped recording.
%
%   traces = getTraces(rec, tspan)
%   traces = getTraces(rec, tspan, channels)
%   traces = getTraces(rec, tspan, channels, unit)
%
% INPUTS
%   rec      - struct from loadRawRecording.
%   tspan    - [t0 t1] window in SECONDS, or [] for the whole recording
%              (careful: the whole thing is huge). Half-open [t0, t1).
%   channels - vector of 1-based channel indices to return (default: all).
%   unit     - 'uV' (default) applies gain/offset and returns double microvolts,
%              or 'raw' returns the native int16 samples untouched.
%
% OUTPUT
%   traces   - [nChannelsRequested x nSamplesInWindow] matrix.
%              Rows follow the order of `channels`; columns are time.
%
% See also loadRawRecording.

if nargin < 3 || isempty(channels)
    channels = 1:rec.nChannels;
end
if nargin < 4 || isempty(unit)
    unit = 'uV';
end

% ---- resolve sample range ---------------------------------------------
if nargin < 2 || isempty(tspan)
    s0 = 1;
    s1 = rec.nSamples;
else
    if numel(tspan) ~= 2
        error('getTraces:tspan', 'tspan must be [t0 t1] in seconds, or [].');
    end
    s0 = floor(tspan(1) * rec.fs) + 1;      % 1-based, inclusive
    s1 = floor(tspan(2) * rec.fs);          % exclusive upper -> last full sample
    s0 = max(s0, 1);
    s1 = min(s1, rec.nSamples);
    if s1 < s0
        traces = zeros(numel(channels), 0);
        return;
    end
end

if any(channels < 1 | channels > rec.nChannels)
    error('getTraces:chan', 'Channel index out of range [1 %d].', rec.nChannels);
end

% ---- read the window from the memmap ----------------------------------
% map.Data.x is [nChannels x nSamples]; slicing pulls only these bytes.
raw = rec.map.Data.x(channels, s0:s1);

switch lower(unit)
    case 'raw'
        traces = raw;
    case {'uv', 'microvolts'}
        g = rec.gain_to_uV(channels);
        o = rec.offset_to_uV(channels);
        traces = double(raw) .* g + o;
    otherwise
        error('getTraces:unit', 'unit must be ''uV'' or ''raw''.');
end
end
