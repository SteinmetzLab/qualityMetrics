function rec = loadRawRecording(rawPath)
% LOADRAWRECORDING  Memory-map a SpikeInterface BinaryRecordingExtractor .raw file.
%
%   rec = loadRawRecording(rawPath)
%
% Loads a raw (un-drift-corrected, unsorted) extracellular recording that was
% cached by SpikeInterface as a flat binary "traces_cached_seg0.raw" file. The
% file is memory-mapped rather than read into RAM, so this works on files far
% larger than physical memory (e.g. the ~143 GB probe00a recording).
%
% INPUT
%   rawPath (char)    = * the .raw file itself, or
%                       * the folder containing binary.json + traces_cached_seg0.raw
%
% OUTPUT (struct rec)
%   .filePath     absolute path to the .raw file
%   .fs           sampling frequency (Hz)
%   .nChannels    number of channels
%   .nSamples     number of time samples
%   .duration     recording duration (s)
%   .dtype        MATLAB class of the samples (e.g. 'int16')
%   .gain_to_uV   [nChannels x 1] multiply raw samples by this to get microvolts
%   .offset_to_uV [nChannels x 1] add after gain to get microvolts
%   .channelIds   cell array of channel id strings (from binary.json)
%   .timeAxis     0 => file laid out (samples x channels) C-order
%   .map          memmapfile handle; map.Data.x is [nChannels x nSamples]
%   .preprocessed logical; true if provenance.json shows the cached data was
%                 already filtered / referenced (see below)
%   .preprocSteps cell array of the SpikeInterface preprocessing step names
%                 found in provenance.json (e.g. {'HighpassFilterRecording',
%                 'CommonReferenceRecording'}); empty when none / no provenance
%
% Metadata is read from the sidecar binary.json so this generalises to other
% recordings (channel count, sample rate, dtype and gain are NOT hardcoded).
%
% GUARD AGAINST DOUBLE-PROCESSING
%   A SpikeInterface "*_pre" cache is byte-for-byte indistinguishable from a raw
%   cache, so loading one and then running preprocessTraces would filter / CMR
%   the signal a second time. This loader reads provenance.json (if present) and
%   sets rec.preprocessed = true (with a warning) when the data already carries
%   filtering / referencing steps. Callers should check rec.preprocessed and
%   skip preprocessTraces when it is true.
%
% See also getTraces, memmapfile.

% ---- resolve paths -----------------------------------------------------
if isfolder(rawPath)
    folder  = rawPath;
    rawFile = fullfile(folder, 'traces_cached_seg0.raw');
else
    [folder, ~, ~] = fileparts(rawPath);
    rawFile = rawPath;
end
if ~isfile(rawFile)
    error('loadRawRecording:noFile', 'Raw file not found: %s', rawFile);
end
jsonFile = fullfile(folder, 'binary.json');
if ~isfile(jsonFile)
    error('loadRawRecording:noJson', ...
        ['binary.json not found next to the raw file (%s).\n' ...
         'It carries the channel count / sample rate / dtype / gain.'], jsonFile);
end

% ---- parse binary.json -------------------------------------------------
j  = jsondecode(fileread(jsonFile));
kw = j.kwargs;

fs        = kw.sampling_frequency;
nChannels = double(kw.num_channels);
timeAxis  = 0;
if isfield(kw, 'time_axis') && ~isempty(kw.time_axis)
    timeAxis = double(kw.time_axis);
end
if timeAxis ~= 0
    error('loadRawRecording:timeAxis', ...
        ['binary.json has time_axis=%d; this loader assumes time_axis=0 ' ...
         '(samples x channels). Handle the transposed layout before mapping.'], timeAxis);
end
fileOffset = 0;
if isfield(kw, 'file_offset') && ~isempty(kw.file_offset)
    fileOffset = double(kw.file_offset);
end

dtype = mapDtype(kw.dtype);          % e.g. '<i2' -> 'int16'
[~, bytesPerSample] = dtypeInfo(dtype);

gain   = columnize(kw.gain_to_uV,   nChannels);
offset = columnize(kw.offset_to_uV, nChannels);

channelIds = {};
if isfield(kw, 'channel_ids')
    channelIds = cellstr(string(kw.channel_ids));
end

% ---- derive sample count from file size & sanity-check -----------------
d = dir(rawFile);
nBytes = d.bytes - fileOffset;
denom  = nChannels * bytesPerSample;
if mod(nBytes, denom) ~= 0
    error('loadRawRecording:badSize', ...
        ['File size (%d bytes after offset) is not divisible by ' ...
         'nChannels*bytesPerSample (%d). Metadata and file disagree.'], nBytes, denom);
end
nSamples = nBytes / denom;

% ---- memory-map --------------------------------------------------------
% Column-major [nChannels x nSamples] matches the C-order (samples x channels)
% on-disk layout: element x(c,s) sits at flat index (s-1)*nChannels + c.
map = memmapfile(rawFile, ...
    'Offset', fileOffset, ...
    'Format', {dtype, [nChannels nSamples], 'x'}, ...
    'Writable', false);

% ---- assemble output ---------------------------------------------------
rec = struct();
rec.filePath     = getAbsPath(rawFile);
rec.fs           = fs;
rec.nChannels    = nChannels;
rec.nSamples     = nSamples;
rec.duration     = nSamples / fs;
rec.dtype        = dtype;
rec.gain_to_uV   = gain;
rec.offset_to_uV = offset;
rec.channelIds   = channelIds;
rec.timeAxis     = timeAxis;
rec.map          = map;

% ---- guard: was this cache already filtered / referenced? --------------
% A "*_pre" folder is structurally identical to a raw one; loading it and then
% running preprocessTraces would double-filter / double-CMR the signal.
[rec.preprocessed, rec.preprocSteps] = detectPreprocessing(folder);
if rec.preprocessed
    warning('loadRawRecording:preprocessed', ...
        ['This recording appears to be ALREADY preprocessed ' ...
         '(provenance.json lists: %s).\n' ...
         '  Downstream filtering/CMR would double-process the signal. ' ...
         'Check rec.preprocessed and skip preprocessTraces when true.'], ...
        strjoin(rec.preprocSteps, ', '));
end

fprintf(['Loaded %s\n  %d channels x %d samples @ %.1f Hz  (%.1f min, %.1f GB)\n' ...
         '  dtype=%s, gain_to_uV=%.6g%s\n'], ...
    rec.filePath, nChannels, nSamples, fs, rec.duration/60, d.bytes/1e9, ...
    dtype, gain(1), ternary(all(gain==gain(1)), ' (uniform)', ' (per-channel)'));
end

% ======================================================================
function out = columnize(v, n)
% Force to an [n x 1] column, broadcasting a scalar.
v = double(v(:));
if isscalar(v)
    out = repmat(v, n, 1);
elseif numel(v) == n
    out = v;
else
    error('loadRawRecording:vecLen', ...
        'Expected 1 or %d values, got %d.', n, numel(v));
end
end

% ----------------------------------------------------------------------
function dt = mapDtype(s)
% Map a numpy-style dtype string ('<i2') to a MATLAB class.
s = char(s);
key = regexprep(s, '^[<>=|]', '');   % strip byte-order char
switch lower(key)
    case 'i2', dt = 'int16';
    case 'i4', dt = 'int32';
    case 'i1', dt = 'int8';
    case 'u2', dt = 'uint16';
    case 'u1', dt = 'uint8';
    case 'f4', dt = 'single';
    case 'f8', dt = 'double';
    otherwise
        % maybe already a MATLAB class name
        if any(strcmp(s, {'int16','int32','int8','uint16','uint8','single','double'}))
            dt = s;
        else
            error('loadRawRecording:dtype', 'Unsupported dtype "%s".', s);
        end
end
% Note: byte-order prefix is not honoured by memmapfile; assumes native/LE.
end

% ----------------------------------------------------------------------
function [cls, nbytes] = dtypeInfo(dt)
cls = dt;
switch dt
    case {'int8','uint8'},   nbytes = 1;
    case {'int16','uint16'}, nbytes = 2;
    case {'int32','single'}, nbytes = 4;
    case 'double',           nbytes = 8;
    otherwise, error('loadRawRecording:dtypeInfo', 'Unknown dtype %s', dt);
end
end

% ----------------------------------------------------------------------
function p = getAbsPath(f)
d = dir(f);
p = fullfile(d.folder, d.name);
end

% ----------------------------------------------------------------------
function out = ternary(cond, a, b)
if cond, out = a; else, out = b; end
end

% ----------------------------------------------------------------------
function [isPre, steps] = detectPreprocessing(folder)
    % Scan provenance.json (if present) for SpikeInterface preprocessing steps.
    % Returns isPre=true and a de-duplicated list of the step class names when the
    % cached data was already filtered / referenced, so callers can avoid
    % double-processing. No provenance.json (older/other recordings) => not flagged.
    isPre = false;
    steps = {};
    provFile = fullfile(folder, 'provenance.json');
    if ~isfile(provFile)
        return;
    end
    txt = fileread(provFile);
    % e.g. "class": "spikeinterface.preprocessing.filter.HighpassFilterRecording"
    tok = regexp(txt, ...
        '"class"\s*:\s*"spikeinterface\.preprocessing\.[^"]*\.(\w+)"', 'tokens');
    if isempty(tok)
        return;
    end
    steps = unique(cellfun(@(c) c{1}, tok, 'UniformOutput', false), 'stable');
    isPre = true;
end
