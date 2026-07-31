function geom = recordingGeometry(recOrFolder)
% recordingGeometry  Per-channel probe geometry aligned to the data channels.
%
%   geom = recordingGeometry(rec)          % rec from loadRawRecording
%   geom = recordingGeometry(folder)       % the probe00a_pre folder
%
% Reads the SpikeInterface property files (properties/location.npy, group.npy)
% that are aligned 1:1 with the recording's channels, plus probe.json for probe
% metadata. Returns geometry and a printed summary (columns, pitch, shanks).
%
% OUTPUT (struct geom)
%   .x, .y        [nCh x 1] site coordinates (µm); y is depth
%   .shank        [nCh x 1] shank id per channel (from group.npy)
%   .nChannels    channel count
%   .model        probe model string (e.g. 'NP2021')
%   .xColumns     unique x positions (µm)
%   .yPitch       median spacing between adjacent depths (µm)
%   .yRange       [min max] depth (µm)
%   .nShanks      number of distinct shanks in the data

if ischar(recOrFolder) || isstring(recOrFolder)
    folder = char(recOrFolder);
    if endsWith(folder, '.raw'), folder = fileparts(folder); end
else
    folder = fileparts(recOrFolder.filePath);
end

propDir = fullfile(folder, 'properties');
loc   = readNPY(fullfile(propDir, 'location.npy'));      % [nCh x 2] (x, y)
group = readNPY(fullfile(propDir, 'group.npy'));         % [nCh x 1] shank/group

geom = struct();
geom.x = double(loc(:,1));
geom.y = double(loc(:,2));
geom.shank = double(group(:));
geom.nChannels = size(loc, 1);

% probe model from probe.json (optional)
geom.model = '';
pj = fullfile(folder, 'probe.json');
if isfile(pj)
    p = jsondecode(fileread(pj));
    try, geom.model = p.probes(1).annotations.model_name; catch, end
end

% summaries
geom.xColumns = unique(round(geom.x, 3));
ys = sort(unique(round(geom.y, 3)));
dy = diff(ys); dy = dy(dy > 0);
geom.yPitch = median(dy);
geom.yRange = [min(geom.y) max(geom.y)];
geom.nShanks = numel(unique(geom.shank));

fprintf('Geometry (%s, %d channels)\n', geom.model, geom.nChannels);
fprintf('  shanks: %d   |   x columns (µm): %s\n', ...
    geom.nShanks, mat2str(geom.xColumns'));
fprintf('  depth y: %.1f – %.1f µm   |   pitch ≈ %.1f µm\n', ...
    geom.yRange(1), geom.yRange(2), geom.yPitch);
end
