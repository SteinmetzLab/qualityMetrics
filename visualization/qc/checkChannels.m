function chk = checkChannels(rec, varargin)
% checkChannels  Per-channel RMS to flag dead / noisy channels.
%
%   chk = checkChannels(rec, Name, Value, ...)
%
% Samples several short windows across the recording, computes each channel's
% RMS in µV (after an optional high-pass so LFP wander doesn't dominate), and
% flags channels that look dead (near-zero RMS) or noisy (RMS far above the
% array median). This is the Step 1 "verify all probe channels are good" check.
%
% NAME-VALUE OPTIONS
%   'NumChunks'   windows sampled across the session   (default 8)
%   'ChunkSec'    length of each window (s)            (default 1)
%   'HighpassHz'  high-pass before RMS, 0 to skip      (default 300)
%   'DeadUV'      RMS below this µV => dead             (default 2)
%   'NoisyMADs'   RMS above median + k*MAD => noisy     (default 6)
%   'Geometry'    geom struct (from recordingGeometry) to order by depth
%   'Plot'        draw RMS-vs-depth/channel figure      (default true)
%
% OUTPUT (struct chk)
%   .rms      [nCh x 1] median RMS across chunks (µV)
%   .status   [nCh x 1] categorical/cellstr: 'good' | 'dead' | 'noisy'
%   .isGood   [nCh x 1] logical
%   .nGood/.nDead/.nNoisy counts
%   .badChannels indices of non-good channels

p = inputParser;
p.addParameter('NumChunks', 8);
p.addParameter('ChunkSec', 1);
p.addParameter('HighpassHz', 300);
p.addParameter('DeadUV', 2);
p.addParameter('NoisyMADs', 6);
p.addParameter('Geometry', []);
p.addParameter('Plot', true);
p.parse(varargin{:});
o = p.Results;

nCh = rec.nChannels;
% chunk start times spread across the session (avoid the very edges)
t0s = linspace(0.02*rec.duration, 0.95*rec.duration, o.NumChunks);

rmsAll = zeros(nCh, o.NumChunks);
for k = 1:o.NumChunks
    X = getTraces(rec, [t0s(k), t0s(k)+o.ChunkSec], 1:nCh, 'uV');
    if ~isempty(o.HighpassHz) && o.HighpassHz > 0
        X = preprocessTraces(X, rec.fs, o.HighpassHz, false);  % HP only, no CMR
    end
    rmsAll(:, k) = sqrt(mean(X.^2, 2));
end
rms = median(rmsAll, 2);

% classify
med = median(rms);
madv = 1.4826 * median(abs(rms - med));
isDead  = rms < o.DeadUV;
isNoisy = ~isDead & (rms > med + o.NoisyMADs*madv);
isGood  = ~isDead & ~isNoisy;

status = repmat({'good'}, nCh, 1);
status(isDead)  = {'dead'};
status(isNoisy) = {'noisy'};

chk = struct();
chk.rms = rms;
chk.status = status;
chk.isGood = isGood;
chk.nGood = nnz(isGood);
chk.nDead = nnz(isDead);
chk.nNoisy = nnz(isNoisy);
chk.badChannels = find(~isGood);

fprintf('Channel check: %d good, %d dead, %d noisy (of %d)\n', ...
    chk.nGood, chk.nDead, chk.nNoisy, nCh);
fprintf('  RMS µV: median %.1f, range %.1f–%.1f\n', med, min(rms), max(rms));
if ~isempty(chk.badChannels)
    fprintf('  flagged channels: %s\n', mat2str(chk.badChannels'));
end

if o.Plot
    figure('Color','w','Position',[100 100 600 800]);
    if ~isempty(o.Geometry)
        y = o.Geometry.y; ylab = 'depth (µm)';
    else
        y = (1:nCh)'; ylab = 'channel index';
    end
    hold on;
    plot(rms(isGood),  y(isGood),  '.', 'Color',[0.2 0.5 0.9], 'MarkerSize',8);
    plot(rms(isNoisy), y(isNoisy), 'o', 'Color',[0.9 0.5 0.1], 'MarkerSize',6, 'LineWidth',1);
    plot(rms(isDead),  y(isDead),  'x', 'Color',[0.85 0.1 0.1], 'MarkerSize',8, 'LineWidth',1.5);
    xline(o.DeadUV, ':', 'dead', 'Color',[0.85 0.1 0.1]);
    xline(med + o.NoisyMADs*madv, ':', 'noisy', 'Color',[0.9 0.5 0.1]);
    xlabel('RMS (µV)'); ylabel(ylab);
    title(sprintf('Per-channel RMS  (%d good / %d)', chk.nGood, nCh));
    box off; set(gca,'TickDir','out');
end
end
