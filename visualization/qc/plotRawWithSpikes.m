function plotRawWithSpikes(ksResults, fs, nChan, z_max, t0, win_ms, varargin)
    %PLOT_RAW_WITH_SPIKES  Grayscale raw-activity image (channels x time) with
    %Kilosort spikes overlaid: GREEN = spike from a 'good' unit, RED = 'mua'.
    %
    % Each channel row is baselined (per-channel mean/median over the window
    % subtracted). Decreases in activity are black, increases in activity are white.
    % Detected spikes are marked at their unit's peak channel.
    %
    % INPUTS
    %   ksResults   (struct)    = kilosort outputs structure. See autoqc_working
    %   fs          (double)    = sampling rate in hertz
    %   nChan       (double)    = number of identified channels
    %   t0          (double)    = window start in seconds
    %   z_max       (double)    = maximum probe depth to analyze in micrometers
    %   win_ms      (double)    = analysis window in milliseconds
    %
    % Name Value Pairs (OPTIONAL):
    %   'clim_uV'       symmetric grayscale saturation, +/- uV (default 30; ~ the "dB" knob)
    %   'gain_uV'       ADC-units -> microvolts scale        (default 3.02734375 for this probe)
    %   'baseline'      'mean' or 'median' per-channel baseline (default 'median')
    %   'traceChannel'  logical for saving figures (default 'false')
    %   'SaveToggle'    logical for saving figures (default 'false')

    %% input parser
    
    p = inputParser;
    p.addParameter('clim_uV',   100,   @isnumeric);
    p.addParameter('gain_uV',   3.02734375, @isnumeric);
    p.addParameter('baseline',  'median', @ischar);
    p.addParameter('traceChannel', [], @isnumeric);
    p.addParameter('SaveToggle',false, @islogical);
    p.parse(varargin{:});
    o = p.Results;

    %% load geometry and masks
    
    z = double(ksResults.chanpos(:,2));
    keepChan = z <= z_max;                  % exclude the top artifact band (z > z_max)
    idxKeep  = find(keepChan);              % channel indices (1-based) we display/plot
    [zsorted, ord] = sort(z(idxKeep));      % ascending depth among kept channels
    order = idxKeep(ord);                   % display-row -> channel index (1-based)
    
    % map ap/lfp trace channels to depth
    
    traceChanStart  = double(ksResults.chanpos(min(o.traceChannel),2));
    traceChanEnd    = double(ksResults.chanpos(max(o.traceChannel),2));
    
    % per-cluster peak channel, from templates (peak-to-peak argmax over channels).
    ptp   = squeeze(max(ksResults.templates,[],2) - min(ksResults.templates,[],2)); % (nClu, nChan)
    [~, peakch] = max(ptp, [], 2);         % 1-based channel index per cluster row
    
    % columns: cluster_id, KSLabel
    isGoodClu = strcmp(ksResults.lbl.KSLabel, 'good');           % logical per cluster row
    goodIds   = ksResults.lbl.cluster_id(isGoodClu);             % 0-based cluster ids
    
    %% read raws in window
    
    s0 = round(t0 * fs);                 % first sample (0-based)
    nS = round(win_ms/1000 * fs);               % samples in window
    fid = fopen(ksResults.recordingPath, 'r');
    if fid < 0, error('cannot open %s', ksResults.recordingPath); end
    fseek(fid, s0 * nChan * 2, 'bof');            % int16 = 2 bytes; sample-major layout
    X = fread(fid, [nChan, nS], '*int16');        % [nChan x nS] (column-major fill = correct)
    fclose(fid);
    X = single(X) * o.gain_uV;                    % -> microvolts
    
    % baseline each channel over the window
    
    switch o.baseline
        case 'mean'   
            X = X - mean(X,2);
        case 'median' 
            X = X - median(X,2);
        otherwise, error('baseline must be mean or median');
    end
    Xdisp = X(order, :);                          % reorder rows to ascending depth
    
    %% load spikes
    
    inWin = ksResults.st >= s0 & ksResults.st < (s0 + nS);
    st_w  = double(ksResults.st(inWin));
    sc_w  = double(ksResults.sc(inWin));
    
    spk_chan  = peakch(sc_w + 1);                 % 1-based peak channel per spike
    keepSpk   = keepChan(spk_chan);               % drop spikes on masked (artifact-band) channels
    st_w = st_w(keepSpk); sc_w = sc_w(keepSpk); spk_chan = spk_chan(keepSpk);
    spk_t     = (st_w - s0) / fs * 1000;          % ms within window (x-axis)
    spk_depth = z(spk_chan);                      % depth of the peak channel (= y-axis)
    spk_good  = ismember(sc_w, goodIds);          % true = good, false = mua
    
    %% plot
    
    f   = figure('Color','w','Position',[100 100 1100 720]);
    ax  = axes(f);
    
    t_ms = (0:nS-1) / fs * 1000;
    imagesc(t_ms, zsorted, Xdisp, [-o.clim_uV o.clim_uV]);
    colormap(gray);
    set(gca, 'YDir', 'normal');                   % depth increases upward
    hold on;
    
    % overlay: mua first (red), good on top (green) so good markers win visually
    
    plot(spk_t(~spk_good), spk_depth(~spk_good), '.', ...
         'Color',[1 0 0], 'MarkerSize', 15);
    plot(spk_t( spk_good), spk_depth( spk_good), '.', ...
         'Color',[0 1 0], 'MarkerSize', 15);
    
    %% cosmetics
    
    % channel inset mapping waveform channels to activity plot
    
    xmax = max(t_ms);
    xr = xmax * 0.012;                            % small inset from left edge
    plot(ax, [xr xr], [traceChanStart traceChanEnd], 'Color',[0 0 1], 'LineWidth', 2);
    
    % labels

    xlabel('time (ms)'); ylabel('depth along probe (\mum)');
    title(sprintf('probe00a raw AP + KS4 spikes  |  T = %g s, %g ms  |  green=good  red=mua', ...
          t0, win_ms));
    axis tight;
    cb = colorbar; cb.Label.String = 'baselined amplitude (\muV)';
    
    fprintf('window %g-%g ms @ %g s: %d spikes (%d good, %d mua) | %d/%d channels kept (z<=%g), %d spikes dropped in artifact band\n', ...
        0, win_ms, t0, numel(sc_w), nnz(spk_good), nnz(~spk_good), ...
        nnz(keepChan), nChan, z_max, nnz(~keepSpk));
    
    %% save

    if o.SaveToggle
        figType = 'ks';
        figName = 'rawWithSpikes';
        qcFigSave(figType, figName);
    end
end