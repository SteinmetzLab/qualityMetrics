function [figs, plotChan] = plotRawTraces(rec, ksResults, recordingType, z_max, t0, varargin)
%PLOTRAWTRACES  Raw-trace waterfalls across several timescales.
%
%   figs = plotRawTraceScales(rec, t0, Name, Value, ...)
%
% Renders a stacked trace figure for the given recording on a window of 40 ms 
% as a multi-scale overview of spiking activity. Plots 32 channels. 
% Optionally saves each figure as its own image.
%
% INPUTS
%   rec             (struct)    = recording file from loadRawRecording
%   ksResults       (struct)    = KS4 outputs
%   recordingType   (char)      = ap or lfp recording
%   z_max           (double)    = out-of-brain boundary, um, if not processed
%   t0              (double)    = window start in seconds
%
% NAME-VALUE OPTIONS
%   'allTimescales'     (logical)   = plot all default timescales in seconds (default [0.01 0.1 0.04 1 10])
%   'manualTimescales'  (double)    = use specified window durations in seconds   
%   'StartChannel'      (double)    = first channel of the span     (default 1)
%   'passHz'            (double)    = manual filter frequency in Hz
%   'SaveToggle'        (logical)   = save each figure to an image  (default false)
%
% OUTPUT
%   figs     (graphics) = array of figure handles, one per timescale
%   plotChan (double)   = plotted channels

    %% input parser

    p = inputParser;
    p.addRequired('rec',            @(x) isstruct(x));
    p.addRequired('ksResults',      @(x) isstruct(x));
    p.addRequired('recordingType',  @ischar);
    p.addRequired('z_max',          @isnumeric);
    p.addRequired('t0',             @isnumeric);

    p.addParameter('allTimescales', false, @islogical);
    p.addParameter('manualTimescales', 0.04, @isnumeric)
    p.addParameter('StartChannel', 1);
    p.addParameter('passHz', [], @isnumeric);
    p.addParameter('SaveToggle', false);
    p.parse(rec, ksResults, recordingType, z_max, t0, varargin{:});
    o = p.Results;

    %% params

    numChannels = 32;
    doCMR       = true;
    
    filePrefix  =  'rawtrace';

    %% load geometry and masks
    
    z = double(ksResults.chanpos(:,2));
    keepChan = z <= z_max;                  % exclude the top artifact band (z > z_max)
    idxKeep  = find(keepChan);              % channel indices (1-based) we display/plot
    [~, ord] = sort(z(idxKeep));            % ascending depth among kept channels
    order = idxKeep(ord);                   % display-row -> channel index (1-based)

    % ready channels
    
    c0 = max(1, round(o.StartChannel));
    c1 = min(rec.nChannels, c0 + round(numChannels) - 1);
    plotChan = order(c0:c1);
    
    timescales = 0.04;
    
    if o.allTimescales
        timescales = [0.01 0.1 0.04 1 10];
    end

    if o.manualTimescales
        timescales = o.manualTimescales;
    end
    
    % ---- guard against running off the end of the recording ----------------
    
    maxDur = max(timescales);
    if t0 + maxDur > rec.duration
        t0 = max(0, rec.duration - maxDur);
        warning('plotRawTraceScales:clamp', ...
            'Window would exceed recording; moved t0 to %.3f s.', t0);
    end

    %% extract, filter, and plot signal
    
    % load padded window for filter
    pad = 0.05;

    figs = gobjects(1, numel(timescales));
    for timescale_i = 1:numel(timescales)
        % extract and filter
        
        dur = timescales(timescale_i);

        % load each window
        a = t0 - pad;               % 260728 this may throw an error if t0 = 0
        b = t0 + dur + pad;
    
        % get channel x amplitude
        X = getTraces(rec, [a b], plotChan, 'uV');   % [nch x nsamp], only these bytes
    
        % absolute time of each returned column 
        idxA  = max(floor(a*rec.fs) + 1, 1);            % index for start of window
        nsamp = size(X, 2);                             % number of timepoints
        tabs  = ((idxA-1) + (0:nsamp-1)) / rec.fs;      % seconds, timepoints in seconds
    
        switch recordingType 
            case 'ap'           
                if ~rec.preprocessed
                    % preprocess
                    [Y, ff] = preprocessTraces(X, rec.fs, 'ap', 'doCMR', doCMR);
                else
                    Y = X;

                    ff.passHz       = 300; 
                    ff.bandpassType = 'highpass';
                    ff.CMR          = true;  
                end
                
                figType = 'ap';
            case 'lfp'
                if ~rec.preprocessed
                    % preprocess
                    [Y, ff] = preprocessTraces(X, rec.fs, 'lfp', 'doCMR', doCMR);
                else
                    [Y, ff] = preprocessTraces(X, rec.fs, 'lfp', 'preprocess', rec.preprocessed, 'passHz', 0.5,'doCMR', doCMR);
                end

                figType = 'lfp';
        end

                    
        mask = tabs >= t0 & tabs < t0 + dur;
        Y  = Y(:, mask);
        td = tabs(mask) - t0;   
        
        if isempty(td)
            error('plotWaterfall:empty', ...
                'No samples in [%.4f, %.4f] s (recording is %.1f s).', t0, t0+dur, rec.duration);
        end

        figs(timescale_i) = plotWaterfall(td, Y, t0, dur, plotChan, ff);
    
        if o.SaveToggle
            figName = sprintf('%s_%s_t%08.3fs_ch%d-%d.%s', ...
                filePrefix, durTag(dur), t0, plotChan(1), plotChan(end), 'png');
            qcFigSave(figType, figName);
        end
    end
end

% ----------------------------------------------------------------------
function s = durTag(dur)
    % Compact filename-safe timescale tag: 0.01 -> '10ms', 10 -> '10s'.
    if dur < 1
        s = sprintf('%gms', dur*1000);
    else
        s = sprintf('%gs', dur);
    end
    s = regexprep(s, '\.', 'p');
end
