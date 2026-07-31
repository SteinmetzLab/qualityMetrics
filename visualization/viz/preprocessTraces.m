function [Y, filtFlags] = preprocessTraces(X, fs, recordingType, varargin)
%PREPROCESSTRACES  Filter and common median reference (CMR) for a recording.
%
%   [Y, filtFlags] = preprocessTraces(X, fs, recordingType, Name, Value, ...)
%
% Filters the traces then applies CMR. AP recordings are high-pass filtered;
% LFP recordings are band-pass filtered (or high-pass filtered when already
% preprocessed). Filtering runs first, then CMR subtracts the across-channel
% median at each time sample.
%
% INPUTS
%   X               (numeric)   = [nChannels x nSamples] traces (uV or raw)
%   fs              (double)    = sampling rate in Hz
%   recordingType   (char)      = ap or lfp recording
%
% NAME-VALUE OPTIONS
%   'preprocess'    (logical)   = traces already preprocessed  (default false)
%   'passHz'        (double)    = filter cutoff(s) in Hz, [] = auto by type
%   'doCMR'         (logical)   = apply common median reference (default true)
%
% OUTPUT
%   Y           (double) = [nChannels x nSamples] filtered + CMR signal
%   filtFlags   (struct) = applied filter settings (passHz, bandpassType, CMR)

    %% input parser

    p = inputParser;
    p.addRequired('X',              @isnumeric);
    p.addRequired('fs',             @isnumeric);
    p.addRequired('recordingType',  @ischar);
    p.addParameter('preprocess', false, @islogical);
    p.addParameter('passHz', [],    @isnumeric);
    p.addParameter('doCMR', true,   @islogical);
    p.parse(X, fs, recordingType, varargin{:});
    o = p.Results;
        
    %% filters

    Y = double(o.X);

    % initialize filter return
    filtFlags.passHz        = [];
    filtFlags.bandpassType  = [];
    filtFlags.CMR           = [];

    switch o.recordingType
        case 'ap'
            % high pass for ap
            if isempty(o.passHz) || o.passHz < 0
                o.passHz = 300;
            end            
            
            Y = highpass(Y, o.passHz, o.fs);

            filtFlags.passHz = o.passHz;
            filtFlags.bandpassType = 'highpass';
        case 'lfp'
            if o.preprocess
                if isempty(o.passHz) || o.passHz < 0 || size(o.passHz, 2) > 1
                    o.passHz = 0.5;
                end

                Y = highpass(Y, o.passHz, o.fs);
                
                % CMR already performed manual set flags
                o.doCMR = false;

                filtFlags.passHz = [o.passHz 300];
                filtFlags.CMR = true;
            else
                % band pass for lfp
                if size(o.passHz, 2) ~= 2 || o.passHz < 0
                    o.passHz = [0.5 300];
                end

                Y = bandpass(Y, o.passHz, o.fs); 

                filtFlags.passHz = o.passHz;
            end
    
            filtFlags.bandpassType = 'bandpass';
    end


    %% CMR

    if o.doCMR
        Y = Y - median(Y, 1);
        filtFlags.CMR = true;
    end
end