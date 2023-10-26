function [neuron_pass, rp_pass, nc_pass, ma_pass] = computeQualityMetrics(spikeTimes, spikeAmps, medAmpThresh, verbose)
% 
%   This code computes the full set of 3 quality metrics below as used in the IBL and Steinmetz Lab.
%   Metrics:
%         1) Sliding Refractory Period metric (slidingRefractory)
%         2) Noise Cutoff metric (noiseCutoff)
%         3) median amplitude threshold (default threshold is 50 uV)
%     Parameters
%     ----------
%     spikeTimes : spike times for the neuron across the entire recording
%     spikeAmps: amplitudes (or template amplitudes) for the neuron across the entire recording
%     Returns
%     -------
%     pass : binary variable which states whether your neuron passes (1) or fails (0)

if nargin<5
    medAmpThresh = 50;
    verbose = true;
    params = struct();
end
    

if verbose
    print("Computing metrics...")
end

%1) Compute RP values:
    [maxConfidenceAt10Cont, minContWith90Confidence, timeOfLowestCont,...
        nSpikesBelow2, confMatrix, cont, rp, nACG, firingRate] ...
        = slidingRP(spikeTimes, params);

    % Also compute whether RP pass:
    if minContWith90Confidence < 10
        rp_pass = true;
    else
        rp_pass = false;
    end

%2) Compute NC values and whether neuron passes NC:
    [nc_pass, cutoff_value, first_low_quantile] = noise_cutoff(spikeAmps);

%3) Compute median amplitude:
    medAmp = nanmedian(spikeAmps) * 1e6;  %median amplitude in microvolts

    if (medAmp > medAmpThresh)
        ma_pass = true;
    else
        ma_pass = false;
    end

%Now compute whether neuron passes all 3 metrics, return as boolean
    neuron_pass = rp_pass & nc_pass & ma_pass;
end





