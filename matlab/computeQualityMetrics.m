function [neuron_pass, rp_pass, nc_pass, ma_pass] = computeQualityMetrics(spikeTimes, spikeAmps)%, medAmpThresh = 50, verbose = True):
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

if nargin<3
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

    % Also compute whether RP pass: #TODO: add nSpikesBelow2 condition
    if minContWith90Confidence < 10
        rp_pass = true;
    else
        rp_pass = false;
    end

%2) Compute NC values and whether neuron passes NC:
    nc_pass, cutoff_value, first_low_quantile = noise_cutoff(spikeAmps)

%3) Compute median amplitude: #TODO: make this optional in params above
    medAmp = nanmedian(spikeAmps) * 1e6  %median amplitude in microvolts

    if (medAmp > medAmpThresh)
        ma_pass = True;
    else
        ma_pass = False;
    end

%Now compute whether neuron passes all 3 metrics, return as boolean
    neuron_pass = rp_pass & nc_pass & ma_pass
end




function [metrics] = computeQualityMetrics_All(spikeTimes, spikeAmps, spikeClusters)
%   
%     This code takes in spike times and amplitudes for a population of recorded neurons and
%      computes the full set of 3 quality metrics below as used in the IBL and Steinmetz Lab
%     (see 'computeQUalityMetrics')
%     Parameters
%     ----------
%     spikeTimes : spike times for the neurons in the population
%     spikeAmps: amplitudes (or template amplitudes) for the neurons in the population
%     spikeClusters: np.array of size (spikeTimes) with the id of the cluster corresponding to that spike
%     Returns
%     -------
%     metrics : struct which includes the outcome of quality metrics
%         (overall_pass is True if a neuron passes, False if a neuron fails)
%     


    if nargin<3
        medAmpThresh = 50
        verbose = true
        params = struct();
    end
    
    cids = np.unique(spikeClusters)

    if verbose
        fprintf(1, 'Computing metrics for %d clusters\n', numel(cids)); 
    end

    for cidx = 1:numel(cids)
        st = spikeTimes(spikeClusters==cids(cidx)); 
        sa = spikeAmps(spikeClusters == cids(cidx));
        verbose = false
        neuron_pass, rp_pass, nc_pass, ma_pass = computeQualityMetrics(st, sa, verbose)


        metrics(cidx).cid = cids(cidx); 
        metrics(cidx).overall_pass = neuron_pass;
        metrics(cidx).rp_pass = rp_pass;
        metrics(cidx).nc_pass = nc_pass;
        metrics(cidx).medAmp_pass = ma_pass;
    end
    
end
