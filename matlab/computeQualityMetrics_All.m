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


    if nargin<4
        medAmpThresh = 50;
        verbose = true;
        params = struct();
    end
    
    cids = unique(spikeClusters);

    if verbose
        fprintf(1, 'Computing metrics for %d clusters\n', numel(cids)); 
        reverseStr = '';
    end

    for cidx = 1:numel(cids)
        st = spikeTimes(spikeClusters==cids(cidx)); 
        sa = spikeAmps(spikeClusters == cids(cidx));
        verboseSingleNeuron = false;
        [neuron_pass, rp_pass, nc_pass, ma_pass] = computeQualityMetrics(st, sa, 50, verboseSingleNeuron);


        %if verbose, display progress bar:
        if verbose
            percentDone = 100 * cidx / numel(cids);
            msg = sprintf('Percent done: %3.1f', percentDone); %Don't forget this semicolon
            fprintf([reverseStr, msg]);
            reverseStr = repmat(sprintf('\b'), 1, length(msg));
        end

        
        metrics(cidx).cid = cids(cidx); 
        metrics(cidx).overall_pass = neuron_pass;
        metrics(cidx).rp_pass = rp_pass;
        metrics(cidx).nc_pass = nc_pass;
        metrics(cidx).medAmp_pass = ma_pass;
    end
    
    if verbose
        fprintf('Number of good clusters is %d out of %d', sum([metrics.overall_pass]), numel(cids))
    end
    
end