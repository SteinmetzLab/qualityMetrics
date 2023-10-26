function [noiseCutoffStruct] = run_noiseCutoff_All(spikeAmps, spikeClusters,verbose,params)
%Inputs:
% spikeAmps : spike amplitudes for the neurons in the population 
%           (in uV or rescaled template amplitudes)
% spikeClusters: np.array of size (spikeAmps) with the id of the
%                cluster corresponding to that spike
% verbose : default is true; set to false to suppress printed output

%Output:
% noiseCutoffStruct which includes:
% - clusterID
% - nc_pass (whether noise cutoff passes)
% - nc_value (test statistic)

if nargin<4
    verbose = true;
    params = struct();
end

cids = unique(spikeClusters);

if verbose
    fprintf(1, 'Computing noise cutoff metric for %d clusters\n', numel(cids)); 
    reverseStr = '';
end

noiseCutoffStruct(numel(cids))=struct();
for cidx = 1:numel(cids)

    %if verbose, display progress bar:
    if verbose
        percentDone = 100 * cidx / numel(cids);
        msg = sprintf('Percent done: %3.1f', percentDone); %Don't forget this semicolon
        fprintf([reverseStr, msg]);
        reverseStr = repmat(sprintf('\b'), 1, length(msg));
    end

    sa = spikeAmps(spikeClusters==cids(cidx)); 

    [nc_pass, cutoff_value, ~] = noise_cutoff(sa);

    noiseCutoffStruct(cidx).cid = cids(cidx); 
    noiseCutoffStruct(cidx).nc_pass = nc_pass;
    noiseCutoffStruct(cidx).cutoff_value = cutoff_value;
end
