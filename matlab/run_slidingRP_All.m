function [slidingRPStruct] = run_slidingRP_All(spikeTimes, spikeClusters,verbose,params)
%Inputs:
% spikeTimes : spike times for the neurons in the population (in s)
% spikeClusters: np.array of size (spikeTimes) with the id of the 
%                cluster corresponding to that spike
% verbose : default is true; set to false to suppress printed output

%Output:
% slidingRPStruct which includes:
% - clusterID
% - rp_pass (whether RP test passes)
% - rp_minContamination (at 90% confidence)
% - rp_maxConfidence (at 10% contamination)
% - rp_rpTime (time of min Contamination)

if nargin<3
    verbose = true;
end
if nargin<4
    params = struct();
end

cids = unique(spikeClusters);

if verbose
    fprintf(1, 'Computing slidingRP metric for %d clusters\n', numel(cids)); 
    reverseStr = '';
end

slidingRPStruct(numel(cids))=struct();
for cidx = 1:numel(cids)

    %if verbose, display progress bar:
    if verbose
        percentDone = 100 * cidx / numel(cids);
        msg = sprintf('Percent done: %3.1f', percentDone); %Don't forget this semicolon
        fprintf([reverseStr, msg]);
        reverseStr = repmat(sprintf('\b'), 1, length(msg));
    end

    st = spikeTimes(spikeClusters==cids(cidx)); 

    [maxConfidenceAt10Cont, minContWith90Confidence, timeOfLowestCont,...
        nSpikesBelow2, confMatrix, cont, rp, nACG, firingRate] ...
        = slidingRP(st, params);

    % Also compute whether RP pass:
    if minContWith90Confidence < 10 %params('contaminationThresh')
        rp_pass = true;
    else
        rp_pass = false;
    end

    slidingRPStruct(cidx).clusterID = cids(cidx); 
    slidingRPStruct(cidx).rp_pass = rp_pass;
    slidingRPStruct(cidx).rp_minContamination = minContWith90Confidence;
    slidingRPStruct(cidx).rp_maxConfidence = maxConfidenceAt10Cont;
    slidingRPStruct(cidx).rp_rpTime = timeOfLowestCont;

end

