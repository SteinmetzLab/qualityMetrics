# qualityMetrics - MATLAB
Code to compute single unit quality metrics on neuropixels data

See script_testQualityMetrics for an example in loading data and computing metrics (you will need to change the filepaths for your github directory and data.)



To compute quality metrics for one neuron, you will need a vector of spike times (spikeTimes), and spike amplitudes (spikeAmps):
[neuron_pass, rp_pass, nc_pass, ma_pass] = computeQualityMetrics(spikeTimes, spikeAmps);

neuron_pass = 1 if your neuron passes all 3 quality metrics, 0 otherwise

rp_pass is 1 if your neuron passes the sliding refractory period metric, 0 otherwise
nc_pass is 1 if your neuron passes the noise cutoff metric, 0 otherwise
ma_pass is 1 if your neuron passes the median amplitude threshold, 0 otherwise



To plot quality plots for an individual neuron: 
plotQualityMetrics(spikeTimes,spikeAmps);



To compute quality metrics for a session (population of neurons), you will need a vector of all spike times (allst), 
all spike amplitudes (allamps), and a vector of the same size that indicates cluster identity. 

(Note: this can take a few minutes to run, by default a progress bar will show)

metrics = computeQualityMetrics_All(allst,allamps,spikeClusters);

This function returns metrics, a struct which contains pass/fail information for each neuron


