function [stability,presence,nspks] = calculate_unit_stability(Cells,params)
    arguments
        Cells (1,1) struct;
        params.ref_event (1,1) string {ref_events_member(params.ref_event,Cells)} = 'cpoke_in';
        params.chunk_size (1,1) {isnumeric,mustBePositive} = 100; % number of trials to average together to locally estimate firing rate
        params.n_repeats (1,1) {isnumeric,mustBePositive} = 1000; % number of times to sample a chunk to obtain a firing rate distribution
    end  
    ncells = length(Cells.spike_time_s.(params.ref_event));
    for k=1:ncells
        matrix(:,k) = cellfun(@numel,Cells.spike_time_s.(params.ref_event){k}); % make a matrix of number of spikes per trial for each cell [ntrials x ncells]
    end
    ntrials = size(matrix,1);
    for i = 1:params.n_repeats
        idx = randperm(ntrials-params.chunk_size+1,1);
        idx_range = idx:idx+params.chunk_size-1;
        idx_shuffled = randperm(ntrials,params.chunk_size);
        samples(i,:) = mean(matrix(idx_range,:)); % distribution of firing rates in time bins of params.chunk_siz trials
        samples_shuffled(i,:) = mean(matrix(idx_shuffled,:)); % time-shuffled distribution for comparison
    end
    stability = std(samples)./std(samples_shuffled); % final metric is the ratio of sd's of this distribution
    %N.B: cutoff of around 4 or 5 seems reasonable for a chunk size of 100 trials
    nspks = mean(matrix);
    presence = mean(matrix>0);
end

function ismem = ref_events_member(ref_event,Cells)
    ismem = ismember(ref_event,fieldnames(Cells.spike_time_s));
end