function [stability,presence,mean_rate] = calculate_unit_stability(Cells,params)
    arguments
        Cells (1,1) struct;
        params.ref_event (1,1) string {ref_events_member(params.ref_event,Cells)} = 'cpoke_in';
        params.chunk_size (1,1) {isnumeric,mustBePositive} = 100; % number of trials to average together to locally estimate firing rate
        params.n_repeats (1,1) {isnumeric,mustBePositive} = 1000; % number of times to sample a chunk to obtain a firing rate distribution
    end  
    ncells = length(Cells.spike_time_s.(params.ref_event));
    for k=ncells:-1:1
        matrix(:,k) = cellfun(@numel,Cells.spike_time_s.(params.ref_event){k}); % make a matrix of number of spikes per trial for each cell [ntrials x ncells]
    end
    ntrials = size(matrix,1);
    if ntrials/params.chunk_size<=1
        error('calculate_unit_stability: number of trials (%g) is less than chunk_size (%g). Cannot proceed.',ntrials,params.chunk_size);
    elseif ntrials/params.chunk_size<1.5
        warning('calculate_unit_stability: only %g trials means a small %1.2g:1 of trials to chunk_size. stability calculation may be inaccurate.',ntrials,ntrials/params.chunk_size);
    end
    for i = params.n_repeats:-1:1
        idx = randperm(ntrials-params.chunk_size+1,1);
        idx_range = idx:idx+params.chunk_size-1;
        idx_shuffled = randperm(ntrials,params.chunk_size);
        samples(i,:) = mean(matrix(idx_range,:)); % distribution of firing rates in time bins of params.chunk_size trials
        samples_shuffled(i,:) = mean(matrix(idx_shuffled,:)); % time-shuffled distribution for comparison
    end
    stability = std(samples)./std(samples_shuffled); % final metric is the ratio of sd's of this distribution
    %N.B: cutoff of around 4 or 5 seems reasonable for a chunk size of 100 trials
    if isfield(Cells,'kSpikeWindowS')
        kspikewindow_range = diff(Cells.kSpikeWindowS.(params.ref_event));
    else
        spike_window = estimate_spike_window(Cells,params.ref_event);
        kspikewindow_range = diff(spike_window);
    end
    mean_rate = mean(matrix) ./ kspikewindow_range;
    presence = mean(matrix>0);
end

function ismem = ref_events_member(ref_event,Cells)
    ismem = ismember(ref_event,fieldnames(Cells.spike_time_s));
end