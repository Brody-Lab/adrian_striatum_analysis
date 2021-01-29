function [corr_mtx,edges,fr_hz] = get_corr_mtx_from_Cells(Cells,eve,bin_lims_s,params)
    arguments
        Cells (1,1) struct
        eve (1,:) char {mustBeNonempty}
        bin_lims_s (1,2) {isnumeric,mustBeFinite}            
        params.bin_size_s (1,1) {isnumeric,mustBePositive,mustBeFinite} = 0.05; 
        params.exclude_trials (:,1) {islogical} = false(numel(Cells.Trials.gamma),1); 
    end
    num_clusters = numel(Cells.spike_time_s.cpoke_in);
    if ~isfield(Cells.spike_time_s,eve)
        warning('Aligned spike times not computed. Trying to do that now.');
        if ~isfield(Cells,'raw_spike_time_s')
            error('');
        end
        if ~isfield(Cells.Trials.stateTimes,eve)
            error('');
        end
        window_s = bin_lims_s + [-1 1]*params.bin_size_s;
        for c=1:num_clusters
            Cells.spike_time_s{c} = group_spike_times(Cells.raw_spike_time_s{c}, Cells.Trials.stateTimes.(eve), window_s);
        end
    end 
    n_bins = round(diff(bin_lims_s)./params.bin_size_s);
    bin_lims_s_used = [1 1]*bin_lims_s(1) + [ 0 (n_bins)*params.bin_size_s];
    if any(abs(bin_lims_s_used-bin_lims_s)>1e-3)
        warning('Due to rounding error, using a time period of [%g to %g].',bin_lims_s_used(1),bin_lims_s_used(2));
    end
    edges=linspace(bin_lims_s_used(1),bin_lims_s_used(2),n_bins+1);
    countfun = @(x)matlab.internal.math.histcounts(x,edges);
    for c=num_clusters:-1:1
        counts = cellfun(countfun,Cells.spike_time_s.(eve){c}(~params.exclude_trials),'uni',0);
        counts=cat(1,counts{:});
        fr_hz(c,1) = nanmean(counts(:));
        corr_mtx{c,1} = corr_columns_fast(counts);
    end
    fr_hz = fr_hz./params.bin_size_s;
end