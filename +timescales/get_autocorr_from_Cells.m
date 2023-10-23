function [autocorr,edges,fr_hz,corr_mtx] = get_autocorr_from_Cells(Cells,eve,bin_lims_s,bin_size_s,params)
    arguments
        Cells (1,1) struct
        eve (1,:) char {mustBeNonempty}
        bin_lims_s (1,2) {isnumeric,mustBeFinite}            
        bin_size_s (1,1) {isnumeric,mustBePositive,mustBeFinite}
        params.exclude_trials (:,1) {islogical} = false(numel(Cells.Trials.gamma),1); 
        params.mask_eve (1,:) char = ''
        params.mask_window_s (1,2) {isnumeric} = [0 0]        
    end
    [corr_mtx,edges,fr_hz] = get_corr_mtx_from_Cells(Cells,eve,bin_lims_s,bin_size_s,params);    
    num_clusters = numel(Cells.raw_spike_time_s);
    for c=num_clusters:-1:1
        autocorr{c} = timescales.get_autocorr_from_corr_mtx(corr_mtx{c});
        if c==1
            edges = (1:numel(autocorr{c})) * bin_size_s;        
        end
    end 
end



function [corr_mtx,edges,fr_hz] = get_corr_mtx_from_Cells(Cells,eve,bin_lims_s,bin_size_s,params)

   
    % set up binning
    n_bins = round(diff(bin_lims_s)./bin_size_s);
    bin_lims_s_used = [1 1]*bin_lims_s(1) + [ 0 (n_bins)*bin_size_s];    
    if any(abs(bin_lims_s_used-bin_lims_s)>1e-3)
        warning('Due to rounding error, using a time period of [%g to %g].',bin_lims_s_used(1),bin_lims_s_used(2));
    end   
    edges=linspace(bin_lims_s_used(1),bin_lims_s_used(2),n_bins+1);
    
    
    % do masking if requested
    out_of_range = false(sum(~params.exclude_trials),n_bins);
    if ~isempty(params.mask_eve)
        if ~isfield(Cells.Trials.stateTimes,params.mask_eve)
            error('No field %s in Cells.Trials.stateTimes.',params.mask_eve);
        end
        if ~isfield(Cells.Trials.stateTimes,eve)
            error('No field %s in Cells.Trials.stateTimes.',eve);
        end        
        mask_times = Cells.Trials.stateTimes.(params.mask_eve)(~params.exclude_trials) - Cells.Trials.stateTimes.(eve)(~params.exclude_trials);
        for t=1:numel(mask_times)
            out_of_range(t,:) = edges(1:end-1) < mask_times(t)+params.mask_window_s(1) | edges(2:end) > mask_times(t)+params.mask_window_s(2);
        end
    end
    
    % get aligned spike counts if not there already
    num_clusters = numel(Cells.raw_spike_time_s);
    if ~isfield(Cells,'spike_time_s') || ~isfield(Cells.spike_time_s,eve)
        warning('Aligned spike times not computed. Trying to do that now.');
        if ~isfield(Cells,'raw_spike_time_s')
            error('No field "raw_spike_time_s" in Cells.');
        end
        if ~isfield(Cells.Trials.stateTimes,eve)
            error('No field %s in Cells.Trials.stateTimes.',eve);
        end
        window_s = bin_lims_s + [-1 1]*bin_size_s;
        Cells.spike_time_s.(eve) = group_spike_times(Cells.raw_spike_time_s, Cells.Trials.stateTimes.(eve)(~params.exclude_trials), window_s);
    else
        for c=1:num_clusters        
            Cells.spike_time_s.(eve){c} = Cells.spike_time_s.(eve){c}(~params.exclude_trials);
        end
    end

    % do binning and make correlation matrices
    countfun = @(x)matlab.internal.math.histcounts(x,edges);
    for c=num_clusters:-1:1
        counts = cellfun(countfun,Cells.spike_time_s.(eve){c},'uni',0);
        counts=cat(1,counts{:});
        counts(out_of_range) = NaN;
        fr_hz(c,1) = nanmean(counts(:));
        corr_mtx{c,1} = corr(counts,'rows','pairwise'); % find non-NaN elements in a pairwise fashion        
    end
    fr_hz = fr_hz./bin_size_s;
end