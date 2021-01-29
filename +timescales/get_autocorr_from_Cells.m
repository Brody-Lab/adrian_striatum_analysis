function [autocorr,edges,fr_hz] = get_autocorr_from_Cells(Cells,eve,bin_lims_s,params)
    arguments
        Cells (1,1) struct
        eve (1,:) char {mustBeNonempty}
        bin_lims_s (1,2) {isnumeric,mustBeFinite}            
        params.bin_size_s (1,1) {isnumeric,mustBePositive,mustBeFinite} = 0.05; 
        params.exclude_trials (:,1) {islogical} = zeros(0,1); % 33ms, making 15 bins over the standard interval [0 0.5]                
    end
    [corr_mtx,edges,fr_hz] = timescales.get_corr_mtx_from_Cells(Cells,eve,bin_lims_s,'bin_size_s',params.bin_size_s,'exclude_trials',params.exclude_trials);    
    num_clusters = numel(Cells.spike_time_s.cpoke_in);
    for c=num_clusters:-1:1
        autocorr{c} = timescales.get_autocorr_from_corr_mtx(corr_mtx{c});
        if c==1
            edges = (1:numel(autocorr{c})) * params.bin_size_s;        
        end
    end 
end