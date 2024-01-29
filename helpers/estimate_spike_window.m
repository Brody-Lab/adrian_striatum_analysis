function cells_window = estimate_spike_window(Cells,ref_event)
    n = numel(Cells.spike_time_s.(ref_event));
    cells_window = [Inf -Inf];
    for c=1:n
        spks = Cells.spike_time_s.(ref_event){c};
        if ~isempty(spks)
            c_spikes = cat(1,spks{:});
        else
            c_spikes=[];
        end
        if ~isempty(c_spikes)
            cells_window(1) = min(min(c_spikes),cells_window(1));
            cells_window(2) = max(max(c_spikes),cells_window(2));     
        end
    end
end