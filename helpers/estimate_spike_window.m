function cells_window = estimate_spike_window(Cells,ref_event)
    n = numel(Cells.raw_spike_time_s);
    cells_window = [Inf -Inf];
    for c=1:n
        c_spikes = cat(1,Cells.spike_time_s.(ref_event){c}{:});
        if ~isempty(c_spikes)
            cells_window(1) = min(min(c_spikes),cells_window(1));
            cells_window(2) = max(max(c_spikes),cells_window(2));     
        end
    end
end