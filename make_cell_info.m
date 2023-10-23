function cell_info = make_cell_info(Cells,save_path)
    % makes a table from the Cells structure with no spike times but just info
    % about each cell which is useful to have alongside, for instance, GLM
    % fits
    % this standardizes the fields so that these tables can be directly
    % concatenated across sessions
    fields = {'recording_name','bank','electrode',...
              'frac_isi_violation','l_ratio','isolation_distance','uVpp','hemisphere',...
              'distance_from_tip','DV','AP','ML','region','ks_good',...
              'like_axon','spatial_spread_um','peak_trough_width_s','peak_width_s','peak_uv'...
              'reliability','signrank','tp',...
              'auc','mi','dp','distance_from_fiber_tip','first_sig_time_s','autocorr_baseline','autocorr_baseline_fr_hz',...
              'autocorr_clicks','autocorr_clicks_fr_hz','presence_ratio','stability','mean_rate_hz','spikes_per_s','num_spikes',...
              'is_in_dorsal_striatum','ap_group'};  
          
% mean_uV and mean_uv are both fields??? these aren't being made properly
% across all Cells generations in format_cells_fiels
          
    missing_vals = num2cell(NaN(numel(fields),1));
    missing_vals{1}=missing; %string recording_name
    missing_vals{8}=missing; %string    hmeisphere
    missing_vals{13}=missing; %string  region      
    missing_vals{28}={[]}; %vector autocorr_baseline
    missing_vals{30}={[]}; %vector autocorr_clicks
    if isfield(Cells,'quality_metrics')
        quality_fields = fieldnames(Cells.quality_metrics);
        for q=1:numel(quality_fields)
            Cells.(quality_fields{q}) = Cells.quality_metrics.(quality_fields{q});
        end
    end
    % keep adding special missings
    cell_info=table();    
    n_cells = numel(Cells.raw_spike_time_s);    
    for f=1:length(fields)
        if iscell(missing_vals{f})
            cell_info.(fields{f}) = cell(n_cells,1);
        else
            cell_info.(fields{f}) = repmat(missing_vals{f},n_cells,1);            
        end
    end
    for f=1:length(fields)
        if isfield(Cells,fields{f})
            if isscalar(Cells.(fields{f}))
                tmp=repmat(Cells.(fields{f}),n_cells,1);
            elseif ischar(Cells.(fields{f}))
                tmp = repmat(string(Cells.(fields{f})(1,:)),n_cells,1);
            elseif isrow(Cells.(fields{f}))
                tmp=Cells.(fields{f})';                
            else
                tmp=Cells.(fields{f});                                
            end
            if ~isvector(tmp)
                warning('Entry %s is not a vector. Cannot insert into table.',fields{f});    
            elseif ~all(size(tmp) == [n_cells,1])
                warning('Size of entry %s (%g,%g) does not match expected size (%g,%g). Cannot insert into table.',fields{f},size(tmp,1),size(tmp,2),n_cells,1);
            else
                cell_info.(fields{f}) = tmp;
            end
        else
           %warning('Fields %s not found.',fields{f}); 
           % happens frequently, warning messages are not useful
        end
    end
    if strncmp(Cells.recording_name,"A249",4)
        cell_info.distance_from_fiber_tip = sqrt((cell_info.AP - 1.6).^2 + (cell_info.DV - 3.9).^2);        
    elseif strncmp(Cells.recording_name,"A256",4)
        cell_info.distance_from_fiber_tip = sqrt((cell_info.AP - 1.6).^2 + (cell_info.DV - 3.9).^2);
    end
    tag_fields = {'reliability','signrank','tp','auc','mi','dp'};
    if Cells.D2Phototagging==1
        cell_info.first_sig_time_s = cellfun(@(x)x.bin_pval2(end),Cells.PPTH.first_sig_time_s); % p=0.05 using the pre-laser times as the null distribution
        % add all statistics from the "prepost" field of
        % Cells.PPTH.combined_pulse_stats (i.e. comparing the 100ms before
        % and after the laser across all pulses)
        for f=1:length(tag_fields)  
            cell_info.(tag_fields{f}) = cellfun(@(x)x.prepost.(tag_fields{f}),Cells.PPTH.combined_pulse_stats);
        end
    end
    
    % save
    save(save_path,'cell_info');
end

