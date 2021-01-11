function cell_info = make_cell_info(Cells,tag_flag)
    % makes a table from the Cells structure with no spike times but just info
    % about each cell which is useful to have alongside, for instance, GLM
    % fits
    % this standardizes the fields so that these tables can be directly
    % concatenated across sessions
    fields = {'last_modified','rat','sess_date','sessid','bank','electrode','unitCount',...
                      'unitISIRatio','unitLRatio','unitIsoDist','unitVppRaw','hemisphere',...
                      'probe_serial','distance_from_tip','DV','AP','ML','regions','ks_good',...
                      'like_axon','mean_uv','peak_trough_width','peak_uv','peak_width_s',...
                      'spike_width_ms','mean_uV','width_ms','reliability','signrank','tp',...
                      'auc','mi','dp','days_implanted','days_since_viral_injection',...
                      'distance_from_fiber_tip','first_sig_time_s','laser_power_mW','mat_file_name',...
                      'session_notes','D2_phototagging'};  
    missing_vals = num2cell(NaN(numel(fields),1));
    missing_vals{1}=NaT;
    missing_vals{3}=NaT;    
    missing_vals{2}=missing; %string
    missing_vals{12}=missing; %string    
    missing_vals{13}=missing; %string        
    missing_vals{18}={};
    missing_vals{39} = missing;
    missing_vals{40} = missing;
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
    region_names = {Cells.penetration.regions.name};
    regions = cell(n_cells,1);
    regions(Cells.regions>0) = region_names(Cells.regions(Cells.regions>0));
    Cells.regions = regions;
    if isfield(Cells,'waveform')
        waveform_fields = fieldnames(Cells.waveform);
        for i=1:length(waveform_fields)
           if ~strcmp(waveform_fields{i},'meanWfGlobalRaw') && ~strcmp(waveform_fields{i},'mean_uv')
               if strcmp(waveform_fields{i},'peak_width_s') % fix bug in peak_width_s calculation in some files
                   Cells.peak_width_s = NaN(n_cells,1);
                    scalar_idx=arrayfun(@(x)numel(x.peak_width_s),Cells.waveform)==1;
                    Cells.peak_width_s(scalar_idx) = cat(1,Cells.waveform(scalar_idx).peak_width_s);
               else
                    Cells.(waveform_fields{i}) = cat(1,Cells.waveform.(waveform_fields{i}));
               end
           end
        end
        Cells = rmfield(Cells,'waveform');
    else
        fprintf('Cells file for sessid %s does not have waveform field.',string(Cells.sessid));
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
    if cell_info.rat(1)=="A249"
        cell_info.days_since_viral_injection = days(datetime(cell_info.sess_date) - datetime('2019-12-17'));
        cell_info.distance_from_fiber_tip = sqrt((cell_info.AP - 1.6).^2 + (cell_info.DV - 3.9).^2);        
    elseif cell_info.rat(1)=="A256"
        cell_info.days_since_viral_injection = days(datetime(cell_info.sess_date) - datetime('2019-12-18'));        
        cell_info.distance_from_fiber_tip = sqrt((cell_info.AP - 1.6).^2 + (cell_info.DV - 3.9).^2);
    end
    tag_fields = {'reliability','signrank','tp','auc','mi','dp'};
    if tag_flag
        cell_info.first_sig_time_s = cellfun(@(x)x.bin_pval2(end),Cells.PPTH.first_sig_time_s); % p=0.05 using the pre-laser times as the null distribution
        % add all statistics from the "prepost" field of
        % Cells.PPTH.combined_pulse_stats (i.e. comparing the 100ms before
        % and after the laser across all pulses)
        for f=1:length(tag_fields)  
            cell_info.(tag_fields{f}) = cellfun(@(x)x.prepost.(tag_fields{f}),Cells.PPTH.combined_pulse_stats);
        end
    end
end

