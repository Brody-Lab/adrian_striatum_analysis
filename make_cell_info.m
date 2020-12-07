function cell_info = make_cell_info(Cells,tag_flag)
    % makes a table from the Cells structure with no spike times but just info
    % about each cell which is useful to have alongside, for instance, GLM
    % fits
    % this standardizes the fields so that these tables can be directly
    % concatenated across sessions
    fields = {'last_modified','rat','sess_date','sessid','bank','electrode','unitCount',...
                      'unitISIRatio','unitLRatio','unitIsoDist','unitVppRaw','hemisphere',...
                      'probe_serial','distance_from_tip','DV','AP','ML','regions','ks_good'};    
    if ~isfield(Cells,'penetration')
        Cells = import_penetration(Cells);
    end
    num_clusters = numel(Cells.spike_time_s.cpoke_in);
    if ~isfield(Cells,'ks_good')
        % calculate refractory period violations a la Kilosort
        for i=1:num_clusters
            Cells.ks_good(i) = is_ks_good(Cells.raw_spike_time_s{i});
        end
    end
    if tag_flag
        %Cells = compute_laser_modulation(Cells);
        % to do: extract relevant fields
    end
    n_cells = numel(Cells.raw_spike_time_s);    
    region_names = {Cells.penetration.regions.name};
    regions = cell(n_cells,1);
    regions(Cells.regions>0) = region_names(Cells.regions(Cells.regions>0));
    Cells.regions = regions;
    cell_info=table();
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
        fields = union(fields,waveform_fields);        
    end
    try
        Cells.probe_serial = Cells.ap_meta.imDatPrb_sn;
    catch
        Cells.probe_serial = Cells.rec.ap_meta.imDatPrb_sn;        
    end
    Cells.probe_serial = string(Cells.probe_serial);    
    for f=1:length(fields)
        tmp = repmat(missing,n_cells,1);                        
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
                tmp = repmat(missing,n_cells,1);                                        
            elseif ~all(size(tmp) == [n_cells,1])
                warning('Size of entry %s (%g,%g) does not match expected size (%g,%g). Cannot insert into table.',fields{f},size(tmp,1),size(tmp,2),n_cells,1);
                tmp = repmat(missing,n_cells,1);                                        
            end
        else
           warning('Fields %s not found.',fields{f}); 
        end
        cell_info.(fields{f}) = tmp;        
    end
end

