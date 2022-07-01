function [sessions_table,cells_table] = create_database(varargin)
    % moves files from bucket into the repo database, and creates index
    % files for easy searching of what's there
    
    %% parse and validate inputs
    P=get_parameters;
    p=inputParser;
    p.addParameter('update',false,@(x)validateattributes(x,{'logical'},{'scalar'})); % add only new sessions?
    p.addParameter('reformat',true);
    p.addParameter('use_local',false);
    p.addParameter('archive',false);
    p.parse(varargin{:});
    params = p.Results;
    
    %% extract relevants rows from recordings log
    recordings_table = get_striatum_glm_recordings_table();    
    n_sessions = height(recordings_table);
    
    %% establish and validate paths
    paths = get_data_paths();
    
    %% loop over sessions to generate local formatted cells files if needed and then make cell info and session info structures
    recording_name = recordings_table.recording_name;
    for i=1:n_sessions
        fprintf('\n\nWorking on session %g of %g: %s\n-----------------\n',i,n_sessions,recording_name(i));
        if paths(i).all_exist && params.update
            fprintf('database does not need to be updated for %s.\n-----------------\n',recording_name(i));                         
        else
            
            % load cells file
            Cells = load_Cells_file(recording_name(i),params.use_local);
            
            if ~params.use_local || params.reformat
                % add extra fields to Cells, correct for database consistency and saves it locally
                Cells = format_Cells_file(Cells,recordings_table(i,:),paths(i).cells_file); 
            end
            
            % make and save reduced cell info table entries that do not include spike times
            make_cell_info(Cells,paths(i).cell_info); 
            
            % make and save very reduced session info table entry with basic information about the recording session
            make_session_info(Cells,paths(i).session_info);  
                     
        end      
    end
    
    %% move unneeded directories to archive folder
    if params.archive
        tmp = dir(fullfile(P.data_path,'cells'));
        folder_names = {tmp(3:end).name};
        to_archive = find(~ismember(folder_names,recordings_table.recording_name));
        for i=1:length(to_archive)
            source = fullfile(P.data_path,'cells',folder_names{to_archive(i)});
            destination = fullfile(P.data_path,'archive','cells',folder_names{to_archive(i)});        
            fprintf('Moving %s to archive location at %s.\n',source,destination);
            movefile(source,destination);
        end
    end
    
    %% make sessions_table and cells_table
    sessions_table = make_sessions_table();
    cells_table = make_cells_table();
    
end

function Cells = format_Cells_file(Cells,recordings_table,save_path)    
    
    P = get_parameters();

    %% remove uninitiated trials
    remove_idx = isnan(Cells.Trials.stateTimes.cpoke_in);
    if any(remove_idx)
        fprintf('Removing %g uninitiated trials.\n',sum(remove_idx));
        Cells.Trials = remove_trials(Cells.Trials,remove_idx);
        fields = fieldnames(Cells.spike_time_s);
        for f=1:length(fields)
            for c = 1:length(Cells.spike_time_s.(fields{f}))
                Cells.spike_time_s.(fields{f}){c} = Cells.spike_time_s.(fields{f}){c}(~remove_idx);
            end
        end
    end

    %% copy some fields from recordings_table
    fields_to_copy = {'laser_power_mW','notes','D2Phototagging','recording_name','cells_file'};
    new_name = {'laser_power_mW','notes','D2Phototagging','recording_name','mat_file_name'};            
    for i=1:length(fields_to_copy)
        Cells.(new_name{i}) = recordings_table.(fields_to_copy{i});
    end

    %% ensure all cells files are up to date, bug free and consistent in structure
    if ~isfield(Cells,'ap_meta') && isfield(Cells,'meta')
        Cells.ap_meta = Cells.meta.ap_meta;
    end    
    Cells = import_penetration(Cells);    
    Cells.n_clusters = numel(Cells.spike_time_s.cpoke_in);                
    if ~isfield(Cells,'ks_good')
        % calculate refractory period violations a la Kilosort
        for k=1:Cells.n_clusters
            Cells.ks_good(k) = is_ks_good(Cells.raw_spike_time_s{k});
        end
    end    
    try
        Cells.probe_serial = Cells.ap_meta.imDatPrb_sn;
    catch
        Cells.probe_serial = Cells.rec.ap_meta.imDatPrb_sn;        
    end
    time_to_clicks = Cells.Trials.stateTimes.clicks_on - Cells.Trials.stateTimes.cpoke_in;
    exclude_trials = Cells.Trials.violated | time_to_clicks<0.5 | Cells.Trials.laser.isOn;
    [Cells.autocorr,~,Cells.autocorr_fr_hz] = timescales.get_autocorr_from_Cells(Cells,'cpoke_in',[0 0.5],'exclude_trials',exclude_trials);                  
    if isfield(Cells,'sess_date')
        try
            Cells.sess_date = datetime(Cells.sess_date);
        catch
            Cells.sess_date = datetime(Cells.sess_date,'InputFormat','uuuu_MM_dd');
        end
    elseif unique(Cells.sessid)==701531 % special case for A242 session
        Cells.sess_date = datetime('2019-06-03');
    elseif unique(Cells.sessid)==702016 % special case for A242 session
        Cells.sess_date = datetime('2019-06-06');      
    elseif unique(Cells.sessid)==703121 % special case for A242 session     
        Cells.sess_date = datetime('2019-06-10');              
    end
    Cells.rat = Cells.rat(1,:);      
    Cells.sessid = Cells.sessid(1);
    Cells.hemisphere = string(Cells.penetration.hemisphere);                              
    Cells.sess_date = Cells.sess_date(1,:);      
    Cells.days_implanted = days(datetime(Cells.sess_date) - datetime(Cells.penetration.date_implanted)); 
    if ~isfield(Cells,'last_modified') && isfield(Cells,'meta')
        Cells.last_modified = Cells.meta.last_modified;
    end
    Cells.last_modified = datetime(Cells.last_modified);
    Cells.probe_serial = string(Cells.probe_serial);
    if isfield(Cells,'PETH')
        Cells = rmfield(Cells,'PETH');
    end
    
    for i=1:numel(P.ap_groups)
        mean_ap_striatal = mean(Cells.AP(Cells.is_in_dorsal_striatum));
        if mean_ap_striatal>P.ap_groups{i}(1) && mean_ap_striatal<P.ap_groups{i}(2)
            Cells.ap_group=i;
            break
        end
    end
    
    % add stateTimes for clicks if not there in Thomas' sessions
    if Cells.rat=="T219" && ~isfield(Cells.Trials.stateTimes,'left_clicks')
        fprintf('Adding click states to this T219 file ...');tic;
        Cells.Trials = add_click_states_to_thomas_files(Cells.Trials);
        fprintf(' took %s.\n-----------------\n',timestr(toc));                   
    end
    
    % computer laser modulation statistics for phototagging sessions
    if Cells.D2Phototagging==1
        fprintf('-----------------');                                                   
        Cells = tagging.compute_laser_modulation(Cells);
        fprintf('-----------------\n');                                   
    end  
    
    % add stability metrics
    [Cells.stability,Cells.presence_ratio,Cells.mean_rate_hz] = calculate_unit_stability(Cells);
    
    % recalculate mean waveforms
    if isfield(Cells,'waveform')
        try
            Cells.waveform.mean_uv = Cells.waveform.mean_uV;
        end
        if isfield(Cells.waveform,'mean_uv')
            ncells = numel(Cells.raw_spike_time_s);
            waveform=Cells.waveform;
            Cells = rmfield(Cells,'waveform');
            if isfield(Cells,'rec')
                samp_rate = Cells.rec.ap_meta.imSampRate;
            else
                try
                    samp_rate = Cells.ap_meta.imSampRate;
                catch
                    samp_rate = 3e4;
                end
            end
            if isscalar(waveform) 
                if size(waveform.mean_uv,1)==ncells                    
                    for c=1:ncells
                        Cells.waveform(c) = calculate_waveform_stats(waveform.mean_uv(c,:),samp_rate) ;   
                    end
                end
            elseif numel(waveform)==ncells
                for c=1:ncells
                    Cells.waveform(c) = calculate_waveform_stats(waveform(c).mean_uv,samp_rate) ;   
                end                
            end
        end
        if ~isempty(Cells.waveform)
            waveform_fields = fieldnames(Cells.waveform);
            for i=1:length(waveform_fields)
               if ~strcmp(waveform_fields{i},'meanWfGlobalRaw') && ~strcmp(waveform_fields{i},'mean_uv')
                   if strcmp(waveform_fields{i},'peak_width_s') % fix bug in peak_width_s calculation in some files
                       Cells.peak_width_s = NaN(Cells.n_clusters,1);
                        scalar_idx=arrayfun(@(x)numel(x.peak_width_s),Cells.waveform)==1;
                        Cells.peak_width_s(scalar_idx) = cat(1,Cells.waveform(scalar_idx).peak_width_s);
                   else
                        Cells.(waveform_fields{i}) = cat(1,Cells.waveform.(waveform_fields{i}));
                   end
               end
            end
            Cells.waveform=[];
        end
    else
        fprintf(' - warning: missing mean waveform data - ');
    end    
    
    
    %% keep a minimal set of spike_time_s fields
    keepfields = {'cpoke_in','iti','cpoke_out','clicks_on'}; % iti needed for computing laser modulation during optotagging sessions
    fields = fieldnames(Cells.spike_time_s);
    for f=1:length(fields)
        if ~ismember(fields{f},keepfields)
            Cells.spike_time_s = rmfield(Cells.spike_time_s,fields{f});
        end
    end    
    
    %% add field specifying where a cell is in dorsal striatum
    Cells.is_in_dorsal_striatum = get_dorsal_striatal_cells(Cells);
    
    %% save cells files 
    fprintf('saving local cells file  ...');tic;      
    lastwarn('');
    try
        save(save_path, '-struct','Cells','-v7','-nocompression'); % these settings create a file that is smaller and MUCH faster to save and load than with the default settings
    catch
        save(save_path, '-struct','Cells','-v7'); % these settings create a file that is smaller and MUCH faster to save and load than with the default settings
    end
    if ~isempty(lastwarn) % warning indicates a variable wasn't saved so try another way
        fprintf('Warning encountered during v7 save. Trying v7.3\n');
        save(save_path, '-struct','Cells','-v7.3');
    end
    fprintf(' took %s.\n-----------------\n',timestr(toc));      
end