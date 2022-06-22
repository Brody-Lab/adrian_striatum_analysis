function session_info = make_session_info(Cells,save_path)
    % cells struct should be a formatted cells struct that's in the local
    % database.
    
    session_info_fields = {'recording_name','sess_date','rat','sessid','probe_serial',...
        'days_implanted','n_clusters','craniotomy_AP','craniotomy_ML','depth_inserted',...
        'hemisphere','region_names','n_trials','violation_rate','percent_correct','nic','stim_dur_range',...
        'laser_power_mW','D2Phototagging','days_since_viral_injection','notes','mat_file_name',...
        'last_modified','ap_group'};
    
    for f=1:length(session_info_fields)
        
        if isfield(Cells,session_info_fields{f})
            session_info.(session_info_fields{f}) = Cells.(session_info_fields{f}); 
        else
            
            if session_info_fields{f}=="region_names"
                session_info.region_names  = {Cells.penetration.regions.name};                        
            
            elseif session_info_fields{f}=="n_trials"
                exclude = validate_trials(Cells.Trials,'mode','agb_glm','quiet',true);
                session_info.n_trials = numel(Cells.Trials.is_hit(~exclude));
            
            elseif session_info_fields{f}=="violation_rate"
                exclude = validate_trials(Cells.Trials,'mode','agb_glm','quiet',true,'exclude_violations',false,...
                    'exclude_no_cpoke_out_trials',false,'exclude_no_spoke_trials',false);                
                session_info.violation_rate = mean(Cells.Trials.violated(~exclude));                
            
            elseif session_info_fields{f}=="percent_correct"
                exclude = validate_trials(Cells.Trials,'mode','agb_glm','quiet',true);                
                session_info.percent_correct = 100*mean(Cells.Trials.is_hit(~exclude));
                
            elseif session_info_fields{f}=="stim_dur_range"
                exclude = validate_trials(Cells.Trials,'mode','agb_glm','quiet',true);                
                session_info.stim_dur_range = [min(Cells.Trials.stim_dur_s(~exclude)) max(Cells.Trials.stim_dur_s(~exclude))];
                
            elseif session_info_fields{f}=="nic"
                session_info.nic = median(Cells.Trials.nic);
            
            elseif session_info_fields{f}=="notes"
                if isfield(Cells,'session_notes')
                    session_info.notes = Cells.session_notes;
                end
                
            elseif session_info_fields{f}=="days_since_viral_injection"
                if session_info.rat=="A249"
                    session_info.days_since_viral_injection = days(datetime(session_info.sess_date) - datetime('2019-12-17'));
                elseif session_info.rat(1)=="A256"
                    session_info.days_since_viral_injection = days(datetime(session_info.sess_date) - datetime('2019-12-18'));        
                else
                    session_info.days_since_viral_injection = NaN;
                end                
                
            else
                session_info.(session_info_fields{f}) = Cells.penetration.(session_info_fields{f}); 
            
            end
        end
        
        if ischar(session_info.(session_info_fields{f}))
            session_info.(session_info_fields{f}) = string(session_info.(session_info_fields{f}));
        end
    end

    save(save_path,'session_info');   
    
end