function session_info = make_session_info(Cells,save_path)
    % cells struct should be a formatted cells struct that's in the local
    % database.
    
    session_info_fields = {'recording_name','sess_date','rat','sessid','probe_serial',...
        'days_implanted','n_clusters','craniotomy_AP','craniotomy_ML','depth_inserted',...
        'hemisphere','region_names','n_trials','violation_rate','percent_correct',...
        'laser_power_mW','D2Phototagging','notes','n_completed_acc_trials_no_stim',...
        'n_completed_acc_trials_no_stim_percent_correct','mat_file_name'};
    
    for f=1:length(session_info_fields)
        
        if isfield(Cells,session_info_fields{f})
            session_info.(session_info_fields{f}) = Cells.(session_info_fields{f}); 
        else
            
            if session_info_fields{f}=="region_names"
                session_info.region_names  = {Cells.penetration.regions.name};                        
            
            elseif session_info_fields{f}=="n_trials"
                session_info.n_trials = numel(Cells.Trials.is_hit);
            
            elseif session_info_fields{f}=="violation_rate"
                session_info.violation_rate = mean(Cells.Trials.violated);                
            
            elseif session_info_fields{f}=="percent_correct"
                session_info.percent_correct = 100*mean(Cells.Trials.is_hit(~Cells.Trials.violated));
            
            elseif session_info_fields{f}=="n_completed_acc_trials_no_stim" 
                % accumulation trials that were completed and had no laser stim (i.e. the ones useful for glm fitting)                                
                session_info.n_completed_acc_trials_no_stim = sum(~Cells.Trials.violated & char(Cells.Trials.trial_type)=='a' & ~Cells.Trials.laser.isOn);  
            
            elseif session_info_fields{f}=="n_completed_acc_trials_no_stim_percent_correct" 
                % accumulation trials that were completed and had no laser stim (i.e. the ones useful for glm fitting)                                
                session_info.n_completed_acc_trials_no_stim_percent_correct = mean(Cells.Trials.is_hit(~Cells.Trials.violated & char(Cells.Trials.trial_type)=='a' & ~Cells.Trials.laser.isOn))*100;                         
            
            elseif session_info_fields{f}=="notes"
                if isfield(Cells,'session_notes')
                    session_info.notes = Cells.session_notes;
                end
                
            else
                session_info.(session_info_fields{f}) = Cells.penetration.(session_info_fields{f}); 
            
            end
        end
    end
    
    save(save_path,'session_info');   
    
end