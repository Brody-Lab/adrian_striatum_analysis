function is_saved = is_glmfit_saved_cell(recording_name,cellno,log)
    if istable(log) && ismember('run',log.Properties.VariableNames) && numel(unique(log.run))==1
    else
        error('Log must be a glmfit_log table with only entries having the same run name, as made by select_glmfit_runs.');
    end     
    [recording_idx,unique_recordings] = findgroups(recording_name);
    n_recordings = numel(unique_recordings);
    cellno=cellno(:);
    is_saved = false(numel(recording_name),1);   
    for i=1:n_recordings       
        if any(log.recording_name == unique_recordings(i))
            saved_cells = log.saved_cells{log.recording_name == unique_recordings(i)};
            is_saved(ismember(cellno,saved_cells) & recording_idx==i) = true;
        end
    end
end