function is_responsive = is_glmfit_responsive_cell(recording_name,cellno,log)
    if istable(log) && ismember('run',log.Properties.VariableNames) && numel(unique(log.run))==1
    else
        error('Log must be a glmfit_log table with only entries having the same run name, as made by select_glmfit_runs.');
    end     
    [recording_idx,unique_recordings] = findgroups(recording_name);
    n_recordings = numel(unique_recordings);
    is_responsive = false(numel(recording_name),1);   
    for i=1:n_recordings       
        responsive_cells = log.responsive_cells{log.recording_name == unique_recordings(i)};
        is_responsive(ismember(cellno,responsive_cells) & recording_idx==i) = true;
    end
    fprintf('is_glmfit_responsive_cell: %g of %g cells were responsive for run %s.\n',sum(is_responsive),numel(is_responsive),log.run(1));
end