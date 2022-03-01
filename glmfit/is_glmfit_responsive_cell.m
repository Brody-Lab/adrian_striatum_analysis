function is_responsive = is_glmfit_responsive_cell(recording_name,cellno,run)
    if istable(run) && ismember('run',run.Properties.VariableNames) && numel(unique(run.run))==1
        glmfit_log=run;
    elseif isscalar(run) && isstring(run)
        glmfit_log = select_glmfit_runs('run',run);        
    else
        error('Run must be a glmfit_log table with only entries having the same run name or it must be a scalar string specifying the runname.');
    end     
    [recording_idx,unique_recordings] = findgroups(recording_name);
    n_recordings = numel(unique_recordings);
    is_responsive = false(numel(recording_name),1);   
    for i=1:n_recordings       
        responsive_cells = glmfit_log.responsive_cells{glmfit_log.recording_name == unique_recordings(i)};
        is_responsive(ismember(cellno,responsive_cells) & recording_idx==i) = true;
    end
    fprintf('is_glmfit_responsive_cell: %g of %g cells were responsive for run %s.\n',sum(is_responsive),numel(is_responsive),glmfit_log.run(1));
end