function [glmfit_log,most_recent_run] = select_glmfit_runs(varargin)
    
    % assumes all sessions and cells have been run successfully for each run and that you always want the most recent given a set of fit params

    P = get_parameters();

    %% read the glmfit catalog (needs to be rebuilt if you added new fits)
    load(P.glmfit_catalog_path);
    
    %% select rows based on selection of catalog fields
    glmfit_log = select_rows(glmfit_log,varargin{:});
    
    %% select most recent run
    sorted_runs = sort(glmfit_log.run);
    most_recent_run = sorted_runs(end);
    glmfit_log = glmfit_log( glmfit_log.run==most_recent_run , :);
    
    %% check that all sessions are represented
    recordings_table = get_striatum_glm_recordings_table();    
    missing_sessions = ~ismember(recordings_table.recording_name,glmfit_log.recording_name);
    if all(missing_sessions)
        error('No glmfit runs matched these criteria.');
    elseif any(missing_sessions)
        error('With these criteria, %g recordings are missing from the glmfit log.',sum(missing_sessions));
    end
    
    %% check that all cells are represented
    if any(glmfit_log.n_missing_cells)
        error('Missing %g cells from this run across all recordings.',sum(glmfit_log.n_missing_cells));
    end

end

function [glmfit_log,include] = select_rows(glmfit_log,varargin)
    p=inputParser;
    p.KeepUnmatched=true;
    p.parse(varargin{:});
    fields = fieldnames(p.Unmatched);
    bad_fields = fields(~ismember(fields,glmfit_log.Properties.VariableNames));
    fields = setdiff(fields,bad_fields);
    if ~isempty(bad_fields)
        fprintf('Some unrecognized fields in the table:\n');
        display(bad_fields)
    end
    include=true(height(glmfit_log),1);
    for f=1:length(fields)
        value = p.Unmatched.(fields{f});
        include = include & glmfit_log.(fields{f})==value;
    end
    glmfit_log=glmfit_log(include,:);
end