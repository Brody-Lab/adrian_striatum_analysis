function [glmfit_log,most_recent_run] = select_glmfit_runs(varargin)

    % assumes glmfit catalog has been validated and you always want the most recent given a set of fit params

    P = get_parameters();
    p=inputParser;
    p.KeepUnmatched=true;
    p.addOptional('glmfit_log',[]);
    p.parse(varargin{:});
    %% select rows based on selection of catalog fields    
    if ismember('glmfit_log',p.UsingDefaults)
        glmfit_log = select_rows(load_glmfit_catalog(),varargin{:});        
    else
        glmfit_log = select_rows(varargin{:});                
    end

    %% select most recent run
    sorted_runs = sort(glmfit_log.run);
    most_recent_run = sorted_runs(end);
    glmfit_log = glmfit_log( glmfit_log.run==most_recent_run , :);
    
    %% check that at least some sessions match criteria
    recordings_table = get_striatum_glm_recordings_table();    
    missing_sessions = ~ismember(recordings_table.recording_name,glmfit_log.recording_name);
    if all(missing_sessions)
        error('No glmfit runs matched these criteria.');
    end
    
    fprintf('select_glmfit_runs: most recent run given criteria is %s.\n',most_recent_run);

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