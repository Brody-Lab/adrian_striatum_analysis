function glmfit_log = catalog_glmfits(varargin)

    %% parse and validate inputs
    p=inputParser;
    p.addParameter('remove_empty_runs',true,@(x)validateattributes(x,{'logical'},{'scalar'}));
    p.addParameter('delete_spikes',false,@(x)validateattributes(x,{'logical'},{'scalar'}));    
    p.addParameter('trim_params',false,@(x)validateattributes(x,{'logical'},{'scalar'}));        
    p.addParameter('remove_missing_mode','sessions',@(x)validatestring(x,{'sessions','cells','none'}));
    p.addParameter('delete_removed',false,@(x)validateattributes(x,{'logical'},{'scalar'})); 
    p.addParameter('max_days_old',Inf);
    p.parse(varargin{:});
    params=p.Results;
    P = get_parameters();    
    switch params.remove_missing_mode
        case 'sessions'
            fprintf('Runs with missing sessions will not be included in the log.\n');
        case 'cells'
            fprintf('Runs with any missing cells will not be included in the log.\n');
        case 'none'
            fprintf('All runs, no matter how imcomplete, will be included in the log.');
    end
    if params.delete_removed && params.remove_missing_mode~="none"
        fprintf('Removed runs will be physically deleted from %s.\n',P.fit_path);
    else
        fprintf('Nothing will be physically deleted from %s.\n',P.fit_path);        
    end
    
    %% load basic params of all glm fit runs in database
    [glmfit_params,saved_cells] = find_glmfit_params(params);
    
    %% convert into a glmfit_catalog table
    glmfit_log = make_catalog_from_params(glmfit_params,saved_cells);
    
    %% validate catalog
    glmfit_log = validate_glmfit_log(glmfit_log,params);
    
    %% save catalog
    fprintf('Saving catalog to %s\n',P.glmfit_catalog_path);tic;
    tic;save(P.glmfit_catalog_path,'glmfit_log','-v7');toc
    fprintf(' ... took %s.\n',timestr(toc));
    
end

function [glmfit_params,saved_cells] = find_glmfit_params(params)
    
    %% get database paths
    paths = get_data_paths();

    %% check all recordings have corresponding fit paths
    for i=length(paths):-1:1
        if ~isfolder(paths(i).fit_path)
            error('Could not find fits path: %s\n',paths(i).fit_path);
        end      
    end
    
    %% loop over recordings and fit runs, loading the params file and checking that files have been made for all cells
    % loop over recordings   
    for i=1:length(paths)
        [glmfit_params{i},saved_cells{i}] = find_glmfit_params_internal(paths(i),params);
    end
    glmfit_params = [glmfit_params{:}];
    saved_cells = [saved_cells{:}];
end

function [glmfit_params,saved_cells] = find_glmfit_params_internal(path,params)
    % get glmfit runs for this recording. (just based on folder existence. run may have failed.)             
    run_list = get_run_list(path.fit_path); 
    % loop over glm fit runs for this recording
    if params.max_days_old<Inf
        days_old = days(datetime(datestr(now)) - cellfun(@run_to_datetime,run_list));
        fprintf('Excluding %g runs more than %g days old.\n',sum(days_old>params.max_days_old),params.max_days_old);
        run_list = run_list(days_old<=params.max_days_old);
    end
    for k=length(run_list):-1:1
        stats_dir = fullfile(run_list{k},'stats');
        if params.delete_spikes
            spikes_dir = fullfile(run_list{k},'spikes');
            if isfolder(spikes_dir)            
                fprintf('Removing spikes folder %s  ...',spikes_dir);tic;
                status=rmdir(spikes_dir,'s');
                if ~status
                    error('Failed to delete spikes folder: %s.\n',spikes_dir);
                end                    
                fprintf('took %s.\n',timestr(toc));
            end                
        end     
        if isfolder(stats_dir)
            tmp=dir(stats_dir);
            size=[tmp.bytes];
            if sum(size)==0
                if params.delete_removed && params.remove_missing_mode~="none"  % empty stats folder 
                    remove_run(run_list{k});
                end
            else
                % load params file
                params_path=fullfile(run_list{k},'glmfit_params.mat');
                if exist(params_path,'file')
                    [~,run_name] = fileparts(run_list{k});
                    fprintf('Loading run %s for recording: %s',run_name,path.recording_name);tic;
                    glmfit_params(k)=load(params_path);
                    if params.trim_params
                        glmfit_params(k).params = trim_params(glmfit_params(k).params); % used to get rid of bases and zscore in older files
                        params = glmfit_params(k).params;
                        save(params_path,'params','-v7');
                    end
                    fprintf('  ... took %s\n',timestr(toc));
                else
                   error('Params file not found: %s.\n',params_path); % something has gone wrong if you've got to this line.
                end
                % check which cells have been saved               
                saved_cells{k} = cellfun(@(x)str2double(x),regexprep({tmp(size>0).name},'.*cell(.*).mat','$1'),'uni',0);
                saved_cells{k}=[saved_cells{k}{:}];
            end
        elseif params.delete_removed && params.remove_missing_mode~="none" % no stats folder made
            remove_run(run_list{k});          
        end
    end 
    
end

function glmfit_log = make_catalog_from_params(glmfit_params,saved_cells)

    P = get_parameters();
    default_params = glmfit_parse_params();    
    default_params_legacy = glmfit_parse_params_legacy();        

    %% make params table
    count=0;
    fields=P.glmfit_catalog_params;
    for t=numel(glmfit_params):-1:1 % scalar indexing across all params
        if ~isempty(glmfit_params(t).params)
            count=count+1;
            T(count).fit_method = "glmnet";            
            for f=1:length(fields)
                this_field=fields{f};
                % special cases
                switch this_field
                    case 'responsive_cells'
                        T(count).responsive_cells = glmfit_params(t).params.cellno(glmfit_params(t).params.responsive_enough);                        
                    case 'saved_cells'
                        T(count).saved_cells = saved_cells{t};                   
                    case 'n_missing_cells'
                        T(count).n_missing_cells = numel(setdiff(T(count).responsive_cells,saved_cells{t}));
                    case 'recording_name'
                        T(count).(this_field) = regexprep(glmfit_params(t).params.save_path,'.*fits(.)(.*)(.)glmfit.*','$2');   
                    case 'run'
                        [~,T(count).(this_field)] = fileparts(glmfit_params(t).params.save_path);                            
                    case 'link'
                        if ~isstruct(glmfit_params(t).params.link)
                            glmfit_params(t).params.link = make_link(glmfit_params(t).params.link,glmfit_params(t).params.distribution);
                        end
                        T(count).(this_field) = func2str(glmfit_params(t).params.link.Link);
                    otherwise
                        if isfield(glmfit_params(t).params,this_field)                
                            T(count).(this_field) = glmfit_params(t).params.(this_field); 
                        elseif this_field=="sess_date"
                           [~, T(count).(this_field)] = extract_info_from_npx_path(glmfit_params(t).params.save_path);
                        elseif this_field=="lambda_correct_zscore" && glmfit_params(t).params.git_commit=="9e53fd2cb4adf6852328e70cb7a072092c1da944"
                            T(count).lambda_correct_zscore=true;
                        elseif this_field=="alpha"
                            T(count).alpha=0;
                            T(count).fit_method="neuroglmfit";                            
                        elseif this_field=="foldseed" || this_field=="max_glmnet_attempts" || this_field=="lambda" || this_field=="glmnet_thresh" || this_field=="maxIter"
                            T(count).(this_field)=NaN;
                        elseif isfield(default_params,this_field)
                            T(count).(this_field) = default_params.(this_field);
                        elseif isfield(default_params_legacy,this_field)
                            T(count).(this_field) = default_params_legacy.(this_field);                            
                        end
                end
                if isfield(T(count),this_field) && ischar(T(count).(this_field))
                    T(count).(this_field) = string(T(count).(this_field));
                end
            end
            if T(count).git_branch=="neuro_glmnet"
                T(count).fit_method = "glmnet";
            end
            T(count).dm.dspec.expt = rmfield(T(count).dm.dspec.expt,'trial');            
            % check that all cells have been saved and throw warnings
            % otherwise
            params_cells = glmfit_params(t).params.cellno;                    
            unknown_cells = setdiff(saved_cells{t},params_cells);
            if ~isempty(unknown_cells)
                error('%g saved cells that are not listed in the params file!',length(unknown_cells));
            end
            unresponsive_cells_fit = setdiff(saved_cells{t},T(count).responsive_cells);
            if ~isempty(unresponsive_cells_fit)
                warning('%g saved cells that are not responsive enough!',length(unresponsive_cells_fit));
            end
            if T(count).n_missing_cells
                warning('%g missing cells for recording name %s run %s!',T(count).n_missing_cells,T(count).recording_name,T(count).run);
            end            
        end
    end
    glmfit_log=struct2table(T);
end

function run_list = get_run_list(this_fit_path)       
    run_list = dir(this_fit_path); 
    run_list = {run_list.name};
    run_list = run_list(contains(run_list,'glmfit_'));
    run_list = cellfun(@(x)fullfile(this_fit_path,x),run_list,'uni',0);
end

function glmfit_log = validate_glmfit_log(glmfit_log,params)
    % checks raw glmfit log (i.e. obtained by a directory search) and
    % optionally removes runs where all cells or sessions are not represented. if
    % "delete_removed" is true, these files are deleted too, not
    % just removed from the log.
    [run_idx,unique_runs] = findgroups(glmfit_log.run);
    n_runs = numel(unique_runs); 
    remove = false(height(glmfit_log),1);
    for i=1:n_runs
        this_idx = run_idx==i;
        %% check that all sessions are represented
        recordings_table = get_striatum_glm_recordings_table();    
        missing_sessions = ~ismember(recordings_table.recording_name,glmfit_log.recording_name(this_idx));
        if any(missing_sessions)
            warning('%g recordings are missing from the glmfit log for run %s.',sum(missing_sessions),unique_runs(i));
            if params.remove_missing_mode~="none"
                remove(i)=true;
            end
        end
        %% check that all cells are represented
        if any(glmfit_log.n_missing_cells(this_idx))
            warning('Missing %g cells from recordings found for run %s.',sum(glmfit_log.n_missing_cells(this_idx)),unique_runs(i));
            if params.remove_missing_mode=="cells"
                remove(i)=true;
            end
        end
        %% optionally delete removed runs    
        if remove(i) && params.delete_removed
            these_recordings = glmfit_log.recording_name(this_idx);
            run_folders = get_run_folder(repmat(unique_runs(i),numel(these_recordings),1),these_recordings);
            for k=1:numel(run_folders)
                remove_run(run_folders(k));
            end
        end
    end
    glmfit_log = glmfit_log(~remove,:);
end

function params = trim_params(params)
    for i=1:numel(params.dm.dspec.covar)
        try
            params.dm.dspec.covar(i).basis = rmfield(params.dm.dspec.covar(i).basis,'B');
        end
        try
            params.dm.dspec.covar(i).basis = rmfield(params.dm.dspec.covar(i).basis,'tr');
        end            
    end
    try
        params = rmfield(params,'zscore');
    end
end

function d = run_to_datetime(run)
    d = regexprep(char(run),'.*glmfit_([0-9][0-9][0-9][0-9]_[0-9][0-9]_[0-9][0-9])_.*','$1');
    d = datetime(d,'InputFormat','yyyy_MM_dd');
end