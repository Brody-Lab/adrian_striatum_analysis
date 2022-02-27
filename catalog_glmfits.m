function glmfit_catalog = catalog_glmfits(varargin)

    %% parse and validate inputs
    p=inputParser;
    p.addParameter('remove_empty_runs',true,@(x)validateattributes(x,{'logical'},{'scalar'}));
    p.parse(varargin{:});
    params=p.Results;
    
    %% load basic params of all glm fit runs in database
    [glmfit_params,saved_cells] = load_glmfit_params(params.remove_empty_runs);
    
    %% convert into a glmfit_catalog table
    glmfit_catalog = make_catalog_from_params(glmfit_params,saved_cells);
    
end

function [glmfit_params,saved_cells] = load_glmfit_params(remove_empty_runs)
    
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
    for i=length(paths):-1:1 
        % get glmfit runs for this recording. (just based on folder existence. run may have failed.)             
        run_list{i} = get_run_list(paths(i).fit_path); 
        % loop over glm fit runs for this recording
        for k=1:length(run_list{i})
            stats_dir = fullfile(run_list{i}{k},'stats');
            if isfolder(stats_dir)
                tmp=dir(stats_dir);
                size=[tmp.bytes];
                if sum(size)==0
                    if remove_empty_runs % empty stats folder 
                        remove_run(run_list{i}{k});
                    end
                else
                    % load params file
                    params_path=fullfile(run_list{i}{k},'glmfit_params.mat');
                    if exist(params_path,'file')
                        [~,run_name] = fileparts(run_list{i}{k});
                        fprintf('Loading run %s for recording %g of %g: %s',run_name,i,length(paths),paths(i).recording_name);tic;
                        glmfit_params(k,i)=load(params_path);
                        fprintf('  ... took %s\n',timestr(toc));
                    else
                       error('Params file not found: %s.\n',params_path); % something has gone wrong if you've got to this line.
                    end
                    % check which scells have been saved               
                    saved_cells{k,i} = cellfun(@(x)str2double(x),regexprep({tmp(size>0).name},'.*cell(.*).mat','$1'),'uni',0);
                    saved_cells{k,i}=[saved_cells{k,i}{:}];
                end
            elseif remove_empty_runs % no stats folder made
                remove_run(run_list{i}{k});          
            end
        end 
    end
    
end

function glmfit_log = make_catalog_from_params(glmfit_params,saved_cells)

    P = get_parameters();

    %% make params table
    count=0;
    fields=P.glmfit_catalog_params;
    for t=1:numel(glmfit_params) % scalar indexing across all params
        if ~isempty(glmfit_params(t).params)
            count=count+1;
            for f=1:length(fields)
                this_field=fields{f};
                % special cases
                switch this_field
                    case 'recording_name'
                        T(count).(this_field) = regexprep(glmfit_params(t).params.save_path,'.*fits/(.*)/glmfit.*','$1');   
                    case 'run'
                        [~,T(count).(this_field)] = fileparts(glmfit_params(t).params.save_path);                            
                    case 'link'
                        T(count).(this_field) = func2str(glmfit_params(t).params.link.Link);
                    otherwise
                        if isfield(glmfit_params(t).params,this_field)                
                            T(count).(this_field) = glmfit_params(t).params.(this_field); 
                            if ischar(T(count).(this_field))
                                T(count).(this_field) = string(T(count).(this_field));
                            end
                            T(count).(this_field) = unique(T(count).(this_field));                            
                        elseif this_field=="sess_date"
                           [~, T(count).(this_field)] = extract_info_from_npx_path(glmfit_params(t).params.save_path);
                        end
                end
                if isfield(T(count),this_field) && ischar(T(count).(this_field))
                    T(count).(this_field) = string(T(count).(this_field));
                end
            end
            % check that all cells have been saved and throw warnings
            % otherwise
            T(count).responsive_cells = glmfit_params(t).params.cellno(glmfit_params(t).params.responsive_enough);
            T(count).saved_cells = saved_cells{t};
            params_cells = glmfit_params(t).params.cellno;                    
            unknown_cells = setdiff(saved_cells{t},params_cells);
            if ~isempty(unknown_cells)
                error('%g saved cells that are not listed in the params file!',length(unknown_cells));
            end
            unresponsive_cells_fit = setdiff(saved_cells{t},T(count).responsive_cells);
            if ~isempty(unresponsive_cells_fit)
                warning('%g saved cells that are not responsive enough!',length(unresponsive_cells_fit));
            end
            T(count).n_missing_cells = ~isempty(setdiff(T(count).responsive_cells,saved_cells{t}));
            if T(count).n_missing_cells
                warning('%g missing cells for recording name %s run %s!',T(count).n_missing_cells,T(count).recording_name,T(count).run);
            end            
        end
    end
    glmfit_log=struct2table(T);
    fprintf('Saving catalog to %s\n',P.glmfit_catalog_path);tic;
    save(P.glmfit_catalog_path,'glmfit_log');
    fprintf(' ... took %s.\n',timestr(toc));
end

function run_list = get_run_list(this_fit_path)       
    run_list = dir(this_fit_path); 
    run_list = {run_list.name};
    run_list = run_list(contains(run_list,'glmfit_'));
    run_list = cellfun(@(x)fullfile(this_fit_path,x),run_list,'uni',0);
end

function remove_run(run_folder)
    status=rmdir(run_folder,'s');
    if ~status
        error('Failed to delete empty folder: %s.\n',run_folder);
    end    
end