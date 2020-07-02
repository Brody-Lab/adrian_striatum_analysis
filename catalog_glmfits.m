function catalog_glmfits(varargin)
    P=get_parameters;
    p=inputParser;
    p.addParameter('delete_empty',true,@(x)validateattributes(x,{'logical'},{'scalar'}));
    p.parse(varargin{:});
    params=p.Results;
    fits_paths = get_data_paths('data_path',fullfile(P.data_path,'fits'),'parent_dir',true,varargin{:});
    for i=1:length(fits_paths)
        if ~isdir(fits_paths{i})
            error('Could not fit fits path: %s\n',fits_paths{i});
        end
        run_list{i} = dir(fits_paths{i});
        run_list{i} = {run_list{i}.name};
        run_list{i} = run_list{i}(contains(run_list{i},'glmfit_'));
        run_list{i} = cellfun(@(x)fullfile(fits_paths{i},x),run_list{i},'uni',0);
        for k=1:length(run_list{i})
            stats_dir = fullfile(run_list{i}{k},'stats');
            if isdir(stats_dir)
                tmp=dir(stats_dir);
                if sum([tmp.bytes])==0
                    if params.delete_empty
                        % this run didn't saved out data
                        status=rmdir(run_list{i}{k},'s');
                        if ~status
                            error('Failed to delete empty folder: %s.\n',run_list{i}{k});
                        end
                    end
                else
                    params_path=fullfile(run_list{i}{k},'glmfit_params.mat');
                    if exist(params_path,'file')
                        glmfit_params(k,i)=load(params_path);
                    else
                       error('Params file not found: %s.\n',params_path); 
                    end
                end
            elseif params.delete_empty
                % this run didn't save out data or it's in an early format
                % and no longer needed
                status=rmdir(run_list{i}{k},'s');
                if ~status
                    error('Failed to delete empty folder: %s.\n',run_list{i}{k});
                end                
            end
        end 
    end
    %% make params table
    count=0;
    fields=P.glmfit_catalog_params;
    for t=1:numel(glmfit_params)
        if ~isempty(glmfit_params(t).params)
            count=count+1;
            for f=1:length(fields)
                this_field=fields{f};
                % special cases
                switch this_field
                    case 'run'
                        [~,T(count).(this_field)] = fileparts(glmfit_params(t).params.save_path);
                    case 'cells_file'
                        T(count).(this_field) = glmfit_params(t).params.dm.dspec.expt.param.mat_file_name;
                    case 'link'
                        T(count).(this_field) = func2str(glmfit_params(t).params.link.Link);
                    otherwise
                        if isfield(glmfit_params(t).params,this_field)                
                            T(count).(this_field) = glmfit_params(t).params.(this_field); 
                            if ischar(T(count).(this_field))
                                T(count).(this_field) = string(T(count).(this_field));
                            end
                            T(count).(this_field) = unique(T(count).(this_field));                            
                        end
                end
            end
        end
    end
    T=struct2table(T);
    writetable(T,P.glmfit_catalog_path);
end

%% what I want this code to do:
% 1. delete empty runs DONE
% 2. get all params for runs and make a table DONE
% 3. warn if there are sessions for which some param sets haven't been run