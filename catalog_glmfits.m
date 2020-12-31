function catalog_glmfits(varargin)
    P=get_parameters;
    p=inputParser;
    p.addParameter('delete_empty',true,@(x)validateattributes(x,{'logical'},{'scalar'}));
    p.addParameter('verbose',false,@(x)validateattributes(x,{'logical'},{'scalar'}));
    p.addParameter('fix_responsive',false,@(x)validateattributes(x,{'logical'},{'scalar'}));
    p.parse(varargin{:});
    params=p.Results;
    fits_paths = get_data_paths('data_path',fullfile(P.data_path,'fits'),'parent_dir',true,varargin{:});
    for i=1:length(fits_paths)
        if ~isdir(fits_paths{i})
            error('Could not fit fits path: %s\n',fits_paths{i});
        end
        run_list{i} = dir(fits_paths{i}); 
        run_list{i} = {run_list{i}.name};
        cells_path = run_list{i}(contains(run_list{i},'Cells'));
        run_list{i} = run_list{i}(contains(run_list{i},'glmfit_'));
        run_list{i} = cellfun(@(x)fullfile(fits_paths{i},x),run_list{i},'uni',0);
        run_list_up_one = dir(fileparts(fits_paths{i})); % search one directory up too because of bug saving some of the fits in the parent folder
        run_list_up_one = {run_list_up_one.name};
        run_list_up_one = run_list_up_one(contains(run_list_up_one,'glmfit_'));
        run_list_up_one = cellfun(@(x)fullfile(fileparts(fits_paths{i}),x),run_list_up_one,'uni',0);   
        run_list{i} = [run_list{i} run_list_up_one];
        if params.fix_responsive
            if length(cells_path)==1
                cells_path = cells_path{1};
            else
                error('Error locating cells file for this session.');
            end                
            cells_path = fullfile(fits_paths{i},cells_path);
            cells_path = strrep(cells_path,'data','cells');
            cells_path = strrep(cells_path,'fits','cells');
            cells_path = correct_file_path(cells_path);
            if exist(cells_path,'file')
                Cells{i} = load(cells_path);
                fields = fieldnames(Cells{i});
                if length(fields)==1
                    Cells{i} = Cells{i}.(fields{1});
                end
            else
                error('Error locating cells file for this session.');
            end
        end
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
                        if params.fix_responsive
                            glmfit_params(k,i).params = fix_responsive(params_path,Cells{i});
                        else
                            glmfit_params(k,i)=load(params_path);
                        end
                    else
                       error('Params file not found: %s.\n',params_path); 
                    end
                    responsive_cells = glmfit_params(k,i).params.cellno(glmfit_params(k,i).params.responsive_enough);
                    params_cells = glmfit_params(k,i).params.cellno;                    
                    saved_cells = cellfun(@(x)str2num(x),regexprep({tmp.name},'.*cell(.*).mat','$1'),'uni',0);
                    saved_cells=[saved_cells{:}];
                    unknown_cells = setdiff(saved_cells,params_cells);
                    if ~isempty(unknown_cells)
                        error('%g saved cells that are not listed in the params file!',length(unknown_cells));
                    end
                    unresponsive_cells_fit = setdiff(saved_cells,responsive_cells);
                    if ~isempty(unresponsive_cells_fit)
                        warning('%g saved cells that are not responsive enough!',length(unresponsive_cells_fit));
                    end
                    missing_cells = setdiff(responsive_cells,saved_cells);
                    if ~isempty(missing_cells)
                        warning('%g missing cells for session %s!',length(missing_cells),run_list{i}{k});
                        n_missing_cells(k,i) = length(missing_cells);
                    else
                        n_missing_cells(k,i) = 0;                        
                        if params.verbose
                            fprintf('All %g cells accounted for in session %s\n',length(responsive_cells),run_list{i}{k});
                        end
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
                        [~,run] = fileparts(glmfit_params(t).params.save_path);                        
                        T(count).(this_field) = glmfit_params(t).params.dm.dspec.expt.param.mat_file_name;
                        T(count).fit_path = fullfile(correct_file_path(fileparts(T(count).(this_field))),run);
                        T(count).fit_path = strrep(T(count).fit_path,'data','fits');
                        T(count).fit_path = strrep(T(count).fit_path,'cells','fits');                        
                        if ~isdir(T(count).fit_path)
                            error('Cannot find putative local fit path: %s',T(count).fit_path);
                        end                        
                    case 'link'
                        T(count).(this_field) = func2str(glmfit_params(t).params.link.Link);
                    case 'save_time'
                        if isfield(glmfit_params(t).params,this_field)
                            if ~isa(glmfit_params(t).params.(this_field),'datetime')
                                T(count).save_time = datetime(glmfit_params(t).params.(this_field),'InputFormat','yyyy_MM_dd_HH_mm_ss');
                            else
                                T(count).save_time = glmfit_params(t).params.save_time;
                            end
                        end
                    otherwise
                        if isfield(glmfit_params(t).params,this_field)                
                            T(count).(this_field) = glmfit_params(t).params.(this_field); 
                            if ischar(T(count).(this_field))
                                T(count).(this_field) = string(T(count).(this_field));
                            end
                            T(count).(this_field) = unique(T(count).(this_field));                            
                        elseif this_field=="sess_date"
                           [~,cells_file_name] = fileparts(glmfit_params(t).params.dm.dspec.expt.param.mat_file_name);
                           [~, T(count).(this_field)] = extract_info_from_npx_path(cells_file_name);
                        end
                end
            end
            T(count).n_missing_cells = n_missing_cells(t);
        end
    end
    T=struct2table(T);
    writetable(T,P.glmfit_catalog_path);    
end