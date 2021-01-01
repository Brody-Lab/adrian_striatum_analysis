function [fits,T] = get_glm_fits(varargin)
    %% parse and validate inputs
    P=get_parameters;
    p=inputParser;
    p.KeepUnmatched=true;
    p.addParameter('DV',[-Inf Inf],@(x)validateattributes(x,{'numeric'},{'increasing','numel',2}));
    p.addParameter('AP',[-Inf Inf],@(x)validateattributes(x,{'numeric'},{'increasing','numel',2}));
    p.addParameter('ML',[0 Inf],@(x)validateattributes(x,{'numeric'},{'increasing','numel',2,'>=',0}));
    p.addParameter('region',{},@(x)validateattributes(x,{'char','cell'},{}));
    p.addParameter('spike_width_ms',[0 Inf],@(x)validateattributes(x,{'numeric'},{'increasing','numel',2,'>=',0}));
    p.addParameter('recency_mode','run',@(x)validateattributes(x,{'char'},{'nonempty'}));
    p.addParameter('reliability',[0 1],@(x)validateattributes(x,{'numeric'},{'increasing','numel',2}));
    p.addParameter('D2Phototagging',NaN,@(x)validateattributes(x,{'numeric'},{'scalar'}));    
    p.parse();
    default_params = p.Results;
    cell_info_fields = {'DV','AP','ML','spike_width_ms','reliability','region','D2Phototagging'};
    p.parse(varargin{:});
    params=p.Results;    
    validatestring(params.recency_mode,{'run','cell'},'get_glm_fits','recency_mode');
    unmatched = fieldnames(p.Unmatched);
    unrecognized = unmatched(~ismember(unmatched,union(P.glmfit_catalog_keys,P.glmfit_catalog_params)));
    for k=1:length(unrecognized)
        error('Unrecognized input: %s !!',unrecognized{k});% inputs must match the var args specified in the parsing lines above or be one of the fields of the catalog table, specified in "get_parameters.m"
    end
    %% read the glmfit catalog (needs to be rebuilt if you added new fits)
    T = read_catalog(P);  
    %% select rows based on selection of catalog fields
    [T,row_idx] = select_rows(T,varargin{:});
    %% select only most recent fit if there are duplicated matching the selected catalog fields
    if strcmp(params.recency_mode,'run')
        T = select_most_recent_run(T);
    end
    %% load all parameters files
    T = add_parameter_struct(T);
    this_count=0;
    %% assemble a fits structure with cells matching criteria but do not load the fits yet.
    fprintf('Determining which cells match criteria from %g sessions ...',height(T));
    for i=1:height(T)
        [saved_cells,file_paths] = get_saved_cellnos(T.fit_path{i});
        is_responsive = T.params{i}.responsive_enough;
        cell_info_path = fullfile(strrep(fileparts(T.fit_path{i}),'fits','cells'),'cell_info.mat'); % b/c the cells path in T is the tiger path. why?
        n_cells = numel(T.params{i}.cellno);
        include = is_responsive(:) & validate_ranges(cell_info_path,params,default_params,cell_info_fields,n_cells);
        include = find(include);
        % comment: below loop can be easily vectorized   
        count=0;
        for k=include'
            % skip load for unresponsive cells. responsive bug meant some
            % of these were saved. would be good to programmatically
            % control this. i want to see what happens with these cells.
            this_count=this_count+1;
            count=count+1;
            fits.stats_path{this_count,1} = file_paths{count};
            fits.fit_path{this_count,1} = T.fit_path{i};
            fits.run_idx(this_count,1) = i; % row index of T
            fits.rat(this_count,1) = T.rat(i);
            fits.sess_date(this_count,1)=T.sess_date(i);
            fits.cellno(this_count,1) = k;
            fits.save_time(this_count,1) = T.save_time(i);
        end        
        % load the fit files for each cell 
    end
    fits = struct2table(fits);
    fprintf(' took %s.\n',timestr(toc));
    %% remove duplicate fits for the same cell
    if strcmp(params.recency_mode,'cell')
        fits = select_most_recent_cell(fits);
    end
    %% now load the fit statistics
    
    fprintf('\n%s: Loading fits for %i cells...\n', mfilename, n_cells);
    tic;    
    stats_path = fits.stats_path;
    parfor i=1:height(fits)
        stats{i,1} = load(stats_path{i});      
    end
    fits.stats=stats;
    fprintf(' took %s.\n',timestr(toc));
end

function [T,include] = select_rows(T,varargin)
    p=inputParser;
    p.KeepUnmatched=true;
    p.parse(varargin{:});
    fields = fieldnames(p.Unmatched);
    bad_fields = fields(~ismember(fields,T.Properties.VariableNames));
    fields = setdiff(fields,bad_fields);
    if ~isempty(bad_fields)
        fprintf('Some unrecognized fields in the table:\n');
        display(bad_fields)
    end
    include=true(height(T),1);
    for f=1:length(fields)
        value = p.Unmatched.(fields{f});
        include = include & T.(fields{f})==value;
    end
    T=T(include,:);
end

function T = add_parameter_struct(T)
    fprintf('Loading fit param files for %g fits ...',height(T));tic
    for t=1:height(T)
        param_path = fullfile(T.fit_path{t},'glmfit_params.mat');
        tmp=load(param_path);
        T.params{t}=tmp.params;
    end
    fprintf(' took %s.\n',timestr(toc));
end

function T = select_most_recent_run(T)
    G = findgroups(T.rat,T.sess_date);
    n_groups = numel(unique(G));
    for i=1:n_groups
       these_sessions = find(G==i);
       dates = T.save_time(these_sessions);
       [~,idx] = max(dates);
       keep_idx(i) = these_sessions(idx);
    end
    T = T(keep_idx,:);
end

function fits = select_most_recent_cell(fits)
    unique_S = unique(fits(:,{'rat','sess_date','cellno'}));
    keep=false(height(fits),1);
    for i=1:height(unique_S)
        [~,idx] = ismember(fits(:,{'rat','sess_date','cellno'}),unique_S(i,:));
        idx=logical(idx);
        if sum(idx)>1
            [~,best] = max(fits.save_time(idx));
            idx=find(idx);
            keep(idx(best))=true;
        else
            keep(idx)=true;
        end
    end        
    fits=fits(keep,:);
end

function include = validate_ranges(cell_info_path,params,default_params,cell_info_fields,n_cells)
    cell_info_needed=false;
    include = true(n_cells,1);    
    check_range=false(numel(cell_info_fields),1);
    for f=1:length(cell_info_fields)
        if iscell(params.(cell_info_fields{f}))
            for i=1:length(params.(cell_info_fields{f}))
                if any(default_params.(cell_info_fields{f}){i}~=params.(cell_info_fields{f}){i})
                    check_range(f)=true;
                    cell_info_needed=true;
                end
            end
        else
            if any(default_params.(cell_info_fields{f})~=params.(cell_info_fields{f})) && ...
                    ~any(isnan(params.(cell_info_fields{f}))) && ~any(isnan(default_params.(cell_info_fields{f})))
                check_range(f)=true;
                cell_info_needed=true;
            end
        end
    end
    if cell_info_needed
        if ~isfile(cell_info_path)
            error('Could not load cell_info file at %s.',cell_info_path);
        else
            load(cell_info_path);
        end
        for i=1:length(check_range)
            if check_range(i)
                if iscell(params.(cell_info_fields{i}))
                    include = include & ismember(cell_info.(cell_info_fields{i}),params.(cell_info_fields{i}));
                elseif isnumeric(params.(cell_info_fields{i}))
                    if numel(params.(cell_info_fields{i}))==2
                        include = include & cell_info.(cell_info_fields{i})>=params.(cell_info_fields{i})(1) & ...
                            cell_info.(cell_info_fields{i})<=params.(cell_info_fields{i})(2);
                    elseif numel(params.(cell_info_fields{i}))==1
                        include = include & cell_info.(cell_info_fields{i})==params.(cell_info_fields{i});
                    end
                else
                    error('');
                end
            end
        end
    end
end