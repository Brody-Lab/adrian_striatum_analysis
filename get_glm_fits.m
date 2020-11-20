function [fits,T] = get_glm_fits(varargin)
    %% parse and validate inputs
    P=get_parameters;
    p=inputParser;
    p.KeepUnmatched=true;
    p.addParameter('dv_range',[-Inf Inf],@(x)validateattributes(x,{'numeric'},{'increasing','numel',2}));
    p.addParameter('ap_range',[-Inf Inf],@(x)validateattributes(x,{'numeric'},{'increasing','numel',2}));
    p.addParameter('ml_range',[0 Inf],@(x)validateattributes(x,{'numeric'},{'increasing','numel',2,'>=',0}));
    p.addParameter('region',{},@(x)validateattributes(x,{'char','cell'},{}));
    p.addParameter('spike_width_ms_range',[0 Inf],@(x)validateattributes(x,{'numeric'},{'increasing','numel',2,'>=',0}));
    p.addParameter('recency_mode','run',@(x)validateattributes(x,{'char'},{'nonempty'}));
    p.addParameter('min_reliability',0,@(x)validateattributes(x,{'numeric'},{'scalar','<=',1,'>=',0}));
    p.addParameter('max_reliability',1,@(x)validateattributes(x,{'numeric'},{'scalar','<=',1,'>=',0}));    
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
    
    % TO DO: figure out which cell_infos to load (based on unique rat and
    % session) and load them.
    count=0;
    %% assemble a fits structure with cells matching criteria but do not load the fits yet.
    fprintf('Determining which cells match criteria from %g sessions ...',height(T));
    for i=1:height(T)
        [saved_cells,file_paths] = get_saved_cellnos(T.fit_path{i});
        responsive_cells = T.params{i}.cellno(T.params{i}.responsive_enough);
        is_responsive = find(ismember(saved_cells,responsive_cells));
        % comment: below loop can be easily vectorized
        for k=is_responsive(:)'
            % skip load for unresponsive cells. responsive bug meant some
            % of these were saved. would be good to programmatically
            % control this. i want to see what happens with these cells.
            count=count+1;
            fits.stats_path{count,1} = file_paths{k};
            fits.fit_path{count,1} = T.fit_path{i};
            fits.run_idx(count,1) = i; % row index of T
            fits.rat(count,1) = T.rat(i);
            fits.sess_date(count,1)=T.sess_date(i);
            fits.cellno(count,1) = saved_cells(k);
            fits.save_time(count,1) = T.save_time(i);
        end
        % TO DO: select list of cells matching anatomical criteria in
        % associated cells files
        % load the fit files for each cell 
    end
    fprintf(' took %s.\n',timestr(toc));
    %% remove duplicate fits for the same cell
    if strcmp(params.recency_mode,'cell')
        fits = select_most_recent_cell(fits);
    end
    %% now load the fit statistics
    fprintf('Loading fit results for %g cells...',height(fits));tic;
    for i=1:height(fits)
        fits.stats{i} = load(fits.stats_path{i});
    end
    fprintf(' took %s.\n',timestr(toc));
end

function [T,include] = select_rows(T,varargin)
    p=inputParser;
    p.KeepUnmatched=true;
    p.parse(varargin{:});
    fields = fieldnames(p.Unmatched);
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
    fits=struct2table(fits);
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