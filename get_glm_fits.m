function [fits,T] = get_glm_fits(varargin)
    P=get_parameters;
    p=inputParser;
    p.KeepUnmatched=true;
    p.addParameter('dv_range',[-Inf Inf],@(x)validateattributes(x,{'numeric'},{'increasing','numel',2}));
    p.addParameter('ap_range',[-Inf Inf],@(x)validateattributes(x,{'numeric'},{'increasing','numel',2}));
    p.addParameter('ml_range',[0 Inf],@(x)validateattributes(x,{'numeric'},{'increasing','numel',2,'>=',0}));
    p.addParameter('region',{},@(x)validateattributes(x,{'char','cell'},{}));
    p.addParameter('spike_width_ms_range',[0 Inf],@(x)validateattributes(x,{'numeric'},{'increasing','numel',2,'>=',0}));
    p.parse(varargin{:});
    params=p.Results;
    unmatched = fieldnames(p.Unmatched);
    unrecognized = unmatched(~ismember(unmatched,union(P.glmfit_catalog_keys,P.glmfit_catalog_params)));
    for k=1:length(unrecognized)
        error('Unrecognized input: %s !!',unrecognized{k});% inputs must match the var args specified in the parsing lines above or be one of the fields of the catalog table, specified in "get_parameters.m"
    end
    T = readtable(P.glmfit_catalog_path);    
    [T,row_idx] = select_rows(T,varargin{:});
    T = add_parameter_struct(T);
    % TO DO: figure out which cell_infos to load (based on unique rat and
    % session) and load them.
    count=0;
    for i=1:height(T)
        [saved_cells,file_paths] = get_saved_cellnos(T.fit_path{i});
        fprintf('Loading %g cells in %s.\n',length(saved_cells),T.fit_path{i});
        for k=1:length(saved_cells)
            count=count+1;
            fits.stats{count}=load(file_paths{k});
            fits.fit_path{count} = T.fit_path{i};
            fits.run_idx(count) = i; % row index of T
        end
        % TO DO: select list of cells matching anatomical criteria in
        % associated cells files
        % load the fit files for each cell 
    end
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
    for t=1:height(T)
        param_path = fullfile(T.fit_path{t},'glmfit_params.mat');
        tmp=load(param_path);
        T.params{t}=tmp.params;
    end
end