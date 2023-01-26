function params = get_pcs(Cells,varargin)
    % function that does PCA on a neural population in a Cells file, aligned to
    % specific trial events. If you select "make_psth", the output field "psth" has the same format as a PETH
    % structure. options include criteria for including units, reference event,
    % time window and time resolution, number of PCs, etc.
    % requires npx_utils repo
    % for optional input arguments see select_units, MakedDataMatrix, 
    % run_PCA, and calc_psth_continuous to which all optional input arguments are passed.
    
    % Adrian Bondy, 2021 -- based partly on code by Wynne Stagnaro    
    
    
    p=inputParser;
    p.KeepUnmatched=true;
    p.addParameter('units',[]);         
    p.parse(varargin{:});
    params=p.Results;
    
    %% select units
    if ismember('units',p.UsingDefaults)
        params.units = select_units(Cells,varargin{:});
    end
    
    %% make data matrix on which to run PCA
    [cells_mat,params] = MakeDataMatrix(Cells,varargin{:},'units',params.units);
    params.times = reshape(params.times,numel(params.time_s),numel(params.trial_idx));    
    params.cells_mat=cells_mat;
    
    %% perform PCA
    [params.pca_output,pca_params] = run_PCA(params.cells_mat,varargin{:});
    % reshape times and scores
    for i=1:numel(params.pca_output)
        params.pca_output(i).score = reshape(params.pca_output(i).score,numel(params.time_s),numel(params.trial_idx),pca_params.npcs);
    end
    
    %% merge params fields
    fields = fieldnames(pca_params);
    for i=1:length(fields)
        params.(fields{i}) = pca_params.(fields{i});
    end

    %% assemble output structure
    params.label = [char(Cells.rat),' - ',char(Cells.sess_date),' - ',char(Cells.probe_serial),' - ',num2str(params.npcs),' PCs'];
    params.repo_path = fileparts(mfilename('fullpath'));
    [params.git_branch,params.commit] = return_git_status(params.repo_path);
    % copy some fields from the cells file for data integrity
    fields={'rec','Trials','sessid','last_modified','rat','sess_date','kSpikeWindowS','penetration','regions','probe_serial'};
    for f=1:length(fields)
        try
            params.(fields{f}) = Cells.(fields{f});
        catch
            params.(fields{f})=[];
        end
    end       
end