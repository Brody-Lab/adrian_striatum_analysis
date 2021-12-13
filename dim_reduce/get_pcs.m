function stats = get_pcs(Cells,varargin)
    % function that does PCA on a neural population in a Cells file, aligned to
    % specific trial events. If you select "make_psth", the output field "psth" has the same format as a PETH
    % structure. options include criteria for including units, reference event,
    % time window and time resolution, number of PCs, etc.
    % requires npx_utils repo
    % for optional input arguments see select_units, MakedDataMatrix, 
    % run_PCA, and calc_psth_continuous to which all optional input arguments are passed.
    
    % Adrian Bondy, 2021 -- based partly on code by Wynne Stagnaro    
    
    
    %% parse and validate inputs
    p=inputParser;
    p.KeepUnmatched=true;
    p.addParameter('make_psth',false,@(x)validateattributes(x,{'logical'},{}));
    p.parse(varargin{:});
    
    %% select units
    units = select_units(Cells,varargin{:});
    
    %% make data matrix on which to run PCA
    [cells_mat,params] = MakeDataMatrix(Cells,varargin{:},'units',units);

    
    %% perform PCA
    stats = run_PCA(cells_mat,varargin{:});
    npcs = size(stats.score,2);
    % reshape times and scores
    stats.score = reshape(stats.score,numel(params.time_s),numel(params.trial_idx),npcs);
    params.times = reshape(params.times,numel(params.time_s),numel(params.trial_idx));    
    
    %% assemble output structure
    fields=fieldnames(params);
    for f=1:length(fields)
        stats.(fields{f}) = params.(fields{f});
    end   
    stats.label = [Cells.rat,' - ',Cells.sess_date,' - ',Cells.probe_serial,' - ',num2str(npcs)',' PCs'];
    stats.repo_path = fileparts(mfilename('fullpath'));
    [stats.git_branch,stats.commit] = return_git_status(stats.repo_path);
    stats.dim = npcs;

    if p.Results.make_psth
        stats.psth = get_pc_psth(Cells,stats,varargin{:});
    end
    
end