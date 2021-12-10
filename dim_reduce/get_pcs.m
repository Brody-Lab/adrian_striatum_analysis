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
    
    %% optionally make PSTHs    
    if p.Results.make_psth
        PB_set_constants;        
        states=kPETH.refEvents;
        states=union(states,{'shuffled_left_clicks','shuffled_right_clicks'});
        kSpikeWindowS=Cells.kSpikeWindowS;
        kSpikeWindowS.shuffled_left_clicks = kSpikeWindowS.left_clicks;
        kSpikeWindowS.shuffled_right_clicks = kSpikeWindowS.right_clicks;   
        kPETH.timeS.shuffled_left_clicks = kPETH.timeS.left_clicks;
        kPETH.timeS.shuffled_right_clicks = kPETH.timeS.right_clicks;  
        kPETH.type.shuffled_left_clicks = kPETH.type.left_clicks;
        kPETH.type.shuffled_right_clicks = kPETH.type.right_clicks;        
        kPETH.stdS.shuffled_left_clicks = kPETH.stdS.left_clicks;
        kPETH.stdS.shuffled_right_clicks = kPETH.stdS.right_clicks;          
        stateTimes = Cells.Trials.stateTimes;
        % uncomment lines to get psths for a particular subset of trials (like
        % hits or misses)
        %trials = ~Cells.Trials.is_hit;    
        %stateTimes = select_trials_from_statetimes(Cells.Trials.stateTimes,trials);
        for i=1:npcs 
            for eves=states(:)'
                eve = eves{:}; % necessary for PETH
                [aligned_times,aligned_values] = group_continuous_values(params.times,squeeze(stats.score(:,:,i)),stateTimes.(eve), kSpikeWindowS.(eve));
                stats.psth(i).(eve) = calc_psth_continuous(aligned_times,aligned_values, kPETH.timeS.(eve), 'kernel','GAUSS', 'tau',kPETH.stdS.(eve),'gpu',true,varargin{:});
            end
        end    
    end
    
end