function aligned_pcs = get_aligned_pcs(Cells,varargin)
    % function that does PCA on a neural population in a Cells file, aligned to
    % specific trial events. the output has the same format as a PETH
    % structure. options include criteria for including units, reference event,
    % time window and time resolution, number of PCs, etc.
    
    % requires npx_utils repo

   
    % for optional input arguments see select_units, MakedDataMatrix and
    % run_PCA to which all optional input arguments are passed.
    
    % make data matrix
    
    % to do -- some of the output times are NaN (presumably 
    
    PB_set_constants;
    
    [cells_mat,params] = MakeDataMatrix(Cells,varargin{:},'units',select_units(Cells,varargin{:}));

    
    
    % perform PCA
    pca_stats = run_PCA(cells_mat,varargin{:});
    npcs = size(pca_stats.score,2);
    % make PSTHs
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
    trials=true(Cells.nTrials,1);
    stateTimes = Cells.Trials.stateTimes;
    % uncomment lines to get psths for a particular subset of trials (like
    % hits or misses)
    %trials = ~Cells.Trials.is_hit;    
    %stateTimes = select_trials_from_statetimes(Cells.Trials.stateTimes,trials);
    for i=1:npcs
        aligned_pcs(i) = import_continuous_variable(params.times,pca_stats.score(:,i),stateTimes,states,kSpikeWindowS,kPETH);
    end
end