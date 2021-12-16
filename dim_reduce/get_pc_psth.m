function psth = get_pc_psth(Cells,stats,varargin)
    %% make PSTHs    
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
    for i=1:stats.dim 
        for eves=states(:)'
            eve = eves{:}; % necessary for PETH
            [aligned_times,aligned_values] = group_continuous_values(stats.times,squeeze(stats.score(:,:,i)),stateTimes.(eve), kSpikeWindowS.(eve));
            psth(i).(eve) = calc_psth_continuous(aligned_times,aligned_values, kPETH.timeS.(eve), 'kernel','LGAUSS', 'tau',kPETH.stdS.(eve),'gpu',true,varargin{:});
        end
    end    
end