function varargout = get_pc_psth(stats,varargin)

    % can resample data (all timepoints across all trials) to estimate
    % error bars. first output is the mean, second output is the s.d.
    
    p = inputParser;
    p.KeepUnmatched=true;
    p.addParameter('trial_idx',~stats.Trials.violated,@(x)validateattributes(x,{'logical','numeric'},{'vector'}));  
    p.addParameter('nresamples',1); 
    p.addParameter('states',{},@(x)validateattributes(x,{'cell'},{}));
    p.addParameter('kPETH',get_PETH_params());
    p.addParameter('gpu',true);
    p.parse(varargin{:});
    params = p.Results;
    kPETH = params.kPETH; 
    % only include desired trials in psth calculation
    if islogical(params.trial_idx)
        params.trial_idx = find(params.trial_idx);
    end
    idx = ismember(stats.trial_idx,params.trial_idx);
    stats.score = stats.score(:,idx,:);
    stats.trial_idx = stats.trial_idx(idx);
    stats.times = stats.times(:,idx);

    if nargout<2
        params.nresamples=1;
    elseif params.nresamples==1
        error('Cannot assign two outputs with only one resample.');
    end
    
    %% make PSTHs    
    if ismember('states',p.UsingDefaults)
        states=kPETH.refEvents;
        states=union(states,{'shuffled_left_clicks','shuffled_right_clicks'});
    else
        states = params.states;
    end
    stats.kSpikeWindowS.shuffled_left_clicks = stats.kSpikeWindowS.left_clicks;
    stats.kSpikeWindowS.shuffled_right_clicks = stats.kSpikeWindowS.right_clicks;        
    for eve=string(states(:)')    
        for i=stats.dim:-1:1 
            if i==stats.dim
                if kPETH.confineToCpoke.(eve)
                    if contains(eve,'left_click')
                        clicks_off = stats.Trials.stateTimes.cpoke_in(stats.Trials.stateTimes.left_click_trial)+1.5;
                    elseif contains(eve,'right_click')
                        clicks_off = stats.Trials.stateTimes.cpoke_in(stats.Trials.stateTimes.right_click_trial)+1.5;                    
                    else
                        clicks_off = stats.Trials.stateTimes.cpoke_in+1.5;
                    end
                    interval = [repmat(stats.kSpikeWindowS.(eve)(1),numel(clicks_off),1) min(clicks_off - stats.Trials.stateTimes.(eve)(:) ,stats.kSpikeWindowS.(eve)(2))];                                
                else
                    interval=stats.kSpikeWindowS.(eve);
                end
                [aligned_times,aligned_values,idx] = group_continuous_values(stats.times,squeeze(stats.score(:,:,i)),stats.Trials.stateTimes.(eve), interval);
                [psth(i).(eve),kernel_weights] = calc_psth_continuous(aligned_times,aligned_values, kPETH.timeS.(eve), 'kernel',kPETH.type.(eve), ...
                    'tau',kPETH.stdS.(eve),'gpu',params.gpu,'nresamples',params.nresamples);                
            else
                v = squeeze(stats.score(:,:,i));
                v=v(:);                
                for k=1:length(idx) 
                    aligned_values{k} = v(idx{k});
                end
                psth(i).(eve) = calc_psth_continuous(aligned_times,aligned_values, kPETH.timeS.(eve), 'kernel',kPETH.type.(eve), ...
                    'tau',kPETH.stdS.(eve),'gpu',params.gpu,'nresamples',params.nresamples,'kernel_weights',kernel_weights);                
            end
            if params.nresamples>1
                psth_std(i).(eve) = std(psth(i).(eve));
                psth(i).(eve) = mean(psth(i).(eve));                
            end
        end
    end    
    if nargout
        varargout{1} = psth;
        if nargout>1
            varargout{2} = psth_std;
        end
    end
end

