function [choice_MI,choice_MI_pval] = get_pref_choice(Cells,cellno,varargin)
    % aligns to cpoke_in and excludes incorrect trials and times within a
    % certain range of cpoke out
    p = inputParser;
    p.KeepUnmatched=true;
    p.addParameter('nresamples',1000); 
    p.addParameter('onlyCorrect',true);
    p.addParameter('cpoke_out_window',0.2); % don't include times within 0.2 seconds of the cpoke out
    p.parse(varargin{:});
    params = p.Results;
    % only include desired trials in psth calculation


     exclude_trials = Cells.Trials.violated;
    exclude_trials = exclude_trials | Cells.Trials.laser.isOn;
    if isfield(Cells.Trials,'stim_dur_s_theoretical')
        stim_dur = Cells.Trials.stim_dur_s_theoretical;
    else
        stim_dur = Cells.Trials.stim_dur_s;
    end    
    if any(isnan(stim_dur))
        warning('Found %g accumulation trials with stimdur of NaN. Removing them. (Probably violations or uninitiated trials).',sum(isnan(stim_dur)));
    end
    exclude_trials = exclude_trials | isnan(stim_dur);    
    trial_idx = ~exclude_trials;
    if islogical(trial_idx)
       trial_idx = find(trial_idx);
    end    
    MIfun = @(x,y)(x-y)./(x+y);
    eve='clicks_on';      
    %% make PSTHs    
    for c=1:numel(cellno)    
        for i=1:2
        if i==1
            these_trials = intersect(trial_idx,Cells.Trials.pokedR==1);
        else
            these_trials = intersect(trial_idx,Cells.Trials.pokedR==0);            
        end
            psth=calc_psth(Cells.spike_time_s.(eve){cellno(c)}(these_trials),0.2:0.1:1,'LGAUSS',0.1);
            nan_inds = kPETH.timeS.clicks_on+params.cpoke_out_window > stim_dur(these_trials);        
            psth(nan_inds)=NaN;
            psth = nanmean(psth,2);
            boots(:,i) = bootstrp(params.nresamples,@nanmean,psth(c).(eve));                 
        end
        MI = MIfun(boots(:,1),boots(:,2));
        choice_MI_pval(c) = empirical_p(0,MI,'both');
        choice_MI(c) = mean(MI);        
    end
end