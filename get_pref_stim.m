function [stim_MI,stim_MI_pval] = get_pref_stim(Cells,cellno,varargin)
    % aligns to cpoke_in and excludes incorrect trials and times within a
    % certain range of cpoke out
    p = inputParser;
    p.KeepUnmatched=true;
    p.addParameter('nresamples',1000); 
    p.addParameter('onlyCorrect',false);
    p.addParameter('cpoke_out_window',0.1); % don't include times within 0.1 seconds of the cpoke out
    p.parse(varargin{:});
    params = p.Results;
    % only include desired trials in psth calculation

    
    exclude_trials = Cells.Trials.violated;
    exclude_trials = exclude_trials | Cells.Trials.laser.isOn;
    if params.onlyCorrect
       exclude_trials = exclude_trials | Cells.Trials.is_hit~=1;
    end
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
    gamma_ranges = [-5 0 5];    
    %% make PSTHs    
    for c=1:numel(cellno)    
        for i=1:2
            these_trials = intersect(trial_idx,find(Cells.Trials.gamma>gamma_ranges(i) & Cells.Trials.gamma<gamma_ranges(i+1)));                    
            psth=calc_psth(Cells.spike_time_s.(eve){cellno(c)}(these_trials),0.1:0.05:0.7,'LGAUSS',0.1);
            nan_inds = ([0.1:0.05:0.7]+params.cpoke_out_window) > stim_dur(these_trials);        
            psth(nan_inds)=NaN;
            psth = nanmean(psth,2);
            boots(:,i) = bootstrp(params.nresamples,@nanmean,psth);                 
        end
        MI = MIfun(boots(:,1),boots(:,2));
        stim_MI_pval(c) = empirical_p(0,MI,'both');
        stim_MI(c) = mean(MI);        
    end
end