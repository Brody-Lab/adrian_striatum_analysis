function [choice_MI,choice_MI_pval] = get_pref_choice(Cells,cellno,varargin)
    % returns a choice modulation index for cells in a recording session.
    % choice modulation index is (R-L)/(R+L) where R and L are avg. spike counts on right- and
    % left-choice trials, respectively. spikes are counted in the period
    % between clicks_on and cpoke_req_end. By default, only correct trials
    % are used.
    
    %% parse inputs
    p = inputParser;
    p.KeepUnmatched=true;
    p.addParameter('nresamples',1000); 
    p.addParameter('onlyCorrect',true);
    p.addParameter('exclude_trials',logical([]),@(x)validateattributes(x,{'logical'},{})); % should be logical
    p.addParameter('cpoke_out_window',0); % exclude this many seconds of the end of the stimulus period (to remove influence of movement prep)
    p.parse(varargin{:});
    params = p.Results;
    exclude_trials = params.exclude_trials;
    if isempty(exclude_trials)
        exclude_trials = validate_trials(Cells.Trials,'mode','agb_glm');
    end
    if params.onlyCorrect
       exclude_trials = exclude_trials | Cells.Trials.is_hit~=1;
    end
    trial_idx = ~exclude_trials;
    trial_idx = find(trial_idx);
    MIfun = @(x,y)(x-y)./(x+y);
    max_time_s = Cells.Trials.stateTimes.cpoke_out - Cells.Trials.stateTimes.first_click - params.cpoke_out_window;    
    these_trials{1} = intersect(trial_idx,find(Cells.Trials.pokedR==1));
    these_trials{2} = intersect(trial_idx,find(Cells.Trials.pokedR==0));            
    n = cellfun(@numel,these_trials);      
    for c=numel(cellno):-1:1
        for i=1:2
            count=0;
            clear spike_counts
            for t=these_trials{i}(:)'
                count=count+1;
                these_times = Cells.spike_time_s.first_click{cellno(c)}{t};
                spike_counts(count) = sum(these_times>0.1 & these_times<max_time_s(t)); % exclude first 100ms of stimulus response
            end
            for k=1:params.nresamples
                boots(k,i) = sum(spike_counts(randi(n(i),n(i),1)));
            end
        end
        boots=boots./n;
        MI = MIfun(boots(:,1),boots(:,2)); 
        choice_MI_pval(c) = empirical_p(0,MI,'both');
        choice_MI(c) = nanmean(MI);                  
    end
end