function gpfa_new = realign_latents(gpfa,ref_event,time_window_s)
    % takes a gpfa struct (output of format_latents) and realigns to a new
    % reference event (without any psth smoothing, just by shifting each
    % trial by approximately the right number of timepoints)
    gpfa_new = gpfa;
    gpfa_new.ref_event = ref_event;
    gpfa_new.time_window_s = time_window_s;
    gpfa_new.repo_path = fileparts(mfilename('fullpath'));
    [gpfa_new.git_branch,gpfa_new.commit] = return_git_status(gpfa_new.repo_path);    
    gpfa_new.average = 0;
    gpfa_new.time_s = gpfa_new.time_window_s(1):gpfa.resolution_s:gpfa_new.time_window_s(2);
    stateTimes = gpfa.Trials.stateTimes.(ref_event)(gpfa.trial_idx);  
    oldstateTimes = gpfa.Trials.stateTimes.(gpfa.ref_event)(gpfa.trial_idx);            
    gpfa_new.old_time_s = bsxfun(@plus,gpfa.time_s',(oldstateTimes-stateTimes)');  % times of old latents in new reference time 
    gpfa_new.times = bsxfun(@plus,gpfa_new.time_s,stateTimes)';    
    ntrials = numel(gpfa.trial_idx);
    % to do: initialize gpfa._new.score as appropriately sized array of NaNs
    gpfa_new.score = NaN(numel(gpfa_new.time_s),ntrials,gpfa.dim);
    for i=1:ntrials
        if isnan(stateTimes(i))
            continue
        end
        idx_old = gpfa_new.old_time_s(:,i)>=time_window_s(1) & gpfa_new.old_time_s(:,i)<time_window_s(2);
        a(i)=sum(idx_old);
        [~,idx_new] = min(abs(gpfa_new.time_s-gpfa_new.old_time_s(find(idx_old,1),i))); % index in time_s whose value is as close as possible to first elements of old_time_s(idx_old)
        gpfa_new.score((1:sum(idx_old)) + idx_new-1,i,:) = gpfa.score(idx_old,i,:);
    end
end