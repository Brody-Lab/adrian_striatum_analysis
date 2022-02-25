function Trials = add_click_states_to_thomas_files(Trials)
    % n.b. in thomas' files clicks_on is the start of the sonud wave, not
    % the time of the first actual click like in my files. 
    nt = numel(Trials.is_hit);
    [Trials.stateTimes.left_clicks,Trials.stateTimes.left_click_times,Trials.stateTimes.left_click_trial,...
        Trials.stateTimes.right_clicks,Trials.stateTimes.right_click_times,Trials.stateTimes.right_click_trial] = deal([]);        
    for i=1:nt
        % left clicks
        nl = numel(Trials.leftBups{i});
        Trials.stateTimes.left_clicks(end+1:end+nl) = Trials.leftBups{i} + Trials.stateTimes.clicks_on(i);
        Trials.stateTimes.left_click_times(end+1:end+nl) = Trials.leftBups{i};        
        Trials.stateTimes.left_click_trial(end+1:end+nl) = i;
        % right clicks
        nr = numel(Trials.rightBups{i});
        Trials.stateTimes.right_clicks(end+1:end+nr) = Trials.rightBups{i} + Trials.stateTimes.clicks_on(i);
        Trials.stateTimes.right_click_times(end+1:end+nr) = Trials.rightBups{i};
        Trials.stateTimes.right_click_trial(end+1:end+nr) = i;        
    end
    Trials = add_shuffled_click_times(Trials);
end
