function stateTimes = select_trials_from_statetimes(stateTimes,correct)
    assert(islogical(correct));
    nt=numel(correct);
    fields = fieldnames(stateTimes);
    for i=1:length(fields)
        nl(i) = numel(stateTimes.(fields{i}));
        if nl(i)==nt
            stateTimes.(fields{i}) = stateTimes.(fields{i})(correct);
        end
    end
    nt_right = nl(strcmp(fields,'right_click_trial'));
    nt_left = nl(strcmp(fields,'left_click_trial'));
    left_fields = fields(nl==nt_left);
    right_fields = fields(nl==nt_right);
    right_idx = ismember(stateTimes.right_click_trial,find(correct));
    left_idx = ismember(stateTimes.left_click_trial,find(correct));
    for i=1:length(left_fields)
        stateTimes.(left_fields{i}) = stateTimes.(left_fields{i})(left_idx);
    end
    for i=1:length(right_fields)
        stateTimes.(right_fields{i}) = stateTimes.(right_fields{i})(right_idx);
    end    
end
