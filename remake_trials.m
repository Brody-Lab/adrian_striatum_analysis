function [mismatch,click_mismatch] = remake_trials(Cells,Trials,varargin)
    % remake all trials structures in database based on new code version of
    % PB_process_data written 6/2022 and redo spikes.
    
    p=inputParser;
    p.addParameter('save',false);
    p.parse(varargin{:});
    params=p.Results;
    
    PB_set_constants;
    paths = get_data_paths();
    
    %% compare old and new trials structs and load data if necessary for the first time load
    for i=1:numel(paths)
        %fprintf('\nFile %g.\n',i);
        Cells{i} = load(paths(i).cells_file);
        %Trials{i} = PB_process_data(paths(i).rat_name,datestr(paths(i).date,'yyyy_mm_dd'));
        %[remove{i},remove2{i}] = compare_old_new_trials(Cells{i}.Trials,Trials{i});  
    end
    
    %% calculate and check new aligned spike times and click times
    for i=1:numel(paths)
        these_trials = remove_trials(Trials{i},remove{i});
        Cells{i}.Trials = remove_trials(Cells{i}.Trials,remove2{i});
        Cells{i}.spike_time_s = remove_trials_from_spikes(Cells{i}.spike_time_s,remove2{i});
        nclust=numel(Cells{i}.raw_spike_time_s);
        nt= numel(these_trials.gamma);
        fprintf('\n\ncells file %g\n---------',i);  
        
        for t=1:nt
            old_left_bups = Cells{i}.Trials.stateTimes.left_clicks(Cells{i}.Trials.stateTimes.left_click_trial==t);        
            old_right_bups = Cells{i}.Trials.stateTimes.right_clicks(Cells{i}.Trials.stateTimes.right_click_trial==t);                    
            left_bups = these_trials.stateTimes.clicks_on(t) + these_trials.leftBups{t};
            right_bups = these_trials.stateTimes.clicks_on(t) + these_trials.rightBups{t};            
            if any(~ismember(left_bups,old_left_bups)) || any(~ismember(right_bups,old_right_bups))
                click_mismatch{i}(t)=1;
            else
                click_mismatch{i}(t)=0;
            end
        end
        % make spike times
        for events = kPETH.refEvents; eve = events{:}; % necessary for PETH
            spike_time_s = cell(nclust,1);
            if ~isfield(Cells{i}.Trials.stateTimes,eve)
                %warning('kPETH refEvent %s is not a field of Cells.Trials.stateTimes. Skipping making spike times.',eve); 
                continue
            end            
            Cells{i}.kSpikeWindowS.(eve) = kSpikeWindowS.(eve);        
            window_s = Cells{i}.kSpikeWindowS.(eve);                            
            event_times = these_trials.stateTimes.(eve);
            times=Cells{i}.raw_spike_time_s;
            if isfield(Cells{i}.spike_time_s,eve)
                aligned_times =  Cells{i}.spike_time_s.(eve);
            else
                aligned_times = cell(numel(Cells{i}.raw_spike_time_s),1);
            end
            fprintf('\n%s',eve);
            m=NaN(nclust,nt);
            parfor c = 1:nclust
                spike_time_s{c,1} = group_spike_times(times{c}, event_times, window_s);
                if ~isempty(aligned_times{c})                    
                    for k=1:nt
                        if isempty(aligned_times{c}{k})
                            if ~isempty(spike_time_s{c,1}{k})
                                m(c,k) = NaN;
                            else
                                m(c,k)=NaN;
                            end
                        else
                            if isempty(spike_time_s{c,1}{k})
                                m(c,k) = NaN;
                            else                            
                                same_frac = mean(ismember(round(aligned_times{c}{k},3),round(spike_time_s{c,1}{k},3)));
                                m(c,k) = same_frac;
                            end
                        end
                    end
                end
            end
            Cells{i}.spike_time_s.(eve) = spike_time_s;
            Cells{i}.Trials = these_trials;
            mismatch.(eve)(i) = nanmean(m(:));
            fprintf(' %g',mismatch.(eve)(i));
        end
    end
    
    if params.save
        for i=1:numel(paths)
            tmp=Cells{i};
            save(paths(i).cells_file,'-struct','tmp');        
        end
    end



end


function [remove_trials,remove_trials2] = compare_old_new_trials(old_trials,new_trials)
    remove_trials = ~ismember(round(new_trials.stateTimes.sending_trialnum,3),round(old_trials.stateTimes.sending_trialnum,3));
    remove_trials2 = ~ismember(old_trials.stateTimes.sending_trialnum,new_trials.stateTimes.sending_trialnum);     
    fprintf('%g trials in new trials structure where rat never cpoked.\n',sum(new_trials.never_cpoked));
    fprintf('%g trials missing from old trials structure. %g trials missing from new trials structure.\n',sum(remove_trials),sum(remove_trials2));
    fields={'gamma','sides','is_hit','reward_loc','trial_type','pokedR','violated'};
    for f=1:numel(fields)
        a=new_trials.(fields{f})(~remove_trials);
        b=old_trials.(fields{f})(~remove_trials2);
        if ~all(a==b | (isnan(a)&isnan(b)))
            nans= isnan(a) | isnan(b);
            if strcmp(fields{f},'pokedR') && any(abs(a(~nans)-b(~nans)))
                error('');
            end
            fprintf('%g%% for %s identical.\n',round(mean(abs(a(~nans)-b(~nans))<5e-4)*100,1),fields{f});
            if strcmp(fields{f},'violated')
                if any(a(~nans)&~b(~nans))
                    fprintf('%g violated in new but not old.\n',sum(a(~nans)&~b(~nans))); % i should have fewer vioations now because of bug in state matrix in 2021
                end
            end
            fprintf('%g%% where only one is NaN.\n',100*mean(xor(isnan(a),isnan(b))));
        end
    end     
    fields={'wait_for_cpoke','cpoke_in','break','error','left_reward','right_reward','clicks_off','cpoke_out','spoke','wait_for_spoke','clicks_on'};
    for f=1:numel(fields)
        a=new_trials.stateTimes.(fields{f})(~remove_trials);
        b=old_trials.stateTimes.(fields{f})(~remove_trials2);
        if ~all(abs(a-b)<5e-4 | (isnan(a)&isnan(b)))
            fprintf('\nrobust mean abs difference for %s is %g seconds.\n',fields{f},trimmean(abs(a-b),5));
            nans= isnan(a) | isnan(b);
            fprintf('%g%% identical.\n',round(mean(abs(a(~nans)-b(~nans))<5e-4)*100,1));
            fprintf('%g%% where only one is NaN.\n',100*mean(xor(isnan(a),isnan(b))));          
            if strcmp(fields{f},'cpoke_out')
                if any(a(~nans)&~b(~nans))
                    warning('%g nan cpoke in new but not old!!\n',sum(a(~nans)&~b(~nans))); % i should have fewer nan cpokes now because they are includes on violations trials now
                end
            end            
        end
    end      
end