function Cells = compute_laser_modulation(Cells, varargin)

    %% parse and validate inputs
    tic;
    p=inputParser;
    p.addParameter('bin_size_s',0.001,@(x)validateattributes(x,{'numeric'},{'scalar','positive'}));
    p.addParameter('exclude_violations',true,@(x)validateattributes(x,{'logical'},{'scalar'}));
    p.addParameter('artifact_removal_period_s',[-0.5 4]/1000,@(x)validateattributes(x,{'numeric'},{'increasing','numel',2})); % do not analyze data in this window
    p.addParameter('nboot',5e3,@(x)validateattributes(x,{'numeric'},{'scalar','>',1}));
    p.addParameter('first_sig_bin_pvals',[1e-3 5e-3 1e-2 5e-2],@(x)validateattributes(x,{'numeric'},{'positive','<',1})); 
    p.parse(varargin{:});
    params=p.Results;
    PB_set_constants;

    %% make TPO, structure with info about the times of the laser pulses
    TPO = find_pulse_times(Cells.Trials);  % TPO = times of pulse on
    
    %% initialize output PPTH structure
    PPTH = table2struct(TPO.parameters,'ToScalar',true);
    PPTH.offtimeMS = 1000/PPTH.freqHz - PPTH.pulseMS;
    PPTH.periodMS = 1000/PPTH.freqHz;
    PPTH.n_pulses = (PPTH.durMS-PPTH.offtimeMS)./PPTH.periodMS; % offtime always comes first.
    PPTH.truncated_pulse = rem(PPTH.n_pulses,1)*PPTH.periodMS<PPTH.pulseMS;
    if any(PPTH.truncated_pulse)
        warning('At least one stimulation condition had a truncated final pulse! Results may not be accurate.');
    end
    PPTH.n_pulses = ceil(PPTH.n_pulses);  %round up because the fraction of a period includes a pulse
    
    %% compute bin properties for each laser condition
    n_laser_cond = length(PPTH.offtimeMS);
    for i=1:n_laser_cond
        PPTH.bin_ctr_s{i} = ((-PPTH.pulseMS(i)/1000):params.bin_size_s:(PPTH.pulseMS(i)+PPTH.offtimeMS(i))/1000)-params.bin_size_s/2; % bin center time in s, using params.bin_size_s, for a SINGLE pulse
        if diff(params.artifact_removal_period_s)
            % remove bins of the ppth that overlap the artifact removal
            % period. set artifact_removal_period_s to [0 0] to forgo this step.
            PPTH.bin_ctr_s{i} = PPTH.bin_ctr_s{i} ((PPTH.bin_ctr_s{i}>=params.artifact_removal_period_s(2)+params.bin_size_s/2-eps ...
                | PPTH.bin_ctr_s{i}<=params.artifact_removal_period_s(1)-params.bin_size_s/2+eps) ...
                & (PPTH.bin_ctr_s{i}>=params.artifact_removal_period_s(2)+params.bin_size_s/2+PPTH.pulseMS(i)./1000-eps ...
                | PPTH.bin_ctr_s{i}<=params.artifact_removal_period_s(1)-params.bin_size_s/2+PPTH.pulseMS(i)./1000+eps));
        end  
        PPTH.laser_on_idx{i} = PPTH.bin_ctr_s{i}>=0 & PPTH.bin_ctr_s{i}<PPTH.pulseMS/1000;   % indices of the laser on for a single pulse 
        PPTH.pre_laser_idx{i} = PPTH.bin_ctr_s{i}<0;        
    end
    n_cells=numel(Cells.raw_spike_time_s) ;
    fprintf('\n%s: Making peri-pulse time histogram for %i cells...\n', mfilename, n_cells);
    tt=tic;
    %% loop over cells
    for c = 1:n_cells   
        %% loop over conditions (i.e. diff. stimulation parameters)
        for i=1:n_laser_cond
            ref_event=char(PPTH.triggerEvent(i));
            %% loop over ctrl v. laser
            for cond = {'ctrl', 'laser'}; cond = cond{:};
                % BELOW SECTION CALCULATES PETHS FOR CONTROL AND LASER TRIALS FOR EACH CELL AND LASER CONDITION
                %% bin spikes relative to pulses
                if params.exclude_violations
                    trial_idx = ~cellfun(@isempty, TPO.(cond){i}) & ~Cells.Trials.violated;  
                else
                    trial_idx = ~cellfun(@isempty, TPO.(cond){i});                  
                end
                if strcmp(cond,'ctrl')
                    trial_idx = trial_idx & ~Cells.Trials.laser.isOn;
                end
                PPTH.nTrials{i}.(cond) = sum(trial_idx);
                tpo_cell = num2cell(cell2mat(TPO.(cond){i}(trial_idx)));
                spike_time_s = Cells.spike_time_s.(ref_event){c}(trial_idx);
                spike_time_s = repmat(spike_time_s, 1, PPTH.n_pulses(i));
                ppth.(cond) = cellfun(@(tpo, t_spk) sum(t_spk >=(PPTH.bin_ctr_s{i} - params.bin_size_s/2) + tpo & ...
                                                 t_spk <  (PPTH.bin_ctr_s{i} + params.bin_size_s/2) + tpo, 1)/params.bin_size_s, tpo_cell, spike_time_s, 'uni', 0);     
                ppth_cat_all.(cond) = cat(1,ppth.(cond){:}); % ppth concatenated for all pulses and trials
                ppth.(cond) = cellfun(@(x) cell2mat(x), num2cell(ppth.(cond),1), 'uni', 0); % separate by pulse      
                
                %% average psth's for visualization of just the pulses
                PPTH.single_pulse_peth_avg.(cond){c,i} = cellfun(@(x) mean(x,1), ppth.(cond), 'uni', 0);
                PPTH.single_pulse_peth_std.(cond){c,i} = cellfun(@(x) std(x,[],1), ppth.(cond), 'uni', 0);
                PPTH.combined_pulse_peth_avg.(cond){c,i} = mean(ppth_cat_all.(cond));
                PPTH.combined_pulse_peth_std.(cond){c,i} =  std(ppth_cat_all.(cond));
                
                %% smoothed PETH of the triggered event, separated by stimulation conditions.
                % the artifact window removal is not done for this purely
                % visual peth of the whole reference event. this peth is
                % not used for any analysis.
                tmp=calc_psth(Cells.spike_time_s.(ref_event){c}(trial_idx), kPETH.timeS.(ref_event), kPETH.type.(ref_event), kPETH.stdS.(ref_event));
                PPTH.whole_trial_peth_avg.(cond){c,i} = nanmean(tmp);                  
                PPTH.whole_trial_peth_std.(cond){c,i} = nanstd(tmp);         
                
                %% mean firing rate during laser on times (combined and separately by pulse) for both laser and control trials
                PPTH.single_pulse_FR.(cond){c,i} = cellfun(@(x)mean(x(:,PPTH.laser_on_idx{i}),2),ppth.(cond),'uni',0);
                PPTH.combined_pulse_FR.(cond){c,i} = mean(ppth_cat_all.(cond)(:,PPTH.laser_on_idx{i}),2);                
            end
            %% comparisons between laser and control for each laser condition                     
            %% calculate the significance of the average FR in each bin relative to baseline
            % precalculate the resampling indices since this is
            % computationally expensive and shared across cells for
            % each laser condition
            if c==1
                randsamples{i}=randsample(sum(PPTH.laser_on_idx{i})*PPTH.nTrials{i}.ctrl,PPTH.nTrials{i}.ctrl*params.nboot,true);
                randsamples{i}=reshape(randsamples{i},[PPTH.nTrials{i}.ctrl params.nboot]);
                randsamples2{i}=randsample(sum(PPTH.pre_laser_idx{i})*PPTH.nTrials{i}.laser,PPTH.nTrials{i}.laser*params.nboot,true);
                randsamples2{i}=reshape(randsamples2{i},[PPTH.nTrials{i}.laser params.nboot]);                        
            end
            baseline_samples=ppth.ctrl{1}(:,PPTH.laser_on_idx{i});
            baseline_samples = baseline_samples(:);
            baseline_samples2 = ppth.laser{1}(:,PPTH.pre_laser_idx{i});
            baseline_samples2 = baseline_samples2(:);
            laser_mean_bins  = mean(ppth_cat_all.laser(:,PPTH.laser_on_idx{i}));
            for pulse=1:PPTH.n_pulses(i)
                pulse_mean_bins{pulse}  = mean(ppth.laser{pulse}(:,PPTH.laser_on_idx{i}));                        
            end
            distribution = mean(baseline_samples(randsamples{i}));
            distribution2 = mean(baseline_samples2(randsamples2{i}));             
            for j=1:sum(PPTH.laser_on_idx{i})
                PPTH.bin_pval{c,i}(j) = empirical_p(laser_mean_bins(j),distribution,'low'); %where baseline distribution is control trials
                PPTH.bin_pval2{c,i}(j) = empirical_p(laser_mean_bins(j),distribution2,'low'); % where baseline distribution is pre-laser time on laser trials                        
                for pulse=1:PPTH.n_pulses(i)
                    PPTH.bin_pval_pulse{c,i}{pulse}(j) = empirical_p(pulse_mean_bins{pulse}(j),distribution,'low'); %where baseline distribution is control trials
                    PPTH.bin_pval2_pulse{c,i}{pulse}(j) = empirical_p(pulse_mean_bins{pulse}(j),distribution2,'low'); % where baseline distribution is pre-laser time on laser trials                                              
                end                
            end
            %% latency to first spike after laser pulse onset (not a very useful statistic)
            latency = cellfun(@(tpo,t_spk)t_spk-tpo,tpo_cell(:,1),spike_time_s(:,1),'uni',0);
            latency=cellfun(@(x)min(x(x>0)),latency,'uni',0);                
            latency = cat(1,latency{:});
            PPTH.mean_latency(c,i) = nanmean(latency);
            PPTH.std_latency(c,i) = nanstd(latency);                
            %% statistical comparisons between firing rates when laser on and off
            % this is calculated two ways: 
                % 1. comparing control v. laser trials during matched times
                % 2. comparing pre- v. post- laser on times on just stimulation trials
            % I had previously broken this analysis down by bins,
            % but this did not add much information and was
            % implemented in a hacky way.
            on_means = mean(ppth_cat_all.laser(:,PPTH.laser_on_idx{i}),2);          
            off_idx = PPTH.bin_ctr_s{i}>-PPTH.pulseMS(i)/1000 ;                    
            off_means = mean(ppth.laser{1}(:,off_idx),2);% always use the pre-laser period from the first pulse, because things get weird between the pulses                    
            PPTH.combined_pulse_stats{c,i}.prepost= compute_paired_stats(mean(reshape(on_means,numel(on_means)./PPTH.n_pulses,PPTH.n_pulses),2),off_means);  
            PPTH.combined_pulse_stats{c,i}.ctrlvlaser = compute_paired_stats(PPTH.combined_pulse_FR.laser{c,i},PPTH.combined_pulse_FR.ctrl{c,i});                                        
            for pulse=1:PPTH.n_pulses(i)
                on_means = mean(ppth.laser{pulse}(:,PPTH.laser_on_idx{i}),2);
                PPTH.single_pulse_stats{c,i}{pulse}.prepost = compute_paired_stats(on_means,off_means);   
                PPTH.single_pulse_stats{c,i}{pulse}.ctrlvlaser = compute_paired_stats(PPTH.single_pulse_FR.laser{c,i}{pulse},PPTH.single_pulse_FR.ctrl{c,i}{pulse});
            end
            fields = {'bin_pval','bin_pval2','bin_pval_pulse','bin_pval2_pulse'};
            for f = 1:length(fields)
                if ~contains(fields{f},'pulse')
                    PPTH.first_sig_time_s{c,i}.(fields{f}) = get_first_sig_time(PPTH.(fields{f}){c,i},params.first_sig_bin_pvals,PPTH.bin_ctr_s{i}(PPTH.laser_on_idx{i}));
                else
                    for k=1:length(PPTH.(fields{f}){c,i})
                        PPTH.first_sig_time_s{c,i}.(fields{f}){k} = get_first_sig_time(PPTH.(fields{f}){c,i}{k},params.first_sig_bin_pvals,PPTH.bin_ctr_s{i}(PPTH.laser_on_idx{i}));            
                    end
                end
            end               
        end
    end  
    Cells.PPTH=PPTH;
    Cells.PPTH.params=params;
    fprintf(' took %.0f seconds.', toc(tt));
end

function TPO = find_pulse_times(Trials)
    %% Make Conditions - structure with information about the different stimulation conditions
    Conditions = PB_make_Conditions(Trials);
    idx_ctrl = ~Conditions.parameters.isOn;
    if sum(idx_ctrl) ~= 1
        error('This code cannot yet accommodate multiple control conditions')
    end
    % set the control condition to be first
    Conditions.parameters = [Conditions.parameters(idx_ctrl,:); Conditions.parameters(~idx_ctrl,:)];
    Conditions.trial_index = [Conditions.trial_index(idx_ctrl,:); Conditions.trial_index(~idx_ctrl,:)];
    %% make TPO, structure with info about the times of the laser pulses
    TPO = struct;   % times of pulse on
    for c = 1:(Conditions.n-1)
        TPO.parameters(c,:) = Conditions.parameters(c+1,:);
        idx = Conditions.trial_index{c+1};
        T_laser_pulses = Trials.stateTimes.laser_pulses;
        T_laser_pulses(~idx) = {[]};
        T_laser_pulses(idx) = cellfun(@(x) x(:,1)', T_laser_pulses(idx), 'uni', 0); % take just the pulse on times
        ref_event = char(TPO.parameters.triggerEvent{c});
        T_ref_event = num2cell(Trials.stateTimes.(ref_event)(idx));
        T_laser_pulses(idx) = cellfun(@(tpo, t_ref) tpo-t_ref, T_laser_pulses(idx), T_ref_event, 'uni', 0);
        TPO.laser{c,1} = T_laser_pulses;
        TPO.ctrl{c,1} = cell(numel(Trials.is_hit),1);
        TPO.ctrl{c,1}(Conditions.trial_index{1}) = {mean(cell2mat(TPO.laser{c}))};
    end
end

function stats = compute_paired_stats(on,off)
% computes statistics for paired numeric vectors (i.e. counts during laser on and off times during matches sets of trials)
    mi_fun = @(x,y)(x-y)./(x+y);
    if numel(on)~=numel(off)
        stats.ranksum = ranksum(on,off,'tail','right');
        [~,stats.tp] = ttest2(on,off,'tail','right');   % p-value of t-test                  
    else % these required paired samples
        stats.reliability = mean(on>off);  
        stats.signrank = signrank(on,off,'tail','right');
        [~,stats.tp] = ttest(on,off,'tail','right');   % p-value of t-test          
    end
    stats.auc = CalcCP(off,on);                    
    stats.mi = mi_fun(mean(on),mean(off));                    
    stats.dp = dprime(off,on);       
end

function first_sig_time = get_first_sig_time(bin_pval,pvals,bin_times)
    for i=1:length(pvals)
        tmp=find(bin_pval<pvals(i),1,'first');        
        if isempty(tmp)
            first_sig_time(i) = NaN;
        else
            first_sig_time(i) = bin_times(tmp);
        end
    end
end