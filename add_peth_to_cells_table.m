function [cells_table,kPETH] = add_peth_to_cells_table(cells_table,varargin)


    % improved version of "get_peths_script" that was used to prepare psths
    % by area for Basal Ganglia GRS in early 2022. this version written
    % 10/2022.

    % to do: make a separate function for adding side preference to
    % cells_table.
    
    %% parse inputs
    p=inputParser;
    p.KeepUnmatched=true;
    p.addParameter('ref_event','first_click',@(x)validateattributes(x,{'char'},{'nonempty'}));    
    p.addParameter('separate_by','gamma',@(x)validateattributes(x,{'char'},{}));
    p.addParameter('column_name','',@(x)validateattributes(x,{'char'},{}));
    p.addParameter('kPETH',get_PETH_params('std_s',0.1,'std_s_clicks',0.1));
    % you can pass mask_states directly to get_psth_from_Cells
    p.parse(varargin{:});
    params=p.Results;
    kPETH=params.kPETH;
    if isempty(params.column_name)
        params.column_name = [params.ref_event,'_peth'];
        warning('Using "%s" as the new column name by default.',params.column_name);
    end
    recording_names = unique(cells_table.recording_name);
    for i=numel(recording_names):-1:1
        these_rows{i} = find(cells_table.recording_name == recording_names{i});
        cellnos{i} = cells_table.cellno(cells_table.recording_name == recording_names{i});
    end
    parfor i=1:numel(recording_names)
        Cells = load_Cells_file(recording_names{i});
        Cells = add_first_click_state(Cells);
        Cells.Trials.signal_strength = calc_signal_strength(Cells.Trials.rightBups,Cells.Trials.leftBups,...
            Cells.Trials.stateTimes.cpoke_req_end - Cells.Trials.stateTimes.first_click);        
        exclude_trials = validate_trials(Cells.Trials,'mode','agb_glm','quiet',true,'require_clicks',true);
        [choice_MI{i},choice_MI_pval{i}] = get_pref_choice(Cells,cellnos{i},'exclude_trials',exclude_trials);
        %[pref_left_stim{i},stim_pval{i}] = get_pref_stim(Cells,cellnos,'onlyCorrect',false);
        psths{i} = get_peth_for_cells_table(Cells,cellnos{i},params,exclude_trials,varargin{:},'use_gpu',false);
    end
    for i=1:numel(recording_names)
        cells_table.(params.column_name)(these_rows{i}) = psths{i};
        cells_table.choice_MI(these_rows{i}) = choice_MI{i};
        cells_table.choice_MI_pval(these_rows{i}) = choice_MI_pval{i};           
        %cells_table.pref_left_stim_MI(these_rows{i}) = pref_left_stim;
        %cells_table.pref_left_stim_MI_pval(these_rows{i}) = stim_pval;                   
    end
end

function psths = get_peth_for_cells_table(Cells,cellnos,params,exclude_trials,varargin)
    P=get_parameters;
    psth = get_psth_from_Cells(Cells,params.ref_event,cellnos,'trial_idx',~exclude_trials,'kPETH',params.kPETH,...
        'return_trials',true,varargin{:});    
    switch params.separate_by
        case 'gamma'
            gamma = Cells.Trials.gamma(~exclude_trials);
            for k=1:(length(P.gamma_ranges)-1)
                trial_idx = gamma>P.gamma_ranges(k) & gamma<P.gamma_ranges(k+1);       
                for c=1:numel(cellnos)
                    psths{c}(k,:) = nanmean(psth{c}(trial_idx,:));
                end
            end
        case 'signal_strength'
            signal_strength = Cells.Trials.signal_strength(~exclude_trials);            
            bin_edges = [-Inf -25 0 25 Inf];%prctile(signal_strength,P.signal_strength_percentiles)
            signal_strength = discretize(signal_strength,bin_edges);
            for k=1:(length(P.signal_strength_percentiles)-1)
                trial_idx = signal_strength==k;    
                for c=1:numel(cellnos)
                    psths{c}(k,:) = nanmean(psth{c}(trial_idx,:));
                end
            end            
        case 'otherwise'
            % this is where you'd add cases for other types of conditions
            % to separate trials
    end
end



