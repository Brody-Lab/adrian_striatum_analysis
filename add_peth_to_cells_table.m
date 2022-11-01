function cells_table = add_peth_to_cells_table(cells_table,varargin)


    % improved version of "get_peths_script" that was used to prepare psths
    % by area for Basal Ganglia GRS in early 2022. this version written
    % 10/2022.

    % to do: varargin ref_event, add programmatic control of whether separating
    % by gamma. make a separate function for adding side preference to
    % cells_table. eventually this should be added to cells table. 
    % programmatic control of new column name
    
    %% parse inputs
    p=inputParser;
    p.KeepUnmatched=true;
    p.addParameter('ref_event','clicks_on',@(x)validateattributes(x,{'char'},{'nonempty'}));    
    p.addParameter('separate_by','',@(x)validateattributes(x,{'char'},{}));
    p.addParameter('column_name','',@(x)validateattributes(x,{'char'},{}));
    p.addParameter('kPETH',get_PETH_params('std_s',0.1,'std_s_clicks',0.1));
    % you can pass mask_states directly to get_psth_from_Cells
    p.parse(varargin{:});
    params=p.Results;
    if isempty(params.column_name)
        params.column_name = params.ref_event;
        warning('Using "%s" (the ref_event) as the new column name by default.',params.ref_event);
    end
    recording_names = unique(cells_table.recording_name);
    for i=numel(recording_names):-1:1
        these_rows{i} = find(cells_table.recording_name == recording_names{i});
        cellnos{i} = cells_table.cellno(cells_table.recording_name == recording_names{i});
    end
    for i=1:numel(recording_names)
        Cells = load_Cells_file(recording_names{i});
        exclude_trials = validate_trials(Cells.Trials,'mode','agb_glm');
        [pref_left{i},pval{i}] = get_pref_choice(Cells,cellnos{i},'exclude_trials',exclude_trials);
        %[pref_left_stim{i},stim_pval{i}] = get_pref_stim(Cells,cellnos,'onlyCorrect',false);
        psths{i} = get_peth_for_cells_table(Cells,cellnos{i},params,exclude_trials,varargin{:},'use_gpu',false);
    end
    for i=1:numel(recording_names)
        cells_table.(params.column_name)(these_rows{i}) = psths{i};
        cells_table.pref_left_MI(these_rows{i}) = pref_left{i};
        cells_table.pref_left_MI_pval(these_rows{i}) = pval{i};           
        %cells_table.pref_left_stim_MI(these_rows{i}) = pref_left_stim;
        %cells_table.pref_left_stim_MI_pval(these_rows{i}) = stim_pval;                   
    end
end

function psths = get_peth_for_cells_table(Cells,cellnos,params,exclude_trials,varargin)
    P=get_parameters;
    switch params.separate_by
        case 'gamma'
            for k=1:(length(P.gamma_ranges)-1)
                trial_idx = find(Cells.Trials.gamma>P.gamma_ranges(k) & Cells.Trials.gamma<P.gamma_ranges(k+1) & ~exclude_trials);        
                psth(:,k) = get_psth_from_Cells(Cells,cellnos,'trial_idx',trial_idx,'kPETH',params.kPETH,...
                    'states',{params.ref_event},'return_trials',false,varargin{:});
            end
            for c=1:numel(cellnos)
                psths{c} = cat(1,psth(c,:).(params.ref_event));
            end
    end
end

