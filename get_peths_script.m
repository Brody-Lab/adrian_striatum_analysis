%% get psths for groups of gammas for all cells in database.
%% generated data for PETH plots for GRS talk

% to do: varargin ref_event, add programmatic control of whether separating
% by gamma. make a separate function for adding side preference to
% cells_table. eventually this should be added to cells table. as should
% explainable variance.

%cells_table = select_cells('is_in_dorsal_striatum',true);

P=get_parameters;

recording_names = unique(cells_table.recording_name);

kPETH = get_PETH_params('std_s',0.1,'std_s_clicks',0.1);    

for i=1:numel(recording_names)
    Cells = load_Cells_file(recording_names{i});
    these_rows = cells_table.recording_name == recording_names{i};
    cellnos = cells_table.cellno(these_rows);
    exclude_trials = validate_trials(Cells.Trials,'mode','agb_glm');
    [pref_left,pval] = get_pref_choice(Cells,cellnos);    
    %[pref_left_stim,stim_pval] = get_pref_stim(Cells,cellnos,'onlyCorrect',false);    
    
    for k=1:(length(P.gamma_ranges)-1)
        trial_idx = find(Cells.Trials.gamma>P.gamma_ranges(k) & Cells.Trials.gamma<P.gamma_ranges(k+1) & ~exclude_trials);        
        psth(:,k) = get_psth_from_Cells(Cells,cellnos,'trial_idx',trial_idx,'kPETH',kPETH,'states',{params.ref_event});
        nan_inds = kPETH.timeS.clicks_on> Cells.Trials.stim_dur_s(trial_idx);     
        if any_nan_inds
            error('why would this ever happen if i am validating trials properly?');
        end
        for c=1:numel(cellnos)
            psth(c,k).clicks_on(nan_inds)=NaN;
            psth(c,k).clicks_on = nanmean(psth(c,k).clicks_on);                        
            if k==4
                this_row = find(these_rows & cells_table.cellno==cellnos(c));
                if numel(this_row)~=1
                    error('');
                end
                cells_table.psth_clicks_on(this_row) = {cat(1,psth(c,:).clicks_on)};
                cells_table.pref_left_MI(this_row) = pref_left(c);
                cells_table.pref_left_MI_pval(this_row) = pval(c);           
                %cells_table.pref_left_stim_MI(this_row) = pref_left_stim(c);
                %cells_table.pref_left_stim_MI_pval(this_row) = stim_pval(c);                   
            end
        end
    end
end

