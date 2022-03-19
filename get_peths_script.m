%% get psths for groups of gammas for all cells in database.

%cells_table = select_cells('is_in_dorsal_striatum',true);

P=get_parameters;

recording_names = unique(cells_table.recording_name);

gamma_ranges = [-5 -2 0 2 5];

kPETH = get_PETH_params('std_s',0.1,'std_s_clicks',0.1);    

for i=1:numel(recording_names)
    Cells = load_Cells_file(recording_names{i});
    these_rows = cells_table.recording_name == recording_names{i};
    cellnos = cells_table.cellno(these_rows);
    
    
    exclude_trials = Cells.Trials.violated;
    exclude_trials = exclude_trials | Cells.Trials.laser.isOn;
    if isfield(Cells.Trials,'stim_dur_s_theoretical')
        stim_dur = Cells.Trials.stim_dur_s_theoretical;
    else
        stim_dur = Cells.Trials.stim_dur_s;
    end    
    if any(isnan(stim_dur))
        warning('Found %g accumulation trials with stimdur of NaN. Removing them. (Probably violations or uninitiated trials).',sum(isnan(stim_dur)));
    end
    exclude_trials = exclude_trials | isnan(stim_dur);    
    %exclude_trials = exclude_trials | Cells.Trials.is_hit~=1;
    
    if isfield(Cells.Trials,'stim_dur_s_actual')
        stim_dur=Cells.Trials.stim_dur_s_actual;
    else
        stim_dur=Cells.Trials.stim_dur_s;
    end
    [pref_left,pval] = get_pref_choice(Cells,cellnos);    
    %[pref_left_stim,stim_pval] = get_pref_stim(Cells,cellnos,'onlyCorrect',false);    
    
    clear psth psths
    for k=1:4
        trial_idx = find(Cells.Trials.gamma>gamma_ranges(k) & Cells.Trials.gamma<gamma_ranges(k+1) & ~exclude_trials);        
        psth(:,k) = get_psth(Cells,cellnos,'trial_idx',trial_idx,'kPETH',kPETH,'states',{'clicks_on'});
        nan_inds = kPETH.timeS.clicks_on> stim_dur(trial_idx);        
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


for i=1:numel(P.ap_groups)
    is_in_group = cells_table.AP>=P.ap_groups{i}(1) & cells_table.AP<P.ap_groups{i}(2);
    cells_table.ap_group(is_in_group) = i;
end

