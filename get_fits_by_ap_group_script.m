%% script for breaking down fits by ap group - 3/1/22

P = get_parameters();
[glmfit_log,run] = select_glmfit_runs('distribution','poisson','lambda',30,'use_trial_history',false);
cells_table = select_cells('is_in_dorsal_striatum',true);
[fits,is_responsive] = get_glm_fits(cells_table.recording_name,cells_table.cellno,glmfit_log);
fits_table = [cells_table(is_responsive,:) fits];
% for i=1:numel(P.ap_groups)
%     is_in_group = fits_table_acausal.AP>=P.ap_groups{i}(1) & fits_table_acausal.AP<P.ap_groups{i}(2);
%     fits_table_acausal.ap_group(is_in_group) = i;
% end

P = get_parameters();
[glmfit_log_no_hist,run] = select_glmfit_runs('choice_time_back_s',0.5,'use_trial_history',false,'lambda',30);
cells_table = select_cells('is_in_dorsal_striatum',true);
[fits,is_responsive] = get_glm_fits(cells_table.recording_name,cells_table.cellno,glmfit_log_no_hist);
fits_table_no_hist = [cells_table(is_responsive,:) fits];
% for i=1:numel(P.ap_groups)
%     is_in_group = fits_table_acausal.AP>=P.ap_groups{i}(1) & fits_table_acausal.AP<P.ap_groups{i}(2);
%     fits_table_acausal.ap_group(is_in_group) = i;
% end

P = get_parameters();
[glmfit_log,run] = select_glmfit_runs('choice_time_back_s',0.5,'use_trial_history',true,'lambda',30);
cells_table = select_cells('is_in_dorsal_striatum',true);
[fits,is_responsive] = get_glm_fits(cells_table.recording_name,cells_table.cellno,glmfit_log);
fits_table = [cells_table(is_responsive,:) fits];
% for i=1:numel(P.ap_groups)
%     is_in_group = fits_table_acausal.AP>=P.ap_groups{i}(1) & fits_table_acausal.AP<P.ap_groups{i}(2);
%     fits_table_acausal.ap_group(is_in_group) = i;
% end