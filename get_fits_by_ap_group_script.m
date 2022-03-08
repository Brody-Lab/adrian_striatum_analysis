%% script for breaking down fits by ap group - 3/1/22

P = get_parameters();
[glmfit_log,run] = select_glmfit_runs('choice_time_back_s',0.75);
cells_table = select_cells('is_in_dorsal_striatum',true);
[fits,is_responsive] = get_glm_fits(cells_table.recording_name,cells_table.cellno,glmfit_log);
fits_table_acausal = [cells_table(is_responsive,:) fits];
for i=1:numel(P.ap_groups)
    is_in_group = fits_table.AP>=P.ap_groups{i}(1) & fits_table_acausal.AP<P.ap_groups{i}(2);
    fits_table_acausal.ap_group(is_in_group) = i;
end

