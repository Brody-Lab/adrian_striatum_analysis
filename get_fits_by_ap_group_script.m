%% script for breaking down fits by ap group - 3/1/22

P = get_parameters();
[glmfit_log,run] = select_glmfit_runs();
cells_table = select_cells('is_in_dorsal_striatum',true);
[fits,is_responsive] = get_glm_fits(cells_table.recording_name,cells_table.cellno,glmfit_log);
fits_table = [cells_table(is_responsive,:) fits];
for i=1:numel(P.ap_groups)
    is_in_group = fits_table.AP>=P.ap_groups{i}(1) & fits_table.AP<P.ap_groups{i}(2);
    fits_table.ap_group(is_in_group) = i;
end


