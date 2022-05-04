%% script for breaking down fits by ap group - 3/1/22


P = get_parameters();
[glmfit_log,run] = select_glmfit_runs('choice_time_back_s',0.5,'use_trial_history',true,'lambda',30);
cells_table = select_cells('is_in_dorsal_striatum',true);
[fits,is_responsive] = get_glm_fits(cells_table.recording_name,cells_table.cellno,glmfit_log);
fits_table = [cells_table(is_responsive,:) fits];

[~,idx] = sort(fits_table.AP);
fits_table = fits_table(idx,:);

group = equalize_groups(double(fits_table.ap_group));

cluster_glmfits(fits_table,glmfit_log.dm{1}.dspec,'max_clust',19,...
    'responsive_cutoff_prctile',60,'reorder_cells',false,'var_cutoff',1e-3,...
    'find_optimal_n_clust',false,'variable_bar_size',true,'reorder_cells',false,'group_by',group,'group_reorder',true,'covariates_to_cluster',P.covariate_order(1:4),'covariates_to_plot',P.covariate_order(1:4));
